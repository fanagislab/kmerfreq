//This program generate kmer frequence table，output as compressed format (*.cz)
//当KmerSize=17时，将占内存4^17 = 2^34 = 16G
//Author: Fan Wei, fanw@genomics.org.cn
//Date: 2012/5/13
//version: 2.3

//To save compuational time, this program treat N as A (sequencing error) automatically
//So it is better to filter N-contained reads before Kmer frequence analysis

//This is the multiple thread version of kmer_freq program, process reads blocks in parallel.
//读文件，相当于是串行的，用pthread_lock控制，各线程排队从文件中读取一块数据
//将reads转换成kmer,以及将kmer频率更新置内存，都是完全并行的。其中，更新内存，用CAS(compare and swap)技术。
//only support one-line fasta and fastq format
//Change the method of chopping reads into kmers, increase two times of the speed.

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <set>
#include <cmath>
#include <ctime>
#include <inttypes.h>
#include <pthread.h>
#include<zlib.h>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include "seqKmer.h"
#include "gzstream.h"

using namespace std;

//command-line options with default values
int KmerSize = 17;   
string prefix = "output";
int Input_file_format = 1;
int threadNum = 1;
int Whether_output_memory = 0;
uint64_t KmerHeadMaskVal = 0;

//global varaibles used in the main routine
uint64_t Kmer_theory_total = 0;
uint64_t Kmer_total_num = 0;
uint64_t Kmer_effect_num = 0;
uint8_t *KmerFreq;
uint64_t *FreqArray;

//variables and routines in parse the read files parallely
int FileEndMark = 0;
int ReadNumInBlock = 10000;
pthread_mutex_t Block_mutex = PTHREAD_RECURSIVE_MUTEX_INITIALIZER_NP;
void *thread_parseBlock(void* currentFilePointer);
void parse_one_reads_file(string &reads_file);

//variables and routines in combine the forward and reverse strands of Kmers parallely
uint64_t *Kmer_effect_num_pthread;
uint64_t *FreqArray_pthread;
void *thread_combineStrand(void* threadId_p);

//variables and routines in compress the kmer frequency table parallely
uint64_t bufBlockId = 0;
uint64_t bufBlockNum = 0; //calculate in the program
uint64_t bufBlockSize = 8*1024*1024; //8M
vector <uint64_t> bufBlockVec;
uLongf *ComprLen;
Bytef *ComprData;
int *ComprStatus;
void *thread_compressKmerFreq(void* threadId_p);


void usage() 
{	cout << "\n1. Function introduction:\
	\nkmerfreq_reads count mer frequency in sequencing reads data . The forward and reverse format of a k-mer are taken as the same k-mer.  The program store kmer frequency in an array, using the k-mer bit-value as index, so the memory usage is 4^K bytes (K is the size of kmer). For example, when K=17, the memory usage is 16G bytes. \
\n\n2.Input and output\
	\nInput is a file that contains the path of all the sequencing reads files in fastq format or one-line fasta format. The program parse each sequencing reads file by the given order, and make a whole statistics. Output is a table file, which contain the k-mer frequency statistics results. Use a parameter (-m) to decide whether or not output the whole kmer frequency table from the computer memory to the hard disk, which can be used for other downstream applications, such as error correction in de novo assembly. The program run in a multiple-threads way.\
\n\n3. Reference papers:\
	\nA more powerful kmer counting program is jellyfish, which balance the kmer size and memory storage:\
		A fast, lock-free approach for efficient parallel counting of occurrences of k-mers. Marcais G, Kingsford. Bioinformatics.\
		A paper describing how to decide the kmer size for various data set and its effects and cautions:\
		Estimation of genomic characteristics by analyzing kmer frequency in de novo genome projects. Binghang Liu et al. arXiv.\n\n";

	cout << "\nkmerfreq_reads  <reads_files_list>" << endl;
	cout << "   Author Fanwei, fanweiagis@126.com"  << endl;
	cout << "   Version 2.3"  << endl;
	cout << "   -k <int>  set the kmer size (<=19), default=" << KmerSize << endl;
	cout << "   -f <int>   set the input file format: 1: fq|gz(one-line), 2: fa|gz(one-line), default=" << Input_file_format << endl;
	cout << "   -t <int>   run the program in multiple thread mode, default=" << threadNum << endl;
	cout << "   -p <str>  set the output prefix, default=" << prefix << endl;
	cout << "   -m <int>  whether to output computer memory data, 1:yes, 0:no, default=" <<  Whether_output_memory << endl;
	cout << "   -h        get help information" << endl;
	cout << "\nNote: support reading fq format (sequencing data), and one-line fa format (sequencing data)" << endl;
	exit(0);
}


int main(int argc, char *argv[])
{	
	//get options from command line
	int c;
	while((c=getopt(argc, argv, "k:f:t:p:m:h")) !=-1) {
		switch(c) {
			case 'k': KmerSize=atoi(optarg); break;
			case 'f': Input_file_format=atoi(optarg); break;
			case 't': threadNum=atoi(optarg); break;
			case 'p': prefix=optarg; break;
			case 'm': Whether_output_memory=atoi(optarg); break;
			case 'h': usage(); break;
			default: usage();
		}
	}
	if (argc < 2) usage();

	clock_t time_start, time_end;
	time_start = clock();

	string reads_file_list = argv[optind++]; //optind, argv[optind++]顺序指向非option的参数
	KmerHeadMaskVal = pow_integer(2, KmerSize*2) - 1;	

	vector <string> reads_files;
	reading_file_list(reads_file_list, reads_files);

/////////////////////////////////////////////////////////////////////////////////
	
	//parse the read files paralelly and create the kmer frequency table
	for(int i=0; i<KmerSize; i++) {
		Kmer_theory_total=(Kmer_theory_total<<2)|0x3;
	}
	Kmer_theory_total ++;
	
	//动态分配的数组必须赋初值
	KmerFreq = new uint8_t[Kmer_theory_total]; 
	memset(KmerFreq, 0, Kmer_theory_total);

	cerr << "Kmer frequency table initialization completed" << endl;
	time_end = clock();
	cerr << "Run time: " << double(time_end - time_start) / CLOCKS_PER_SEC << endl;

	//parse each read file sequentially, and parse reads blocks parallely
	for (int i=0; i<reads_files.size(); i++)
	{	
		parse_one_reads_file(reads_files[i]);
		cerr << "Finished to parse reads file: " << reads_files[i] << endl;
	}

	cerr << "Parsed all the reads files completed" << endl;
	time_end = clock();
	cerr << "Run time: " << double(time_end - time_start) / CLOCKS_PER_SEC << endl;

/////////////////////////////////////////////////////////////////////////////////


	//This is the new multiple-thread codes to combine forward and reverse strand Kmers parallely
	FreqArray = new uint64_t[256];
	for (int i=0; i<256; i++)
	{	FreqArray[i] = 0;
	}

	pthread_t *pthread = new pthread_t[threadNum];
	int *pthreadId = new int[threadNum];
	Kmer_effect_num_pthread = new uint64_t[threadNum];
	FreqArray_pthread = new uint64_t[256*threadNum];

	if (pthread == NULL || pthreadId == NULL)
	{	cerr << "Out of memory!" << endl;
		exit(0);
	}
	
	for (int i=0; i<threadNum; i++)
	{ /* create children threads */
		Kmer_effect_num_pthread[i] = 0;
		for (int j=0; j<256; j++)
		{	FreqArray_pthread[i*256+j] = 0;
		}
		
		pthreadId[i] = i;
		pthread_create(pthread+i, NULL, thread_combineStrand, (void*)(pthreadId+i));
	}
	cerr << threadNum << " threads were created!" << endl;

	/* wait threads to exit and free resource after threads exit */
	for (int i=0; i<threadNum; i++)
	{ 
		pthread_join(pthread[i], NULL);
	}

	//combine the statistic results
	for (int i=0; i<threadNum; i++)
	{	Kmer_effect_num += Kmer_effect_num_pthread[i];
	}

	for (int i=0; i<256; i++)
	{	for (int j=0; j<threadNum; j++)
		{	FreqArray[i] += FreqArray_pthread[j*256+i];
		}
	}
	
	delete [] Kmer_effect_num_pthread;
	delete [] FreqArray_pthread;

	cerr << "Complete to combine forward and reverse strands of Kmers" << endl;
	time_end = clock();
	cerr << "Run time: " << double(time_end - time_start) / CLOCKS_PER_SEC << endl;
	

/////////////////////////////////////////////////////////////////////////////////

	if (Whether_output_memory == 1)
	{
		//compress parallely with zlib and pthread
		string freq_compress_file = prefix + ".freq.cz";
		ofstream CompDataFile ( freq_compress_file.c_str(), ofstream::binary );
		if ( ! CompDataFile )
		{	cerr << "fail to create freq stat file" << freq_compress_file << endl;
		}
		string freq_comprlen_file = prefix + ".freq.cz.len";
		ofstream  CompLenInfoFile ( freq_comprlen_file.c_str() );
		if ( ! CompLenInfoFile )
		{	cerr << "fail to create freq stat file" << freq_comprlen_file << endl;
		}
		
		ComprLen = new uLongf[threadNum];
		ComprData = new Bytef[bufBlockSize * threadNum];
		ComprStatus = new int[threadNum];
		bufBlockNum = Kmer_theory_total / bufBlockSize;
		for (uint64_t i=0; i<bufBlockNum; i ++)
		{	bufBlockVec.push_back(bufBlockSize);
		}
		uint64_t tail_remain_size = Kmer_theory_total % bufBlockSize;
		if(tail_remain_size > 0)
		{	bufBlockVec.push_back(tail_remain_size);
		}
		bufBlockNum = bufBlockVec.size();

		pthread = new pthread_t[threadNum];
		pthreadId = new int[threadNum];
		
		//creat children threads with all comprstatus set to 0
		for (int i=0; i<threadNum; i++)
		{
			pthreadId[i] = i;
			pthread_create(pthread+i, NULL, thread_compressKmerFreq, (void*)(pthreadId+i));
		}
		cerr << threadNum << " threads were created!" << endl;
		
		bufBlockId = 0;
		while (1)
		{	//initialize variables
			for (int i=0; i<threadNum; i++)
			{	ComprStatus[i] = 1;
				ComprLen[i] = bufBlockSize;
			}
			
			//等子线程，直到ComprStatus都是0，说明当前一批压缩任务已全部完成
			while (1)
			{
				usleep(1);
				int i=0;
				for (; i<threadNum; i++)
				{	if (ComprStatus[i] == 1)
					{	break;
					}
				}
				if (i == threadNum)
				{	break;
				}
			}

			//将压缩结果写入到结果文件中
			for (int i=0; i<threadNum; i++)
			{	if (bufBlockId + i < bufBlockNum)
				{	CompLenInfoFile << ComprLen[i] << "\n";
					CompDataFile.write((const char*)(ComprData+bufBlockSize*i), ComprLen[i]);
				}
			}

			bufBlockId += threadNum;
			
			if (bufBlockId >= bufBlockNum)
			{	for (int i=0; i<threadNum; i++)
				{	ComprStatus[i] = 2;
				}
				break;
			}
		}
		
		//等待全部子线程结束
		for (int i=0; i<threadNum; i++)
		{	pthread_join(pthread[i], NULL);
		}
		
		delete[] ComprLen;
		delete[] ComprData;
		delete[] ComprStatus;
		
		cerr << "Generate the kmer frequency compressed file" << endl;
		time_end = clock();
		cerr << "Run time: " << double(time_end - time_start) / CLOCKS_PER_SEC << endl;
	}

/////////////////////////////////////////////////////////////////////////////////

	//generate the kmer frequency statistics file
	delete[] KmerFreq;
	string freq_stat_file = prefix + ".freq.stat";
	ofstream outfile ( freq_stat_file.c_str() );
	if ( ! outfile )
	{	cerr << "fail to create freq stat file" << freq_stat_file << endl;
	}
	
	outfile << "#KmerSize: " << KmerSize << endl;
	outfile << "#Kmer_theory_max: " << Kmer_theory_total << endl;
	outfile << "#Kmer_indivdual_num: " << Kmer_total_num << endl;
	outfile << "#Kmer_species_num: " << Kmer_effect_num << endl << endl;
	
	outfile << "##Kmer_Frequence\tKmer_Species_Number\tKmer_Species_Ratio\tKmer_Species_accumulate_Ratio\tKmer_Individual_Number\tKmer_Individual_Ratio\tKmer_Individual_accumulate_ratio" << endl;
	
	uint64_t accum_kmer_individual = 0;
	uint64_t accum_kmer_species = 0;
	double accum_kmer_individual_rate = 0.0;
	double accum_kmer_species_rate = 0.0;
	for (uint64_t i=1; i<=254; i++)
	{	
		accum_kmer_species += FreqArray[i];
		accum_kmer_species_rate = (double)accum_kmer_species / (double)Kmer_effect_num;
		accum_kmer_individual += i * FreqArray[i];
		accum_kmer_individual_rate = (double)accum_kmer_individual / (double)Kmer_total_num;
		
		outfile << i << "\t" << FreqArray[i] << "\t" << (double)FreqArray[i] / (double)Kmer_effect_num << "\t" << accum_kmer_species_rate << "\t" << i * FreqArray[i] << "\t" << (double)i * (double)FreqArray[i] / (double)Kmer_total_num << "\t" << accum_kmer_individual_rate << endl;
	}
	outfile << ">=255" << "\t" << FreqArray[255] << "\t" << (double)FreqArray[255] / (double)Kmer_effect_num << "\t" << 1.0 << "\t" << Kmer_total_num - accum_kmer_individual << "\t" << (Kmer_total_num - accum_kmer_individual) / (double)Kmer_total_num << "\t" << 1.0 << endl;
	
	delete [] FreqArray;

	cerr << "Generate the kmer frequency statistics file" << endl;
	time_end = clock();
	cerr << "Run time: " << double(time_end - time_start) / CLOCKS_PER_SEC << endl;
}



//this is the thread routine to parse files parallely
void *thread_parseBlock(void* currentFilePointer)
{
	igzstream *aCurrentFile = (igzstream *)currentFilePointer;

	while (1)
	{
		pthread_mutex_lock(&Block_mutex);
	
		if (FileEndMark == 0)
		{ 
			int readNum = 0;
			vector <string> ReadStrVec;
			string LineStr;

			//step 1: read one block of the input file, 针对fa和fq两种文件格式
			if (Input_file_format == 1)
			{	//读取fq格式文件, support one-line fastq format
				while ( readNum < ReadNumInBlock && getline( *aCurrentFile, LineStr, '\n' ) )
				{	
					if (LineStr[0] == '@') 
					{	
						getline( *aCurrentFile, LineStr, '\n');
						ReadStrVec.push_back(LineStr);
						getline( *aCurrentFile, LineStr, '\n');
						getline( *aCurrentFile, LineStr, '\n');
						readNum ++;
					}	
				}
			}
			else
			{	//读取fa格式文件, support one-line fasta format
				while ( readNum < ReadNumInBlock && getline( *aCurrentFile, LineStr, '\n' ) )
				{	
					if (LineStr[0] == '>') 
					{	
						getline( *aCurrentFile, LineStr, '\n');
						ReadStrVec.push_back(LineStr);
						readNum ++;
					}	
				}
			}
				
			if (readNum < ReadNumInBlock)
			{	FileEndMark = 1;
			}

			pthread_mutex_unlock(&Block_mutex);

////////////////////////////////////////////////////////////////////////////////////////////////			
			
			//step 2: do the further analysis, get the kmers in read
			vector <uint64_t> kmerVec;
			for (uint64_t i=0; i<readNum; i++)
			{	
				string &readi = ReadStrVec[i];
				string kseq_pre = readi.substr(0,KmerSize-1);
				uint64_t kbit = seq2bit(kseq_pre);
				for (int j=KmerSize-1; j<readi.size(); j++)
				{	
					kbit = ((kbit<<2)|alphabet[readi[j]]) & KmerHeadMaskVal;
					kmerVec.push_back(kbit);
				}
			}
			uint64_t block_kmer_num = kmerVec.size();


			//step 3: add the kmers into KmerFreq table, use the gcc built-in sync operations
			 __sync_add_and_fetch(&Kmer_total_num,block_kmer_num);
			for (uint64_t i=0; i<block_kmer_num; i++)
			{	
				int done = 0;
				uint8_t oldVal = 0;
				while ( (oldVal = KmerFreq[kmerVec[i]]) < 255 && !done)
				{	done = __sync_bool_compare_and_swap(KmerFreq+kmerVec[i], oldVal, oldVal + 1);
				}
			}

		}
		else
		{
			pthread_mutex_unlock(&Block_mutex);
			break;
		}
	}

	return NULL;
}

//parse one reads file, invoked by the main function
void parse_one_reads_file(string &reads_file)
{
	FileEndMark = 0;
	
	igzstream currentFile;
	currentFile.open(reads_file.c_str());
	
	//read and parse each block parallely
	pthread_t *pthread = new pthread_t[threadNum];
	int *pthreadId = new int[threadNum];
	
	for (int i=0; i<threadNum; i++)
	{	pthreadId[i] = i;
		pthread_create(pthread+i, NULL, thread_parseBlock, (void*)(&currentFile));
	}
	cerr << threadNum << " threads were created!" << endl;

	/* wait threads to exit and free resource after threads exit */
	for (int i=0; i<threadNum; i++)
	{	pthread_join(pthread[i], NULL);
	}

	//currentFile.close();

	cerr << "##*******" << reads_file << "******" << Kmer_total_num << endl;
}




//this is the thread routine to compress kmer frequency table in parallel
void *thread_compressKmerFreq(void* threadId_p)
{	
	int threadId = *((int*)threadId_p);
		
	while (1)
	{	
		usleep(1);
		
		if (bufBlockId + threadId >= bufBlockNum)
		{	ComprStatus[threadId] = 2;
			return NULL;
		}
		
		if (ComprStatus[threadId] == 1)
		{
			compress(ComprData+bufBlockSize*threadId, &ComprLen[threadId], (const Bytef*)(KmerFreq+(bufBlockId+threadId)*bufBlockSize), bufBlockVec[bufBlockId+threadId]);
			
			ComprStatus[threadId] = 0;
		}

		if (ComprStatus[threadId] == 2)
		{	return NULL;
		}
	}
	
}



//this is the thread routine to combine forward and reverse strand Kmers parallely
//仅当正负链kmer均有frequency时，才合并到kbit小的位置上，否则不管
//多个线程同时操作，互不影响，内存不需加锁
void *thread_combineStrand(void* threadId_p)
{	
	int threadId = *((int*)threadId_p);

	for (uint64_t i = threadId; i<Kmer_theory_total; i += threadNum)
	{	
		if (KmerFreq[i] == 0) { continue; }  //freq为０的跳过
		int add_freq = KmerFreq[i]; //default value
		uint64_t rc_i = get_rev_com_kbit(i, KmerSize);
		
		//仅当正负链kmer均有frequency时，才合并到kbit小的位置上，否则不管
		if ( i < rc_i && KmerFreq[rc_i] > 0 )
		{	add_freq = KmerFreq[i] + KmerFreq[rc_i];
			if(add_freq > 255) { add_freq = 255; }
			KmerFreq[i] = add_freq; 
			KmerFreq[rc_i] = 0; 
		}
		if (i > rc_i && KmerFreq[rc_i] > 0)
		{	continue;  //不处理，不计数
		}
		
		Kmer_effect_num_pthread[threadId] ++;
		FreqArray_pthread[threadId*256+add_freq]++;

	}
	
	return NULL;
}
