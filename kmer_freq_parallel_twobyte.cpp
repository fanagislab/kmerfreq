//当KmerSize=17时，将占内存2* 4^17 = 2^34 = 32G
//Author: Fan Wei, fanweiagis@126.com
//Date: 2015/10/19
//version: 3.0

//To save compuational time, this program treat N as A (sequencing error) automatically
//So it is better to filter N-contained reads before Kmer frequence analysis

//Combine the forward and reverse strand of a kmer, and use the kmer strand with smaller bit-value
//将reads转换成kmer,以及将kmer频率更新置内存，都是完全并行的。其中，更新内存，用CAS(compare and swap)技术。
//only support one-line fasta and fastq format
//Change the method of chopping reads into kmers, increase two times of the speed.
//Store a kmer frequency using two bytes, ranging from 0 to 65535.
//the father thread load reads sequence into memory, while the children threads chop reads into kmers and count frequency

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
int threadNum = 10;
int Whether_output_kmerfreq = 0;
int Output_kmerfreq_cutoff = 10000;
uint64_t KmerHeadMaskVal = 0;
uint64_t KmerRCOrVal[4];
int FreqMax = 65535;

//global varaibles used in the main routine
uint64_t Kmer_theory_total = 0;
uint64_t Kmer_total_num = 0;
uint64_t Kmer_effect_num = 0;
uint64_t Kmer_hifreq_num = 0;

uint16_t *KmerFreq;
uint64_t *FreqArray;

//define entity structure to store the data passed to thread routines
typedef struct {
        void *pointer;
		uint64_t memsize;
		int value;

} THREAD;
void *thread_memset(void* paras);
void *memset_parallel( void *pointer, int value, uint64_t memsize, int threadNum );


//variables and routines in parse the read files parallely
const int BufferNum = 10000;  //max number of reads can be stored in the memeory buffer 
string *RawReads;  //store reads sequence
uint8_t *Signal; //store the status: 0, empty; 1, filled; 2, end of file;
int ReadsNum = 0;   //real number of reads stored in the buffer block
void *thread_parseBlock(void* threadId_p);
void parse_one_reads_file(string &reads_file);
//this is the thread routine to parse files parallely

//variables for compress the memory data in parallel
int Whether_output_memory = 0;
int Low_high_freq_cutoff = 5;
uint64_t bufBlockId = 0;
uint64_t bufBlockNum = 0; //calculate in the program
uint64_t bufBlockSize = 8*1024*1024; //8M kmer number for each compressing block
uLongf *ComprLen;
Bytef *ComprData;   //zlib processed memory block each time must be smaller than 1G bytes
uint8_t *BitData;   //zlib processed memory block each time must be smaller than 1G bytes
void *thread_compressKmerFreq(void* threadId_p);



void usage() 
{	cout << "\n1. Function introduction:\
	\nkmerfreq_16bit count Kmer frequency in sequencing reads data. The forward and reverse strand of a k-mer are taken as the same k-mer, and the kmer strand with smaller bit-value is used to represent the kmer.  The program store kmer frequency in an array of unit 16bit(two bytes, max value 65535), using the k-mer bit-value as index, so the memory usage is 2* 4^K bytes (K is the size of kmer). For example, when K=17, the memory usage is 32G bytes. \
\n\n2.Input and output\
	\nInput is a file that contains the path of all the sequencing reads files in fastq format or one-line fasta format. The program parse each sequencing reads file by the given order, and make a whole statistics. Output is a table file, which contain the k-mer frequency statistics results. Use a parameter (-m) to decide whether or not output the whole kmer frequency table from the computer memory to the hard disk, which can be used for other downstream applications, such as error correction in de novo assembly. The program run in a multiple-threads way.\
\n\n3. Reference papers:\
	\nA more powerful kmer counting program is jellyfish, which balance the kmer size and memory storage:\
		A fast, lock-free approach for efficient parallel counting of occurrences of k-mers. Marcais G, Kingsford. Bioinformatics.\
		A paper describing how to decide the kmer size for various data set and its effects and cautions:\
		Estimation of genomic characteristics by analyzing kmer frequency in de novo genome projects. Binghang Liu et al. arXiv.\n\n";

	cout << "\nkmerfreq_16bit  <reads_files.lib>" << endl;
	cout << "   Author Fanwei, fanweiagis@126.com"  << endl;
	cout << "   Version 2.4"  << endl;
	cout << "   -k <int>  kmer size, range from 12 to 19, default=" << KmerSize << endl;
	cout << "   -f <int>  input file format: 1: fq|gz(one-line), 2: fa|gz(one-line), default=" << Input_file_format << endl;
	cout << "   -t <int>  thread number to use, default=" << threadNum << endl;
	cout << "   -p <str>  output file prefix, default=" << prefix << endl;
	cout << "   -w <int>  whether output kmer sequence and frequency value, , 1:yes, 0:no, default=" <<  Whether_output_kmerfreq << endl;
	cout << "   -c <int>  kmer frequency larger than cutoff will be output, co-used with -w, default=" <<  Output_kmerfreq_cutoff << endl;
	cout << "   -m <int>  whether to output computer memory data(high and low frequency), 1:yes, 0:no, default=" <<  Whether_output_memory << endl;
	cout << "   -q <int>  frequency cutoff to differentiate low and high frequency kmers, co-used with -m,  default=" <<  Low_high_freq_cutoff << endl;
	cout << "   -h        get help information" << endl;
	cout << "\nNote: support reading fq format (sequencing data), and one-line fa format (sequencing data)" << endl;
	cout << "\nExample: kmerfreq_16bit  -k 17 -t 10 -p Ecoli  reads_files.lib\n" << endl;
	exit(0);
}


int main(int argc, char *argv[])
{	
	//get options from command line
	int c;
	while((c=getopt(argc, argv, "k:f:t:p:w:c:m:q:h")) !=-1) {
		switch(c) {
			case 'k': KmerSize=atoi(optarg); break;
			case 'f': Input_file_format=atoi(optarg); break;
			case 't': threadNum=atoi(optarg); break;
			case 'p': prefix=optarg; break;
			case 'w': Whether_output_kmerfreq=atoi(optarg); break;
			case 'c': Output_kmerfreq_cutoff=atoi(optarg); break;
			case 'm': Whether_output_memory=atoi(optarg); break;
			case 'q': Low_high_freq_cutoff=atoi(optarg); break;
			case 'h': usage(); break;
			default: usage();
		}
	}
	if (argc < 2) usage();

	clock_t time_start, time_end;
	time_start = clock();

	if ( KmerSize < 12 || KmerSize > 19)
	{	cerr << "\nError: Kmer size is too small or too large, please set kmer size between 12 to 19\n" << endl;
		exit(0);
	}

	string reads_file_list = argv[optind++]; //optind, argv[optind++]顺序指向非option的参数
	KmerHeadMaskVal = pow_integer(2, KmerSize*2) - 1;	
	
	KmerRCOrVal[3] = 0;
	KmerRCOrVal[1] = pow_integer(2,KmerSize*2-1);
	KmerRCOrVal[2] = pow_integer(2,KmerSize*2-1-1);;
	KmerRCOrVal[0] = KmerRCOrVal[1] + KmerRCOrVal[2];


	vector <string> reads_files;
	reading_file_list(reads_file_list, reads_files);

/////////////////////////////////////////////////////////////////////////////////
	
	//parse the read files paralelly and create the kmer frequency table
	for(int i=0; i<KmerSize; i++) {
		Kmer_theory_total=(Kmer_theory_total<<2)|0x3;
	}
	Kmer_theory_total ++;

	
	//动态分配的数组必须赋初值
	KmerFreq = new uint16_t[Kmer_theory_total]; 
	//memset(KmerFreq, 0, Kmer_theory_total*2);  //memset processing unit is byte
	memset_parallel( KmerFreq, 0, Kmer_theory_total*2, threadNum );
	

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


	//Make statistics of kmer frequency and output kmer and frequency values whose frequency are larger than the frequency cutoff
	FreqArray = new uint64_t[FreqMax+1];
	for (int i=0; i<FreqMax+1; i++)
	{	FreqArray[i] = 0;
	}
	
	string highfreq_kmer_file = prefix + ".highfreq.kmers.gz";
	ogzstream hifile;
	if (Whether_output_kmerfreq == 1){
		hifile.open(highfreq_kmer_file.c_str());
		if ( ! hifile )
		{	cerr << "fail to create file" << highfreq_kmer_file << endl;
		}
	}

	for (uint64_t i = 0; i < Kmer_theory_total; i ++) {
		if ( KmerFreq[i] > 0) {
			Kmer_effect_num ++;
			FreqArray[KmerFreq[i]] ++;
			if ( KmerFreq[i] >= Output_kmerfreq_cutoff && Whether_output_kmerfreq == 1) {
				hifile << bit2seq(i, KmerSize) << "\t" << KmerFreq[i] << "\n";
			}
		}
	}
	hifile.close();

	cerr << "Make statistics of kmer frequency and output kmer and frequency values whose frequency are larger than the frequency cutoff finished" << endl;
	time_end = clock();
	cerr << "Run time: " << double(time_end - time_start) / CLOCKS_PER_SEC << endl;



/////////////////////////////////////////////////////////////////////////////////

	//generate the kmer frequency statistics file
	
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
	for (uint64_t i=1; i<=FreqMax-1; i++)
	{	
		accum_kmer_species += FreqArray[i];
		accum_kmer_species_rate = (double)accum_kmer_species / (double)Kmer_effect_num;
		accum_kmer_individual += i * FreqArray[i];
		accum_kmer_individual_rate = (double)accum_kmer_individual / (double)Kmer_total_num;
		
		outfile << i << "\t" << FreqArray[i] << "\t" << (double)FreqArray[i] / (double)Kmer_effect_num << "\t" << accum_kmer_species_rate << "\t" << i * FreqArray[i] << "\t" << (double)i * (double)FreqArray[i] / (double)Kmer_total_num << "\t" << accum_kmer_individual_rate << endl;
	}
	outfile << ">=65535" << "\t" << FreqArray[FreqMax] << "\t" << (double)FreqArray[FreqMax] / (double)Kmer_effect_num << "\t" << 1.0 << "\t" << Kmer_total_num - accum_kmer_individual << "\t" << (Kmer_total_num - accum_kmer_individual) / (double)Kmer_total_num << "\t" << 1.0 << endl;
	
	delete [] FreqArray;

	cerr << "Generate the kmer frequency statistics file" << endl;
	time_end = clock();
	cerr << "Run time: " << double(time_end - time_start) / CLOCKS_PER_SEC << endl;


/////////////////////////////////////////////////////////////////////////////////
	if (Whether_output_memory == 0)
	{	return 0;
	}
	
	cerr << "begin output memory" << endl;


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
	ComprData = new Bytef[bufBlockSize / 8 * threadNum];
	BitData  = new uint8_t[bufBlockSize / 8 * threadNum];
	bufBlockNum = Kmer_theory_total / bufBlockSize;
	
	
	for ( bufBlockId = 0; bufBlockId < bufBlockNum; bufBlockId += threadNum )
	{	
		//assign 0 values to the buffer memory
		for (int i=0; i<threadNum; i++ )
		{	ComprLen[i] = bufBlockSize / 8;
		}
		memset(ComprData, 0, bufBlockSize / 8 * threadNum);
		memset(BitData, 0, bufBlockSize / 8 * threadNum);

		
		//create and run children threads to parse the reads that have loaded into the memory buffer
		pthread_t *pthread = new pthread_t[threadNum];
		int *pthreadId = new int[threadNum];
		for (int i=0; i<threadNum; i++)
		{
			if (bufBlockId + i < bufBlockNum)
			{
				pthreadId[i] = i;
				pthread_create((pthread+i), NULL, thread_compressKmerFreq, (void*)(pthreadId+i));
			}

		}

		
		//等待全部子线程分析数据结束
		for (int i=0; i<threadNum; i++)
		{
			if (bufBlockId + i < bufBlockNum)
			{
				pthread_join(pthread[i], NULL);
			}
		}

		//output result data to file
		for (int i=0; i<threadNum; i++)
		{	if (bufBlockId + i < bufBlockNum)
			{	CompLenInfoFile << ComprLen[i] << "\n";
				CompDataFile.write((const char*)(ComprData+bufBlockSize/8*i), ComprLen[i]);
			}
		}

	}

	cerr << "finished output memory" << endl;
	cerr << "High frequency kmer number[Kmer_hifreq_num] " << Kmer_hifreq_num << endl;
	
	CompLenInfoFile.close();
	CompDataFile.close();
	delete [] ComprLen;
	delete [] ComprData;
	delete [] BitData;

	delete[] KmerFreq;

}





//this is the thread routine to parse files parallely
void *thread_compressKmerFreq(void* threadId_p)
{
	int threadId = *((int*)threadId_p);

	uint16_t *KmerFreqStart = KmerFreq + (bufBlockId + threadId) * bufBlockSize;
	Bytef *ComprDataStart = ComprData + threadId * bufBlockSize / 8 ;
	uint8_t *BitDataStart = BitData + threadId * bufBlockSize / 8 ;
	
	
	//Convert data from KmerFreq to BitData, [16 bits -> 1 bit]
	for (uint64_t i = 0; i < bufBlockSize; i++ )
	{	if ( KmerFreqStart[i] >= Low_high_freq_cutoff )
		{
			BitDataStart[i/8] |= bitAll[i%8];

			__sync_add_and_fetch(&Kmer_hifreq_num, 1);
			//cout << i << "\t" << get_freq(BitDataStart, i) << "\t" << KmerFreqStart[i] << endl;
		}
	}
	
	//compress the memory in BitData to ComprData by zlib routine,
	//note that the second parameter must be pre-assigned length value, and when the routine finished it change to the new length value after compressed
	compress(ComprDataStart, &ComprLen[threadId], (const Bytef*)(BitDataStart), bufBlockSize / 8);

}



//this is the thread routine to parse files parallely
void *thread_parseBlock(void* threadId_p)
{
	int threadId = *((int*)threadId_p);

	for (uint64_t i = threadId; i < BufferNum; i += threadNum)
	{		
		while (1)
		{
			if (Signal[i] == 1)
			{  
				//get the kmers from a read sequence
				string &readi = RawReads[i];
				vector<uint64_t> kmerVec; 
				
				if (readi.size() < KmerSize)
				{	break;
				}
				
				//for (int j=0; j<10000000; j++) { int m = j + j; int yy = m * j * j * m; }

				string kseq = readi.substr(0,KmerSize);
				uint64_t kbit = seq2bit(kseq);
				uint64_t kbit_rc =  get_rev_com_kbit(kbit, KmerSize);
				uint64_t kbit_final = (kbit <= kbit_rc) ? kbit : kbit_rc;
				kmerVec.push_back(kbit_final);


				for (int j=KmerSize; j<readi.size(); j++)
				{	
					int base_bit = alphabet[readi[j]]; 
					kbit = ((kbit<<2)|base_bit) & KmerHeadMaskVal;
					kbit_rc = (kbit_rc>>2)|KmerRCOrVal[base_bit]; 
					kbit_final = (kbit <= kbit_rc) ? kbit : kbit_rc;
					kmerVec.push_back(kbit_final);
				}
			

				//add the kmers into KmerFreq table, use the gcc built-in sync operations
				uint64_t block_kmer_num = kmerVec.size();
				__sync_add_and_fetch(&Kmer_total_num,block_kmer_num);
				for (uint64_t i=0; i<block_kmer_num; i++)
				{	
					int done = 0;
					uint16_t oldVal = 0;
					while ( (oldVal = KmerFreq[kmerVec[i]]) < FreqMax && !done)
					{	done = __sync_bool_compare_and_swap(KmerFreq+kmerVec[i], oldVal, oldVal + 1);
					}
				}
				
				
				break;
		
			}else
			{
				if (Signal[i] == 0)
				{	usleep(1); // 相当于0.01s, 休息0.01秒，不占计算资源
				}else  
				{	break;   //when signal[i] equals 2, reach the end of file
				}
			}
		}


	}		
			
	return NULL;
}

//parse one reads file, invoked by the main function
void parse_one_reads_file(string &reads_file)
{
	string LineStr;
	igzstream currentFile;  //Caution: gzstream.cpp has a bug, which can be overcomed by gcc/g++ optimation (-O1, -O2, -O3).
	currentFile.open(reads_file.c_str());
	
	RawReads = new string[BufferNum];
	Signal = new uint8_t[BufferNum];

	while (1)
	{
		ReadsNum = 0;
		memset(Signal, 0, BufferNum);  
		
		//create and run children threads to parse the reads that have loaded into the memory buffer
		pthread_t *pthread = new pthread_t[threadNum];
		int *pthreadId = new int[threadNum];
		for (int i=0; i<threadNum; i++)
		{
			pthreadId[i] = i;
			pthread_create((pthread+i), NULL, thread_parseBlock, (void*)(pthreadId+i));

		}
		cerr << threadNum << " threads creation done!" << endl;
		
		//load the reads data into memory buffer by the father thread
		if (Input_file_format == 1)
		{	//读取fq格式文件, support one-line fastq format
			while ( ReadsNum < BufferNum && getline( currentFile, LineStr, '\n' ) )
			{	
				if (LineStr[0] == '@') 
				{	
					getline( currentFile, RawReads[ReadsNum], '\n');			
					getline( currentFile, LineStr, '\n');
					getline( currentFile, LineStr, '\n');
					Signal[ReadsNum] = 1;
					ReadsNum ++;
				}	
			}
		}
		else
		{	//读取fa格式文件, support one-line fasta format
			while ( ReadsNum < BufferNum && getline( currentFile, LineStr, '\n' )  )
			{	
				if (LineStr[0] == '>') 
				{	getline( currentFile, RawReads[ReadsNum], '\n');	
					Signal[ReadsNum] = 1;	
					ReadsNum ++;
				}	
			}
		}

		cerr << "Loading " << ReadsNum << " reads into memory\n" << endl;
		
		//judge the end of file, and make sigle for the children threads
		if ( ReadsNum < BufferNum )
		{
			for (int i = ReadsNum; i < BufferNum; i ++) 
			{	
				Signal[i] = 2;
			}
		}

		//当父线程读取数据结束后，等待全部子线程分析数据结束
		for (int i=0; i<threadNum; i++)
		{
			pthread_join(pthread[i], NULL);
		}
		
		//when reach the end of file, and all children threads finished
		if ( ReadsNum < BufferNum )
		{	break;
		}
		
	
	}


	//currentFile.close();

	cerr << "##*******" << reads_file << "******" << Kmer_total_num << endl;
}



//this is the thread routine to do memset parallely
void *thread_memset(void* paras)
{

	THREAD threpara = *((THREAD*)paras);
	
	memset(threpara.pointer, threpara.value, threpara.memsize); 

}

//the parallel version of memset routine, using pthread functions inside
//void pointer is equivalent to char pointer, the unit is one byte, 以字节为计算单位，void指针就是按照字节计算的
void *memset_parallel( void *pointer, int value, uint64_t memsize, int threadNum )
{
	//Note the memsize uinte is a byte
	uint64_t block_size = memsize / threadNum;
	uint64_t tail_size = memsize % threadNum;

	threadNum += 1;  //needed threads to parse all the divided blocks
	

	pthread_t *pthread = new pthread_t[threadNum];
	for (int i=0; i<threadNum; i++)
	{
		THREAD *threpara = new THREAD;
		threpara->pointer = (char*)pointer + i * block_size;
		threpara->memsize = (i != threadNum - 1) ? block_size : tail_size;
		threpara->value = 0;

		pthread_create( (pthread+i), NULL, thread_memset, (void*)threpara );
	}
	
	//当父线程读取数据结束后，等待全部子线程分析数据结束
	for (int i=0; i<threadNum; i++)
	{
		pthread_join(pthread[i], NULL);
	}

	delete [] pthread;

}


