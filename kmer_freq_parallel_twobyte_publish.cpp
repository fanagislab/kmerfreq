
//Author: Fan Wei, fanweiagis@126.com
//Date: 2019/1/2
//version: 4.0

//Memory usage: Store a kmer frequency using two bytes, ranging from 0 to 65535. 当KmerSize=17时，将占内存2* 4^17 = 2^34 = 32G

//Parallel style: the disk reading (main thread) and kmer computation (children threads) are performed simutaneously, use CAS(compare and swap) lock-free approach in the children threads, the loading data from disk and Kmer computation are doing simutaneously.

//Data pre-filter: To save compuational time, this program treat N as A (sequencing error) automatically, So it is better to filter N-contained reads before Kmer frequence analysis

//Strand combination: Combine the forward and reverse strand of a kmer, and use the kmer strand with smaller bit-value

//The advanced method of chopping reads into kmers, increase two times of the speed comparing to the original method.

//Note that the gzipped memory data can be ungzipped parallely by routine make_kmerFreq_1bit_table_from_1BitGz_pthread() firstly written for correct_error/main_parallel_senior.cpp


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

//global varaibles for commond-line options
int KmerSize = 17;   
string prefix;
int Input_file_format = 1;
int threadNum = 10;

int Whether_output_kmerfreq = 0;
int Output_kmerfreq_cutoff = 5;
uint64_t nonLow_Kmer_species_num = 0;

//global varaibles for the major compution 
uint64_t KmerHeadMaskVal = 0;  //used in generating a new kmer bit data from previous kmer bit data
uint64_t KmerRCOrVal[4];  //used in generating a new kmer bit data from previous kmer bit data
int FreqMax = 65535;    //the max value of kmer frequency for using 16-bit storage.

uint64_t Kmer_theory_total = 0;  //the max theoretic number of kmer species with given KmerSize
uint64_t Kmer_total_num = 0;  //the total number of kmer individuals in the input data
uint64_t Kmer_effect_num = 0;  //the total number of kmer species in the input data


uint16_t *KmerFreq;   //store frequency values of each kmer
uint64_t *FreqArray;  //store the statistis distribution of kmer frequency

//define entity structure to store the data passed to thread routines
typedef struct {
        void *pointer;
		uint64_t memsize;
		int value;

} THREAD;
void *thread_memset(void* paras);  //this is the thread routine to do memset parallely
void *memset_parallel( void *pointer, int value, uint64_t memsize, int threadNum );  //the parallel version of memset routine, using pthread functions inside


//variables and routines in parsing the read files parallely
int BufferNum = 10000;  //max number of reads can be stored in the memeory buffer 
string *RawReads;  //store reads sequence
uint8_t *Signal; //store the status: 0, empty; 1, filled; 2, end of file;
int ReadsNum = 0;   //real number of reads stored in the buffer block
void *thread_parseBlock(void* threadId_p);
void parse_one_reads_file(string &reads_file);
//this is the thread routine to parse files parallely

//variables for compress the memory data in parallel
int Whether_output_memory = 0;
uint64_t Kmer_hifreq_num = 0;  //the total number of high-frequency kmer species in the input data
int Low_high_freq_cutoff = 5;   //cutoff to define low and non-low frequency kmer

uint64_t bufBlockId = 0;
uint64_t bufBlockNum = 0; //calculate in the program
uint64_t bufBlockSize = 8*1024*1024; //8M kmer species space number for each compressing block
uLongf *ComprLen;
Bytef *ComprData;   //zlib processed memory block each time must be smaller than 1G bytes
uint8_t *BitData;   //zlib processed memory block each time must be smaller than 1G bytes
void *thread_compressKmerFreq(void* threadId_p);



void usage() 
{	cout << "\nFunction introduction:\
	\nkmerfreq count K-mer (with size K) frequency from the input sequence data, typically sequencing reads data, and reference genome data is also applicable. The forward and reverse strand of a k-mer are taken as the same k-mer, and only the kmer strand with smaller bit-value is used to represent the kmer. It adopts a 16-bit integer with max value 65535 to store the frequency value of a unique K-mer, and any K-mer with frequency larger than 65535 will be recorded as 65535. The program store all kmer frequency values in a 4^K size array of 16-bit integer (2 bytes), using the k-mer bit-value as index, so the total memory usage is 2* 4^K bytes. For K-mer size 15, 16, 17, 18, 19, it will consume constant 2G, 8G 32G 128G 512G memory, respectively. kmerfreq works in a highly simple and parallel style, to achieve as fast speed as possible.  \
\n\nInput and output:\
	\nInput is a library file that contains the path of all the input sequence files in fastq format or one-line fasta format, with each line represents the path of a single sequence file. The program parse each sequencing reads file by the given order, and make a whole statistics. Output by default is a table file for the distribution of the kmer frequency, which can be analyzed for estimating sequencing quality and genomic characteristics. Two other optional output files are also avialable. By setting parameter \"-w 1\", the Kmer sequence with corresponding frequency values will be output for those K-mers with frequency value >= cutoff set by parameter -c; This is a human readable file that can facilate other kmer associated analysis. By setting parameter \"-m 1\", the whole kmer frequency table from the computer memory will be compressed and output, using 1-bit for each unique kmer, with 0 for kmer with frequency lower than cutoff set by parameter -q, and 1 for other kmers. This is a binary file, which can be reloaded to memory by other program that uses the kmer frequency data as input.\
\n\nReference papers:\
	\nA paper describing how to decide the kmer size for various data set and the potential applications for kmers:\
	\nBinghang Liu, Yujian Shi, Jianying Yuan, et al. and Wei Fan*. Estimation of genomic characteristics by analyzing k-mer frequency in de novo genome project. arXiv.org arXiv: 1308.2012. (2013)\n\n";

	cout << "\nkmerfreq  [options] <reads_files.lib>" << endl;
	cout << "   Author Wei Fan, fanweiagis@126.com"  << endl;
	cout << "   Version 4.0"  << endl;
	cout << "   -k <int>  kmer size, recommand value 13 to 19, default=" << KmerSize << endl;
	cout << "   -f <int>  input file format: 1: fq|gz(one-line), 2: fa|gz(one-line), default=" << Input_file_format << endl;
	cout << "   -p <str>  output file prefix, default=" << "reads_files.lib" << endl;
	cout << "   -r <int>  number of reads stored in buffer memory, default=" << BufferNum << endl;
	cout << "   -t <int>  thread number to use in parallel, default=" << threadNum << endl;
	cout << "   -w <int>  whether output kmer sequence and frequency value, , 1:yes, 0:no, default=" <<  Whether_output_kmerfreq << endl;
	cout << "   -c <int>  kmer frequency cutoff, equal or larger will be output, co-used with -w, default=" <<  Output_kmerfreq_cutoff << endl;
	cout << "   -m <int>  whether output computer memory data, 1:yes, 0:no, default=" <<  Whether_output_memory << endl;
	cout << "   -q <int>  kmer frequency cutoff, 0 for lower, 1 for equal and larger, co-used with -m,  default=" <<  Low_high_freq_cutoff << endl;
	cout << "   -h        get help information\n" << endl;
	
	cout << "Example: kmerfreq  reads_files.lib" << endl;
	cout << "         kmerfreq  -k 17 -t 10 -p Ecoli_K17 reads_files.lib" << endl;
	cout << "         kmerfreq  -k 17 -t 10 -p Ecoli_K17 -w 1 -c 5 reads_files.lib" << endl;
	cout << "         kmerfreq  -k 17 -t 10 -p Ecoli_K17 -m 1 -q 5 reads_files.lib\n" << endl;
	exit(0);
}


int main(int argc, char *argv[])
{	
	//get options from command line
	int c;
	while((c=getopt(argc, argv, "k:f:t:r:p:w:c:m:q:h")) !=-1) {
		switch(c) {
			case 'k': KmerSize=atoi(optarg); break;
			case 'f': Input_file_format=atoi(optarg); break;
			case 't': threadNum=atoi(optarg); break;
			case 'r': BufferNum=atoi(optarg); break;
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

	string reads_file_list = argv[optind++]; //optind, argv[optind++]顺序指向非option的参数
	
	if (prefix.size() == 0)
	{	prefix = reads_file_list;
	}

	//used in generating a new kmer bit data from previous kmer bit data
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
	

	cerr << "\nKmer frequency table initialization completed" << endl;
	time_end = clock();
	cerr << "CPU time: " << double(time_end - time_start) / CLOCKS_PER_SEC << " seconds" << endl << endl;

	//parse each read file sequentially, and parse reads blocks parallely
	for (int i=0; i<reads_files.size(); i++)
	{	
		parse_one_reads_file(reads_files[i]);
		cerr << "Finished to parse reads file: " << reads_files[i] << endl;
	}

	cerr << "\nParsing all the reads files completed" << endl;
	time_end = clock();
	cerr << "CPU time: " << double(time_end - time_start) / CLOCKS_PER_SEC << " seconds" << endl  << endl;

/////////////////////////////////////////////////////////////////////////////////


	//Make statistics of kmer frequency and output kmer and frequency values whose frequency are larger than the frequency cutoff
	FreqArray = new uint64_t[FreqMax+1];
	for (int i=0; i<FreqMax+1; i++)
	{	FreqArray[i] = 0;
	}
	
	string nonLowFreq_kmer_file = prefix + ".nonLow.kmer.freq.gz";
	ogzstream hifile;
	if (Whether_output_kmerfreq == 1){
		hifile.open(nonLowFreq_kmer_file.c_str());
		if ( ! hifile )
		{	cerr << "fail to create file" << nonLowFreq_kmer_file << endl;
		}
	}
	
	
	for (uint64_t i = 0; i < Kmer_theory_total; i ++) {
		if ( KmerFreq[i] > 0) {
			Kmer_effect_num ++;
			FreqArray[KmerFreq[i]] ++;
			if ( KmerFreq[i] >= Output_kmerfreq_cutoff && Whether_output_kmerfreq == 1) {
				hifile << bit2seq(i, KmerSize) << "\t" << KmerFreq[i] << "\n";
				nonLow_Kmer_species_num ++;
			}
		}
	}
	hifile.close();

	cerr << "Make statistics of kmer frequency done" << endl;

	if (Whether_output_kmerfreq == 1)
	{	cerr << "Output nonLow kmer sequence and frequency values into " << nonLowFreq_kmer_file << endl;
		cerr << "Number of Kmer species with frequency >= " << Output_kmerfreq_cutoff << ": " << nonLow_Kmer_species_num << endl;
	}

	time_end = clock();
	cerr << "CPU time: " << double(time_end - time_start) / CLOCKS_PER_SEC << " seconds" << endl << endl;



/////////////////////////////////////////////////////////////////////////////////

	//generate the kmer frequency statistics file
	
	string freq_stat_file = prefix + ".kmer.freq.stat";
	ofstream outfile ( freq_stat_file.c_str() );
	if ( ! outfile )
	{	cerr << "fail to create freq stat file" << freq_stat_file << endl;
	}
	
	outfile << "#Kmer size: " << KmerSize << endl;
	outfile << "#Maximum Kmer frequency: " << FreqMax << endl;
	outfile << "#Kmer indivdual number: " << Kmer_total_num << endl;
	outfile << "#Kmer species number: " << Kmer_effect_num << endl;
	outfile << "#Theoretic space of Kmer species: " << Kmer_theory_total << "  occupied ratio: " << (double)Kmer_effect_num / (double)Kmer_theory_total << endl << endl;

	outfile << "#Kmer_Frequency\tKmer_Species_Number\tKmer_Species_Ratio\tKmer_Species_accumulate_Ratio\tKmer_Individual_Number\tKmer_Individual_Ratio\tKmer_Individual_accumulate_ratio" << endl;
	
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
	outfile << "65535" << "\t" << FreqArray[FreqMax] << "\t" << (double)FreqArray[FreqMax] / (double)Kmer_effect_num << "\t" << 1.0 << "\t" << Kmer_total_num - accum_kmer_individual << "\t" << (Kmer_total_num - accum_kmer_individual) / (double)Kmer_total_num << "\t" << 1.0 << endl;
	
	delete [] FreqArray;

	cerr << "Generate the kmer frequency statistics file" << endl;
	time_end = clock();
	cerr << "CPU time: " << double(time_end - time_start) / CLOCKS_PER_SEC << " seconds" << endl  << endl;


/////////////////////////////////////////////////////////////////////////////////
	if (Whether_output_memory == 0)
	{	delete[] KmerFreq;
		return 0;
	}
	
	//first convert 16-bit frequency value to 1-bit value, and second compress the 1-bit value data structure by zlib and pthread parallely
	//in the resulting 1-bit data structure, 0 represents low-frequency kmer, and 1 represent high-frequency kmer
	//Note that the gzipped data can be ungzipped parallely by routine make_kmerFreq_1bit_table_from_1BitGz_pthread() firstly written for correct_error/main_parallel_senior.cpp
	string freq_compress_file = prefix + ".kmer.freq.cz";
	ofstream CompDataFile ( freq_compress_file.c_str(), ofstream::binary );
	if ( ! CompDataFile )
	{	cerr << "fail to create freq stat file" << freq_compress_file << endl;
	}
	string freq_comprlen_file = prefix + ".kmer.freq.cz.len";
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

		
		//create and run children threads 
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

	CompLenInfoFile.close();
	CompDataFile.close();
	delete [] ComprLen;
	delete [] ComprData;
	delete [] BitData;

	delete[] KmerFreq;

	cerr << "Finished output memory into " << freq_compress_file << endl;
	cerr << "Number of Kmer species with frequency >= " << Low_high_freq_cutoff << ": " << Kmer_hifreq_num << endl;
	time_end = clock();
	cerr << "CPU time: " << double(time_end - time_start) / CLOCKS_PER_SEC << " seconds" << endl  << endl;


}





//first convert 16-bit frequency value to 1-bit value, and second compress the 1-bit value data structure by zlib and pthread parallely
//in the resulting 1-bit data structure, 0 represents low-frequency kmer, and 1 represent high-frequency kmer
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



//this is the thread routine to parse reads data stored in the buffer-memory
//use sigal to control the children thread working: 0 means reads not loaded and should wait; 1 means reads loaded and should do; 2 means nothing need to do and should break out
//with the help of sigal for each reads in the buffer, the children threads don't need the whole buffer to be fully loaded, and will start to work immediately when some reads being loaded
//the children threads works with the main thread (reading file from disk to memory buffer) simutaneously to achieve the fastest speed
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
		//cerr << threadNum << " threads creation done!" << endl;
		
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

		//cerr << "Loading " << ReadsNum << " reads into memory\n" << endl;
		
		//judge the end of file, and make sigle for the children threads
		//For the last block of reading an input file, the buffer region will not be fully loaded
		//So we should give signal value 2 to the left empty positions in order to tell the children threads to do nothing for these positions
		if ( ReadsNum < BufferNum )
		{
			//将没有填充reads数据的buffer区域的signal设置为2，这样子线程可以什么也不用做了
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

	//cerr << "##*******" << reads_file << "******" << Kmer_total_num << endl;
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


