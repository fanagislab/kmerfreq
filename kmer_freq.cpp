//This program generate kmer frequence table，output as compressed format (*.cz)
//当KmerSize=17时，将占内存4^17 = 2^34 = 16G
//Author: Fan Wei, fanw@genomics.org.cn
//Date: 2012/5/13
//version: 2.3

//To save compuational time, this program treat N as A (sequencing error) automatically
//So it is better to filter N-contained reads before Kmer frequence analysis
//support multi-line fasta and one-line fastq format
//Change the method of chopping reads into kmers, increase two times of the speed.

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <set>
#include <cmath>
#include <inttypes.h>
#include<zlib.h>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include "seqKmer.h"
#include "gzstream.h"

using namespace std;

int KmerSize = 17;   
string prefix = "output";
int Input_file_format = 1; // 1: one-line fastq;  2: multiple-line fasta; 
int Whether_output_memory = 0;
uint64_t bufBlockSize = 8*1024*1024;  //for compression of kmerfreq table
uint64_t KmerHeadMaskVal = 0;

//read the kmers from reads_file_list into memory, store frequency value in 8 bit, range 0-255
uint8_t* make_kmerFreq_8bit_table_from_FaFqFiles(string &seq_file_list, int Ksize, uint64_t &total, uint64_t &num_total_kmers, uint64_t &num_effect_kmers, uint64_t *freqArray);

void usage() 
{	
	cout << "\n1. Function introduction:\
	\nkmerfreq_single count kmer frequency in sequencing data or a reference genome sequence. The forward and reverse format of a k-mer are taken as the same k-mer.  The program store kmer frequency in an array, using the k-mer bit-value as index, so the memory usage is 4^K bytes (K is the size of kmer). For example, when K=17, the memory usage is 16G bytes. \
\n\n2.Input and output\
	\nInput is a file that contains the path of all the genome sequence files in multiple-line fasta format. The program parse each sequence file by the given order, and make a whole statistics. Output is a table file, which contain the k-mer frequency statistics results. Use a parameter (-m) to decide whether or not output the whole kmer frequency table from the computer memory to the hard disk, which can be used for other downstream applications, such as error correction in de novo assembly. The program run in a single thread way.\
\n\n3. Reference papers:\
	\nA more powerful kmer counting program is jellyfish, which balance the kmer size and memory storage:\
		A fast, lock-free approach for efficient parallel counting of occurrences of k-mers. Marcais G, Kingsford. Bioinformatics.\
		A paper describing how to decide the kmer size for various data set and its effects and cautions:\
		Estimation of genomic characteristics by analyzing kmer frequency in de novo genome projects. Binghang Liu et al. arXiv.\n\n";


	cout << "kmerfreq_single  <sequence_files_list>" << endl;
	cout << "   Author Fanwei, fanweiagis@126.com"  << endl;
	cout << "   Version 2.3"  << endl;
	cout << "   -k <int>  set the kmer size (<=19), default=" << KmerSize << endl;
	cout << "   -f <int>   set the input file format: 1: fq|gz(one-line), 2: fa|gz(multi-line), default=" << Input_file_format << endl;
	cout << "   -p <str>  set the output prefix, default=" << prefix << endl;
	cout << "   -m <int>  whether to output computer memory data, 1:yes, 0:no, default=" <<  Whether_output_memory << endl;
	cout << "   -h        get help information" << endl;
	cout << "\nNote: support reading one-line fq[.gz] or multiple-lines fa[.gz] format (genome sequence)\n" << endl;
	exit(0);
}


int main(int argc, char *argv[])
{	


	//get options from command line
	int c;
	while((c=getopt(argc, argv, "k:f:p:m:h")) !=-1) {
		switch(c) {
			case 'k': KmerSize=atoi(optarg); break;
			case 'f': Input_file_format=atoi(optarg); break;
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
	uint8_t *KmerFreq;
	uint64_t *FreqArray = new uint64_t[256];
	uint64_t Kmer_theory_total = 0;
	uint64_t Kmer_total_num = 0;
	uint64_t Kmer_effect_num = 0;

	KmerFreq = make_kmerFreq_8bit_table_from_FaFqFiles(reads_file_list, KmerSize, Kmer_theory_total, Kmer_total_num, Kmer_effect_num, FreqArray);

	cerr << "\nConstruction of Kmer frequency table completed" << endl;
	time_end = clock();
	cerr << "Run time: " << double(time_end - time_start) / CLOCKS_PER_SEC << endl;

///////////////////////////////////////////////////////////////////////////////////
	
	if (Whether_output_memory == 1)
	{
		//use the compress() function in zlib to compress the kmer frequency table 
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
		
		Bytef *ComprData = new Bytef[bufBlockSize];
		uint64_t bufBlockNum = Kmer_theory_total / bufBlockSize;
		vector<uint64_t> bufBlockVec;
		for (uint64_t i=0; i<bufBlockNum; i ++)
		{	bufBlockVec.push_back(bufBlockSize);
		}
		uint64_t tail_remain_size = Kmer_theory_total % bufBlockSize;
		if(tail_remain_size > 0)
		{	bufBlockVec.push_back(tail_remain_size);
		}
		bufBlockNum = bufBlockVec.size();

		for (uint64_t i=0; i<bufBlockNum; i++)
		{	uLongf ComprLen = bufBlockSize;
			compress(ComprData, &ComprLen, (const Bytef*)(KmerFreq+i*bufBlockSize), bufBlockVec[i]);
			CompLenInfoFile << ComprLen << "\n";
			CompDataFile.write((const char*)(ComprData), ComprLen);
		}
		
		delete [] ComprData;

		cerr << "Generate the compressed kmer frequency file" << endl;
		time_end = clock();
		cerr << "Run time: " << double(time_end - time_start) / CLOCKS_PER_SEC << endl;
		
	}

///////////////////////////////////////////////////////////////////////////////////

	//generate the kmer frequency statistics file
	delete [] KmerFreq;

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


//read the kmers from reads_file_list into memory, store frequency value in 8 bit, range 0-255
uint8_t* make_kmerFreq_8bit_table_from_FaFqFiles(string &seq_file_list, int Ksize, uint64_t &total, uint64_t &num_total_kmers, uint64_t &num_effect_kmers, uint64_t *freqArray)
{	
	//得到理论上最大的kmer个数
	total = 0;
	num_total_kmers = 0;
	num_effect_kmers = 0;

	for(int i=0; i<Ksize; i++) {
		total=(total<<2)|0x3;
	}
	total++;
	
	//动态分配的数组必须赋初值
	uint8_t *freq=new uint8_t[total]; 
	memset(freq, 0, total);

	cerr << "KmerFreq initializtion completed" << endl;
	
	vector <string> seq_files;
	reading_file_list(seq_file_list, seq_files);

	for (int i=0; i<seq_files.size(); i++)
	{	
		cerr << "Parse read_file: " << seq_files[i] << endl;

		igzstream infile;
		infile.open(seq_files[i].c_str());
		if ( ! infile )
		{	cerr << "fail to open input file" << seq_files[i] << endl;
		}

		string read, read_head, no_used_line;
		if (Input_file_format == 2)
		{	
			//support multi-lines fasta format
			getline( infile, read_head, '>' );
			while (getline( infile, read_head, '\n' ))
			{	getline( infile, read, '>' );
				boost::erase_all(read, "\n");
				string kseq_pre = read.substr(0,Ksize-1);
				uint64_t kbit = seq2bit(kseq_pre);
				for (int i=Ksize-1; i<read.size(); i++)
				{	
					kbit = ((kbit<<2)|alphabet[read[i]]) & KmerHeadMaskVal;
					if (freq[kbit] < 255)
					{	freq[kbit] ++;
					}
					num_total_kmers++;
				}
			}

		}else
		{	
			//support one-line fastaq format
			while ( getline( infile, read_head, '\n' ) )
			{	
				if (read_head[0] == '@') 
				{	getline( infile, read, '\n' );
					getline( infile, no_used_line, '\n' );
					getline( infile, no_used_line, '\n' );
				}
				
				//if (! check_seq(read)) { continue; }  //pass reads with N

				string kseq_pre = read.substr(0,Ksize-1);
				uint64_t kbit = seq2bit(kseq_pre);
				for (int i=Ksize-1; i<read.size(); i++)
				{	
					kbit = ((kbit<<2)|alphabet[read[i]]) & KmerHeadMaskVal;
					if (freq[kbit] < 255)
					{	freq[kbit] ++;
					}
					num_total_kmers++;
				}
			}
		}
		infile.close();
	}
	
	cerr << "Finished to parse all the read files" << endl;

	//正反链处理，仅当正负链kmer均有frequency时，才合并到kbit小的位置上，否则不管
	//freq非0的位置之和，即为effect_kmer_num, 这样做可以提高压缩率
	for (uint64_t i=0; i<total; i++)
	{	
		if (freq[i] == 0) { continue; }  //freq为０的跳过
		int add_freq = freq[i]; //default value
		uint64_t rc_i = get_rev_com_kbit(i, Ksize);
		
		//注意：仅当正负链kmer均有frequency时，才合并到kbit小的位置上
		if (i < rc_i && freq[rc_i] > 0)
		{	add_freq = freq[i] + freq[rc_i];
			if(add_freq > 255) { add_freq = 255; }
			freq[i] = add_freq; 
			freq[rc_i] = 0; 
		}
		
		//注意： i == rc_i，当KmerSize为偶数时可能发生，这时不做任何处理

		num_effect_kmers ++;
		freqArray[add_freq] ++;
	}

	cerr << "Finished to combine the forward and reverse strand of Kmers" << endl;

	return freq;
}

