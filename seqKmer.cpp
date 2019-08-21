//author: Fan Wei, email: fanw@genomics.org.cn, date: 2010-2-1
//collect useful C/C++ subroutines related with Kmers
//注意：本文件内的函数和全局变量均可以供外部程序调用
//Date: 2011-7-18
//version: 2.0

#include "seqKmer.h"
#include <inttypes.h>


//the masking values used  by   128bit-reverse-complementary-function
 const __uint128_t mask11 = ((mask11 | 0x3333333333333333LLU) << 64) | 0x3333333333333333LLU;
 const __uint128_t mask12 = ((mask12 | 0xCCCCCCCCCCCCCCCCLLU) << 64) | 0xCCCCCCCCCCCCCCCCLLU;
 const __uint128_t mask21 = ((mask21 | 0x0F0F0F0F0F0F0F0FLLU) << 64) | 0x0F0F0F0F0F0F0F0FLLU;
 const __uint128_t mask22 = ((mask22 | 0xF0F0F0F0F0F0F0F0LLU) << 64) | 0xF0F0F0F0F0F0F0F0LLU;
 const __uint128_t mask31 = ((mask31 | 0x00FF00FF00FF00FFLLU) << 64) | 0x00FF00FF00FF00FFLLU;
 const __uint128_t mask32 = ((mask32 | 0xFF00FF00FF00FF00LLU) << 64) | 0xFF00FF00FF00FF00LLU;
 const __uint128_t mask41 = ((mask41 | 0x0000FFFF0000FFFFLLU) << 64) | 0x0000FFFF0000FFFFLLU;
 const __uint128_t mask42 = ((mask42 | 0xFFFF0000FFFF0000LLU) << 64) | 0xFFFF0000FFFF0000LLU;
 const __uint128_t mask51 = ((mask51 | 0x00000000FFFFFFFFLLU) << 64) | 0x00000000FFFFFFFFLLU;
 const __uint128_t mask52 = ((mask52 | 0xFFFFFFFF00000000LLU) << 64) | 0xFFFFFFFF00000000LLU;


//由ＡＣＧＴ到ASCII码到０１２３，能自动处理大小写
//256个字母表alphabet数组,用8bit的char型存储,A=a=N=n=0,C=c=1,G=g=2,T=t=3,其他的字母都为4
char alphabet[128] =
{
 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
 4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 0, 4, 
 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
 4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 0, 4,
 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4
};

//由０１２３ 4到ＡＣＧＴ N
char bases[5] ={
		'A', 'C', 'G', 'T', 'N'
};

//由０１２３ 4到ＡＣＧＴ N
char c_bases[5] ={
		'T', 'G', 'C', 'A', 'N'
};

//used for bit-operation
uint8_t bitAll[8] = {128,64,32,16,8,4,2,1};

uint64_t bitLeft[4] = {
	0x0000000000000000,
	0x4000000000000000,
	0x8000000000000000,
	0xc000000000000000
};


//convert kmer-seq to kmer-bit，64bit最多能装32bp
//需要事先确定序列中只含有ACGT(或acgt)4个碱基
uint64_t seq2bit(string &kseq)
{
	uint64_t kbit=0;
	for(int i=0; i<kseq.size(); i++) {
			kbit=(kbit<<2)|alphabet[kseq[i]];
	}
	return kbit;
}


__uint128_t seq2bit128 (string &kseq)
{	__uint128_t kbit=0; 
	for(int i=0; i<kseq.size(); i++) 
	{
		kbit=(kbit<<2)|alphabet[kseq[i]];
	}
	return kbit;
}



//convert kmer-bit to kmer-seq
//此处必须给定kmer的长度值
string bit2seq(uint64_t kbit, int kmerSize)
{
	string kseq;
	for(int i=0; i<kmerSize; i++) {
			kseq.push_back(bases[(kbit>>(kmerSize-1-i)*2)&0x3]);
	}
	return kseq;
}

string bit2seq128 (__uint128_t kbit, int kmerSize)
{
	string kseq;
	for(int i=0; i<kmerSize; i++) {	
		kseq.push_back(bases[(kbit>>(kmerSize-1-i)*2)&0x3]);
	}
	return kseq;
}



//check whether a sequence contain non base characters, such as "N"
int check_seq (string &seq)
{       
	int is_good = 1;
	for (int i = 0; i < seq.size(); i++)
	{       if (alphabet[seq[i]] == 4)
			{   is_good = 0;
				break;
			}
	}
	return is_good;
}

//get the reverse and complement sequence
void reverse_complement (string &in_str, string &out_str)
{	
	for (int64_t i=in_str.size()-1; i>=0; i--)
	{	
		out_str.push_back(c_bases[alphabet[in_str[i]]]);
	}
}

//get the reverse and complement kbit, independent of the sequence 
uint64_t get_rev_com_kbit(uint64_t kbit, uint8_t ksize){
	kbit = ~kbit;
	kbit = ((kbit & 0x3333333333333333LLU)<<  2) | ((kbit & 0xCCCCCCCCCCCCCCCCLLU)>>  2);
	kbit = ((kbit & 0x0F0F0F0F0F0F0F0FLLU)<<  4) | ((kbit & 0xF0F0F0F0F0F0F0F0LLU)>>  4);
	kbit = ((kbit & 0x00FF00FF00FF00FFLLU)<<  8) | ((kbit & 0xFF00FF00FF00FF00LLU)>>  8);
	kbit = ((kbit & 0x0000FFFF0000FFFFLLU)<< 16) | ((kbit & 0xFFFF0000FFFF0000LLU)>> 16);
	kbit = ((kbit & 0x00000000FFFFFFFFLLU)<< 32) | ((kbit & 0xFFFFFFFF00000000LLU)>> 32);
	return kbit >> (64 - (ksize<<1));
}


//get the reverse and complement kbit, independent of the sequence 
__uint128_t get_rev_com_kbit128(__uint128_t kbit, uint8_t ksize){
	kbit = ~kbit;
	kbit = ((kbit & mask11)<<  2) | ((kbit & mask12)>>  2);
	kbit = ((kbit & mask21)<<  4) | ((kbit & mask22)>>  4);
	kbit = ((kbit & mask31)<<  8) | ((kbit & mask32)>>  8);
	kbit = ((kbit & mask41)<< 16) | ((kbit & mask42)>> 16);
	kbit = ((kbit & mask51)<< 32) | ((kbit & mask52)>> 32);
	kbit = (kbit << 64) | (kbit>> 64);
	return kbit >> (128 - (ksize<<1));
}


//used in routine: make_kmerFreq_2bit_table
//get the frequence value for a given index
int get_freq(uint8_t *freq, uint64_t idx) 
{
	int value = ( freq[idx/8] >> (7-idx%8) ) & 0x1u;
	return value;
}


//read file_list into a vector
void reading_file_list(string &file_list, vector<string> &files)
{	
	ifstream infile ( file_list.c_str() );
	if ( ! infile )
	{	cerr << "fail to open input file" << file_list << endl;
	}

	string line_str;
	while ( getline( infile, line_str, '\n' ) )
	{	string file_name;
		for(int i=0; i<line_str.size();i++)
		{	if(line_str[i] != ' ' && line_str[i] != '\t' && line_str[i] != '\n')
			{	file_name.push_back(line_str[i]);
			}
		}
		if (file_name.size())
		{	files.push_back(file_name);
		}
	}
}

//chengfang calculation for small integers
uint64_t pow_integer(int base, int exponent)
{       uint64_t result = 1;
        for (int i = 1; i<=exponent; i++)
        {       result *= base;
        }
        return result;
}

		
			
