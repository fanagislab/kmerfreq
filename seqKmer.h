#ifndef __SEQ_KMER_H_
#define __SEQ_KMER_H_

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <set>
#include <cmath>
#include <zlib.h>
#include <inttypes.h>
#include <memory.h>
using namespace std;

extern char alphabet[128];
extern char bases[5];
extern char c_bases[5];
extern uint8_t bitAll[8];
extern uint64_t bitLeft[4];

//convert kmer-seq to kmer-bit，64bit最多能装32bp
//需要事先确定序列中只含有ACGT(或acgt)4个碱基
uint64_t seq2bit(string &kseq);

__uint128_t seq2bit128 (string &kseq);
string bit2seq128 (__uint128_t kbit, int kmerSize);

__uint128_t get_rev_com_kbit128(__uint128_t kbit, uint8_t ksize);

//convert kmer-bit to kmer-seq
//此处必须给定kmer的长度值
string bit2seq(uint64_t kbit, int kmerSize);

//check whether a sequence contain non base characters, such as "N"
int check_seq (string &seq);

//get the reverse and complement sequence
void reverse_complement (string &in_str, string &out_str);

//used in routine: make_kmerFreq_2bit_table
//get the frequence value for a given index
int get_freq(uint8_t *freq, uint64_t idx);

//read file_list into a vector
void reading_file_list(string &file_list, vector<string> &files);

//get the reverse and complement kbit, independent of the sequence 
uint64_t get_rev_com_kbit(uint64_t kbit, uint8_t ksize);

//calculate pow for uint64_t integers
uint64_t pow_integer(int base, int exponent);


#endif
