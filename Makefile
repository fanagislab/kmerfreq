#all: kmerfreq_single  kmerfreq_parallel kmerfreq_16bit  kmerfreq
all: kmerfreq

kmerfreq_single: kmer_freq.cpp seqKmer.cpp gzstream.cpp
	g++ -O3  -o $@  $^ -lz

kmerfreq_parallel: kmer_freq_parallel.cpp seqKmer.cpp gzstream.cpp
	g++  -O3  -o $@  $^ -lz -lpthread

##For this program, g++ optimiztion is essential to overcome the bug in igzstream
kmerfreq_16bit: kmer_freq_parallel_twobyte.cpp seqKmer.cpp gzstream.cpp
	g++ -O3   -o $@  $^ -lz -lpthread

##For this program, g++ optimiztion is essential to overcome the bug in igzstream
kmerfreq: kmer_freq_parallel_twobyte_publish.cpp seqKmer.cpp gzstream.cpp
	g++ -O3   -o $@  $^ -lz -lpthread



clean:
	#rm kmerfreq_single kmerfreq_parallel kmerfreq_16bit kmerfreq
	rm kmerfreq	
