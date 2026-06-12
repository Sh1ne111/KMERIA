#ifndef UTIL_H
#define UTIL_H

#define CALLOC(ptr, len) ((ptr) = (__typeof__(ptr))calloc((len), sizeof(*(ptr))))
#define MALLOC(ptr, len) ((ptr) = (__typeof__(ptr))malloc((len) * sizeof(*(ptr))))
#define REALLOC(ptr, len) ((ptr) = (__typeof__(ptr))realloc((ptr), (len) * sizeof(*(ptr))))

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <stdint.h>
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>


typedef unsigned __int128 uint128_t;

static unsigned char seq_nt4_table[256] = { // translate ACGT to 0123
	0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};


static inline uint64_t seq2kmer(char *seq, int k)
{
    uint64_t kmer = 0, krev = 0;
    const uint64_t shift = (k - 1) * 2;
    for(int i = 0; i < k; i++) {
        int c = seq_nt4_table[(uint8_t)seq[i]] & 0x3 ; 
        kmer = (kmer << 2) | c;
        krev = krev >> 2 | (uint64_t)(~c) << shift;
    }
    return kmer < krev ? kmer : krev;
}

static inline void kmer2seq(char * seq, int ksize, uint64_t kmer)
{
    for(int i = 0; i < ksize; i++) {
        seq[i] = "ACGTN"[(kmer >> ((ksize - 1 - i) * 2)) & 0x03];
    }
}

static inline void kmer2seq_64(char *seq, int kmer_size, uint64_t kmer_val) {
    const char *nt4_chars = "ACGT";
    for (int i = 0; i < kmer_size; ++i) {
        seq[kmer_size - 1 - i] = nt4_chars[kmer_val & 0x03];
        kmer_val >>= 2;
    }
    seq[kmer_size] = '\0';
}

static inline void kmer2seq_128(char *seq, int kmer_size, uint128_t kmer_val) {
    const char *nt4_chars = "ACGT";
    for (int i = 0; i < kmer_size; ++i) {
        seq[kmer_size - 1 - i] = nt4_chars[kmer_val & 0x03];
        kmer_val >>= 2;
    }
    seq[kmer_size] = '\0';
}


static inline uint64_t hash(uint64_t key, uint64_t mask) // invertible integer hash function
{
	//return key&mask;
	key = (~key + (key << 21)) & mask; // key = (key << 21) - key - 1;
	key = key ^ key >> 24;
	key = ((key + (key << 3)) + (key << 8)) & mask; // key * 265
	key = key ^ key >> 14;
	key = ((key + (key << 2)) + (key << 4)) & mask; // key * 21
	key = key ^ key >> 28;
	key = (key + (key << 31)) & mask;
	return key;
}

void print_usage();


#endif
