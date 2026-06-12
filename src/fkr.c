/**
 * @file fetch_reads_with_bgzf.c
 * @brief Fetch reads associated with k-mers from FASTQ data with BGZF compression
 * @author Chen Shuai
 * @date 2024/09/30
 * 
 * Modified version based on:
 * https://github.com/voichek/fetch_reads_with_kmers/blob/master/fetch_reads.cpp
 * 
 */

#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

/* HTSlib BGZF support for compression */
#ifdef USE_BGZF
#include "../include/bgzf.h"


static inline BGZF* bgzf_open_with_level(const char *fn, int level) {
    char mode[8];
    snprintf(mode, sizeof(mode), "w%d", level);
    return bgzf_open(fn, mode);
}

#define OUTPUT_HANDLE BGZF*
#define OUTPUT_OPEN(fn, level) bgzf_open_with_level(fn, level)
#define OUTPUT_CLOSE(fp) bgzf_close(fp)
#define OUTPUT_WRITE_FMT(fp, fmt, ...) \
    do { \
        char _buf[8192]; \
        int _len = snprintf(_buf, sizeof(_buf), fmt, __VA_ARGS__); \
        if (_len > 0 && _len < sizeof(_buf)) { \
            ssize_t _ret = bgzf_write(fp, _buf, _len); \
            if (_ret < 0) { \
                fprintf(stderr, "Warning: BGZF write error\n"); \
            } \
        } \
    } while(0)
#define OUTPUT_SET_THREADS(fp, n) bgzf_mt(fp, n, 256)
#define OUTPUT_SET_COMPRESSION(fp, level) ((void)0) 
#else

static inline gzFile gzopen_with_level(const char *fn, int level) {
    char mode[8];
    snprintf(mode, sizeof(mode), "wb%d", level);
    return gzopen(fn, mode);
}

#define OUTPUT_HANDLE gzFile
#define OUTPUT_OPEN(fn, level) gzopen_with_level(fn, level)
#define OUTPUT_CLOSE(fp) gzclose(fp)
#define OUTPUT_WRITE_FMT(fp, fmt, ...) gzprintf(fp, fmt, __VA_ARGS__)
#define OUTPUT_SET_THREADS(fp, n) ((void)0)
#define OUTPUT_SET_COMPRESSION(fp, level) ((void)0)
#endif



#include "khash.h"
#include "kseq.h"


#define MAX_KMER_LENGTH 32
#define MIN_KMER_LENGTH 1
#define PROGRESS_INTERVAL 1000000
#define INITIAL_HASH_SIZE (1 << 16)
#define MAX_FILENAME_LEN 256
#define COMPRESSION_LEVEL 6
#define NUM_THREADS 4

/* =============================================================================
 * Nucleotide Encoding Table
 * ========================================================================== */

/**
 * @brief Lookup table to convert ASCII nucleotides to 2-bit encoding
 * A=0, C=1, G=2, T=3, N=4 (invalid)
 */
static const unsigned char seq_nt4_table[128] = {
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,  // A=65, C=67, G=71
    4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  // T=84
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,  // a=97, c=99, g=103
    4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4   // t=116
};

/* =============================================================================
 * Hash Function and Hash Table Initialization
 * ========================================================================== */

/**
 * @brief 64-bit hash function to reduce collisions
 * @param key The 64-bit key to hash
 * @return Hashed 64-bit value
 */
static inline uint64_t hash_64(uint64_t key) {
    key = (~key) + (key << 21);
    key = key ^ (key >> 24);
    key = (key + (key << 3)) + (key << 8);
    key = key ^ (key >> 14);
    key = (key + (key << 2)) + (key << 4);
    key = key ^ (key >> 28);
    key = key + (key << 31);
    return key;
}

/* Hash table: k-mer (64-bit integer) -> k-mer name (string) */
KHASH_INIT(kmer_hash, khint64_t, char*, 1, hash_64, kh_int64_hash_equal)

/* Initialize kseq for reading gzipped FASTQ files */
KSEQ_INIT(gzFile, gzread)

/* =============================================================================
 * K-mer Processing Functions
 * ========================================================================== */

/**
 * @brief Insert all k-mers from a sequence into the hash table
 * 
 * @param h Hash table to store k-mers
 * @param k K-mer length
 * @param seq_len Length of the sequence
 * @param seq Nucleotide sequence
 * @param kmer_name Name/identifier for this k-mer
 */
static void insert_kmers_from_sequence(
    khash_t(kmer_hash) *h,
    int k,
    int seq_len,
    const char *seq,
    const char *kmer_name
) {
    uint64_t forward_kmer = 0, reverse_kmer = 0;
    uint64_t mask = (1ULL << (k * 2)) - 1;
    uint64_t shift = (k - 1) * 2;
    int valid_len = 0;
    
    for (int i = 0; i < seq_len; i++) {
        int nucleotide = (uint8_t)seq[i] < 128 ? seq_nt4_table[(uint8_t)seq[i]] : 4;
        
        if (nucleotide < 4) {  // Valid nucleotide (not N)
            forward_kmer = ((forward_kmer << 2) | nucleotide) & mask;
            
            reverse_kmer = (reverse_kmer >> 2) | ((uint64_t)(3 - nucleotide) << shift);
            
            if (++valid_len >= k) {  // Complete k-mer found
                int absent;
                khint_t iter;
                
                iter = kh_put(kmer_hash, h, forward_kmer, &absent);
                if (absent) {
                    kh_val(h, iter) = strdup(kmer_name);
                }
                
                iter = kh_put(kmer_hash, h, reverse_kmer, &absent);
                if (absent) {
                    kh_val(h, iter) = strdup(kmer_name);
                }
            }
        } else {
            valid_len = 0;
            forward_kmer = 0;
            reverse_kmer = 0;
        }
    }
}

/**
 * @brief Look up k-mers in a sequence and return the first match
 * 
 * @param h Hash table containing k-mers
 * @param k K-mer length
 * @param seq_len Length of the sequence
 * @param seq Nucleotide sequence to search
 * @return K-mer name if found, NULL otherwise
 */
static const char* lookup_kmer_in_sequence(
    khash_t(kmer_hash) *h,
    int k,
    int seq_len,
    const char *seq
) {
    uint64_t kmer = 0;
    uint64_t mask = (1ULL << (k * 2)) - 1;
    int valid_len = 0;
    
    for (int i = 0; i < seq_len; i++) {
        int nucleotide = (uint8_t)seq[i] < 128 ? seq_nt4_table[(uint8_t)seq[i]] : 4;
        
        if (nucleotide < 4) {  // Valid nucleotide
            kmer = ((kmer << 2) | nucleotide) & mask;
            
            if (++valid_len >= k) {  // Complete k-mer
                khint_t iter = kh_get(kmer_hash, h, kmer);
                if (iter != kh_end(h)) {
                    return kh_val(h, iter);
                }
            }
        } else {
            valid_len = 0;
            kmer = 0;
        }
    }
    
    return NULL;
}

/* =============================================================================
 * File I/O Functions
 * ========================================================================== */

/**
 * @brief Load k-mers from FASTA file into hash table
 * 
 * @param fasta_file Path to FASTA file containing k-mers
 * @param k K-mer length
 * @return Populated hash table, or NULL on error
 */
static khash_t(kmer_hash)* load_kmers_from_fasta(const char *fasta_file, int k) {
    khash_t(kmer_hash) *h = kh_init(kmer_hash);
    kh_resize(kmer_hash, h, INITIAL_HASH_SIZE);
    
    gzFile fp = gzopen(fasta_file, "r");
    if (!fp) {
        fprintf(stderr, "Error: Cannot open k-mer file: %s\n", fasta_file);
        kh_destroy(kmer_hash, h);
        return NULL;
    }
    
    kseq_t *ks = kseq_init(fp);
    uint64_t kmer_count = 0;
    
    while (kseq_read(ks) >= 0) {
        insert_kmers_from_sequence(h, k, ks->seq.l, ks->seq.s, ks->name.s);
        kmer_count++;
    }
    
    fprintf(stderr, "[K-mer Loading] Loaded %llu k-mers from %s\n", 
            (unsigned long long)kmer_count, fasta_file);
    
    kseq_destroy(ks);
    gzclose(fp);
    
    return h;
}

/**
 * @brief Write a FASTQ record to compressed output
 * 
 * @param out Output file handle (BGZF or gzFile)
 * @param kmer_name K-mer name (prefix), can be NULL
 * @param read_name Read name/identifier
 * @param seq Sequence string
 * @param qual Quality string
 */
static void write_fastq_record(
    OUTPUT_HANDLE out,
    const char *kmer_name,
    const char *read_name,
    const char *seq,
    const char *qual
) {
    if (kmer_name != NULL) {
        OUTPUT_WRITE_FMT(out, "@%s_%s\n%s\n+\n%s\n", kmer_name, read_name, seq, qual);
    } else {
        OUTPUT_WRITE_FMT(out, "@%s\n%s\n+\n%s\n", read_name, seq, qual);
    }
}


/**
 * @brief Process paired-end FASTQ files and extract reads containing k-mers
 * 
 * @param h Hash table containing k-mers
 * @param k K-mer length
 * @param r1_file Path to R1 FASTQ file
 * @param r2_file Path to R2 FASTQ file
 * @param output_prefix Output file prefix
 * @return Number of read pairs found
 */
static uint64_t process_fastq_files(
    khash_t(kmer_hash) *h,
    int k,
    const char *r1_file,
    const char *r2_file,
    const char *output_prefix
) {
    gzFile fp_r1 = gzopen(r1_file, "r");
    gzFile fp_r2 = gzopen(r2_file, "r");
    
    if (!fp_r1 || !fp_r2) {
        fprintf(stderr, "Error: Cannot open FASTQ files\n");
        if (fp_r1) gzclose(fp_r1);
        if (fp_r2) gzclose(fp_r2);
        return 0;
    }
    
    kseq_t *ks_r1 = kseq_init(fp_r1);
    kseq_t *ks_r2 = kseq_init(fp_r2);
    
    // Open compressed output files
    char r1_out_fn[MAX_FILENAME_LEN];
    char r2_out_fn[MAX_FILENAME_LEN];
    snprintf(r1_out_fn, MAX_FILENAME_LEN, "%s_R1.fastq.gz", output_prefix);
    snprintf(r2_out_fn, MAX_FILENAME_LEN, "%s_R2.fastq.gz", output_prefix);
    
    //OUTPUT_HANDLE out_r1 = OUTPUT_OPEN(r1_out_fn);
    //OUTPUT_HANDLE out_r2 = OUTPUT_OPEN(r2_out_fn);
    OUTPUT_HANDLE out_r1 = OUTPUT_OPEN(r1_out_fn, COMPRESSION_LEVEL);
    OUTPUT_HANDLE out_r2 = OUTPUT_OPEN(r2_out_fn, COMPRESSION_LEVEL);    
    if (!out_r1 || !out_r2) {
        fprintf(stderr, "Error: Cannot create output files\n");
        kseq_destroy(ks_r1);
        kseq_destroy(ks_r2);
        gzclose(fp_r1);
        gzclose(fp_r2);
        if (out_r1) OUTPUT_CLOSE(out_r1);
        if (out_r2) OUTPUT_CLOSE(out_r2);
        return 0;
    }
    
    // Set compression parameters
    OUTPUT_SET_COMPRESSION(out_r1, COMPRESSION_LEVEL);
    OUTPUT_SET_COMPRESSION(out_r2, COMPRESSION_LEVEL);
    
    OUTPUT_SET_THREADS(out_r1, NUM_THREADS);
    OUTPUT_SET_THREADS(out_r2, NUM_THREADS);
    
    uint64_t total_reads = 0;
    uint64_t matched_reads = 0;
    
#ifdef USE_BGZF
    fprintf(stderr, "[Read Processing] Starting with BGZF parallel compression (%d threads)...\n", NUM_THREADS);
#else
    fprintf(stderr, "[Read Processing] Starting with standard gzip compression...\n");
#endif
    
    while (kseq_read(ks_r1) >= 0 && kseq_read(ks_r2) >= 0) {
        total_reads++;
        
        if (total_reads % PROGRESS_INTERVAL == 0) {
            fprintf(stderr, "[Progress] Processed %llu reads, found %llu matches (%.2f%%)\n",
                    (unsigned long long)total_reads, 
                    (unsigned long long)matched_reads,
                    total_reads > 0 ? (100.0 * matched_reads / total_reads) : 0.0);
        }
        
        const char *kmer_name_r1 = lookup_kmer_in_sequence(h, k, ks_r1->seq.l, ks_r1->seq.s);
        const char *kmer_name_r2 = lookup_kmer_in_sequence(h, k, ks_r2->seq.l, ks_r2->seq.s);
        
        if (kmer_name_r1 != NULL || kmer_name_r2 != NULL) {
            matched_reads++;
            write_fastq_record(out_r1, kmer_name_r1, ks_r1->name.s, ks_r1->seq.s, ks_r1->qual.s);
            write_fastq_record(out_r2, kmer_name_r2, ks_r2->name.s, ks_r2->seq.s, ks_r2->qual.s);
        }
    }
    
    kseq_destroy(ks_r1);
    kseq_destroy(ks_r2);
    gzclose(fp_r1);
    gzclose(fp_r2);
    OUTPUT_CLOSE(out_r1);
    OUTPUT_CLOSE(out_r2);
    
    fprintf(stderr, "[Summary] Total reads: %llu, Matched reads: %llu (%.2f%%)\n",
            (unsigned long long)total_reads,
            (unsigned long long)matched_reads,
            total_reads > 0 ? (100.0 * matched_reads / total_reads) : 0.0);
    
    return matched_reads;
}


/**
 * @brief Main function for fetching reads with k-mers
 */
int fetch_reads(int argc, char *argv[]) {
    // Print usage
    if (argc != 6) {
        fprintf(stderr, "\n");
        fprintf(stderr, "======================================================\n");
#ifdef USE_BGZF
        fprintf(stderr, "  K-mer Based Read Fetcher\n");
#else
        fprintf(stderr, "  K-mer Based Read Fetcher\n");
#endif
        fprintf(stderr, "======================================================\n");
        fprintf(stderr, "Function: Fetch reads associated with k-mers from FASTQ data\n");
#ifdef USE_BGZF
        fprintf(stderr, "          Outputs BGZF compressed files (parallel gzip)\n\n");
#else
        fprintf(stderr, "          Outputs gzip compressed files\n\n");
#endif
        fprintf(stderr, "Based on: https://github.com/voichek/fetch_reads_with_kmers\n");
        fprintf(stderr, "Author:   Chen Shuai\n");
        fprintf(stderr, "Date:     2024/09/30\n\n");
        fprintf(stderr, "Usage:\n");
        fprintf(stderr, "  %s <R1.fq.gz> <R2.fq.gz> <kmers.fasta> <k> <output_prefix>\n\n", argv[0]);
        fprintf(stderr, "Arguments:\n");
        fprintf(stderr, "  R1.fq.gz       - Forward reads (gzipped FASTQ)\n");
        fprintf(stderr, "  R2.fq.gz       - Reverse reads (gzipped FASTQ)\n");
        fprintf(stderr, "  kmers.fasta    - K-mers in FASTA format\n");
        fprintf(stderr, "  kmer_size      - K-mer length (2-31)\n");
        fprintf(stderr, "  output_prefix  - Output file prefix\n\n");
        fprintf(stderr, "Output:\n");
        fprintf(stderr, "  <prefix>_R1.fastq.gz\n");
        fprintf(stderr, "  <prefix>_R2.fastq.gz\n\n");
//        fprintf(stderr, "Compilation:\n");
//#ifdef USE_BGZF
//        fprintf(stderr, "  gcc -O3 -DUSE_BGZF -o fetch_reads %s -lz -lhts -lpthread\n", argv[0]);
//#else
//        fprintf(stderr, "  gcc -O3 -o fetch_reads %s -lz\n", argv[0]);
//#endif
        fprintf(stderr, "======================================================\n\n");
        return 1;
    }
    
    const char *r1_file = argv[1];
    const char *r2_file = argv[2];
    const char *kmer_file = argv[3];
    int k = atoi(argv[4]);
    const char *output_prefix = argv[5];
    
    if (k < MIN_KMER_LENGTH || k > MAX_KMER_LENGTH) {
        fprintf(stderr, "Error: K-mer length must be between %d and %d\n",
                MIN_KMER_LENGTH, MAX_KMER_LENGTH);
        return 1;
    }
    
    fprintf(stderr, "[Configuration] K-mer length: %d\n", k);
#ifdef USE_BGZF
    fprintf(stderr, "[Configuration] Compression: BGZF (level %d, %d threads)\n", 
            COMPRESSION_LEVEL, NUM_THREADS);
#else
    fprintf(stderr, "[Configuration] Compression: gzip (level %d)\n", COMPRESSION_LEVEL);
#endif
    
    khash_t(kmer_hash) *h = load_kmers_from_fasta(kmer_file, k);
    if (!h) {
        return 1;
    }
    
    uint64_t matched_reads = process_fastq_files(h, k, r1_file, r2_file, output_prefix);
    
    for (khint_t i = 0; i < kh_end(h); i++) {
        if (kh_exist(h, i)) {
            free(kh_val(h, i));
        }
    }
    kh_destroy(kmer_hash, h);
    
    fprintf(stderr, "[Complete] Successfully processed %llu read pairs\n",
            (unsigned long long)matched_reads);
    
    return 0;
}
