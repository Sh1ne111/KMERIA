/*
   The MIT License

   Copyright (c) Chen Shuai, 2024- <chensss1209@gmail.com>

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#include "khash.h"
#include "kseq.h"

#define PROG_NAME    "Extract_kmer_read_from_tgs"
#define PROG_VERSION "1.0.0"
#define PROG_AUTHOR  "Chen Shuai"
#define PROG_DATE    "2025-06-12"

#define MIN_KMER_LEN 17
#define MAX_KMER_LEN 31
#define INITIAL_MATCH_SIZE 4
#define MAX_NAME_LEN 1024
#define MAX_PREFIX_LEN 800
#define PROGRESS_INTERVAL 100000
#define HASH_INIT_SIZE (1<<16)

static unsigned char seq_nt4_table[128] = { 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

static inline uint64_t hash_64(uint64_t key)
{ 
    key = (~key + (key << 21));
    key = key ^ key >> 24;
    key = ((key + (key << 3)) + (key << 8));
    key = key ^ key >> 14;
    key = ((key + (key << 2)) + (key << 4));
    key = key ^ key >> 28;
    key = (key + (key << 31));
    return key;
}

KHASH_MAP_INIT_STR(kmer_name, char*)
KHASH_INIT(64, khint64_t, char*, 1, hash_64, kh_int64_hash_equal)

KSEQ_INIT(gzFile, gzread)

static void print_usage(const char *prog)
{
    fprintf(stderr, "\n");
    fprintf(stderr, "====================================================================\n");
    fprintf(stderr, " %s v%s\n", PROG_NAME, PROG_VERSION);
    fprintf(stderr, "====================================================================\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Author: %s\n", PROG_AUTHOR);
    fprintf(stderr, "Date:   %s\n", PROG_DATE);
    fprintf(stderr, "Description:\n");
    fprintf(stderr, "  Extract reads containing specific k-mers from TGS (Third Generation\n");
    fprintf(stderr, "  Sequencing) long-read data. Supports FASTQ format.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage:\n");
    fprintf(stderr, "  %s <long_reads.fq.gz> <asso_kmers.fasta> <kmer_length> <output_prefix>\n", prog);
    fprintf(stderr, "\n");
    fprintf(stderr, "Required Arguments:\n");
    fprintf(stderr, "  <long_reads.fq.gz>    Long-read FASTQ file (can be gzipped)\n");
    fprintf(stderr, "  <asso_kmers.fasta>    Target k-mer sequences in FASTA format\n");
    fprintf(stderr, "  <kmer_length>         K-mer length [%d-%d]\n", MIN_KMER_LEN, MAX_KMER_LEN);
    fprintf(stderr, "  <output_prefix>       Output file prefix\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Output:\n");
    fprintf(stderr, "  <output_prefix>_filtered.fastq\n");
    fprintf(stderr, "    - Reads containing target k-mers\n");
    fprintf(stderr, "    - Read name format: <kmer_name>_<original_read_name>\n");
    fprintf(stderr, "    - Multiple k-mer matches: <kmer1,kmer2,...>_<original_read_name>\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Example:\n");
    fprintf(stderr, "  %s tgs_reads.fq.gz kmers.fa 21 output\n", prog);
    fprintf(stderr, "\n");
    fprintf(stderr, "====================================================================\n");
    fprintf(stderr, "\n");
}

static void insert_seq(khash_t(64) *h, int k, int len, char *seq, char *name)
{
    int i, l;
    uint64_t x[2], mask = (1ULL << k*2) - 1, shift = (k - 1) * 2;
    
    for (i = l = 0, x[0] = x[1] = 0; i < len; ++i) {
        int c = (uint8_t)seq[i] < 128 ? seq_nt4_table[(uint8_t)seq[i]] : 4;
        
        if (c < 4) { 
            // Forward k-mer encoding
            x[0] = (x[0] << 2 | c) & mask;
            // Reverse complement k-mer encoding
            x[1] = x[1] >> 2 | (uint64_t)(3 - c) << shift;
            
            if (++l >= k) { 
                khint_t itr;
                int absent;
                
                // Insert forward k-mer
                itr = kh_put(64, h, x[0], &absent);
                if (absent) kh_val(h, itr) = strdup(name);
                
                itr = kh_put(64, h, x[1], &absent);
                if (absent) kh_val(h, itr) = strdup(name);
            }
        } else {
            l = 0;
            x[0] = x[1] = 0;
        }
    }
}

static char** look_up_all_matches(khash_t(64) *h, int k, int len, char *seq, int *match_count)
{
    int i, l;
    uint64_t x, mask = (1ULL << k*2) - 1;
    char **matches = NULL;
    int matches_size = 0;
    *match_count = 0;
    
    for (i = l = 0, x = 0; i < len; ++i) {
        int c = (uint8_t)seq[i] < 128 ? seq_nt4_table[(uint8_t)seq[i]] : 4;
        
        if (c < 4) { 
            x = (x << 2 | c) & mask;
            
            if (++l >= k) { 
                khint_t itr = kh_get(64, h, x);
                
                if (itr != kh_end(h)) {
                    char *kmer_name = kh_val(h, itr);
                    int already_recorded = 0;
                    
                    // Check if this k-mer is already recorded
                    for (int j = 0; j < *match_count; j++) {
                        if (strcmp(matches[j], kmer_name) == 0) {
                            already_recorded = 1;
                            break;
                        }
                    }
                    
                    if (!already_recorded) {
                        // Dynamically expand array
                        if (*match_count >= matches_size) {
                            matches_size = matches_size == 0 ? INITIAL_MATCH_SIZE : matches_size * 2;
                            matches = realloc(matches, matches_size * sizeof(char*));
                        }
                        matches[*match_count] = strdup(kmer_name);
                        (*match_count)++;
                    }
                }
            }
        } else {
            l = 0;
            x = 0;
        }
    }
    
    return matches;
}

static int load_kmers(const char *kmer_file, khash_t(64) *h, int k)
{
    gzFile fp = gzopen(kmer_file, "r");
    if (!fp) {
        fprintf(stderr, "[ERROR] Open k-mer file failed: %s\n", kmer_file);
        return -1;
    }
    
    kseq_t *ks = kseq_init(fp);
    kh_resize(64, h, HASH_INIT_SIZE);
    
    int cnt = 0;
    while (kseq_read(ks) >= 0) {
        insert_seq(h, k, ks->seq.l, ks->seq.s, ks->name.s);
        cnt++;
    }
    
    kseq_destroy(ks);
    gzclose(fp);
    
    fprintf(stderr, "[INFO] Loaded %d k-mers\n", cnt);
    return cnt;
}

static int process_reads(const char *reads_file, const char *output_file, 
                        khash_t(64) *h, int k)
{
    gzFile fp = gzopen(reads_file, "r");
    if (!fp) {
        fprintf(stderr, "[Error] Open read file failed: %s\n", reads_file);
        return -1;
    }
    
    FILE *out = fopen(output_file, "w");
    if (!out) {
        fprintf(stderr, "[Error] Can not create output file: %s\n", output_file);
        gzclose(fp);
        return -1;
    }
    
    kseq_t *ks = kseq_init(fp);
    uint64_t cnt_total = 0, cnt_found = 0;
    
    fprintf(stderr, "[INFO] Start processing reads...\n");
    
    while (kseq_read(ks) >= 0) {
        cnt_total++;
        
        if ((cnt_total % PROGRESS_INTERVAL) == 0) {
            fprintf(stderr, "[INFO] Processed %lu reads\n", cnt_total);
        }
        
        int match_count = 0;
        char **matched_kmers = look_up_all_matches(h, k, ks->seq.l, ks->seq.s, &match_count);
        
        if (match_count > 0) {
            cnt_found++;
            
            char new_name[MAX_NAME_LEN];
            if (match_count == 1) {
                snprintf(new_name, sizeof(new_name), "%s_%s", 
                        matched_kmers[0], ks->name.s);
            } else {
                strcpy(new_name, matched_kmers[0]);
                for (int i = 1; i < match_count && strlen(new_name) < MAX_PREFIX_LEN; i++) {
                    strcat(new_name, ",");
                    strcat(new_name, matched_kmers[i]);
                }
                strcat(new_name, "_");
                strcat(new_name, ks->name.s);
            }
            
            fprintf(out, "@%s\n%s\n+\n%s\n", 
                   new_name, ks->seq.s, ks->qual.s);
            
            for (int i = 0; i < match_count; i++) {
                free(matched_kmers[i]);
            }
            free(matched_kmers);
        }
    }
    
    kseq_destroy(ks);
    gzclose(fp);
    fclose(out);
    
    fprintf(stderr, "[INFO] Found %lu/%lu reads include the target k-mer (%.2f%%)\n", 
           cnt_found, cnt_total, (double)cnt_found / cnt_total * 100.0);
    
    return 0;
}

int fetch_reads_tgs(int argc, char *argv[])
{
    if (argc != 5) {
        print_usage(argv[0]);
        return 1;
    }
    
    const char *reads_file = argv[1];
    const char *kmer_file = argv[2];
    int k = atoi(argv[3]);
    const char *output_prefix = argv[4];
    
    if (k < MIN_KMER_LEN || k > MAX_KMER_LEN) {
        fprintf(stderr, "[ERROR] k-mer length must within the range of %d-%d\n", 
                MIN_KMER_LEN, MAX_KMER_LEN);
        return 1;
    }
    
    fprintf(stderr, "\n");
    fprintf(stderr, "====================================================================\n");
    fprintf(stderr, " %s v%s\n", PROG_NAME, PROG_VERSION);
    fprintf(stderr, "====================================================================\n");
    fprintf(stderr, "[OPTION] k-mer length: %d\n", k);
    fprintf(stderr, "[OPTION] reads file: %s\n", reads_file);
    fprintf(stderr, "[OPTION] k-mer file: %s\n", kmer_file);
    fprintf(stderr, "[OPTION] output prefix: %s\n", output_prefix);
    fprintf(stderr, "====================================================================\n");
    fprintf(stderr, "\n");
    
    khash_t(64) *h = kh_init(64);
    
    if (load_kmers(kmer_file, h, k) < 0) {
        kh_destroy(64, h);
        return 1;
    }
    
    char output_file[256];
    snprintf(output_file, sizeof(output_file), "%s_filtered.fastq", output_prefix);
    
    int ret = process_reads(reads_file, output_file, h, k);
    
    kh_destroy(64, h);
    
    if (ret == 0) {
        fprintf(stderr, "\n[FINISH] Result reserved to: %s\n\n", output_file);
    }
    
    return ret;
}
