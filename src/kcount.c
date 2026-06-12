#include "ketopt.h"
#include "kstring.h"
#include "kthread.h"
#include "kseq.h"
#include "kvec.h"
#include "ksort.h"
#include "util.h"

KSEQ_INIT(gzFile, gzread)
#include "khashl.h"

#define HASH_BITS 16
#define HASH_MAX ((1 << HASH_BITS) - 1)

#define KMERIA_MAGIC "KMER"
static const uint8_t KMERIA_VERSION = 2;


static inline khint_t hash_func_64(uint64_t x) {
    return (khint_t)hash(x >> HASH_BITS, ~(uint64_t)0);
}
static inline int hash_eq_64(uint64_t x, uint64_t y) { return (x >> HASH_BITS) == (y >> HASH_BITS); }
KHASHL_SET_INIT(static inline, hash_set64_t, hash_set64, uint64_t, hash_func_64, hash_eq_64)
typedef kvec_t(uint64_t) uint64_vec_t;

static inline khint_t hash_func_128(uint128_t x) {
    uint64_t mask = ~(uint64_t)0;
    uint64_t lo = (uint64_t)(x >> HASH_BITS);
    uint64_t hi = (uint64_t)(x >> (HASH_BITS + 64));
    return (khint_t)(hash(lo, mask) ^ hash(hi, mask));
}
static inline int hash_eq_128(uint128_t x, uint128_t y) { return (x >> HASH_BITS) == (y >> HASH_BITS); }
KHASHL_SET_INIT(static inline, hash_set128_t, hash_set128, uint128_t, hash_func_128, hash_eq_128)
typedef kvec_t(uint128_t) uint128_vec_t;

// ---------------------------------------------------------
// Partitioning Structures
// ---------------------------------------------------------
typedef struct {
    void *table; // Dynamically stores hash_set64_t* or hash_set128_t*
    uint64_t partition_key;
} partition_hash_t;

#define partition_cmp(a, b) ((a).partition_key < (b).partition_key)
KSORT_INIT(partition_sort, partition_hash_t, partition_cmp)
typedef kvec_t(partition_hash_t) partition_vec_t;

// Filter parameters structure
typedef struct {
    int min_count;
    int max_count;
} filter_params_t;

static partition_vec_t *init_partition_vec(int part_bits, int is_128) {
    partition_vec_t *part_vec = calloc(1, sizeof(partition_vec_t));
    kv_init(*part_vec);
    kv_resize(partition_hash_t, *part_vec, 1 << part_bits);
    kv_size(*part_vec) = 1 << part_bits;
    for (int idx = 0; idx < 1 << part_bits; ++idx) {
        kv_A(*part_vec, idx).partition_key = idx;
        if (is_128) {
            kv_A(*part_vec, idx).table = hash_set128_init();
        } else {
            kv_A(*part_vec, idx).table = hash_set64_init();
        }
    }
    return part_vec;
}

static void release_partition_vec(partition_vec_t *part_vec, int is_128) {
    for (int idx = 0; idx < kv_size(*part_vec); ++idx) {
        if (kv_A(*part_vec, idx).table) {
            if (is_128) {
                hash_set128_destroy((hash_set128_t*)kv_A(*part_vec, idx).table);
            } else {
                hash_set64_destroy((hash_set64_t*)kv_A(*part_vec, idx).table);
            }
        }
    }
    kv_destroy(*part_vec);
}


// Extraction & Buffering
static void add_kmers_to_buf_64(const char *sequence, int seq_len, int kmer_size, int part_bits, int compress_homo,
                                int canonical_mode, uint64_vec_t *kmer_buf) {
    int pos, valid_len;
    uint64_t fwd[2], kmer_mask = (1ULL << kmer_size * 2) - 1, shift_val = (kmer_size - 1) * 2;
    for (pos = valid_len = 0, fwd[0] = fwd[1] = 0; pos < seq_len; ++pos) {
        if (compress_homo && pos && sequence[pos] == sequence[pos - 1]) continue;
        int base = seq_nt4_table[(uint8_t)sequence[pos]];
        if (base < 4) {
            fwd[0] = (fwd[0] << 2 | base) & kmer_mask;
            fwd[1] = fwd[1] >> 2 | (uint64_t)(3 - base) << shift_val;
            if (++valid_len >= kmer_size) {
                uint64_t current_kmer = canonical_mode ? (fwd[0] < fwd[1] ? fwd[0] : fwd[1]) : fwd[0];
                kv_push(uint64_t, kmer_buf[current_kmer << (64 - kmer_size * 2) >> (64 - part_bits)], current_kmer & kmer_mask);
            }
        } else {
            valid_len = 0, fwd[0] = fwd[1] = 0;
        }
    }
}

static void add_kmers_to_buf_128(const char *sequence, int seq_len, int kmer_size, int part_bits, int compress_homo,
                                 int canonical_mode, uint128_vec_t *kmer_buf) {
    int pos, valid_len;
    uint128_t fwd[2];
    uint128_t kmer_mask = (kmer_size == 64) ? ~(uint128_t)0 : ((((uint128_t)1) << (kmer_size * 2)) - 1);
    int shift_val = (kmer_size - 1) * 2;
    
    for (pos = valid_len = 0, fwd[0] = fwd[1] = 0; pos < seq_len; ++pos) {
        if (compress_homo && pos && sequence[pos] == sequence[pos - 1]) continue;
        int base = seq_nt4_table[(uint8_t)sequence[pos]];
        if (base < 4) {
            fwd[0] = (fwd[0] << 2 | base) & kmer_mask;
            fwd[1] = fwd[1] >> 2 | ((uint128_t)(3 - base)) << shift_val;
            if (++valid_len >= kmer_size) {
                uint128_t current_kmer = canonical_mode ? (fwd[0] < fwd[1] ? fwd[0] : fwd[1]) : fwd[0];
                int shift_up = 128 - kmer_size * 2;
                int shift_down = 128 - part_bits;
                uint64_t part_idx = (uint64_t)((current_kmer << shift_up) >> shift_down);
                kv_push(uint128_t, kmer_buf[part_idx], current_kmer & kmer_mask);
            }
        } else {
            valid_len = 0, fwd[0] = fwd[1] = 0;
        }
    }
}

typedef struct {
    int kmer_size, chunk_size, thread_count, part_bits, compress_homo;
    int canonical_mode;
    int is_128;
    kseq_t *seq_reader;
    partition_vec_t *partitions;
    uint64_t total_seq_len;
    uint64_t total_kmer_count;
} pipeline_data_t;

typedef struct {
    pipeline_data_t *pl_data;
    int seq_count, seq_capacity, chunk_seq_len, chunk_kmer_est;
    int *seq_lengths;
    char **sequences;
    void *kmer_buffers; // Points to uint64_vec_t* or uint128_vec_t*
} chunk_data_t;

// ---------------------------------------------------------
// Hash Hitting Processing
// ---------------------------------------------------------
static void process_chunk_64(void *chunk_ptr, long part_idx, int thread_id) {
    chunk_data_t *chunk = (chunk_data_t *)chunk_ptr;
    uint64_vec_t *buf_part = &((uint64_vec_t*)chunk->kmer_buffers)[part_idx];
    hash_set64_t *hash_table = (hash_set64_t*)kv_A(*chunk->pl_data->partitions, part_idx).table;
    int part_bits = chunk->pl_data->part_bits;
    for (int buf_idx = 0; buf_idx < kv_size(*buf_part); ++buf_idx) {
        khint_t hash_iter; int is_new;
        uint64_t val = kv_A(*buf_part, buf_idx);
        uint64_t key = (val & ((1ULL << (chunk->pl_data->kmer_size * 2 - part_bits)) - 1)) << HASH_BITS;
        hash_iter = hash_set64_put(hash_table, key, &is_new);
        if ((kh_key(hash_table, hash_iter) & HASH_MAX) < HASH_MAX) {
            ++kh_key(hash_table, hash_iter);
        }
    }
}

static void process_chunk_128(void *chunk_ptr, long part_idx, int thread_id) {
    chunk_data_t *chunk = (chunk_data_t *)chunk_ptr;
    uint128_vec_t *buf_part = &((uint128_vec_t*)chunk->kmer_buffers)[part_idx];
    hash_set128_t *hash_table = (hash_set128_t*)kv_A(*chunk->pl_data->partitions, part_idx).table;
    int part_bits = chunk->pl_data->part_bits;
    for (int buf_idx = 0; buf_idx < kv_size(*buf_part); ++buf_idx) {
        khint_t hash_iter; int is_new;
        uint128_t val = kv_A(*buf_part, buf_idx);
        uint128_t mask = (chunk->pl_data->kmer_size * 2 - part_bits >= 128) ? ~(uint128_t)0 : ((((uint128_t)1) << (chunk->pl_data->kmer_size * 2 - part_bits)) - 1);
        uint128_t key = (val & mask) << HASH_BITS;
        hash_iter = hash_set128_put(hash_table, key, &is_new);
        if ((kh_key(hash_table, hash_iter) & HASH_MAX) < HASH_MAX) {
            ++kh_key(hash_table, hash_iter);
        }
    }
}

static void *pipeline_processor(void *pl_ptr, int stage, void *input) {
    pipeline_data_t *pl = (pipeline_data_t *)pl_ptr;
    if (stage == 0) {
        int read_status;
        chunk_data_t *chunk;
        CALLOC(chunk, 1);
        chunk->pl_data = pl;
        while ((read_status = kseq_read(pl->seq_reader)) >= 0) {
            int curr_len = pl->seq_reader->seq.l;
            if (curr_len < pl->kmer_size) continue;
            
            if (chunk->seq_count == chunk->seq_capacity) {
                chunk->seq_capacity = chunk->seq_capacity < 16 ? 16 : chunk->seq_capacity + (chunk->seq_count >> 1);
                REALLOC(chunk->seq_lengths, chunk->seq_capacity);
                REALLOC(chunk->sequences, chunk->seq_capacity);
            }
            MALLOC(chunk->sequences[chunk->seq_count], curr_len);
            memcpy(chunk->sequences[chunk->seq_count], pl->seq_reader->seq.s, curr_len);
            chunk->seq_lengths[chunk->seq_count++] = curr_len;
            chunk->chunk_seq_len += curr_len;
            chunk->chunk_kmer_est += curr_len - pl->kmer_size + 1;
            if (chunk->chunk_seq_len >= pl->chunk_size) break;
        }
        if (chunk->chunk_seq_len == 0) {
            free(chunk);
        } else {
            return chunk;
        }
    } else if (stage == 1) {
        chunk_data_t *chunk = (chunk_data_t *)input;
        int part_count = 1 << pl->part_bits, est_size;
        est_size = (int)(chunk->chunk_kmer_est * 1.2 / part_count) + 1;
        
        if (pl->is_128) {
            uint128_vec_t *bufs; CALLOC(bufs, part_count); chunk->kmer_buffers = bufs;
            for (int p = 0; p < part_count; ++p) kv_resize(uint128_t, bufs[p], est_size > 8 ? est_size : 8);
            for (int seq_idx = 0; seq_idx < chunk->seq_count; ++seq_idx) {
                add_kmers_to_buf_128(chunk->sequences[seq_idx], chunk->seq_lengths[seq_idx], pl->kmer_size, 
                                     pl->part_bits, pl->compress_homo, pl->canonical_mode, bufs);
                free(chunk->sequences[seq_idx]);
            }
        } else {
            uint64_vec_t *bufs; CALLOC(bufs, part_count); chunk->kmer_buffers = bufs;
            for (int p = 0; p < part_count; ++p) kv_resize(uint64_t, bufs[p], est_size > 8 ? est_size : 8);
            for (int seq_idx = 0; seq_idx < chunk->seq_count; ++seq_idx) {
                add_kmers_to_buf_64(chunk->sequences[seq_idx], chunk->seq_lengths[seq_idx], pl->kmer_size, 
                                    pl->part_bits, pl->compress_homo, pl->canonical_mode, bufs);
                free(chunk->sequences[seq_idx]);
            }
        }
        
        pl->total_seq_len += chunk->chunk_seq_len;
        for (int p = 0; p < part_count; ++p) {
            pl->total_kmer_count += pl->is_128 ? kv_size(((uint128_vec_t*)chunk->kmer_buffers)[p]) : kv_size(((uint64_vec_t*)chunk->kmer_buffers)[p]);
        }
        free(chunk->sequences);
        free(chunk->seq_lengths);
        return chunk;
    } else if (stage == 2) {
        chunk_data_t *chunk = (chunk_data_t *)input;
        int part_count = 1 << pl->part_bits;
        
        if (pl->is_128) {
            kt_for(pl->thread_count, process_chunk_128, chunk, part_count);
            for (int p = 0; p < part_count; ++p) kv_destroy(((uint128_vec_t*)chunk->kmer_buffers)[p]);
        } else {
            kt_for(pl->thread_count, process_chunk_64, chunk, part_count);
            for (int p = 0; p < part_count; ++p) kv_destroy(((uint64_vec_t*)chunk->kmer_buffers)[p]);
        }
        free(chunk->kmer_buffers);
        free(chunk);
    }
    return 0;
}

static partition_vec_t *process_input_file(const char *filename, int kmer_size, int part_bits, int chunk_size,
                                           int thread_count, int compress_homo, int canonical_mode,
                                           partition_vec_t *existing_parts, uint64_t* total_len, uint64_t* total_kmers) {
    pipeline_data_t pl;
    gzFile file_handle;
    file_handle = filename && strcmp(filename, "-") ? gzopen(filename, "r") : gzdopen(fileno(stdin), "r");
    if (file_handle == 0) return 0;
    
    pl.seq_reader = kseq_init(file_handle);
    pl.kmer_size = kmer_size;
    pl.thread_count = thread_count;
    pl.part_bits = part_bits;
    pl.is_128 = (kmer_size >= 32);
    pl.partitions = existing_parts ? existing_parts : init_partition_vec(part_bits, pl.is_128);
    pl.chunk_size = chunk_size;
    pl.compress_homo = compress_homo;
    pl.canonical_mode = canonical_mode;
    pl.total_seq_len = 0;
    pl.total_kmer_count = 0;
    
    kt_pipeline(3, pipeline_processor, &pl, 3);
    *total_len = pl.total_seq_len;
    *total_kmers = pl.total_kmer_count;
    kseq_destroy(pl.seq_reader);
    gzclose(file_handle);
    return pl.partitions;
}

typedef struct {
    uint32_t counts[HASH_MAX + 1]; // HASH_MAX+1
} count_buffer_t;

typedef struct {
    const partition_vec_t *partitions;
    count_buffer_t *thread_counts;
    int is_128;
} hist_helper_t;

static void compute_hist_part(void *helper_ptr, long part_idx, int thread_id) {
    hist_helper_t *helper = (hist_helper_t *)helper_ptr;
    uint32_t *thread_cnt = helper->thread_counts[thread_id].counts;
    
    if (helper->is_128) {
        hash_set128_t *hash_table = (hash_set128_t*)kv_A(*helper->partitions, part_idx).table;
        for (khint_t hash_iter = 0; hash_iter < kh_end(hash_table); ++hash_iter) {
            if (kh_exist(hash_table, hash_iter)) {
                int freq = kh_key(hash_table, hash_iter) & HASH_MAX;
                ++thread_cnt[freq]; // freq
            }
        }
    } else {
        hash_set64_t *hash_table = (hash_set64_t*)kv_A(*helper->partitions, part_idx).table;
        for (khint_t hash_iter = 0; hash_iter < kh_end(hash_table); ++hash_iter) {
            if (kh_exist(hash_table, hash_iter)) {
                int freq = kh_key(hash_table, hash_iter) & HASH_MAX;
                ++thread_cnt[freq]; // freq: 0..HASH_MAX
            }
        }
    }
}

static int output_histogram(const partition_vec_t *partitions, int thread_count, const char *filename, int is_128) {
    FILE *output_file;
    hist_helper_t helper;
    uint32_t aggregated_counts[HASH_MAX + 1]; 
    int idx, thread_idx;
    helper.partitions = partitions;
    helper.is_128 = is_128;
    if ((output_file = strcmp(filename, "-") ? fopen(filename, "wb") : stdout) == 0) return -1;
    CALLOC(helper.thread_counts, thread_count);
    kt_for(thread_count, compute_hist_part, &helper, kv_size(*partitions));
    
    for (idx = 0; idx <= HASH_MAX; ++idx) aggregated_counts[idx] = 0;
    for (thread_idx = 0; thread_idx < thread_count; ++thread_idx) {
        for (idx = 0; idx <= HASH_MAX; ++idx)
            aggregated_counts[idx] += helper.thread_counts[thread_idx].counts[idx];
    }
    free(helper.thread_counts);

    // freq=1..HASH_MAX-1
    for (idx = 1; idx < HASH_MAX; ++idx) {
        if (aggregated_counts[idx]) fprintf(output_file, "%d\t%u\n", idx, aggregated_counts[idx]);
    }
    if (aggregated_counts[HASH_MAX])
        fprintf(output_file, ">=%d\t%u\n", HASH_MAX, aggregated_counts[HASH_MAX]);

    fclose(output_file);
    return 0;
}

// ---------------------------------------------------------
// Key-Value Outputs & Structs
// ---------------------------------------------------------
#pragma pack(2)
typedef struct { uint64_t kmer_val; uint16_t freq; } kmer_freq_pair64_t;
typedef struct { uint128_t kmer_val; uint16_t freq; } kmer_freq_pair128_t;
#pragma pack()

typedef kvec_t(kmer_freq_pair64_t) kmer_pair_vec64_t;
typedef kvec_t(kmer_freq_pair128_t) kmer_pair_vec128_t;

#define kmer_pair_cmp64(a, b) ((a).kmer_val < (b).kmer_val)
KSORT_INIT(kmer_pair_sort64, kmer_freq_pair64_t, kmer_pair_cmp64)

#define kmer_pair_cmp128(a, b) ((a).kmer_val < (b).kmer_val)
KSORT_INIT(kmer_pair_sort128, kmer_freq_pair128_t, kmer_pair_cmp128)

typedef struct {
    const partition_vec_t *partitions;
    void *kmer_vectors; // kmer_pair_vec64_t* or kmer_pair_vec128_t*
    int kmer_len, part_bits, is_128;
    filter_params_t *filter; 
} kmer_output_helper_t;

static void extract_kmers_part(void *helper_ptr, long part_idx, int thread_id) {
    kmer_output_helper_t *helper = (kmer_output_helper_t *)helper_ptr;
    int k_len = helper->kmer_len, p_bits = helper->part_bits;

    if (helper->is_128) {
        hash_set128_t *hash_table = (hash_set128_t*)kv_A(*helper->partitions, part_idx).table;
        if (!kh_size(hash_table)) return;
        kmer_pair_vec128_t *vec = &((kmer_pair_vec128_t*)helper->kmer_vectors)[part_idx];
        kv_init(*vec);
        kv_resize(kmer_freq_pair128_t, *vec, kh_size(hash_table));
        
        uint128_t kmer_mask = (k_len == 64) ? ~(uint128_t)0 : ((((uint128_t)1) << (k_len * 2)) - 1);
        for (khint_t hash_iter = 0; hash_iter < kh_end(hash_table); ++hash_iter) {
            if (kh_exist(hash_table, hash_iter)) {
                const uint128_t hash_key = kh_key(hash_table, hash_iter);
                const uint16_t freq_val = hash_key & HASH_MAX;
                if (helper->filter) {
                    if (freq_val < helper->filter->min_count) continue;
                    if (helper->filter->max_count >= 0 && freq_val > helper->filter->max_count) continue;
                }
                const uint128_t reconstructed_kmer = ((hash_key >> HASH_BITS) | (((uint128_t)kv_A(*helper->partitions, part_idx).partition_key) << (k_len * 2 - p_bits))) & kmer_mask;
                const kmer_freq_pair128_t pair = {reconstructed_kmer, freq_val};
                kv_push(kmer_freq_pair128_t, *vec, pair);
            }
        }
        ks_introsort(kmer_pair_sort128, kv_size(*vec), vec->a);
        hash_set128_destroy(hash_table);
    } else {
        hash_set64_t *hash_table = (hash_set64_t*)kv_A(*helper->partitions, part_idx).table;
        if (!kh_size(hash_table)) return;
        kmer_pair_vec64_t *vec = &((kmer_pair_vec64_t*)helper->kmer_vectors)[part_idx];
        kv_init(*vec);
        kv_resize(kmer_freq_pair64_t, *vec, kh_size(hash_table));
        
        uint64_t kmer_mask = (1ULL << (k_len * 2)) - 1;
        for (khint_t hash_iter = 0; hash_iter < kh_end(hash_table); ++hash_iter) {
            if (kh_exist(hash_table, hash_iter)) {
                const uint64_t hash_key = kh_key(hash_table, hash_iter);
                const uint16_t freq_val = hash_key & HASH_MAX;
                if (helper->filter) {
                    if (freq_val < helper->filter->min_count) continue;
                    if (helper->filter->max_count >= 0 && freq_val > helper->filter->max_count) continue;
                }
                const uint64_t reconstructed_kmer = ((hash_key >> HASH_BITS) | (((uint64_t)kv_A(*helper->partitions, part_idx).partition_key) << (k_len * 2 - p_bits))) & kmer_mask;
                const kmer_freq_pair64_t pair = {reconstructed_kmer, freq_val};
                kv_push(kmer_freq_pair64_t, *vec, pair);
            }
        }
        ks_introsort(kmer_pair_sort64, kv_size(*vec), vec->a);
        hash_set64_destroy(hash_table);
    }
    kv_A(*helper->partitions, part_idx).table = NULL;
}


static long write_binary_header(FILE *fp, int kmer_len, int is_128) {
    uint8_t padding[3] = {0, 0, 0};
    fwrite(KMERIA_MAGIC, 1, 4, fp);               // 4字节 magic "KMER"
    fwrite(&KMERIA_VERSION, 1, 1, fp);             // 1字节版本
    fwrite(padding, 1, 3, fp);                     // 3字节 padding，共8字节前缀
    long header_data_offset = ftell(fp);           // 记录 header_data 起始偏移
    uint32_t header_data[3] = {(uint32_t)kmer_len, (uint32_t)(is_128 ? 2 : 1), 0};
    fwrite(header_data, 4, 3, fp);                 // 12字节，其中 [2]=total_pairs 待回写
    return header_data_offset;
}

static int output_kmers_in_memory(const partition_vec_t *partitions, int thread_count,
                                  int kmer_len, int part_bits, const char *filename, int is_text,
                                  filter_params_t *filter, int is_128) {
    FILE *output_file;
    if ((output_file = (filename && strcmp(filename, "-")) ? fopen(filename, "wb") : stdout) == 0) return -1;
    ks_introsort(partition_sort, kv_size(*partitions), partitions->a);
    
    kmer_output_helper_t helper;
    helper.partitions = partitions;
    helper.kmer_len = kmer_len;
    helper.part_bits = part_bits;
    helper.filter = filter; 
    helper.is_128 = is_128;
    
    if (is_128) helper.kmer_vectors = calloc(kv_size(*partitions), sizeof(kmer_pair_vec128_t));
    else helper.kmer_vectors = calloc(kv_size(*partitions), sizeof(kmer_pair_vec64_t));
    
    kt_for(thread_count, extract_kmers_part, &helper, kv_size(*partitions));
    
    if (!is_text) {
        // Bug5修复：统计 total_pairs 后用 write_binary_header 写完整头部
        uint32_t total_pairs = 0;
        for (int part_idx = 0; part_idx < kv_size(*partitions); ++part_idx) {
            total_pairs += is_128
                ? kv_size(((kmer_pair_vec128_t*)helper.kmer_vectors)[part_idx])
                : kv_size(((kmer_pair_vec64_t*)helper.kmer_vectors)[part_idx]);
        }
        long header_data_offset = write_binary_header(output_file, kmer_len, is_128);
        fseek(output_file, header_data_offset + 8, SEEK_SET); // +8 跳过 kmer_len(4) + type(4)
        fwrite(&total_pairs, 4, 1, output_file);
        fseek(output_file, 0, SEEK_END); 
    }
    
    for (int part_idx = 0; part_idx < kv_size(*partitions); ++part_idx) {
        if (!is_text) {
            if (is_128) {
                kmer_pair_vec128_t *vec = &((kmer_pair_vec128_t*)helper.kmer_vectors)[part_idx];
                if (kv_size(*vec)) fwrite(vec->a, sizeof(kmer_freq_pair128_t), kv_size(*vec), output_file);
            } else {
                kmer_pair_vec64_t *vec = &((kmer_pair_vec64_t*)helper.kmer_vectors)[part_idx];
                if (kv_size(*vec)) fwrite(vec->a, sizeof(kmer_freq_pair64_t), kv_size(*vec), output_file);
            }
        } else {
            char kmer_str[kmer_len + 1];
            if (is_128) {
                kmer_pair_vec128_t *vec = &((kmer_pair_vec128_t*)helper.kmer_vectors)[part_idx];
                for (int vec_idx = 0; vec_idx < kv_size(*vec); ++vec_idx) {
                    kmer2seq_128(kmer_str, kmer_len, vec->a[vec_idx].kmer_val);
                    fprintf(output_file, "%s\t%u\n", kmer_str, vec->a[vec_idx].freq);
                }
            } else {
                kmer_pair_vec64_t *vec = &((kmer_pair_vec64_t*)helper.kmer_vectors)[part_idx];
                for (int vec_idx = 0; vec_idx < kv_size(*vec); ++vec_idx) {
                    kmer2seq_64(kmer_str, kmer_len, vec->a[vec_idx].kmer_val);
                    fprintf(output_file, "%s\t%u\n", kmer_str, vec->a[vec_idx].freq);
                }
            }
        }
    }
    
    for (int part_idx = 0; part_idx < kv_size(*partitions); ++part_idx) {
        if (is_128) kv_destroy(((kmer_pair_vec128_t*)helper.kmer_vectors)[part_idx]);
        else kv_destroy(((kmer_pair_vec64_t*)helper.kmer_vectors)[part_idx]);
    }
    free(helper.kmer_vectors);
    fclose(output_file);
    return 0;
}

// ---------------------------------------------------------
// External Merge Sort
// ---------------------------------------------------------
#define BUFFER_CAPACITY 8192

typedef struct {
    FILE *file_ptr;
    void *data_buf;
    int buf_capacity;
    int buf_index;
    int buf_valid;   // 当前buffer中有效元素数（最后一次read可能不足BUFFER_CAPACITY）
    int is_end;
    char temp_filename[256];
    int is_128;
} merge_stream_t;

static merge_stream_t* create_merge_stream(const char *temp_file, int is_128) {
    merge_stream_t *stream = calloc(1, sizeof(merge_stream_t));
    if (!stream) return NULL;
    stream->file_ptr = fopen(temp_file, "rb");
    if (!stream->file_ptr) { free(stream); return NULL; }
    
    stream->is_128 = is_128;
    stream->buf_capacity = BUFFER_CAPACITY;
    size_t elem_sz = is_128 ? sizeof(kmer_freq_pair128_t) : sizeof(kmer_freq_pair64_t);
    stream->data_buf = malloc(stream->buf_capacity * elem_sz);
    
    int read_count = fread(stream->data_buf, elem_sz, stream->buf_capacity, stream->file_ptr);
    if (read_count == 0) { fclose(stream->file_ptr); free(stream->data_buf); free(stream); return NULL; }
    
    stream->buf_index = 0;
    stream->buf_valid = read_count;
    stream->is_end = (read_count < stream->buf_capacity);
    strcpy(stream->temp_filename, temp_file);
    return stream;
}

static inline uint128_t get_current_kmer(merge_stream_t *stream) {
    if (stream->is_128) return ((kmer_freq_pair128_t*)stream->data_buf)[stream->buf_index].kmer_val;
    return ((kmer_freq_pair64_t*)stream->data_buf)[stream->buf_index].kmer_val;
}

static inline uint16_t get_current_freq(merge_stream_t *stream) {
    if (stream->is_128) return ((kmer_freq_pair128_t*)stream->data_buf)[stream->buf_index].freq;
    return ((kmer_freq_pair64_t*)stream->data_buf)[stream->buf_index].freq;
}

static int advance_stream(merge_stream_t *stream) {
    stream->buf_index++;
    if (stream->buf_index >= stream->buf_valid) {
        if (stream->is_end) return 0;
        size_t elem_sz = stream->is_128 ? sizeof(kmer_freq_pair128_t) : sizeof(kmer_freq_pair64_t);
        int read_count = fread(stream->data_buf, elem_sz, stream->buf_capacity, stream->file_ptr);
        if (read_count == 0) { stream->is_end = 1; return 0; }
        stream->buf_index = 0;
        stream->buf_valid = read_count;
        stream->is_end = (read_count < stream->buf_capacity);
    }
    return 1;
}

static void destroy_merge_stream(merge_stream_t *stream) {
    if (stream->file_ptr) fclose(stream->file_ptr);
    if (stream->data_buf) free(stream->data_buf);
    unlink(stream->temp_filename);
    free(stream);
}

// Unified heap node using uint128_t to handle both types cleanly
typedef struct {
    merge_stream_t *stream;
    uint128_t kmer_val;
    int index;
} heap_node_t;

#define heap_cmp(a, b) ((a).kmer_val > (b).kmer_val)

static void heapify_down(heap_node_t *heap, int size, int pos) {
    int smallest = pos;
    int left = 2 * pos + 1;
    int right = 2 * pos + 2;
    if (left < size && heap_cmp(heap[smallest], heap[left])) smallest = left;
    if (right < size && heap_cmp(heap[smallest], heap[right])) smallest = right;
    if (smallest != pos) {
        heap_node_t temp = heap[pos];
        heap[pos] = heap[smallest];
        heap[smallest] = temp;
        heapify_down(heap, size, smallest);
    }
}

static void heap_extract_min(heap_node_t *heap, int *size) {
    heap[0] = heap[--(*size)];
    heapify_down(heap, *size, 0);
}

static void heap_insert(heap_node_t *heap, int *size, heap_node_t node) {
    int pos = (*size)++;
    heap[pos] = node;
    while (pos > 0) {
        int parent = (pos - 1) / 2;
        if (heap_cmp(heap[parent], heap[pos])) {
            heap_node_t temp = heap[pos];
            heap[pos] = heap[parent];
            heap[parent] = temp;
            pos = parent;
        } else break;
    }
}

typedef struct {
    const partition_vec_t *partitions;
    int kmer_len;
    int part_bits;
    char temp_dir[256];
    int *temp_files_exist;
    filter_params_t *filter;
    int is_128;
} output_helper_t;

static void process_partition_bucket(void *helper_ptr, long part_idx, int thread_id) {
    output_helper_t *helper = (output_helper_t*)helper_ptr;
    char temp_file[512];
    snprintf(temp_file, sizeof(temp_file), "%s/part_%06ld.tmp", helper->temp_dir, part_idx);

    if (helper->is_128) {
        hash_set128_t *hash_table = (hash_set128_t*)kv_A(*helper->partitions, part_idx).table;
        if (!kh_size(hash_table)) { helper->temp_files_exist[part_idx] = 0; return; }
        kmer_pair_vec128_t pair_list; kv_init(pair_list);
        kv_resize(kmer_freq_pair128_t, pair_list, kh_size(hash_table));
        
        uint128_t kmer_mask = (helper->kmer_len == 64) ? ~(uint128_t)0 : ((((uint128_t)1) << (helper->kmer_len * 2)) - 1);
        for (khint_t hash_iter = 0; hash_iter < kh_end(hash_table); ++hash_iter) {
            if (kh_exist(hash_table, hash_iter)) {
                const uint128_t hash_key = kh_key(hash_table, hash_iter);
                const uint16_t freq_val = hash_key & HASH_MAX;
                if (helper->filter) {
                    if (freq_val < helper->filter->min_count) continue;
                    if (helper->filter->max_count >= 0 && freq_val > helper->filter->max_count) continue;
                }
                const uint128_t reconstructed_kmer = ((hash_key >> HASH_BITS) | (((uint128_t)kv_A(*helper->partitions, part_idx).partition_key) << (helper->kmer_len * 2 - helper->part_bits))) & kmer_mask;
                const kmer_freq_pair128_t pair = {reconstructed_kmer, freq_val};
                kv_push(kmer_freq_pair128_t, pair_list, pair);
            }
        }
        ks_introsort(kmer_pair_sort128, kv_size(pair_list), pair_list.a);
        FILE *temp_fp = fopen(temp_file, "wb");
        if (temp_fp) { fwrite(pair_list.a, sizeof(kmer_freq_pair128_t), kv_size(pair_list), temp_fp); fclose(temp_fp); helper->temp_files_exist[part_idx] = 1; }
        kv_destroy(pair_list); hash_set128_destroy(hash_table);
    } else {
        hash_set64_t *hash_table = (hash_set64_t*)kv_A(*helper->partitions, part_idx).table;
        if (!kh_size(hash_table)) { helper->temp_files_exist[part_idx] = 0; return; }
        kmer_pair_vec64_t pair_list; kv_init(pair_list);
        kv_resize(kmer_freq_pair64_t, pair_list, kh_size(hash_table));
        
        uint64_t kmer_mask = (1ULL << (helper->kmer_len * 2)) - 1;
        for (khint_t hash_iter = 0; hash_iter < kh_end(hash_table); ++hash_iter) {
            if (kh_exist(hash_table, hash_iter)) {
                const uint64_t hash_key = kh_key(hash_table, hash_iter);
                const uint16_t freq_val = hash_key & HASH_MAX;
                if (helper->filter) {
                    if (freq_val < helper->filter->min_count) continue;
                    if (helper->filter->max_count >= 0 && freq_val > helper->filter->max_count) continue;
                }
                const uint64_t reconstructed_kmer = ((hash_key >> HASH_BITS) | (((uint64_t)kv_A(*helper->partitions, part_idx).partition_key) << (helper->kmer_len * 2 - helper->part_bits))) & kmer_mask;
                const kmer_freq_pair64_t pair = {reconstructed_kmer, freq_val};
                kv_push(kmer_freq_pair64_t, pair_list, pair);
            }
        }
        ks_introsort(kmer_pair_sort64, kv_size(pair_list), pair_list.a);
        FILE *temp_fp = fopen(temp_file, "wb");
        if (temp_fp) { fwrite(pair_list.a, sizeof(kmer_freq_pair64_t), kv_size(pair_list), temp_fp); fclose(temp_fp); helper->temp_files_exist[part_idx] = 1; }
        kv_destroy(pair_list); hash_set64_destroy(hash_table);
    }
    kv_A(*helper->partitions, part_idx).table = NULL;
}

static int output_kmers_external(const partition_vec_t *partitions, int thread_count,
                                 int kmer_len, int part_bits, const char *filename, int is_text,
                                 filter_params_t *filter, int is_128) {
    char temp_dir[256];
    snprintf(temp_dir, sizeof(temp_dir), "/tmp/kmer_count_%d", getpid());
    if (mkdir(temp_dir, 0755) != 0 && errno != EEXIST) {
        fprintf(stderr, "[ERROR] Failed to create temporary directory: %s\n", temp_dir);
        return -1;
    }
    
    fprintf(stderr, "[INFO] Performing external sort (temp dir: %s)\n", temp_dir);
    output_helper_t helper;
    helper.partitions = partitions;
    helper.kmer_len = kmer_len;
    helper.part_bits = part_bits;
    helper.filter = filter;
    helper.is_128 = is_128;
    strcpy(helper.temp_dir, temp_dir);
    helper.temp_files_exist = calloc(kv_size(*partitions), sizeof(int));
    
    kt_for(thread_count, process_partition_bucket, &helper, kv_size(*partitions));
    
    int file_count = 0;
    for (int idx = 0; idx < kv_size(*partitions); ++idx) if (helper.temp_files_exist[idx]) file_count++;
    fprintf(stderr, "[INFO] Generated %d temp files\n", file_count);
    
    if (file_count == 0) {
        free(helper.temp_files_exist);
        rmdir(temp_dir);
        return 0;
    }
    FILE *output_fp = (filename && strcmp(filename, "-")) ? fopen(filename, "wb") : stdout;
    if (!output_fp) { free(helper.temp_files_exist); return -1; }
    
    long header_data_offset = 0;
    if (!is_text) {
        header_data_offset = write_binary_header(output_fp, kmer_len, is_128);
    }
    
    heap_node_t *merge_heap = calloc(file_count, sizeof(heap_node_t));
    int heap_size = 0, stream_idx = 0;
    for (int part_idx = 0; part_idx < kv_size(*partitions); ++part_idx) {
        if (!helper.temp_files_exist[part_idx]) continue;
        char temp_file[512];
        snprintf(temp_file, sizeof(temp_file), "%s/part_%06d.tmp", temp_dir, part_idx);
        merge_stream_t *new_stream = create_merge_stream(temp_file, is_128);
        if (new_stream) {
            heap_node_t node = {new_stream, get_current_kmer(new_stream), stream_idx++};
            heap_insert(merge_heap, &heap_size, node);
        }
    }
    
    fprintf(stderr, "[INFO] Merging %d streams using heap...\n", heap_size);
    uint64_t output_count = 0;
    char kmer_buffer[kmer_len + 1];
    kmer_buffer[kmer_len] = '\0';
    
    size_t elem_sz = is_128 ? sizeof(kmer_freq_pair128_t) : sizeof(kmer_freq_pair64_t);
    void *output_buf = malloc(BUFFER_CAPACITY * elem_sz);
    int buf_pos = 0;
    
    while (heap_size > 0) {
        heap_node_t min_node = merge_heap[0];
        uint128_t current_kmer = get_current_kmer(min_node.stream);
        uint16_t current_freq = get_current_freq(min_node.stream);

        if (!is_text) {
            if (is_128) {
                ((kmer_freq_pair128_t*)output_buf)[buf_pos].kmer_val = current_kmer;
                ((kmer_freq_pair128_t*)output_buf)[buf_pos].freq = current_freq;
            } else {
                ((kmer_freq_pair64_t*)output_buf)[buf_pos].kmer_val = (uint64_t)current_kmer;
                ((kmer_freq_pair64_t*)output_buf)[buf_pos].freq = current_freq;
            }
            buf_pos++;
            if (buf_pos >= BUFFER_CAPACITY) {
                fwrite(output_buf, elem_sz, buf_pos, output_fp);
                buf_pos = 0;
            }
        } else {
            if (is_128) kmer2seq_128(kmer_buffer, kmer_len, current_kmer);
            else kmer2seq_64(kmer_buffer, kmer_len, (uint64_t)current_kmer);
            fprintf(output_fp, "%s\t%u\n", kmer_buffer, current_freq);
        }
        
        output_count++;
        if (!advance_stream(min_node.stream)) {
            destroy_merge_stream(min_node.stream);
            heap_extract_min(merge_heap, &heap_size);
        } else {
            merge_heap[0].kmer_val = get_current_kmer(min_node.stream);
            heapify_down(merge_heap, heap_size, 0);
        }
        if (output_count % 10000000 == 0) fprintf(stderr, "\r[INFO] Processed %lu k-mers...", output_count);
    }
    
    if (!is_text && buf_pos > 0) fwrite(output_buf, elem_sz, buf_pos, output_fp);
    free(output_buf);
    free(merge_heap);
    
    fprintf(stderr, "\r[INFO] Merge finished. Total k-mers: %lu \n", output_count);

    if (!is_text && output_fp != stdout) {
        fseek(output_fp, header_data_offset + 8, SEEK_SET);
        uint32_t total_count = (uint32_t)output_count;
        fwrite(&total_count, 4, 1, output_fp);
    }
    if (output_fp != stdout) fclose(output_fp);
    free(helper.temp_files_exist);
    rmdir(temp_dir);
    return 0;
}

static int output_kmers_smart(const partition_vec_t *partitions, int thread_count,
                              int kmer_len, int part_bits, const char *filename, int is_text,
                              filter_params_t *filter, int is_128) {
    uint64_t total_unique_kmers = 0;
    for (int idx = 0; idx < kv_size(*partitions); ++idx) {
        if (is_128) total_unique_kmers += kh_size((hash_set128_t*)kv_A(*partitions, idx).table);
        else total_unique_kmers += kh_size((hash_set64_t*)kv_A(*partitions, idx).table);
    }
    uint64_t mem_estimate = total_unique_kmers * (is_128 ? sizeof(kmer_freq_pair128_t) : sizeof(kmer_freq_pair64_t));
    uint64_t mem_limit = 1ULL * 1024 * 1024 * 1024 * 64; 
    
    fprintf(stderr, "[INFO] Unique k-mers: %lu\n", total_unique_kmers);
    fprintf(stderr, "[INFO] Memory estimate: %.2f GB\n", mem_estimate / 1024.0 / 1024.0 / 1024.0);
    
    if (mem_estimate < mem_limit) {
        fprintf(stderr, "[INFO] Using memory-based sort\n");
        return output_kmers_in_memory(partitions, thread_count, kmer_len, part_bits, filename, is_text, filter, is_128);
    } else {
        fprintf(stderr, "[INFO] Using disk-based merge sort\n");
        return output_kmers_external(partitions, thread_count, kmer_len, part_bits, filename, is_text, filter, is_128);
    }
}

int kcount(int argc, char *argv[]) {
    partition_vec_t *partitions = 0;
    int opt_char, kmer_size = 21, part_bits = HASH_BITS, chunk_size = 100000000, thread_count = 4;
    int compress_homo = 0;
    int canonical_mode = 1;
    int text_output = 0;
    char *output_file = "-";
    char *hist_output = NULL;
    
    filter_params_t filter;
    filter.min_count = 1; 
    filter.max_count = -1; 
    int use_filter = 0; 
    
    ketopt_t opt_parser = KETOPT_INIT;
    while ((opt_char = ketopt(&opt_parser, argc, argv, 1, "k:p:Tt:H:o:cCm:M:h", 0)) >= 0) {
        switch (opt_char) {
            case 'k': kmer_size = atoi(opt_parser.arg); break;
            case 'p': part_bits = atoi(opt_parser.arg); break;
            case 'T': text_output = 1; break;
            case 't': thread_count = atoi(opt_parser.arg); break;
            case 'H': hist_output = opt_parser.arg; break;
            case 'o': output_file = opt_parser.arg; break;
            case 'c': compress_homo = 1; break;
            case 'C': canonical_mode = 0; break;
            case 'm': filter.min_count = atoi(opt_parser.arg); use_filter = 1; break;
            case 'M': filter.max_count = atoi(opt_parser.arg); use_filter = 1; break;
            case 'h':
                fprintf(stderr, "Usage: kmeria count [options] <input.fa|fa.gz|fq|fq.gz>\n\n");
                fprintf(stderr, "Options:\n");
                fprintf(stderr, " -k INT        Length of k-mer [2-63] [default: %d]\n", kmer_size);
                fprintf(stderr, " -t INT        Thread number [default: %d]\n", thread_count);
                fprintf(stderr, " -C            Count strands separately (no canonical)\n");
                fprintf(stderr, " -H FILE       Histogram output file\n");
                fprintf(stderr, " -o FILE       Result output [default: stdout]\n");
                fprintf(stderr, " -T            Output in plain text format\n");
                fprintf(stderr, " -p INT        Partitioning bits [default: %d]\n", HASH_BITS);
                fprintf(stderr, "\nFiltering Options:\n");
                fprintf(stderr, " -m INT        Minimum k-mer abundance to output [default: 1]\n");
                fprintf(stderr, " -M INT        Maximum k-mer abundance to output\n");
                fprintf(stderr, "\nAdvanced:\n");
                fprintf(stderr, " -c            Compress homopolymers (experimental)\n");
                fprintf(stderr, "\nExamples:\n");
                fprintf(stderr, "  kmeria count -k 21 -t 8 input.fq.gz -o results.bin\n");
                fprintf(stderr, "  kmeria count -k 21 -m 5 input.fq.gz -o results.bin\n");
                fprintf(stderr, "  kmeria count -k 21 -m 5 -M 100 input.fq.gz -o results.bin\n");
                fprintf(stderr, "  kmeria count -k 21 -C input.fq.gz -o results.bin\n");
                fprintf(stderr, "  kmeria count -k 21 -H hist.txt input.fq.gz -o results.bin\n");
                fprintf(stderr, "  kmeria count -k 21 -T input.fq.gz -o results.txt\n");
		return 0;
            default:
                fprintf(stderr, "[ERROR] Invalid option: %c\n", opt_char);
                return 1;
        }
    }
    if (argc - opt_parser.ind < 1) {
        fprintf(stderr, "[ERROR] Input file required\n\n");
	fprintf(stderr, "Usage: kmeria count [options] <input.fa|fa.gz|fq|fq.gz>\n");
	fprintf(stderr, "Run 'kmeria count -h' for details\n");
        return 1;
    }
    
    if (kmer_size < 2 || kmer_size > 63) {
        fprintf(stderr, "[ERROR] k-mer length must be 2-63\n");
        return 1;
    }
    
    // Core decision branch logic based on k-mer parameter 
    int is_128 = (kmer_size >= 32); 

    partitions = 0;
    uint64_t file_seq_len, file_kmer_cnt, overall_seq_len = 0, overall_kmer_cnt = 0;
    
    fprintf(stderr, "[INFO] Initiating k-mer count...\n");
    fprintf(stderr, "[INFO] k-mer length: %d\n", kmer_size);
    fprintf(stderr, "[INFO] Internal engine(bit): %s\n", is_128 ? "128-bit" : "64-bit");
    fprintf(stderr, "[INFO] Threads: %d\n", thread_count);
    fprintf(stderr, "[INFO] Canonical: %s\n", canonical_mode ? "on" : "off");
    fprintf(stderr, "[INFO] Format: %s\n", text_output ? "text" : "binary");
    
    if (use_filter) {
        fprintf(stderr, "[INFO] Filter: min_count=%d, max_count=%d\n", filter.min_count, filter.max_count >= 0 ? filter.max_count : -1);
    }
    
    for (int file_idx = opt_parser.ind; file_idx < argc; ++file_idx) {
        fprintf(stderr, "[INFO] Handling file: %s\n", argv[file_idx]);
        partitions = process_input_file(argv[file_idx], kmer_size, part_bits, chunk_size, thread_count, compress_homo,
                                        canonical_mode, partitions, &file_seq_len, &file_kmer_cnt);
        if (!partitions) {
            fprintf(stderr, "[ERROR] Error handling file: %s\n", argv[file_idx]);
            return 1;
        }
        overall_seq_len += file_seq_len;
        overall_kmer_cnt += file_kmer_cnt;
        fprintf(stderr, "[INFO] File summary - bases: %lu, k-mers: %lu\n", file_seq_len, file_kmer_cnt);
    }
    
    fprintf(stderr, "\n[STATS] Overall bases: %lu\n", overall_seq_len);
    fprintf(stderr, "[STATS] Overall k-mers: %lu\n", overall_kmer_cnt);
    
    if (hist_output) {
        fprintf(stderr, "[INFO] Creating histogram: %s\n", hist_output);
        output_histogram(partitions, thread_count, hist_output, is_128);
    }
    
    fprintf(stderr, "[INFO] Outputting to: %s\n", strcmp(output_file, "-") == 0 ? "stdout" : output_file);
    output_kmers_smart(partitions, thread_count, kmer_size, part_bits, output_file, text_output,
                       use_filter ? &filter : NULL, is_128);
                       
    fprintf(stderr, "[INFO] Releasing resources...\n");
    release_partition_vec(partitions, is_128);
    free(partitions);
    fprintf(stderr, "[INFO] Complete!\n");
    return 0;
}
