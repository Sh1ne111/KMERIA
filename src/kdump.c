#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "ketopt.h"
#include "kstring.h"
#include "kdump.h"
#include "util.h"

static int read_kmeria_header(FILE *instream, uint32_t *kmer_len,
                              uint32_t *kmer_type, uint32_t *total_kmers) {
    char magic[4];
    if (fread(magic, 1, 4, instream) != 4 || strncmp(magic, "KMER", 4) != 0) {
        fprintf(stderr, "[ERROR] Invalid file format. Expected 'KMER' magic header.\n");
        return -1;
    }

    uint8_t version;
    uint8_t padding[3];
    if (fread(&version, 1, 1, instream) != 1 ||
        fread(padding, 1, 3, instream) != 3) {
        fprintf(stderr, "[ERROR] Could not read file header prefix.\n");
        return -1;
    }
    if (version != 2) {
        fprintf(stderr, "[WARNING] Unexpected KMERIA version %u, expected 2.\n", version);
    }

    // kmer_len(4) + type(4) + total_kmers(4)
    uint32_t metadata[3];
    if (fread(metadata, sizeof(uint32_t), 3, instream) != 3) {
        fprintf(stderr, "[ERROR] Could not read file metadata.\n");
        return -1;
    }

    *kmer_len   = metadata[0];
    *kmer_type  = metadata[1]; // 1=64-bit(k<32), 2=128-bit(k≥32)
    *total_kmers = metadata[2];

    if (*kmer_len < 2 || *kmer_len > 63) {
        fprintf(stderr, "[ERROR] Invalid k-mer length %u (must be 2-63).\n", *kmer_len);
        return -1;
    }

    int expected_type = (*kmer_len >= 32) ? 2 : 1;
    if (*kmer_type != (uint32_t)expected_type) {
        fprintf(stderr, "[WARNING] type field %u inconsistent with kmer_len %u"
                        " (expected type=%d).\n",
                *kmer_type, *kmer_len, expected_type);
    }

    return 0;
}

static void sync_filter_vectors(int_vector_t *vec, size_t target_size) {
    if (kv_size(*vec) < target_size) {
        int last_val = kv_A(*vec, kv_size(*vec) - 1);
        for (size_t i = kv_size(*vec); i < target_size; ++i) {
            kv_push(int, *vec, last_val);
        }
    } else {
        kv_size(*vec) = target_size;
    }
}

int convert_kmer_binary_to_text(FILE *instream, FILE *outstream,
                                int min_sample_pass, int max_sample_pass,
                                int_vector_t *min_count_filters,
                                int_vector_t *max_count_filters,
                                int output_format)
{
    uint32_t kmer_len, kmer_type, total_kmers;

    if (read_kmeria_header(instream, &kmer_len, &kmer_type, &total_kmers) != 0) {
        return -1;
    }

    int is_128 = (kmer_type == 2);

    uint32_t sample_count = 1;
    if (sample_count == 0) {
        output_format = FORMAT_KMER_ONLY;
    }

    sync_filter_vectors(min_count_filters, sample_count > 0 ? sample_count : 1);
    sync_filter_vectors(max_count_filters, sample_count > 0 ? sample_count : 1);

    fprintf(stderr, "[INFO] K-mer Length : %u\n", kmer_len);
    fprintf(stderr, "[INFO] Engine       : %s\n", is_128 ? "128-bit" : "64-bit");
    fprintf(stderr, "[INFO] Total K-mers : %u\n", total_kmers);

    kstring_t output_buffer = {0, 0, 0};
    ks_resize(&output_buffer, kmer_len + 64);

    uint64_t kmers_processed = 0;
    uint64_t kmers_written   = 0;

    if (is_128) {
#pragma pack(push, 2)
        typedef struct { uint128_t kmer_val; uint16_t freq; } pair128_t;
#pragma pack(pop)

        pair128_t record;
        char kmer_str[kmer_len + 1];
        kmer_str[kmer_len] = '\0';

        while (fread(&record, sizeof(pair128_t), 1, instream) == 1) {
            kmers_processed++;

            uint16_t cnt = record.freq;
            int min_thresh = kv_A(*min_count_filters, 0);
            int max_thresh = kv_A(*max_count_filters, 0);
            if (cnt < min_thresh || (max_thresh >= 0 && cnt > max_thresh)) {
                if (kmers_processed % 10000000 == 0) {
                    fprintf(stderr, "\r[PROGRESS] Scanned: %lu, Emitted: %lu",
                            kmers_processed, kmers_written);
                    fflush(stderr);
                }
                continue;
            }

            kmer2seq_128(kmer_str, kmer_len, record.kmer_val);

            output_buffer.l = 0;
            kputs(kmer_str, &output_buffer);

            if (output_format != FORMAT_KMER_ONLY) {
                kputc('\t', &output_buffer);
                kputw(cnt, &output_buffer);
            }
            kputc('\n', &output_buffer);
            fwrite(output_buffer.s, 1, output_buffer.l, outstream);
            kmers_written++;

            if (kmers_processed % 10000000 == 0) {
                fprintf(stderr, "\r[PROGRESS] Scanned: %lu, Emitted: %lu",
                        kmers_processed, kmers_written);
                fflush(stderr);
            }
        }

    } else {
#pragma pack(push, 2)
        typedef struct { uint64_t kmer_val; uint16_t freq; } pair64_t;
#pragma pack(pop)

        pair64_t record;
        char kmer_str[kmer_len + 1];
        kmer_str[kmer_len] = '\0';

        while (fread(&record, sizeof(pair64_t), 1, instream) == 1) {
            kmers_processed++;

            uint16_t cnt = record.freq;
            int min_thresh = kv_A(*min_count_filters, 0);
            int max_thresh = kv_A(*max_count_filters, 0);
            if (cnt < min_thresh || (max_thresh >= 0 && cnt > max_thresh)) {
                if (kmers_processed % 10000000 == 0) {
                    fprintf(stderr, "\r[PROGRESS] Scanned: %lu, Emitted: %lu",
                            kmers_processed, kmers_written);
                    fflush(stderr);
                }
                continue;
            }

            kmer2seq(kmer_str, kmer_len, record.kmer_val); 

            output_buffer.l = 0;
            kputs(kmer_str, &output_buffer);

            if (output_format != FORMAT_KMER_ONLY) {
                kputc('\t', &output_buffer);
                kputw(cnt, &output_buffer);
            }
            kputc('\n', &output_buffer);
            fwrite(output_buffer.s, 1, output_buffer.l, outstream);
            kmers_written++;

            if (kmers_processed % 10000000 == 0) {
                fprintf(stderr, "\r[PROGRESS] Scanned: %lu, Emitted: %lu",
                        kmers_processed, kmers_written);
                fflush(stderr);
            }
        }
    }

    fprintf(stderr, "\r[PROGRESS] Scanned: %lu, Emitted: %lu\n",
            kmers_processed, kmers_written);
    free(output_buffer.s);

    if (kmers_processed != total_kmers) {
        fprintf(stderr, "[WARNING] Header expected %u k-mers, file contained %lu.\n",
                total_kmers, kmers_processed);
    }

    return 0;
}


static void print_kdump_usage(FILE* stream) {
    fprintf(stream, "\n");
    fprintf(stream, "Program: kdump - A tool to convert binary k-mer files to text.\n");
    fprintf(stream, "Version: 2.0.1\n");
    fprintf(stream, "Usage:   kdump [options] <input.bin>\n\n");
    fprintf(stream, "Basic Options:\n");
    fprintf(stream, "  -o FILE       Output file path [default: standard output]\n");
    fprintf(stream, "  -O STR        Output format [default: ALL_COUNTS]\n");
    fprintf(stream, "                  ALL_COUNTS - k-mer followed by all counts\n");
    fprintf(stream, "                  SUM        - k-mer followed by sum of counts\n");
    fprintf(stream, "                  DIFFERENCE - k-mer followed by (cnt1 - cnt2 - ...)\n");
    fprintf(stream, "                  MIN        - k-mer followed by minimum count\n");
    fprintf(stream, "                  MAX        - k-mer followed by maximum count\n");
    fprintf(stream, "                  N          - k-mer followed by Nth sample count (1-based)\n\n");
    fprintf(stream, "Filtering Options:\n");
    fprintf(stream, "  -c INT[,...]  Minimum count threshold [default: 1]\n");
    fprintf(stream, "  -C INT[,...]  Maximum count threshold [default: no limit]\n");
    fprintf(stream, "  -n INT        Minimum samples passing filters [default: 1]\n");
    fprintf(stream, "  -N INT        Maximum samples passing filters [default: all]\n\n");
    fprintf(stream, "Other Options:\n");
    fprintf(stream, "  -h            Display this help message and exit\n");
    fprintf(stream, "  -v            Enable verbose output\n\n");
    fprintf(stream, "Supported k-mer lengths: 2-63\n\n");
}

static void parse_threshold_arg(const char* arg_str, int_vector_t* thresholds) {
    kstring_t s = {0};
    kputs(arg_str, &s);
    int n_tokens;
    int* tokens = ksplit(&s, ',', &n_tokens);
    for (int i = 0; i < n_tokens; ++i) {
        kv_push(int, *thresholds, atoi(s.s + tokens[i]));
    }
    free(s.s);
    free(tokens);
}

int kdump_main(int argc, char *argv[]) {
    int options_valid = 1;
    int min_req_samples = 1, max_req_samples = -1;
    int format_mode = FORMAT_ALL_COUNTS;
    int_vector_t min_thresholds, max_thresholds;
    kv_init(min_thresholds);
    kv_init(max_thresholds);

    FILE *input_stream = NULL, *output_stream = stdout;
    ketopt_t opt_parser = KETOPT_INIT;
    int c;

    while ((c = ketopt(&opt_parser, argc, argv, 1, "c:C:n:N:O:o:hv", 0)) >= 0) {
        switch (c) {
            case 'c': parse_threshold_arg(opt_parser.arg, &min_thresholds); break;
            case 'C': parse_threshold_arg(opt_parser.arg, &max_thresholds); break;
            case 'n': min_req_samples = atoi(opt_parser.arg); break;
            case 'N': max_req_samples = atoi(opt_parser.arg); break;
            case 'o':
                if ((output_stream = (strcmp(opt_parser.arg, "-") == 0)
                        ? stdout : fopen(opt_parser.arg, "w")) == NULL) {
                    fprintf(stderr, "[ERROR] Could not open output file: %s\n", opt_parser.arg);
                    options_valid = 0;
                }
                break;
            case 'O':
                for (char *p = opt_parser.arg; *p; ++p) *p = toupper(*p);
                if      (strcmp(opt_parser.arg, "ALL_COUNTS") == 0) format_mode = FORMAT_ALL_COUNTS;
                else if (strcmp(opt_parser.arg, "SUM")        == 0) format_mode = FORMAT_SUM;
                else if (strcmp(opt_parser.arg, "DIFFERENCE") == 0) format_mode = FORMAT_DIFFERENCE;
                else if (strcmp(opt_parser.arg, "MIN")        == 0) format_mode = FORMAT_MIN;
                else if (strcmp(opt_parser.arg, "MAX")        == 0) format_mode = FORMAT_MAX;
                else if (atoi(opt_parser.arg) > 0) format_mode = FORMAT_TYPE_BOUNDARY + atoi(opt_parser.arg);
                else {
                    fprintf(stderr, "[ERROR] Unrecognized output format: %s\n", opt_parser.arg);
                    options_valid = 0;
                }
                break;
            case 'h': options_valid = 0; break;
            case 'v': break;
            default:
                fprintf(stderr, "[ERROR] Unknown option: -%c\n",
                        opt_parser.opt ? opt_parser.opt : ':');
                options_valid = 0;
                break;
        }
    }

    if (kv_size(min_thresholds) == 0) kv_push(int, min_thresholds, 1);
    if (kv_size(max_thresholds) == 0) kv_push(int, max_thresholds, -1);

    if (argc - opt_parser.ind < 1 || !options_valid) {
        print_kdump_usage(stderr);
        kv_destroy(min_thresholds);
        kv_destroy(max_thresholds);
        return options_valid ? 0 : 1;
    }

    if (strcmp(argv[opt_parser.ind], "-") == 0) {
        input_stream = stdin;
        fprintf(stderr, "[INFO] Reading from standard input.\n");
    } else {
        input_stream = fopen(argv[opt_parser.ind], "rb");
        if (!input_stream) {
            fprintf(stderr, "[ERROR] Failed to open input file: %s\n", argv[opt_parser.ind]);
            options_valid = 0;
        }
    }

    if (!options_valid) {
        kv_destroy(min_thresholds);
        kv_destroy(max_thresholds);
        if (output_stream != stdout) fclose(output_stream);
        return 1;
    }

    fprintf(stderr, "[INFO] Processing: %s\n",
            (input_stream == stdin) ? "stdin" : argv[opt_parser.ind]);

    int result = convert_kmer_binary_to_text(
        input_stream, output_stream,
        min_req_samples, max_req_samples,
        &min_thresholds, &max_thresholds, format_mode);

    kv_destroy(min_thresholds);
    kv_destroy(max_thresholds);
    if (input_stream  && input_stream  != stdin)  fclose(input_stream);
    if (output_stream && output_stream != stdout) fclose(output_stream);

    fprintf(stderr, result == 0
            ? "[INFO] Conversion finished successfully.\n"
            : "[ERROR] An error occurred during conversion.\n");

    return result == 0 ? 0 : 1;
}

int kdump(int argc, char *argv[]) {
    return kdump_main(argc, argv);
}
