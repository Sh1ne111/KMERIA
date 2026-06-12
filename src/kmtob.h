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

#ifndef KMER_BIMBAM_H
#define KMER_BIMBAM_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <limits>
#include <algorithm>
#include <getopt.h>
#include <thread>
#include <future>
#include <stdexcept>
#include <cstdlib>
#include <dirent.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <memory>
#include <atomic>
#include <mutex>
#include <chrono>
#include <iomanip>
#include <cmath>
#include <numeric>

extern "C" {
#include "../include/bgzf.h"
#include "../include/kstring.h"
}

class BgzfOutputStream {
private:
    BGZF* bgzf_file;
    std::string filename;
    bool is_open;
    std::string buffer;
    size_t buffer_size;
    static constexpr size_t DEFAULT_BUFFER_SIZE = 131072; // 128KB buffer
    
    void flush_buffer();
    
public:
    explicit BgzfOutputStream(const std::string& fname, int compression_level = 6, int threads = 4);
    ~BgzfOutputStream();
    
    bool open(const std::string& fname, int compression_level = 6, int threads = 4);
    void close();
    bool is_valid() const { return is_open && bgzf_file != nullptr; }
    
    bool write_string(const std::string& str);
    void set_buffer_size(size_t size) { buffer_size = size; }
    
    // Stream operators
    BgzfOutputStream& operator<<(const std::string& str);
    BgzfOutputStream& operator<<(const char* str);
    BgzfOutputStream& operator<<(char c);
    BgzfOutputStream& operator<<(int val);
    BgzfOutputStream& operator<<(long val);
    BgzfOutputStream& operator<<(long long val);
    BgzfOutputStream& operator<<(unsigned int val);
    BgzfOutputStream& operator<<(unsigned long val);
    BgzfOutputStream& operator<<(unsigned long long val);
    BgzfOutputStream& operator<<(float val);
    BgzfOutputStream& operator<<(double val);
    BgzfOutputStream& operator<<(std::ostream& (*manip)(std::ostream&));
};

struct BimbamConfig {
    std::string input_dir = "filtered_kmatrices";
    std::string output_dir = "bimbam_kmatrices";
    std::string input_prefix = "filtered_";
    int max_threads = 8;
    bool verbose = false;
    std::string input_format = "auto";
    bool normalize_values = true;
    double min_range = 0.0;
    double max_range = 2.0;
    // Quantile norm config
    bool use_quantile_normalization = false;
    double lower_quantile = 0.05;            
    double upper_quantile = 0.95;   

    bool output_stats = false;
    std::string allele_info = "X, Y";
    
    // BGZF compression options
    bool compress_output = true;       
    int compression_level = 6;          
    int bgzf_threads = 4;                
    size_t write_buffer_size = 131072;
    
    void display_config() const;
    bool validate() const;
};

struct BimbamKmerData {
    std::string kmer_sequence;
    uint32_t kmer_length;
    std::vector<uint64_t> kmer_binary;
    std::vector<double> frequencies;
    //uint32_t kmer_length;
    
    BimbamKmerData() : kmer_length(0) {}
    BimbamKmerData(uint32_t k_len, size_t num_samples) 
        : kmer_length(k_len), frequencies(num_samples, 0.0) {
        uint32_t words = (k_len + 31) >> 5;
        kmer_binary.resize(words, 0);
    }
    
    std::string binary_to_string() const;
    void normalize_frequencies(double min_val, double max_val);

    void quantile_normalize_frequencies(double lower_quantile, 
                                       double upper_quantile,
                                       double min_val, 
                                       double max_val);    
    struct KmerStats {
        double min_freq, max_freq, mean_freq, std_dev;
        size_t zero_count, non_zero_count;
    };
    KmerStats get_statistics() const;
};

class FrequencyDecompressor {
public:
    enum EncodingType : uint8_t {
        FULL_32BIT = 0,
        FULL_16BIT = 1,
        FULL_8BIT = 2,
        SPARSE = 3,
        RLE = 4
    };
    
    struct SparseEncoding {
        std::vector<uint16_t> indices;
        std::vector<uint32_t> values;
    };
    
    static std::vector<uint32_t> decode_8bit(const std::vector<uint8_t>& encoded);
    static std::vector<uint16_t> decode_16bit(const std::vector<uint16_t>& encoded);
    static std::vector<uint32_t> decode_sparse(const SparseEncoding& encoded, size_t total_samples);
};

class CompressedBinaryKmerReader {
private:
    std::ifstream file_stream;
    uint32_t kmer_size;
    uint32_t sample_count;
    uint64_t total_kmers;
    uint64_t current_position;
    uint32_t words_per_kmer;
    bool is_valid;
    bool is_compressed;
    
    bool read_uncompressed_kmer(BimbamKmerData& kmer_data);
    bool read_compressed_kmer(BimbamKmerData& kmer_data);

public:
    explicit CompressedBinaryKmerReader(const std::string& filename);
    ~CompressedBinaryKmerReader();
    
    bool is_open() const { return is_valid; }
    uint32_t get_kmer_size() const { return kmer_size; }
    uint32_t get_sample_count() const { return sample_count; }
    uint64_t get_total_kmers() const { return total_kmers; }
    bool has_more_data() const { return current_position < total_kmers; }
    bool compressed() const { return is_compressed; }
    
    bool read_next_kmer(BimbamKmerData& kmer_data);
    void reset_to_beginning();
};

class BinaryKmerReader {
private:
    std::ifstream file_stream;
    uint32_t kmer_size;
    uint32_t sample_count;
    uint64_t total_kmers;
    uint64_t current_position;
    uint32_t words_per_kmer;
    bool is_valid;

public:
    explicit BinaryKmerReader(const std::string& filename);
    ~BinaryKmerReader();
    
    bool is_open() const { return is_valid; }
    uint32_t get_kmer_size() const { return kmer_size; }
    uint32_t get_sample_count() const { return sample_count; }
    uint64_t get_total_kmers() const { return total_kmers; }
    bool has_more_data() const { return current_position < total_kmers; }
    
    bool read_next_kmer(BimbamKmerData& kmer_data);
    void reset_to_beginning();
};

class BimbamUtils {
public:
    static std::vector<double> normalize_to_range(const std::vector<double>& numbers,
                                                 double min_range, double max_range);

    static std::vector<double> quantile_normalize(const std::vector<double>& numbers,
                                                  double lower_quantile,
                                                  double upper_quantile,
                                                  double min_range,
                                                  double max_range);

    static double calculate_quantile(const std::vector<double>& sorted_data, double quantile);    
    static double calculate_mean(const std::vector<double>& numbers);
    static double calculate_std_dev(const std::vector<double>& numbers, double mean);
    static std::pair<double, double> get_min_max(const std::vector<double>& numbers);
};

class TextBimbamProcessor {
public:
    static bool process_text_file(const std::string& input_file,
                                 const std::string& output_file,
                                 const BimbamConfig& config,
                                 uint64_t& processed_count);
};

class BinaryBimbamProcessor {
public:
    static bool process_binary_file(const std::string& input_file,
                                   const std::string& output_file,
                                   const BimbamConfig& config,
                                   uint64_t& processed_count);
};

class EnhancedBimbamConverter {
private:
    BimbamConfig config;
    std::atomic<uint64_t> processed_files{0};
    std::atomic<uint64_t> total_kmers_processed{0};
    std::atomic<uint64_t> failed_files{0};
    std::atomic<uint64_t> total_bytes_written{0};
    std::mutex progress_mutex;

    std::vector<std::string> discover_input_files();
    std::string determine_file_format(const std::string& filename);
    std::string generate_output_filename(const std::string& input_file);
    void process_single_file(const std::string& input_file, const std::string& output_file);
    void update_progress(const std::string& filename, uint64_t processed, bool success);

public:
    explicit EnhancedBimbamConverter(const BimbamConfig& cfg);
    
    bool execute_conversion();
    void display_statistics() const;
};

void display_bimbam_help_message(const char* program_name);
BimbamConfig parse_bimbam_arguments(int argc, char** argv);
bool create_output_directory(const std::string& dir_path);
std::vector<std::string> get_files_with_prefix(const std::string& dir, 
                                              const std::string& prefix,
                                              const std::vector<std::string>& extensions);
std::string detect_file_format(const std::string& filename);
bool validate_bimbam_output(const std::string& filename);

struct FileProcessingStats {
    uint64_t total_files;
    uint64_t successful_files;
    uint64_t failed_files;
    uint64_t total_kmers;
    double processing_time_seconds;
    
    void display() const;
};

extern "C" int kmtob(int argc, char* argv[]);

#endif
