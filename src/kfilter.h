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


#ifndef KMER_FILTER_H
#define KMER_FILTER_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <limits>
#include <algorithm>
#include <getopt.h>
#include <cmath>
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

struct FilterConfig {
    std::string input_dir = "kmatrices";
    std::string output_dir = "filtered_kmatrices";
    std::string depth_file;
    int max_processes = 8;
    int max_coverage = 1000;
    double missing_rate = 0.8;
    int ploidy = 4;
    bool verbose = false;
    std::string input_format = "auto"; // auto, text, binary, compressed
    std::string output_format = "text"; // text, binary, compressed
    
    void display_config() const;
    bool validate() const;
};

struct KmerData {
    std::string kmer_sequence;
    std::vector<uint64_t> kmer_binary;
    std::vector<uint32_t> frequencies;
    uint32_t kmer_length;
    
    KmerData() : kmer_length(0) {}
    KmerData(uint32_t k_len, size_t num_samples) 
        : kmer_length(k_len), frequencies(num_samples, 0) {
        uint32_t words = (k_len + 31) >> 5;
        kmer_binary.resize(words, 0);
    }
    
    std::string binary_to_string() const;
    void string_to_binary(const std::string& seq);
    bool apply_filtering(const std::vector<double>& depths, 
                        const std::vector<int>& ploidies,
                        const FilterConfig& config, 
                        double min_haploid_depth);
};

class FrequencyCompressor {
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
    
    static uint8_t determine_optimal_bits(const std::vector<uint32_t>& frequencies);
    static uint8_t determine_optimal_bits(uint32_t max_value);
    static std::vector<uint8_t> encode_frequencies_8bit(const std::vector<uint32_t>& frequencies);
    static std::vector<uint16_t> encode_frequencies_16bit(const std::vector<uint32_t>& frequencies);
    static SparseEncoding encode_sparse(const std::vector<uint32_t>& frequencies);
    static std::vector<uint32_t> decode_sparse(const SparseEncoding& encoded, size_t total_samples);
    static EncodingType choose_best_encoding(const std::vector<uint32_t>& frequencies, size_t& estimated_size);
    
    // Decoding functions
    static std::vector<uint32_t> decode_8bit(const std::vector<uint8_t>& encoded);
    static std::vector<uint32_t> decode_16bit(const std::vector<uint16_t>& encoded);
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
    
    bool read_uncompressed_kmer(KmerData& kmer_data);
    bool read_compressed_kmer(KmerData& kmer_data);

public:
    explicit CompressedBinaryKmerReader(const std::string& filename);
    ~CompressedBinaryKmerReader();
    
    bool is_open() const { return is_valid; }
    uint32_t get_kmer_size() const { return kmer_size; }
    uint32_t get_sample_count() const { return sample_count; }
    uint64_t get_total_kmers() const { return total_kmers; }
    bool has_more_data() const { return current_position < total_kmers; }
    bool compressed() const { return is_compressed; }
    
    bool read_next_kmer(KmerData& kmer_data);
    void reset_to_beginning();
};

class CompressedBinaryKmerWriter {
private:
    std::ofstream file_stream;
    uint32_t kmer_size;
    uint32_t sample_count;
    uint64_t written_kmers;
    uint32_t words_per_kmer;
    std::streampos kmer_count_position;
    bool header_written;
    bool enable_compression;
    
    // Compression statistics
    uint64_t total_uncompressed_bytes;
    uint64_t total_compressed_bytes;
    
    void write_uncompressed_kmer(const KmerData& kmer_data);
    void write_compressed_kmer(const KmerData& kmer_data);

public:
    CompressedBinaryKmerWriter(const std::string& filename, uint32_t k_size, 
                              uint32_t num_samples, bool compress = false);
    ~CompressedBinaryKmerWriter();
    
    bool is_open() const { return file_stream.is_open(); }
    bool write_kmer(const KmerData& kmer_data);
    void finalize();
    uint64_t get_written_count() const { return written_kmers; }
    double get_compression_ratio() const;
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
    
    bool read_next_kmer(KmerData& kmer_data);
    void reset_to_beginning();
};

class BinaryKmerWriter {
private:
    std::ofstream file_stream;
    uint32_t kmer_size;
    uint32_t sample_count;
    uint64_t written_kmers;
    uint32_t words_per_kmer;
    std::streampos kmer_count_position;
    bool header_written;

public:
    BinaryKmerWriter(const std::string& filename, uint32_t k_size, uint32_t num_samples);
    ~BinaryKmerWriter();
    
    bool is_open() const { return file_stream.is_open(); }
    bool write_kmer(const KmerData& kmer_data);
    void finalize();
    uint64_t get_written_count() const { return written_kmers; }
};

class TextKmerProcessor {
public:
    static bool process_text_file(const std::string& input_file,
                                 const std::string& output_file,
                                 const std::vector<double>& depths,
                                 const std::vector<int>& ploidies,
                                 double min_haploid_depth,
                                 const FilterConfig& config);
};

class EnhancedKmerProcessor {
private:
    FilterConfig config;
    std::vector<double> sample_depths;
    std::vector<int> sample_ploidies;
    double min_haploid_depth;
    std::atomic<uint64_t> processed_files{0};
    std::atomic<uint64_t> total_kmers_processed{0};
    std::atomic<uint64_t> total_kmers_kept{0};
    std::mutex progress_mutex;

    void load_depth_information();
    std::vector<std::string> discover_input_files();
    std::string determine_file_format(const std::string& filename);
    std::string generate_output_filename(const std::string& input_file);
    void process_single_file(const std::string& input_file, const std::string& output_file);
    void update_progress(const std::string& filename, uint64_t processed, uint64_t kept);

public:
    explicit EnhancedKmerProcessor(const FilterConfig& cfg);
    
    bool execute_filtering();
    void display_statistics() const;
};

void display_help_message(const char* program_name);
FilterConfig parse_filter_arguments(int argc, char** argv);
bool create_output_directory(const std::string& dir_path);
std::vector<std::string> get_files_with_extensions(const std::string& dir, 
                                                  const std::vector<std::string>& extensions);
std::string detect_file_format(const std::string& filename);

extern "C" int kflt(int argc, char* argv[]);

#endif
