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

#ifndef KMATRIX_H
#define KMATRIX_H

#include <cinttypes>
#include <limits>
#include <memory>
#include <cstdio>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <queue>
#include <atomic>
#include <future>
#include <algorithm>
#include <getopt.h>
#include <iomanip>
#include <sstream>
#include <chrono>
#include <dirent.h>
#include <sys/stat.h>
#include "kmer_api.h"
#include "kmc_file.h"

struct ProgramConfig {
    std::string input_file_list;
    std::string output_prefix;
    bool text_output_mode = false;
    size_t batch_processing_size = 1000;
    size_t entries_per_block = 10000000;
    size_t write_buffer_threshold = 100000;
    int worker_thread_count = 0;
    bool verbose_logging = false;
    bool compress_output = true;
    std::string delimiter = "\t";
    bool include_headers = true;
    bool zero_based_indexing = false;
    int minimum_frequency_threshold = 1;
    bool memory_efficient_mode = true;
    size_t max_queue_depth = 5;
    bool monitor_memory = false;
    
    void display_configuration() const;
};

class KmerDataAccessor : public CKmerAPI {
public:
    explicit KmerDataAccessor(uint32_t sequence_length = 0);
    const uint64_t* extract_raw_data() const;
    uint32_t get_kmer_length() const;
    std::string to_string_representation() const;
};

using sequence_kmer = KmerDataAccessor;

// Enum for database format types
enum class DatabaseFormat {
    KMC,      // KMC database format
    KMERIA    // KMERIA binary format (kdump)
};

// Base class for database readers
class IDatabaseReader {
public:
    virtual ~IDatabaseReader() = default;
    virtual const std::string& get_database_path() const = 0;
    virtual uint32_t get_kmer_length() const = 0;
    virtual uint64_t get_total_kmers() const = 0;
    virtual bool is_exhausted() const = 0;
    virtual const sequence_kmer& get_current_kmer() const = 0;
    virtual uint32_t get_current_frequency() const = 0;
    virtual void advance_to_next() = 0;
    virtual DatabaseFormat get_format() const = 0;
};

// KMERIA format reader for kdump binary files
class KmeriaReader : public IDatabaseReader {
private:
    std::string database_path;
    FILE* file_handle;
    uint32_t kmer_length;
    uint32_t sample_count;
    uint64_t total_kmers;
    uint64_t current_position;
    sequence_kmer current_kmer;
    uint32_t current_frequency;
    bool exhausted;
    
    bool read_header();
    bool read_next_kmer();
    void kmer_from_binary(uint64_t binary_kmer);

public:
    explicit KmeriaReader(const std::string& db_path);
    ~KmeriaReader() noexcept override;
    
    const std::string& get_database_path() const override { return database_path; }
    uint32_t get_kmer_length() const override { return kmer_length; }
    uint64_t get_total_kmers() const override { return total_kmers; }
    bool is_exhausted() const override { return exhausted; }
    const sequence_kmer& get_current_kmer() const override { return current_kmer; }
    uint32_t get_current_frequency() const override { return current_frequency; }
    void advance_to_next() override;
    DatabaseFormat get_format() const override { return DatabaseFormat::KMERIA; }
};

// KMC format reader
class KMCReader : public IDatabaseReader {
private:
    std::unique_ptr<CKMCFile> database_handle;
    uint64_t total_sequence_count;
    uint64_t current_position;
    sequence_kmer current_sequence;
    uint32_t current_frequency;
    uint32_t kmer_length;
    std::string database_path;

public:
    explicit KMCReader(const std::string& db_path);
    ~KMCReader() noexcept override;
    
    const std::string& get_database_path() const override { return database_path; }
    uint32_t get_kmer_length() const override { return kmer_length; }
    uint64_t get_total_kmers() const override { return total_sequence_count; }
    bool is_exhausted() const override { return current_position >= total_sequence_count; }
    const sequence_kmer& get_current_kmer() const override { return current_sequence; }
    uint32_t get_current_frequency() const override { return current_frequency; }
    void advance_to_next() override;
    DatabaseFormat get_format() const override { return DatabaseFormat::KMC; }
};

// Unified database reader wrapper
class DatabaseReader {
private:
    std::unique_ptr<IDatabaseReader> reader;

public:
    DatabaseReader(DatabaseReader&&) = default;
    DatabaseReader& operator=(DatabaseReader&&) = default;
    explicit DatabaseReader(const std::string& db_path);
    
    const std::string& get_database_path() const { return reader->get_database_path(); }
    uint32_t get_kmer_length() const { return reader->get_kmer_length(); }
    uint64_t get_total_kmers() const { return reader->get_total_kmers(); }
    bool is_exhausted() const { return reader->is_exhausted(); }
    const sequence_kmer& get_current_kmer() const { return reader->get_current_kmer(); }
    uint32_t get_current_frequency() const { return reader->get_current_frequency(); }
    void advance_to_next() { reader->advance_to_next(); }
    DatabaseFormat get_format() const { return reader->get_format(); }
};

struct KmerMatrixRow {
    std::vector<uint64_t> kmer_representation;
    std::vector<uint32_t> frequency_counts;
    std::string kmer_sequence;
    
    KmerMatrixRow(size_t num_words, size_t num_samples, bool store_sequence = false);
    bool meets_frequency_threshold(int threshold) const;
    size_t count_non_zero_samples() const;
    uint32_t get_max_frequency() const;
};

class FrequencyCompressor {
public:
    static uint8_t determine_optimal_bits(const std::vector<uint32_t>& frequencies);
    static uint8_t determine_optimal_bits(uint32_t max_value);
    static std::vector<uint8_t> encode_frequencies_8bit(const std::vector<uint32_t>& frequencies);
    static std::vector<uint16_t> encode_frequencies_16bit(const std::vector<uint32_t>& frequencies);
    
    struct SparseEncoding {
        std::vector<uint16_t> indices;
        std::vector<uint32_t> values;
        size_t encoded_size() const {
            return indices.size() * sizeof(uint16_t) + values.size() * sizeof(uint32_t);
        }
    };
    
    static SparseEncoding encode_sparse(const std::vector<uint32_t>& frequencies);
    static std::vector<uint32_t> decode_sparse(const SparseEncoding& encoded, size_t total_samples);
    
    enum EncodingType : uint8_t {
        FULL_32BIT = 0,
        FULL_16BIT = 1,
        FULL_8BIT = 2,
        SPARSE = 3,
        RLE = 4
    };
    
    static EncodingType choose_best_encoding(const std::vector<uint32_t>& frequencies,
                                            size_t& estimated_size);
};

template<typename DataType>
class ThreadSafeQueue {
private:
    std::queue<DataType> data_queue;
    mutable std::mutex queue_mutex;
    std::condition_variable data_available;
    std::condition_variable space_available;
    std::atomic<bool> processing_complete{false};
    size_t max_queue_size;

public:
    explicit ThreadSafeQueue(size_t max_size = 10) : max_queue_size(max_size) {}
    void enqueue(DataType item);
    bool dequeue(DataType& item);
    void signal_completion();
    size_t size() const;
};

class OutputManager {
private:
    const ProgramConfig& config;
    uint64_t current_block_id;
    std::mutex write_mutex;
    uint64_t total_rows_written;
    
    std::unique_ptr<std::ofstream> current_file_stream;
    uint64_t current_file_row_count;
    std::string current_filename;
    bool file_is_open;
    std::streampos row_count_position;
    
    uint64_t total_uncompressed_bytes;
    uint64_t total_compressed_bytes;
    
    void write_text_rows_to_stream(const std::vector<KmerMatrixRow>& rows);
    void write_binary_rows_to_stream(const std::vector<KmerMatrixRow>& rows,
                                    uint32_t kmer_size, uint32_t sample_count);
    void write_compressed_rows_to_stream(const std::vector<KmerMatrixRow>& rows,
                                        uint32_t kmer_size, uint32_t sample_count);
    
    void open_new_output_file(uint32_t kmer_size, uint32_t sample_count, 
                             const std::vector<std::string>& sample_names);
    void close_current_file();
    void append_to_current_file(const std::vector<KmerMatrixRow>& data_block,
                               uint32_t kmer_size, uint32_t sample_count,
                               const std::vector<std::string>& sample_names);
    void finalize_current_file();

public:
    explicit OutputManager(const ProgramConfig& cfg);
    ~OutputManager();
    
    void write_data_block(const std::vector<KmerMatrixRow>& data_block, 
                         uint32_t kmer_size, uint32_t sample_count,
                         const std::vector<std::string>& sample_names);
    
    void finalize();
    uint64_t get_total_rows_written() const;
    double get_compression_ratio() const;
};

class ParallelKmerProcessor {
private:
    std::vector<DatabaseReader> database_readers;
    ThreadSafeQueue<std::vector<KmerMatrixRow>> processing_queue;
    OutputManager file_manager;
    const ProgramConfig& config;
    uint32_t kmer_size;
    uint32_t sample_count;
    uint32_t words_per_kmer;
    std::vector<std::string> sample_names;
    std::atomic<uint64_t> processed_kmers{0};
    
    static constexpr size_t DEFAULT_BATCH_SIZE = 1000;
    static constexpr size_t WRITE_BUFFER_THRESHOLD = 100000;

    void generate_sample_names();
    void print_database_statistics();
    void producer_thread();
    void writer_thread();
    size_t find_minimum_kmer_source();
    KmerMatrixRow build_matrix_row(const sequence_kmer& target_kmer);

public:
    ParallelKmerProcessor(std::vector<DatabaseReader>&& readers, const ProgramConfig& cfg);
    void execute_processing();
};

bool validate_kmer_consistency(const std::vector<DatabaseReader>& database_list);
void display_help_message(const char* program_name);
ProgramConfig parse_command_arguments(int argc, char** argv);
DatabaseFormat detect_database_format(const std::string& path);

class BinaryToTextConverter {
private:
    std::string input_file;
    std::string output_file;
    std::string delimiter;
    bool include_header;
    bool verbose;
    
    struct BinaryFileHeader {
        uint32_t kmer_size;
        uint32_t sample_count;
        uint64_t total_kmers;
        uint8_t compression_enabled; // compression flag
    };
    
    BinaryFileHeader read_binary_header(std::ifstream& file);
    std::string decode_kmer(const std::vector<uint64_t>& binary_data, uint32_t kmer_length);

public:
    BinaryToTextConverter(const std::string& input, const std::string& output,
                         const std::string& delim = "\t", bool header = true, bool verb = false);
    bool convert();
    void display_conversion_stats(uint64_t kmers_converted, double time_seconds);
};

class BatchBinaryConverter {
private:
    std::string input_dir;
    std::string output_dir;
    std::string delimiter;
    bool include_header;
    bool verbose;
    int max_threads;
    
    std::vector<std::string> find_binary_files();
    void convert_single_file(const std::string& input_file, const std::string& output_file);

public:
    BatchBinaryConverter(const std::string& in_dir, const std::string& out_dir,
                        const std::string& delim = "\t", bool header = true,
                        bool verb = false, int threads = 1);
    bool convert_all();
};

int run_conversion_mode(int argc, char** argv);

// Template implementations
template<typename DataType>
void ThreadSafeQueue<DataType>::enqueue(DataType item) {
    std::unique_lock<std::mutex> lock(queue_mutex);
    space_available.wait(lock, [this]{ 
        return data_queue.size() < max_queue_size || processing_complete.load(); 
    });
    
    if (!processing_complete.load()) {
        data_queue.push(std::move(item));
        data_available.notify_one();
    }
}

template<typename DataType>
bool ThreadSafeQueue<DataType>::dequeue(DataType& item) {
    std::unique_lock<std::mutex> lock(queue_mutex);
    data_available.wait(lock, [this]{ 
        return !data_queue.empty() || processing_complete.load(); 
    });
    
    if (data_queue.empty()) {
        return false;
    }
    
    item = std::move(data_queue.front());
    data_queue.pop();
    space_available.notify_one();
    return true;
}

template<typename DataType>
void ThreadSafeQueue<DataType>::signal_completion() {
    processing_complete.store(true);
    data_available.notify_all();
    space_available.notify_all();
}

template<typename DataType>
size_t ThreadSafeQueue<DataType>::size() const {
    std::lock_guard<std::mutex> lock(queue_mutex);
    return data_queue.size();
}

#endif // KMATRIX_H
