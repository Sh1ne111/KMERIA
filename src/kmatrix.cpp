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


#include "kmatrix.h"
#include <fstream>
#include <iostream>
#include <cstring>

void ProgramConfig::display_configuration() const {
    std::cout << "=== Program Configuration ===\n";
    std::cout << "Input file list: " << input_file_list << "\n";
    std::cout << "Output prefix: " << output_prefix << "\n";
    std::cout << "Output format: " << (text_output_mode ? "Text format" : "Binary format") << "\n";
    if (!text_output_mode) {
        std::cout << "Compression: " << (compress_output ? "Enabled (auto)" : "Disabled") << "\n";
    }
    std::cout << "Batch size: " << batch_processing_size << "\n";
    std::cout << "Entries per block: " << entries_per_block << "\n";
    std::cout << "Buffer threshold: " << write_buffer_threshold << "\n";
    std::cout << "Worker threads: " << (worker_thread_count == 0 ? "Auto-detect" : std::to_string(worker_thread_count)) << "\n";
    std::cout << "Verbose logging: " << (verbose_logging ? "Enabled" : "Disabled") << "\n";
    std::cout << "Min frequency threshold: " << minimum_frequency_threshold << "\n";
    std::cout << "Memory efficient mode: " << (memory_efficient_mode ? "Enabled" : "Disabled") << "\n";
    std::cout << "Max queue depth: " << max_queue_depth << "\n";
    if (text_output_mode) {
        std::cout << "Field delimiter: '" << delimiter << "'\n";
        std::cout << "Include headers: " << (include_headers ? "Yes" : "No") << "\n";
    }
    std::cout << "==============================\n\n";
}

KmerDataAccessor::KmerDataAccessor(uint32_t sequence_length) : CKmerAPI(sequence_length) {}

const uint64_t* KmerDataAccessor::extract_raw_data() const { 
    return kmer_data; 
}

uint32_t KmerDataAccessor::get_kmer_length() const {
    return kmer_length;
}

/*
std::string KmerDataAccessor::to_string_representation() const {
    std::string result;
    const_cast<KmerDataAccessor*>(this)->to_string(result);
    return result;
}
*/

static constexpr char NUCLEOTIDE_LOOKUP[4] = {'A', 'C', 'G', 'T'};

std::string KmerDataAccessor::to_string_representation() const {
    std::string result;
    uint32_t k_len = this->kmer_length;
    result.reserve(k_len);
    for(int32_t i = k_len - 1; i >= 0; --i){
//    for (uint32_t i = 0; i < k_len; ++i) {
        uint32_t word_idx = i >> 5; 
        uint32_t bit_offset = (i & 31) << 1;
        uint32_t symbol = (kmer_data[word_idx] >> bit_offset) & 3;
        result += NUCLEOTIDE_LOOKUP[symbol];
    }
    
    return result;
}


DatabaseFormat detect_database_format(const std::string& path) {
    // Check for KMERIA binary format
    FILE* test_file = fopen(path.c_str(), "rb");
    if (test_file) {
        char magic[8];
//        if (fread(magic, 1, 7, test_file) == 7) {
	if (fread(magic, 1, 4, test_file) == 4) {
            fclose(test_file);
            if (strncmp(magic, "KMER", 4) == 0) {
                return DatabaseFormat::KMERIA;
            }
        } else {
            fclose(test_file);
        }
    }
    
    // Check for KMC format (has .kmc_pre and .kmc_suf files)
    std::string pre_file = path + ".kmc_pre";
    std::string suf_file = path + ".kmc_suf";
    
    std::ifstream pre(pre_file);
    std::ifstream suf(suf_file);
    
    if (pre.good() && suf.good()) {
        return DatabaseFormat::KMC;
    }
    
    // Check if path already ends with .kmc_pre
    if (path.size() > 8 && path.substr(path.size() - 8) == ".kmc_pre") {
        std::string base = path.substr(0, path.size() - 8);
        suf_file = base + ".kmc_suf";
        std::ifstream suf_check(suf_file);
        if (suf_check.good()) {
            return DatabaseFormat::KMC;
        }
    }
    
    // Default to KMERIA if file exists
    std::ifstream test(path);
    if (test.good()) {
        return DatabaseFormat::KMERIA;
    }
    
    return DatabaseFormat::KMC;
}

KmeriaReader::KmeriaReader(const std::string& db_path)
    : database_path(db_path), file_handle(nullptr), kmer_length(0),
      sample_count(0), total_kmers(0), current_position(0),
      current_frequency(0), exhausted(false) {
    
    file_handle = fopen(database_path.c_str(), "rb");
    if (!file_handle) {
        std::cerr << "Error: Cannot open KMERIA database " << database_path << "\n";
        exit(EXIT_FAILURE);
    }
    
    if (!read_header()) {
        std::cerr << "Error: Invalid KMERIA format in " << database_path << "\n";
        fclose(file_handle);
        exit(EXIT_FAILURE);
    }
    
    current_kmer = sequence_kmer(kmer_length);
    
    // Read first k-mer
    if (!read_next_kmer()) {
        exhausted = true;
    }
}

KmeriaReader::~KmeriaReader() noexcept {
    if (file_handle) {
        fclose(file_handle);
    }
}

/*
bool KmeriaReader::read_header() {
    char magic[4];
    if (fread(magic, 1, 4, file_handle) != 4) {
        return false;
    }
    
    if (strncmp(magic, "KMER", 4) != 0) {
        std::cerr << "Error: Invalid KMERIA magic header\n";
        return false;
    }
    
    uint32_t metadata[3];
    if (fread(metadata, sizeof(uint32_t), 3, file_handle) != 3) {
        std::cerr << "Error: Could not read KMERIA metadata\n";
        return false;
    }
    
    kmer_length = metadata[0];
    sample_count = metadata[1];
    total_kmers = metadata[2];
    
    if (kmer_length > 32) {
        std::cerr << "Error: K-mer length " << kmer_length << " exceeds maximum of 32\n";
        return false;
    }
    
    return true;
}
*/

bool KmeriaReader::read_header() {
    char magic[4];
    if (fread(magic, 1, 4, file_handle) != 4) {
        return false;
    }
    if (strncmp(magic, "KMER", 4) != 0) {
        std::cerr << "Error: Invalid KMERIA magic header\n";
        return false;
    }
    
    uint8_t version;
    uint8_t padding[3];
    if (fread(&version, 1, 1, file_handle) != 1 ||
        fread(padding, 1, 3, file_handle) != 3) {
        std::cerr << "Error: Could not read KMERIA file header prefix\n";
        return false;
    }
    if (version != 2) {
        std::cerr << "Warning: Unexpected KMERIA version " << (int)version 
                  << ", expected 2\n";
    }
    
    uint32_t metadata[3];
    if (fread(metadata, sizeof(uint32_t), 3, file_handle) != 3) {
        std::cerr << "Error: Could not read KMERIA metadata\n";
        return false;
    }
    
    kmer_length  = metadata[0];  // kmer_len，偏移8
    sample_count = metadata[1];  // type (1=64bit/2=128bit)，偏移12
    total_kmers  = metadata[2];  // total_pairs，偏移16
    
    if (kmer_length > 32) {
        std::cerr << "Error: K-mer length " << kmer_length << " exceeds maximum of 32\n";
        return false;
    }
    
    return true;
}



void KmeriaReader::kmer_from_binary(uint64_t binary_kmer) {
    // Convert binary k-mer to CKmerAPI format
    uint64_t* kmer_data = const_cast<uint64_t*>(current_kmer.extract_raw_data());
    
    // KMERIA stores k-mer in single 64-bit word with 2 bits per base
    kmer_data[0] = binary_kmer;
    
    // Clear remaining words if k-mer length requires multiple words
    uint32_t num_words = (kmer_length + 31) >> 5;
    for (uint32_t i = 1; i < num_words; ++i) {
        kmer_data[i] = 0;
    }
}

bool KmeriaReader::read_next_kmer() {
    if (current_position >= total_kmers) {
        exhausted = true;
        return false;
    }
    
    uint64_t binary_kmer;
    if (fread(&binary_kmer, sizeof(uint64_t), 1, file_handle) != 1) {
        exhausted = true;
        return false;
    }
    
    if (sample_count > 0) {
        std::vector<uint16_t> counts(sample_count);
        if (fread(counts.data(), sizeof(uint16_t), sample_count, file_handle) != sample_count) {
            exhausted = true;
            return false;
        }
        
        // Use first sample's count (or could sum/max them)
        current_frequency = counts[0];
    } else {
        current_frequency = 1;
    }
    
    kmer_from_binary(binary_kmer);
    current_position++;
    
    return true;
}

void KmeriaReader::advance_to_next() {
    if (!read_next_kmer()) {
        exhausted = true;
    }
}

KMCReader::KMCReader(const std::string& db_path) 
    : database_path(db_path), current_position(0) {
    database_handle = std::make_unique<CKMCFile>();
    
    if (!database_handle->OpenForListing(database_path)) {
        std::cerr << "Error: Cannot open KMC database " << database_path << "\n";
        exit(EXIT_FAILURE);
    }
    
    if (database_handle->IsKMC2()) {
        std::cerr << "Error: Need sorted KMC database " << database_path << "\n";
        std::cerr << "DEBUG: Use kmc_tools transform <input_database> sort <output_database> for sorting\n";
        exit(EXIT_FAILURE);
    }
    
    CKMCFileInfo database_info;
    database_handle->Info(database_info);
    kmer_length = database_info.kmer_length;
    total_sequence_count = database_info.total_kmers;
    
    current_sequence = sequence_kmer(kmer_length);
    
    if (!is_exhausted()) {
        advance_to_next();
    }
}

KMCReader::~KMCReader() noexcept {
    if (database_handle) {
        database_handle->Close();
    }
}

void KMCReader::advance_to_next() {
    ++current_position;
    if (!database_handle->ReadNextKmer(current_sequence, current_frequency)) {
        std::cerr << "Critical error: k-mer reading failed " << __FILE__ << ":" << __LINE__ << "\n";
        exit(EXIT_FAILURE);
    }
}

DatabaseReader::DatabaseReader(const std::string& db_path) {
    DatabaseFormat format = detect_database_format(db_path);
    
    if (format == DatabaseFormat::KMERIA) {
        reader = std::make_unique<KmeriaReader>(db_path);
    } else {
        reader = std::make_unique<KMCReader>(db_path);
    }
}

KmerMatrixRow::KmerMatrixRow(size_t num_words, size_t num_samples, bool store_sequence)
    : kmer_representation(num_words), frequency_counts(num_samples, 0) {
    if (store_sequence) {
        kmer_sequence.reserve(32);
    }
}

bool KmerMatrixRow::meets_frequency_threshold(int threshold) const {
    return std::any_of(frequency_counts.begin(), frequency_counts.end(),
                      [threshold](uint32_t count) { return static_cast<int>(count) >= threshold; });
}

size_t KmerMatrixRow::count_non_zero_samples() const {
    return std::count_if(frequency_counts.begin(), frequency_counts.end(),
                       [](uint32_t count) { return count > 0; });
}

uint32_t KmerMatrixRow::get_max_frequency() const {
    if (frequency_counts.empty()) return 0;
    return *std::max_element(frequency_counts.begin(), frequency_counts.end());
}

// FrequencyCompressor implementation
uint8_t FrequencyCompressor::determine_optimal_bits(const std::vector<uint32_t>& frequencies) {
    if (frequencies.empty()) return 32;
    uint32_t max_val = *std::max_element(frequencies.begin(), frequencies.end());
    return determine_optimal_bits(max_val);
}

uint8_t FrequencyCompressor::determine_optimal_bits(uint32_t max_value) {
    if (max_value <= 255) return 8;
    if (max_value <= 65535) return 16;
    return 32;
}

std::vector<uint8_t> FrequencyCompressor::encode_frequencies_8bit(const std::vector<uint32_t>& frequencies) {
    std::vector<uint8_t> encoded;
    encoded.reserve(frequencies.size());
    for (uint32_t freq : frequencies) {
        encoded.push_back(static_cast<uint8_t>(std::min(freq, 255u)));
    }
    return encoded;
}

std::vector<uint16_t> FrequencyCompressor::encode_frequencies_16bit(const std::vector<uint32_t>& frequencies) {
    std::vector<uint16_t> encoded;
    encoded.reserve(frequencies.size());
    for (uint32_t freq : frequencies) {
        encoded.push_back(static_cast<uint16_t>(std::min(freq, 65535u)));
    }
    return encoded;
}

FrequencyCompressor::SparseEncoding FrequencyCompressor::encode_sparse(const std::vector<uint32_t>& frequencies) {
    SparseEncoding result;
    
    for (size_t i = 0; i < frequencies.size(); ++i) {
        if (frequencies[i] > 0) {
            result.indices.push_back(static_cast<uint16_t>(i));
            result.values.push_back(frequencies[i]);
        }
    }
    
    return result;
}

std::vector<uint32_t> FrequencyCompressor::decode_sparse(const SparseEncoding& encoded, size_t total_samples) {
    std::vector<uint32_t> decoded(total_samples, 0);
    
    for (size_t i = 0; i < encoded.indices.size(); ++i) {
        if (encoded.indices[i] < total_samples) {
            decoded[encoded.indices[i]] = encoded.values[i];
        }
    }
    
    return decoded;
}

FrequencyCompressor::EncodingType FrequencyCompressor::choose_best_encoding(
    const std::vector<uint32_t>& frequencies, size_t& estimated_size) {
    
    if (frequencies.empty()) {
        estimated_size = 0;
        return FULL_32BIT;
    }
    
    size_t sample_count = frequencies.size();
    uint32_t max_val = *std::max_element(frequencies.begin(), frequencies.end());
    size_t non_zero_count = std::count_if(frequencies.begin(), frequencies.end(),
                                         [](uint32_t v) { return v > 0; });
    
    size_t size_32bit = sample_count * 4 + 1;
    size_t size_16bit = sample_count * 2 + 1;
    size_t size_8bit = sample_count * 1 + 1;
    size_t size_sparse = non_zero_count * 2 + non_zero_count * 4 + 2 + 1;
    
    EncodingType best_type = FULL_32BIT;
    size_t min_size = size_32bit;

    if (max_val <= 65535 && size_16bit < min_size) {
        min_size = size_16bit;
        best_type = FULL_16BIT;
    }

    if (max_val <= 255 && size_8bit < min_size) {
        min_size = size_8bit;
        best_type = FULL_8BIT;
    }
    
    if (non_zero_count < sample_count * 0.3 && size_sparse < min_size) {
        min_size = size_sparse;
        best_type = SPARSE;
    }
    
    estimated_size = min_size;
    return best_type;
}

OutputManager::OutputManager(const ProgramConfig& cfg) 
    : config(cfg), current_block_id(1), total_rows_written(0),
      current_file_row_count(0), file_is_open(false),
      total_uncompressed_bytes(0), total_compressed_bytes(0) {}

OutputManager::~OutputManager() {
    if (file_is_open) {
        finalize_current_file();
        close_current_file();
    }
}

void OutputManager::write_data_block(const std::vector<KmerMatrixRow>& data_block, 
                     uint32_t kmer_size, uint32_t sample_count,
                     const std::vector<std::string>& sample_names) {
    std::lock_guard<std::mutex> lock(write_mutex);
    append_to_current_file(data_block, kmer_size, sample_count, sample_names);
    total_rows_written += data_block.size();
}

void OutputManager::open_new_output_file(uint32_t kmer_size, uint32_t sample_count,
                                        const std::vector<std::string>& sample_names) {
    if (file_is_open) {
        finalize_current_file();
        close_current_file();
    }
    
    char filename[1024];
    if (config.text_output_mode) {
        snprintf(filename, sizeof(filename), 
                "%s.%04" PRIu64 ".txt", config.output_prefix.c_str(), current_block_id);
    } else {
        snprintf(filename, sizeof(filename), 
                "%s.%04" PRIu64 ".bin", config.output_prefix.c_str(), current_block_id);
    }
    
    current_filename = filename;
    current_file_stream = std::make_unique<std::ofstream>();
    
    if (config.text_output_mode) {
        current_file_stream->open(current_filename);
        if (!current_file_stream->is_open()) {
            std::cerr << "File creation error: " << current_filename << "\n";
            exit(EXIT_FAILURE);
        }
        
        if (config.include_headers && current_block_id == 1) {
            *current_file_stream << "kmer";
            for (const auto& name : sample_names) {
                *current_file_stream << config.delimiter << name;
            }
            *current_file_stream << "\n";
        }
    } else {
        current_file_stream->open(current_filename, std::ios::binary);
        if (!current_file_stream->is_open()) {
            std::cerr << "File creation error: " << current_filename << "\n";
            exit(EXIT_FAILURE);
        }
        
        current_file_stream->write(reinterpret_cast<const char*>(&kmer_size), sizeof(uint32_t));
        current_file_stream->write(reinterpret_cast<const char*>(&sample_count), sizeof(uint32_t));
        
        row_count_position = current_file_stream->tellp();
        uint64_t placeholder = 0;
        current_file_stream->write(reinterpret_cast<const char*>(&placeholder), sizeof(uint64_t));
        
        uint8_t compression_flag = config.compress_output ? 1 : 0;
        current_file_stream->write(reinterpret_cast<const char*>(&compression_flag), sizeof(uint8_t));
    }
    
    file_is_open = true;
    current_file_row_count = 0;
    
    if (config.verbose_logging) {
        std::cout << "Opened new output file: " << current_filename 
                  << (config.compress_output && !config.text_output_mode ? " (compressed)" : "") << "\n";
    }
}

void OutputManager::append_to_current_file(const std::vector<KmerMatrixRow>& data_block,
                                          uint32_t kmer_size, uint32_t sample_count,
                                          const std::vector<std::string>& sample_names) {
    if (!file_is_open || current_file_row_count >= config.entries_per_block) {
        open_new_output_file(kmer_size, sample_count, sample_names);
    }
    
    size_t rows_to_write = data_block.size();
    size_t space_in_current_file = config.entries_per_block - current_file_row_count;
    
    if (rows_to_write > space_in_current_file) {
        std::vector<KmerMatrixRow> first_part(data_block.begin(), 
                                              data_block.begin() + space_in_current_file);
        std::vector<KmerMatrixRow> second_part(data_block.begin() + space_in_current_file,
                                               data_block.end());
        
        if (config.text_output_mode) {
            write_text_rows_to_stream(first_part);
        } else if (config.compress_output) {
            write_compressed_rows_to_stream(first_part, kmer_size, sample_count);
        } else {
            write_binary_rows_to_stream(first_part, kmer_size, sample_count);
        }
        current_file_row_count += first_part.size();
        
        append_to_current_file(second_part, kmer_size, sample_count, sample_names);
    } else {
        if (config.text_output_mode) {
            write_text_rows_to_stream(data_block);
        } else if (config.compress_output) {
            write_compressed_rows_to_stream(data_block, kmer_size, sample_count);
        } else {
            write_binary_rows_to_stream(data_block, kmer_size, sample_count);
        }
        current_file_row_count += rows_to_write;
    }
}

void OutputManager::write_text_rows_to_stream(const std::vector<KmerMatrixRow>& rows) {
    for (const auto& row : rows) {
        if (!row.meets_frequency_threshold(config.minimum_frequency_threshold)) {
            continue;
        }
        
        *current_file_stream << row.kmer_sequence;
        for (uint32_t count : row.frequency_counts) {
            *current_file_stream << config.delimiter << count;
        }
        *current_file_stream << "\n";
    }
}

void OutputManager::write_binary_rows_to_stream(const std::vector<KmerMatrixRow>& rows,
                                               uint32_t kmer_size, uint32_t sample_count) {
    uint32_t words_per_kmer = (kmer_size + 31) >> 5;
    
    for (const auto& row : rows) {
        if (!row.meets_frequency_threshold(config.minimum_frequency_threshold)) {
            continue;
        }
        
        current_file_stream->write(reinterpret_cast<const char*>(row.kmer_representation.data()), 
                                  sizeof(uint64_t) * words_per_kmer);
        current_file_stream->write(reinterpret_cast<const char*>(row.frequency_counts.data()), 
                                  sizeof(uint32_t) * sample_count);
        
        total_uncompressed_bytes += sizeof(uint64_t) * words_per_kmer + sizeof(uint32_t) * sample_count;
        total_compressed_bytes += sizeof(uint64_t) * words_per_kmer + sizeof(uint32_t) * sample_count;
    }
}

void OutputManager::write_compressed_rows_to_stream(const std::vector<KmerMatrixRow>& rows,
                                                   uint32_t kmer_size, uint32_t sample_count) {
    uint32_t words_per_kmer = (kmer_size + 31) >> 5;
    
    for (const auto& row : rows) {
        if (!row.meets_frequency_threshold(config.minimum_frequency_threshold)) {
            continue;
        }
        
        current_file_stream->write(reinterpret_cast<const char*>(row.kmer_representation.data()), 
                                  sizeof(uint64_t) * words_per_kmer);
        
        size_t estimated_size;
        auto encoding_type = FrequencyCompressor::choose_best_encoding(row.frequency_counts, estimated_size);
        
        current_file_stream->write(reinterpret_cast<const char*>(&encoding_type), sizeof(uint8_t));
        
        size_t uncompressed_size = sizeof(uint64_t) * words_per_kmer + sizeof(uint32_t) * sample_count;
        size_t compressed_size = sizeof(uint64_t) * words_per_kmer + 1;
        
        switch (encoding_type) {
            case FrequencyCompressor::FULL_8BIT: {
                auto encoded = FrequencyCompressor::encode_frequencies_8bit(row.frequency_counts);
                current_file_stream->write(reinterpret_cast<const char*>(encoded.data()), 
                                          encoded.size());
                compressed_size += encoded.size();
                break;
            }
            case FrequencyCompressor::FULL_16BIT: {
                auto encoded = FrequencyCompressor::encode_frequencies_16bit(row.frequency_counts);
                current_file_stream->write(reinterpret_cast<const char*>(encoded.data()), 
                                          encoded.size() * sizeof(uint16_t));
                compressed_size += encoded.size() * sizeof(uint16_t);
                break;
            }
            case FrequencyCompressor::SPARSE: {
                auto sparse = FrequencyCompressor::encode_sparse(row.frequency_counts);
                uint16_t count = static_cast<uint16_t>(sparse.indices.size());
                current_file_stream->write(reinterpret_cast<const char*>(&count), sizeof(uint16_t));
                current_file_stream->write(reinterpret_cast<const char*>(sparse.indices.data()), 
                                          sparse.indices.size() * sizeof(uint16_t));
                current_file_stream->write(reinterpret_cast<const char*>(sparse.values.data()), 
                                          sparse.values.size() * sizeof(uint32_t));
                compressed_size += sizeof(uint16_t) + sparse.indices.size() * sizeof(uint16_t) 
                                 + sparse.values.size() * sizeof(uint32_t);
                break;
            }
            case FrequencyCompressor::FULL_32BIT:
            default:
                current_file_stream->write(reinterpret_cast<const char*>(row.frequency_counts.data()), 
                                          sizeof(uint32_t) * sample_count);
                compressed_size += sizeof(uint32_t) * sample_count;
                break;
        }
        
        total_uncompressed_bytes += uncompressed_size;
        total_compressed_bytes += compressed_size;
    }
}

void OutputManager::finalize_current_file() {
    if (!file_is_open) return;
    
    if (!config.text_output_mode && current_file_stream->is_open()) {
        auto current_pos = current_file_stream->tellp();
        current_file_stream->seekp(row_count_position);
        current_file_stream->write(reinterpret_cast<const char*>(&current_file_row_count), 
                                  sizeof(uint64_t));
        current_file_stream->seekp(current_pos);
    }
    
    if (config.verbose_logging) {
        std::cout << "Finalized file " << current_filename << " with " 
                  << current_file_row_count << " k-mers\n";
    }
}

void OutputManager::close_current_file() {
    if (current_file_stream && current_file_stream->is_open()) {
        current_file_stream->close();
        current_file_stream.reset();
        file_is_open = false;
        ++current_block_id;
    }
}

void OutputManager::finalize() {
    std::lock_guard<std::mutex> lock(write_mutex);
    if (file_is_open) {
        finalize_current_file();
        close_current_file();
    }
}

uint64_t OutputManager::get_total_rows_written() const { 
    return total_rows_written; 
}

double OutputManager::get_compression_ratio() const {
    if (total_uncompressed_bytes == 0) return 0.0;
    double saved = (double)(total_uncompressed_bytes - total_compressed_bytes);
    return (saved / total_uncompressed_bytes) * 100.0;
}

ParallelKmerProcessor::ParallelKmerProcessor(std::vector<DatabaseReader>&& readers, const ProgramConfig& cfg)
    : database_readers(std::move(readers)), 
      processing_queue(cfg.max_queue_depth),
      file_manager(cfg), 
      config(cfg) {
    
    if (database_readers.empty()) {
        std::cerr << "Error: No available database files\n";
        exit(EXIT_FAILURE);
    }
    
    kmer_size = database_readers.front().get_kmer_length();
    sample_count = static_cast<uint32_t>(database_readers.size());
    words_per_kmer = (kmer_size + 31) >> 5;
    
    generate_sample_names();
    
    if (config.verbose_logging) {
        print_database_statistics();
        
        size_t kmer_row_size = sizeof(KmerMatrixRow) + 
                              words_per_kmer * sizeof(uint64_t) + 
                              sample_count * sizeof(uint32_t);
        size_t batch_size = config.batch_processing_size > 0 ? 
                          config.batch_processing_size : DEFAULT_BATCH_SIZE;
        size_t batch_memory = batch_size * kmer_row_size;
        size_t write_threshold = std::min(config.entries_per_block, WRITE_BUFFER_THRESHOLD);
        size_t max_memory = batch_memory * config.max_queue_depth + 
                          write_threshold * kmer_row_size;
        
        std::cout << "Memory estimation:\n";
        std::cout << "  K-mer row size: ~" << kmer_row_size << " bytes\n";
        std::cout << "  Batch size: " << batch_size << " k-mers\n";
        std::cout << "  Queue depth: " << config.max_queue_depth << " batches\n";
        std::cout << "  Estimated peak memory: ~" << (max_memory / 1024 / 1024) << " MB\n\n";
    }
}

void ParallelKmerProcessor::execute_processing() {
    auto start_time = std::chrono::high_resolution_clock::now();
    
    auto writer_future = std::async(std::launch::async, 
                                   &ParallelKmerProcessor::writer_thread, this);
    
    producer_thread();
    
    writer_future.wait();
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
    
    std::cout << "\n=== Processing Complete ===\n";
    std::cout << "Total processing time: " << duration.count() << " seconds\n";
    std::cout << "k-mers processed: " << processed_kmers.load() << "\n";
    std::cout << "Rows written: " << file_manager.get_total_rows_written() << "\n";
    if (duration.count() > 0) {
        std::cout << "Average processing speed: " << (processed_kmers.load() / duration.count()) << " k-mers/second\n";
    }
    
    if (config.compress_output && !config.text_output_mode) {
        double compression_ratio = file_manager.get_compression_ratio();
        if (compression_ratio > 0) {
            std::cout << "\n=== Compression Statistics ===\n";
            std::cout << "Compression ratio: " << std::fixed << std::setprecision(1) 
                      << compression_ratio << "%\n";
            std::cout << "Space saved: ~" << std::fixed << std::setprecision(1)
                      << compression_ratio << "% of original size\n";
            std::cout << "============================\n";
        }
    }
}

void ParallelKmerProcessor::generate_sample_names() {
    for (size_t i = 0; i < database_readers.size(); ++i) {
        std::string path = database_readers[i].get_database_path();
        size_t last_slash = path.find_last_of("/\\");
        if (last_slash != std::string::npos) {
            path = path.substr(last_slash + 1);
        }
        
        // Remove .kmc_ extension for KMC files
        size_t dot_pos = path.find(".kmc_");
        if (dot_pos != std::string::npos) {
            path = path.substr(0, dot_pos);
        }
        
        // Remove .bin extension for KMERIA files
        if (path.size() > 4 && path.substr(path.size() - 4) == ".bin") {
            path = path.substr(0, path.size() - 4);
        }
        
        sample_names.push_back(path);
    }
}

void ParallelKmerProcessor::print_database_statistics() {
    std::cout << "\n=== Database Statistics ===\n";
    uint64_t total_kmers = 0;
    for (size_t i = 0; i < database_readers.size(); ++i) {
        uint64_t db_kmers = database_readers[i].get_total_kmers();
        total_kmers += db_kmers;
        std::string format_str = (database_readers[i].get_format() == DatabaseFormat::KMC) ? "KMC" : "KMERIA";
        std::cout << "Sample " << (i + 1) << " (" << sample_names[i] << ", " << format_str << "): " 
                  << db_kmers << " k-mers\n";
    }
    std::cout << "Total k-mer count: " << total_kmers << "\n";
    std::cout << "k-mer length: " << kmer_size << "\n";
    std::cout << "===========================\n\n";
}

void ParallelKmerProcessor::producer_thread() {
    std::vector<KmerMatrixRow> current_batch;
    size_t batch_size = config.batch_processing_size > 0 ? 
                      config.batch_processing_size : DEFAULT_BATCH_SIZE;
    current_batch.reserve(batch_size);
    
    while (true) {
        size_t min_database_index = find_minimum_kmer_source();
        
        if (min_database_index == std::numeric_limits<size_t>::max()) {
            break;
        }
        
        auto minimum_kmer = database_readers[min_database_index].get_current_kmer();
        KmerMatrixRow matrix_row = build_matrix_row(minimum_kmer);
        
        current_batch.push_back(std::move(matrix_row));
        processed_kmers.fetch_add(1);
        
        if (current_batch.size() >= batch_size) {
            processing_queue.enqueue(std::move(current_batch));
            current_batch.clear();
            current_batch.reserve(batch_size);
        }
        
        if (config.verbose_logging && processed_kmers.load() % 100000 == 0) {
            std::cout << "Processed " << processed_kmers.load() << " k-mers, "
                      << "Queue size: " << processing_queue.size() << "\n";
        }
    }
    
    if (!current_batch.empty()) {
        processing_queue.enqueue(std::move(current_batch));
    }
    
    processing_queue.signal_completion();
}

void ParallelKmerProcessor::writer_thread() {
    std::vector<KmerMatrixRow> write_buffer;
    std::vector<KmerMatrixRow> batch_data;
    
    size_t buffer_threshold = config.write_buffer_threshold;
    write_buffer.reserve(buffer_threshold);
    
    while (processing_queue.dequeue(batch_data)) {
        write_buffer.insert(write_buffer.end(), 
                           std::make_move_iterator(batch_data.begin()),
                           std::make_move_iterator(batch_data.end()));
        
        batch_data.clear();
        batch_data.shrink_to_fit();
        
        if (write_buffer.size() >= buffer_threshold) {
            file_manager.write_data_block(write_buffer, kmer_size, sample_count, sample_names);
            write_buffer.clear();
            write_buffer.shrink_to_fit();
            write_buffer.reserve(buffer_threshold);
        }
    }
    
    if (!write_buffer.empty()) {
        file_manager.write_data_block(write_buffer, kmer_size, sample_count, sample_names);
        write_buffer.clear();
        write_buffer.shrink_to_fit();
    }
    
    file_manager.finalize();
}

size_t ParallelKmerProcessor::find_minimum_kmer_source() {
    size_t minimum_index = std::numeric_limits<size_t>::max();
    
    for (size_t i = 0; i < database_readers.size(); ++i) {
        if (!database_readers[i].is_exhausted()) {
            if (minimum_index == std::numeric_limits<size_t>::max() || 
                database_readers[i].get_current_kmer() < database_readers[minimum_index].get_current_kmer()) {
                minimum_index = i;
            }
        }
    }
    
    return minimum_index;
}

KmerMatrixRow ParallelKmerProcessor::build_matrix_row(const sequence_kmer& target_kmer) {
    KmerMatrixRow row(words_per_kmer, sample_count, config.text_output_mode);
    
    const uint64_t* kmer_data = target_kmer.extract_raw_data();
    std::copy(kmer_data, kmer_data + words_per_kmer, row.kmer_representation.begin());
    
    if (config.text_output_mode) {
        row.kmer_sequence = target_kmer.to_string_representation();
    }
    
    for (size_t i = 0; i < database_readers.size(); ++i) {
        if (!database_readers[i].is_exhausted() && 
            database_readers[i].get_current_kmer() == target_kmer) {
            row.frequency_counts[i] = database_readers[i].get_current_frequency();
            database_readers[i].advance_to_next();
        }
    }
    
    return row;
}

bool validate_kmer_consistency(const std::vector<DatabaseReader>& database_list) {
    if (database_list.empty()) {
        return true;
    }
    
    uint32_t reference_length = database_list.front().get_kmer_length();
    return std::all_of(database_list.begin(), database_list.end(),
                      [reference_length](const DatabaseReader& db) {
                          return db.get_kmer_length() == reference_length;
                      });
}

void display_help_message(const char* program_name) {
    std::cout << "K-mer Matrix Builder v2.1\n\n";
    std::cout << "=== MATRIX BUILDING MODE ===\n";
    std::cout << "Usage: " << program_name << " [options] -i <input_file> -o <output_prefix>\n\n";
    std::cout << "Supported formats:\n";
    std::cout << "  - KMC databases (sorted KMC3 format)\n";
    std::cout << "  - KMERIA binary files (kdump format)\n";
    std::cout << "  - Mixed: can combine both formats in one matrix\n\n";
    std::cout << "Required arguments:\n";
    std::cout << "  -i, --input FILE        File containing paths to databases (one per line)\n";
    std::cout << "  -o, --output PREFIX     Output file prefix\n\n";
    std::cout << "Output format options:\n";
    std::cout << "  -t, --text              Output text format (default: binary format)\n";
    std::cout << "  -d, --delimiter STR     Field delimiter for text format (default: tab)\n";
    std::cout << "  --no-header             Do not include header row in text format\n\n";
    std::cout << "Performance and memory options:\n";
    std::cout << "  -j, --threads N         Number of worker threads (default: auto-detect)\n";
    std::cout << "  -b, --batch-size N      Batch processing size (default: 1000)\n";
    std::cout << "  -s, --block-size N      K-mers per output file (default: 10000000)\n";
    std::cout << "  --buffer-size N         Write buffer threshold (default: 100000)\n";
    std::cout << "  --queue-depth N         Max queue depth for backpressure (default: 5)\n";
    std::cout << "  --memory-efficient      Enable memory-efficient mode (default: on)\n\n";
    std::cout << "Compression options (binary format only):\n";
    std::cout << "  --no-compression        Disable compression (default: enabled)\n";
    std::cout << "                          Note: Compression can reduce file size by 30-90%\n\n";
    std::cout << "Filtering options:\n";
    std::cout << "  -m, --min-freq N        Minimum frequency threshold (default: 1)\n\n";
    std::cout << "Other options:\n";
    std::cout << "  -v, --verbose           Verbose output mode with memory stats\n";
    std::cout << "  -h, --help              Show this help information\n\n";
    std::cout << "=== CONVERSION MODE ===\n";
    std::cout << "Convert binary k-mer matrix files to text format\n\n";
    std::cout << "Usage: " << program_name << " --convert [options]\n\n";
    std::cout << "Single file conversion:\n";
    std::cout << "  -i, --input FILE        Input binary file\n";
    std::cout << "  -o, --output FILE       Output text file\n\n";
    std::cout << "Batch conversion:\n";
    std::cout << "  --input-dir DIR         Input directory with binary files\n";
    std::cout << "  --output-dir DIR        Output directory for text files\n";
    std::cout << "  -t, --threads N         Parallel threads for batch mode (default: 1)\n\n";
    std::cout << "Examples:\n";
    std::cout << "  # Build matrix from mixed formats\n";
    std::cout << "  echo \"sample1.kmc_pre\" > databases.txt\n";
    std::cout << "  echo \"sample2.bin\" >> databases.txt\n";
    std::cout << "  " << program_name << " -i databases.txt -o output_matrix\n\n";
    std::cout << "  # Convert binary to text\n";
    std::cout << "  " << program_name << " --convert -i matrix.0001.bin -o matrix.0001.txt -v\n\n";
}

ProgramConfig parse_command_arguments(int argc, char** argv) {
    ProgramConfig config;
    
    static struct option long_options[] = {
        {"input", required_argument, 0, 'i'},
        {"output", required_argument, 0, 'o'},
        {"text", no_argument, 0, 't'},
        {"delimiter", required_argument, 0, 'd'},
        {"no-header", no_argument, 0, 1001},
        {"threads", required_argument, 0, 'j'},
        {"batch-size", required_argument, 0, 'b'},
        {"block-size", required_argument, 0, 's'},
        {"buffer-size", required_argument, 0, 1004},
        {"min-freq", required_argument, 0, 'm'},
        {"queue-depth", required_argument, 0, 1002},
        {"memory-efficient", no_argument, 0, 1003},
        {"no-compression", no_argument, 0, 1005},
        {"verbose", no_argument, 0, 'v'},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };
    
    int option_index = 0;
    int c;
    
    while ((c = getopt_long(argc, argv, "i:o:td:j:b:s:m:vh", long_options, &option_index)) != -1) {
        switch (c) {
            case 'i':
                config.input_file_list = optarg;
                break;
            case 'o':
                config.output_prefix = optarg;
                break;
            case 't':
                config.text_output_mode = true;
                break;
            case 'd':
                config.delimiter = optarg;
                break;
            case 1001:
                config.include_headers = false;
                break;
            case 'j':
                config.worker_thread_count = std::atoi(optarg);
                break;
            case 'b':
                config.batch_processing_size = std::atol(optarg);
                break;
            case 's':
                config.entries_per_block = std::atol(optarg);
                break;
            case 1004:
                config.write_buffer_threshold = std::atol(optarg);
                break;
            case 'm':
                config.minimum_frequency_threshold = std::atoi(optarg);
                break;
            case 1002:
                config.max_queue_depth = std::atol(optarg);
                break;
            case 1003:
                config.memory_efficient_mode = true;
                break;
            case 1005:
                config.compress_output = false;
                break;
            case 'v':
                config.verbose_logging = true;
                break;
            case 'h':
                display_help_message(argv[0]);
                exit(EXIT_SUCCESS);
            default:
                std::cerr << "Unknown option. Use -h for help.\n";
                exit(EXIT_FAILURE);
        }
    }
    
    if (config.input_file_list.empty() || config.output_prefix.empty()) {
        std::cerr << "Error: Must specify input file list (-i) and output prefix (-o)\n";
        std::cerr << "Use -h for complete usage instructions\n";
        exit(EXIT_FAILURE);
    }
    
    if (config.batch_processing_size == 0) {
        config.batch_processing_size = 1000;
    }
    if (config.entries_per_block == 0) {
        config.entries_per_block = 100000;
    }
    if (config.max_queue_depth == 0) {
        config.max_queue_depth = 5;
    }
    
    return config;
}

BinaryToTextConverter::BinaryToTextConverter(const std::string& input, const std::string& output,
                                           const std::string& delim, bool header, bool verb)
    : input_file(input), output_file(output), delimiter(delim), 
      include_header(header), verbose(verb) {}

BinaryToTextConverter::BinaryFileHeader BinaryToTextConverter::read_binary_header(std::ifstream& file) {
    BinaryFileHeader header;
    file.read(reinterpret_cast<char*>(&header.kmer_size), sizeof(uint32_t));
    file.read(reinterpret_cast<char*>(&header.sample_count), sizeof(uint32_t));
    file.read(reinterpret_cast<char*>(&header.total_kmers), sizeof(uint64_t));
    file.read(reinterpret_cast<char*>(&header.compression_enabled), sizeof(uint8_t)); // compression flag
    if (file.fail()) {
        throw std::runtime_error("Failed to read binary file header");
    }
    
    return header;
}


static constexpr char NUCLEOTIDE_TABLE[4] = {'A', 'C', 'G', 'T'};
std::string BinaryToTextConverter::decode_kmer(const std::vector<uint64_t>& binary_data, uint32_t kmer_length) {
    std::string result;
    result.reserve(kmer_length);
    
    for (int32_t i = kmer_length - 1; i >= 0; --i) {
        uint32_t word_idx = i >> 5;
        uint32_t bit_offset = (i & 31) << 1;
        uint32_t symbol = (binary_data[word_idx] >> bit_offset) & 3;
        result += NUCLEOTIDE_TABLE[symbol];
    }
    
    return result;
}


bool BinaryToTextConverter::convert() {
    auto start_time = std::chrono::high_resolution_clock::now();
    
    std::ifstream input(input_file, std::ios::binary);
    if (!input.is_open()) {
        std::cerr << "Error: Cannot open input file " << input_file << "\n";
        return false;
    }
    
    std::ofstream output(output_file);
    if (!output.is_open()) {
        std::cerr << "Error: Cannot create output file " << output_file << "\n";
        return false;
    }
    
    BinaryFileHeader header;
    try {
        header = read_binary_header(input);
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return false;
    }
    
    if (verbose) {
        std::cout << "Converting binary file: " << input_file << "\n";
        std::cout << "  K-mer length: " << header.kmer_size << "\n";
        std::cout << "  Sample count: " << header.sample_count << "\n";
        std::cout << "  Total k-mers: " << header.total_kmers << "\n";
        std::cout << "  Compression: " << (header.compression_enabled ? "Enabled" : "Disabled") << "\n";
    }
    
    if (include_header) {
        output << "kmer";
        for (uint32_t i = 0; i < header.sample_count; ++i) {
            output << delimiter << "sample" << (i + 1);
        }
        output << "\n";
    }
    
    uint32_t words_per_kmer = (header.kmer_size + 31) >> 5;
    std::vector<uint64_t> kmer_binary(words_per_kmer);
    std::vector<uint32_t> frequencies(header.sample_count);
    
    uint64_t converted_count = 0;
    while (converted_count < header.total_kmers) {
        
        input.read(reinterpret_cast<char*>(kmer_binary.data()), 
                  sizeof(uint64_t) * words_per_kmer);
        
        if (input.fail()) {
            std::cerr << "Warning: Unexpected end of file at k-mer " << converted_count << "\n";
            break;
        }
        
        
        if (header.compression_enabled) {
           
            uint8_t encoding_type;
            input.read(reinterpret_cast<char*>(&encoding_type), sizeof(uint8_t));
            
            if (input.fail()) {
                std::cerr << "Warning: Failed to read encoding type at k-mer " << converted_count << "\n";
                break;
            }
            
            switch (encoding_type) {
                case 0: // FULL_32BIT
                    input.read(reinterpret_cast<char*>(frequencies.data()), 
                              sizeof(uint32_t) * header.sample_count);
                    break;
                    
                case 1: { // FULL_16BIT
                    std::vector<uint16_t> freq_16(header.sample_count);
                    input.read(reinterpret_cast<char*>(freq_16.data()), 
                              sizeof(uint16_t) * header.sample_count);
                    for (uint32_t i = 0; i < header.sample_count; ++i) {
                        frequencies[i] = freq_16[i];
                    }
                    break;
                }
                
                case 2: { // FULL_8BIT
                    std::vector<uint8_t> freq_8(header.sample_count);
                    input.read(reinterpret_cast<char*>(freq_8.data()), 
                              header.sample_count);
                    for (uint32_t i = 0; i < header.sample_count; ++i) {
                        frequencies[i] = freq_8[i];
                    }
                    break;
                }
                
                case 3: { // SPARSE
                    
                    std::fill(frequencies.begin(), frequencies.end(), 0);
                    
                    
                    uint16_t non_zero_count;
                    input.read(reinterpret_cast<char*>(&non_zero_count), sizeof(uint16_t));
                    
                    
                    std::vector<uint16_t> indices(non_zero_count);
                    input.read(reinterpret_cast<char*>(indices.data()), 
                              sizeof(uint16_t) * non_zero_count);
                    
                    
                    std::vector<uint32_t> values(non_zero_count);
                    input.read(reinterpret_cast<char*>(values.data()), 
                              sizeof(uint32_t) * non_zero_count);
                    
                    for (uint16_t i = 0; i < non_zero_count; ++i) {
                        if (indices[i] < header.sample_count) {
                            frequencies[indices[i]] = values[i];
                        }
                    }
                    break;
                }
                
                default:
                    std::cerr << "Error: Unknown encoding type " << (int)encoding_type 
                              << " at k-mer " << converted_count << "\n";
                    return false;
            }
        } else {
            input.read(reinterpret_cast<char*>(frequencies.data()), 
                      sizeof(uint32_t) * header.sample_count);
        }
        
        if (input.fail()) {
            std::cerr << "Warning: Failed to read frequency data at k-mer " << converted_count << "\n";
            break;
        }
        
        std::string kmer_seq = decode_kmer(kmer_binary, header.kmer_size);
        output << kmer_seq;
        
        for (uint32_t freq : frequencies) {
            output << delimiter << freq;
        }
        output << "\n";
        
        ++converted_count;
        
        if (verbose && converted_count % 1000000 == 0) {
            std::cout << "  Converted " << converted_count << " / " << header.total_kmers << " k-mers...\n";
        }
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    double seconds = duration.count() / 1000.0;
    
    if (verbose || converted_count != header.total_kmers) {
        display_conversion_stats(converted_count, seconds);
    }
    
    return converted_count == header.total_kmers;
}



void BinaryToTextConverter::display_conversion_stats(uint64_t kmers_converted, double time_seconds) {
    std::cout << "\n=== Conversion Complete ===\n";
    std::cout << "K-mers converted: " << kmers_converted << "\n";
    std::cout << "Time: " << std::fixed << std::setprecision(2) << time_seconds << " seconds\n";
    if (time_seconds > 0) {
        std::cout << "Speed: " << static_cast<uint64_t>(kmers_converted / time_seconds) << " k-mers/second\n";
    }
    std::cout << "Output file: " << output_file << "\n";
    std::cout << "==========================\n";
}

BatchBinaryConverter::BatchBinaryConverter(const std::string& in_dir, const std::string& out_dir,
                                         const std::string& delim, bool header, bool verb, int threads)
    : input_dir(in_dir), output_dir(out_dir), delimiter(delim),
      include_header(header), verbose(verb), max_threads(threads) {}

std::vector<std::string> BatchBinaryConverter::find_binary_files() {
    std::vector<std::string> files;
    DIR* dir = opendir(input_dir.c_str());
    if (!dir) {
        std::cerr << "Error: Cannot open directory " << input_dir << "\n";
        return files;
    }
    
    struct dirent* entry;
    while ((entry = readdir(dir)) != nullptr) {
        std::string filename = entry->d_name;
        if (filename.find(".bin") != std::string::npos) {
            std::string full_path = input_dir + "/" + filename;
            
            struct stat file_stat;
            if (lstat(full_path.c_str(), &file_stat) == 0) {
                if (S_ISREG(file_stat.st_mode) || S_ISLNK(file_stat.st_mode)) {
                    files.push_back(full_path);
                }
            }
        }
    }
    closedir(dir);
    
    std::sort(files.begin(), files.end());
    return files;
}

void BatchBinaryConverter::convert_single_file(const std::string& input_file, const std::string& output_file) {
    BinaryToTextConverter converter(input_file, output_file, delimiter, include_header, verbose);
    converter.convert();
}

bool BatchBinaryConverter::convert_all() {
    auto files = find_binary_files();
    
    if (files.empty()) {
        std::cerr << "No binary files found in " << input_dir << "\n";
        return false;
    }
    
    std::cout << "Found " << files.size() << " binary files to convert\n";
    
    mkdir(output_dir.c_str(), 0755);
    
    std::vector<std::future<void>> futures;
    
    for (const auto& input_file : files) {
        std::string base_name = input_file.substr(input_file.find_last_of('/') + 1);
        size_t dot_pos = base_name.find_last_of('.');
        if (dot_pos != std::string::npos) {
            base_name = base_name.substr(0, dot_pos);
        }
        std::string output_file = output_dir + "/" + base_name + ".txt";
        
        if (futures.size() >= static_cast<size_t>(max_threads)) {
            futures[0].wait();
            futures.erase(futures.begin());
        }
        
        futures.emplace_back(std::async(std::launch::async,
                                       &BatchBinaryConverter::convert_single_file,
                                       this, input_file, output_file));
    }
    
    for (auto& future : futures) {
        future.wait();
    }
    
    std::cout << "All conversions completed!\n";
    return true;
}

int run_conversion_mode(int argc, char** argv) {
    std::string input_file;
    std::string output_file;
    std::string input_dir;
    std::string output_dir;
    std::string delimiter = "\t";
    bool no_header = false;
    bool verbose = false;
    bool batch_mode = false;
    int threads = 1;
    
    static struct option long_options[] = {
        {"input", required_argument, 0, 'i'},
        {"output", required_argument, 0, 'o'},
        {"input-dir", required_argument, 0, 1001},
        {"output-dir", required_argument, 0, 1002},
        {"delimiter", required_argument, 0, 'd'},
        {"no-header", no_argument, 0, 1003},
        {"threads", required_argument, 0, 't'},
        {"verbose", no_argument, 0, 'v'},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };
    
    int opt, option_index = 0;
    while ((opt = getopt_long(argc, argv, "i:o:d:t:vh", long_options, &option_index)) != -1) {
        switch (opt) {
            case 'i':
                input_file = optarg;
                break;
            case 'o':
                output_file = optarg;
                break;
            case 1001:
                input_dir = optarg;
                batch_mode = true;
                break;
            case 1002:
                output_dir = optarg;
                break;
            case 'd':
                delimiter = optarg;
                break;
            case 1003:
                no_header = true;
                break;
            case 't':
                threads = std::atoi(optarg);
                break;
            case 'v':
                verbose = true;
                break;
            case 'h':
                std::cout << "Binary to Text Conversion Mode\n\n";
                std::cout << "Usage: " << argv[0] << " --convert [options]\n\n";
                return EXIT_SUCCESS;
            default:
                return EXIT_FAILURE;
        }
    }
    
    if (batch_mode) {
        if (input_dir.empty() || output_dir.empty()) {
            std::cerr << "Error: Both input-dir and output-dir are required for batch mode\n";
            return EXIT_FAILURE;
        }
        
        BatchBinaryConverter converter(input_dir, output_dir, delimiter, !no_header, verbose, threads);
        return converter.convert_all() ? EXIT_SUCCESS : EXIT_FAILURE;
    } else {
        if (input_file.empty() || output_file.empty()) {
            std::cerr << "Error: Both input and output files are required\n";
            return EXIT_FAILURE;
        }
        
        BinaryToTextConverter converter(input_file, output_file, delimiter, !no_header, verbose);
        return converter.convert() ? EXIT_SUCCESS : EXIT_FAILURE;
    }
}


//int main(int argc, char* argv[]) {
extern "C" int kctm(int argc, char* argv[]) {
    if (argc > 1 && std::string(argv[1]) == "--convert") {
        for (int i = 1; i < argc - 1; ++i) {
            argv[i] = argv[i + 1];
        }
        argc--;
        return run_conversion_mode(argc, argv);
    }
    
    ProgramConfig config = parse_command_arguments(argc, argv);
    
    if (config.verbose_logging) {
        config.display_configuration();
    }
    
    std::ifstream database_list_file(config.input_file_list);
    if (!database_list_file) {
        std::cerr << "File opening error: " << config.input_file_list << "\n";
        return EXIT_FAILURE;
    }
    
    std::string database_path;
    std::vector<DatabaseReader> readers;
    
    std::cout << "Loading databases...\n";
    while (std::getline(database_list_file, database_path)) {
        if (database_path.empty() || database_path[0] == '#') {
            continue;
        }
        
        database_path.erase(database_path.find_last_not_of(" \t\r\n") + 1);
        
        if (!database_path.empty()) {
            if (config.verbose_logging) {
                DatabaseFormat fmt = detect_database_format(database_path);
                std::string fmt_str = (fmt == DatabaseFormat::KMC) ? "KMC" : "KMERIA";
                std::cout << "  Loading " << fmt_str << ": " << database_path << "\n";
            }
            readers.emplace_back(database_path);
        }
    }
    
    if (readers.empty()) {
        std::cout << "No valid database files found\n";
        return EXIT_SUCCESS;
    }
    
    if (!validate_kmer_consistency(readers)) {
        std::cerr << "Error: All databases must have the same k-mer length\n";
        return EXIT_FAILURE;
    }
    
    std::cout << "Successfully loaded " << readers.size() << " database files\n";
    std::cout << "k-mer length: " << readers.front().get_kmer_length() << "\n";
    
    if (config.worker_thread_count == 0) {
        config.worker_thread_count = std::thread::hardware_concurrency();
        if (config.worker_thread_count == 0) {
            config.worker_thread_count = 2;
        }
    }
    
    if (config.verbose_logging) {
        std::cout << "Using " << config.worker_thread_count << " worker threads\n";
    }
    
    std::cout << "\nStarting k-mer matrix construction...\n";
    
    ParallelKmerProcessor processor(std::move(readers), config);
    processor.execute_processing();
    
    return EXIT_SUCCESS;
}
