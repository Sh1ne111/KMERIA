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

#include "kmtob.h"


BgzfOutputStream::BgzfOutputStream(const std::string& fname, int compression_level, int threads) 
    : bgzf_file(nullptr), filename(fname), is_open(false), buffer_size(DEFAULT_BUFFER_SIZE) {
    open(fname, compression_level, threads);
}

BgzfOutputStream::~BgzfOutputStream() {
    close();
}

bool BgzfOutputStream::open(const std::string& fname, int compression_level, int threads) {
    if (is_open) {
        close();
    }
    
    filename = fname;
    
    if (compression_level < 0 || compression_level > 9) {
        compression_level = 6;
    }
    
    if (threads < 1) threads = 1;
    if (threads > 256) threads = 256;
    
    std::string mode = "w" + std::to_string(compression_level);
    
    bgzf_file = bgzf_open(filename.c_str(), mode.c_str());
    
    if (bgzf_file == nullptr) {
        std::cerr << "Error: Cannot create BGZF file " << filename << "\n";
        is_open = false;
        return false;
    }
    
    if (threads > 1) {
        if (bgzf_mt(bgzf_file, threads, 256) < 0) {
            std::cerr << "Warning: Failed to set BGZF threads to " << threads 
                      << ", using single thread\n";
        }
    }
    
    // Pre-alloc
    buffer.reserve(buffer_size);
    is_open = true;
    
    return true;
}

void BgzfOutputStream::flush_buffer() {
    if (!buffer.empty() && is_valid()) {
        ssize_t written = bgzf_write(bgzf_file, buffer.c_str(), buffer.length());
        if (written < 0 || static_cast<size_t>(written) != buffer.length()) {
            std::cerr << "Error: BGZF write failed for " << filename << "\n";
        }
        buffer.clear();
    }
}

void BgzfOutputStream::close() {
    if (bgzf_file != nullptr) {
        flush_buffer();
        
        // flush & close BGZF
        if (bgzf_close(bgzf_file) < 0) {
            std::cerr << "Error: Failed to close BGZF file " << filename << "\n";
        }
        bgzf_file = nullptr;
    }
    is_open = false;
}

bool BgzfOutputStream::write_string(const std::string& str) {
    if (!is_valid()) {
        return false;
    }
    
    buffer += str;
    
    if (buffer.size() >= buffer_size) {
        flush_buffer();
    }
    
    return true;
}

BgzfOutputStream& BgzfOutputStream::operator<<(const std::string& str) {
    write_string(str);
    return *this;
}

BgzfOutputStream& BgzfOutputStream::operator<<(const char* str) {
    write_string(std::string(str));
    return *this;
}

BgzfOutputStream& BgzfOutputStream::operator<<(char c) {
    buffer += c;
    if (buffer.size() >= buffer_size) {
        flush_buffer();
    }
    return *this;
}

BgzfOutputStream& BgzfOutputStream::operator<<(int val) {
    write_string(std::to_string(val));
    return *this;
}

BgzfOutputStream& BgzfOutputStream::operator<<(long val) {
    write_string(std::to_string(val));
    return *this;
}

BgzfOutputStream& BgzfOutputStream::operator<<(long long val) {
    write_string(std::to_string(val));
    return *this;
}

BgzfOutputStream& BgzfOutputStream::operator<<(unsigned int val) {
    write_string(std::to_string(val));
    return *this;
}

BgzfOutputStream& BgzfOutputStream::operator<<(unsigned long val) {
    write_string(std::to_string(val));
    return *this;
}

BgzfOutputStream& BgzfOutputStream::operator<<(unsigned long long val) {
    write_string(std::to_string(val));
    return *this;
}

BgzfOutputStream& BgzfOutputStream::operator<<(float val) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(6) << val;
    write_string(oss.str());
    return *this;
}

BgzfOutputStream& BgzfOutputStream::operator<<(double val) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(6) << val;
    write_string(oss.str());
    return *this;
}

BgzfOutputStream& BgzfOutputStream::operator<<(std::ostream& (*manip)(std::ostream&)) {
    if (manip == static_cast<std::ostream& (*)(std::ostream&)>(std::endl)) {
        write_string("\n");
    }
    return *this;
}

// ============================================================================
// BimbamConfig implementation
// ============================================================================

void BimbamConfig::display_config() const {
    std::cout << "=== BIMBAM Converter Configuration (BGZF) ===\n";
    std::cout << "Input directory: " << input_dir << "\n";
    std::cout << "Output directory: " << output_dir << "\n";
    std::cout << "Input prefix filter: " << input_prefix << "\n";
    std::cout << "File processing threads: " << max_threads << "\n";
    std::cout << "Input format: " << input_format << "\n";
    if (use_quantile_normalization) {
        std::cout << "Normalization: Quantile-based (Enabled)\n";
        std::cout << "Lower quantile threshold: " << (lower_quantile * 100) << "%\n";
        std::cout << "Upper quantile threshold: " << (upper_quantile * 100) << "%\n";
        std::cout << "Output range: [" << min_range << ", " << max_range << "]\n";
    } else if (normalize_values) {
        std::cout << "Normalization: Linear (Enabled)\n";
        std::cout << "Normalization range: [" << min_range << ", " << max_range << "]\n";
    } else {
        std::cout << "Normalization: Disabled\n";
    }
    
    std::cout << "Allele information: " << allele_info << "\n"; 

/*
   std::cout << "Normalize values: " << (normalize_values ? "Enabled" : "Disabled") << "\n";
    if (normalize_values) {
        std::cout << "Normalization range: [" << min_range << ", " << max_range << "]\n";
    }
    std::cout << "Allele information: " << allele_info << "\n";
*/    
    if (compress_output) {
        std::cout << "Compression: BGZF (Level " << compression_level << ")\n";
        std::cout << "BGZF threads: " << bgzf_threads << "\n";
        std::cout << "Write buffer: " << (write_buffer_size / 1024) << " KB\n";
    } else {
        std::cout << "Compression: Disabled\n";
    }
    
    std::cout << "Verbose logging: " << (verbose ? "Enabled" : "Disabled") << "\n";
    std::cout << "============================================\n\n";
}

bool BimbamConfig::validate() const {
    if (input_dir.empty()) {
        std::cerr << "Error: Input directory not specified\n";
        return false;
    }
    if (output_dir.empty()) {
        std::cerr << "Error: Output directory not specified\n";
        return false;
    }
    if (max_threads <= 0) {
        std::cerr << "Error: Number of threads must be positive\n";
        return false;
    }
    if (normalize_values && min_range >= max_range) {
        std::cerr << "Error: Invalid normalization range\n";
        return false;
    }
 
    if (use_quantile_normalization) {
        if (lower_quantile < 0.0 || lower_quantile >= 1.0) {
            std::cerr << "Error: Lower quantile must be in [0, 1)\n";
            return false;
        }
        if (upper_quantile <= 0.0 || upper_quantile > 1.0) {
            std::cerr << "Error: Upper quantile must be in (0, 1]\n";
            return false;
        }
        if (lower_quantile >= upper_quantile) {
            std::cerr << "Error: Lower quantile must be less than upper quantile\n";
            return false;
        }
    }

   if (compression_level < 0 || compression_level > 9) {
        std::cerr << "Error: Compression level must be between 0 and 9\n";
        return false;
    }
    if (bgzf_threads < 1 || bgzf_threads > 256) {
        std::cerr << "Error: BGZF threads must be between 1 and 256\n";
        return false;
    }
    if (write_buffer_size < 1024 || write_buffer_size > 10485760) {
        std::cerr << "Error: Write buffer size must be between 1KB and 10MB\n";
        return false;
    }
    return true;
}


static constexpr char NUCLEOTIDE_LOOKUP[4] = {'A', 'C', 'G', 'T'};

std::string BimbamKmerData::binary_to_string() const {
    std::string result;
    uint32_t k_len = this->kmer_length;
    result.reserve(k_len);
    
    for (int32_t i = k_len - 1; i >= 0; --i) {
        uint32_t word_idx = i >> 5; 
        uint32_t bit_offset = (i & 31) << 1;
        uint32_t symbol = (kmer_binary[word_idx] >> bit_offset) & 0x03;
        result += NUCLEOTIDE_LOOKUP[symbol];
    }
    
    return result;
}

void BimbamKmerData::normalize_frequencies(double min_val, double max_val) {
    if (frequencies.empty()) return;
    
    auto minmax = std::minmax_element(frequencies.begin(), frequencies.end());
    double data_min = *minmax.first;
    double data_max = *minmax.second;
    double data_range = data_max - data_min;
    
    if (data_range == 0.0) {
        std::fill(frequencies.begin(), frequencies.end(), max_val);
        return;
    }
    
    double target_range = max_val - min_val;
    for (double& freq : frequencies) {
        freq = min_val + (freq - data_min) / data_range * target_range;
    }
}

void BimbamKmerData::quantile_normalize_frequencies(double lower_quantile, 
                                                   double upper_quantile,
                                                   double min_val, 
                                                   double max_val) {
    if (frequencies.empty()) return;

    std::vector<double> sorted_freqs = frequencies;
    std::sort(sorted_freqs.begin(), sorted_freqs.end());
    double q_lower = BimbamUtils::calculate_quantile(sorted_freqs, lower_quantile);
    double q_upper = BimbamUtils::calculate_quantile(sorted_freqs, upper_quantile);

    if (q_lower == q_upper) {
        std::fill(frequencies.begin(), frequencies.end(), max_val);
        return;
    }
    
    double quantile_range = q_upper - q_lower;
    double target_range = max_val - min_val;
    for (double& freq : frequencies) {
        if (freq <= q_lower) {
            freq = min_val;
        } else if (freq >= q_upper) {
            freq = max_val;
        } else {
            
            freq = min_val + (freq - q_lower) / quantile_range * target_range;
        }
    }
}


BimbamKmerData::KmerStats BimbamKmerData::get_statistics() const {
    KmerStats stats;
    
    if (frequencies.empty()) {
        stats = {0.0, 0.0, 0.0, 0.0, 0, 0};
        return stats;
    }
    
    auto minmax = std::minmax_element(frequencies.begin(), frequencies.end());
    stats.min_freq = *minmax.first;
    stats.max_freq = *minmax.second;
    
    double sum = std::accumulate(frequencies.begin(), frequencies.end(), 0.0);
    stats.mean_freq = sum / frequencies.size();
    
    double sq_sum = std::accumulate(frequencies.begin(), frequencies.end(), 0.0,
                                   [&stats](double acc, double freq) {
                                       double diff = freq - stats.mean_freq;
                                       return acc + diff * diff;
                                   });
    stats.std_dev = std::sqrt(sq_sum / frequencies.size());
    
    stats.zero_count = std::count_if(frequencies.begin(), frequencies.end(),
                                    [](double f) { return f == 0.0; });
    stats.non_zero_count = frequencies.size() - stats.zero_count;
    
    return stats;
}

// ============================================================================
// FrequencyDecompressor implementation
// ============================================================================

std::vector<uint32_t> FrequencyDecompressor::decode_8bit(const std::vector<uint8_t>& encoded) {
    std::vector<uint32_t> decoded;
    decoded.reserve(encoded.size());
    for (uint8_t val : encoded) {
        decoded.push_back(static_cast<uint32_t>(val));
    }
    return decoded;
}

std::vector<uint16_t> FrequencyDecompressor::decode_16bit(const std::vector<uint16_t>& encoded) {
    return encoded;
}

std::vector<uint32_t> FrequencyDecompressor::decode_sparse(const SparseEncoding& encoded, 
                                                          size_t total_samples) {
    std::vector<uint32_t> decoded(total_samples, 0);
    
    for (size_t i = 0; i < encoded.indices.size(); ++i) {
        if (encoded.indices[i] < total_samples) {
            decoded[encoded.indices[i]] = encoded.values[i];
        }
    }
    
    return decoded;
}

// ============================================================================
// CompressedBinaryKmerReader implementation
// ============================================================================

CompressedBinaryKmerReader::CompressedBinaryKmerReader(const std::string& filename) 
    : current_position(0), is_valid(false), is_compressed(false) {
    
    file_stream.open(filename, std::ios::binary);
    if (!file_stream.is_open()) {
        std::cerr << "Error: Cannot open binary file " << filename << "\n";
        return;
    }
    
    file_stream.read(reinterpret_cast<char*>(&kmer_size), sizeof(uint32_t));
    file_stream.read(reinterpret_cast<char*>(&sample_count), sizeof(uint32_t));
    file_stream.read(reinterpret_cast<char*>(&total_kmers), sizeof(uint64_t));
    
    if (file_stream.fail()) {
        std::cerr << "Error: Failed to read header from " << filename << "\n";
        return;
    }
    
    uint8_t compression_flag = 0;
    file_stream.read(reinterpret_cast<char*>(&compression_flag), sizeof(uint8_t));
    
    if (file_stream.fail()) {
        file_stream.clear();
        file_stream.seekg(sizeof(uint32_t) * 2 + sizeof(uint64_t), std::ios::beg);
        is_compressed = false;
    } else {
        is_compressed = (compression_flag != 0);
    }
    
    words_per_kmer = (kmer_size + 31) >> 5;
    is_valid = true;
}

CompressedBinaryKmerReader::~CompressedBinaryKmerReader() {
    if (file_stream.is_open()) {
        file_stream.close();
    }
}

bool CompressedBinaryKmerReader::read_uncompressed_kmer(BimbamKmerData& kmer_data) {
    file_stream.read(reinterpret_cast<char*>(kmer_data.kmer_binary.data()),
                    sizeof(uint64_t) * words_per_kmer);
    
    std::vector<uint32_t> temp_frequencies(sample_count);
    file_stream.read(reinterpret_cast<char*>(temp_frequencies.data()),
                    sizeof(uint32_t) * sample_count);
    
    if (file_stream.fail()) {
        return false;
    }
    
    for (size_t i = 0; i < sample_count; ++i) {
        kmer_data.frequencies[i] = static_cast<double>(temp_frequencies[i]);
    }
    
    return true;
}

bool CompressedBinaryKmerReader::read_compressed_kmer(BimbamKmerData& kmer_data) {
    file_stream.read(reinterpret_cast<char*>(kmer_data.kmer_binary.data()),
                    sizeof(uint64_t) * words_per_kmer);
    
    if (file_stream.fail()) {
        return false;
    }
    
    uint8_t encoding_type;
    file_stream.read(reinterpret_cast<char*>(&encoding_type), sizeof(uint8_t));
    
    if (file_stream.fail()) {
        return false;
    }
    
    std::vector<uint32_t> temp_frequencies;
    
    switch (static_cast<FrequencyDecompressor::EncodingType>(encoding_type)) {
        case FrequencyDecompressor::FULL_8BIT: {
            std::vector<uint8_t> encoded(sample_count);
            file_stream.read(reinterpret_cast<char*>(encoded.data()), sample_count);
            if (file_stream.fail()) return false;
            temp_frequencies = FrequencyDecompressor::decode_8bit(encoded);
            break;
        }
        
        case FrequencyDecompressor::FULL_16BIT: {
            std::vector<uint16_t> encoded(sample_count);
            file_stream.read(reinterpret_cast<char*>(encoded.data()), 
                           sample_count * sizeof(uint16_t));
            if (file_stream.fail()) return false;
            temp_frequencies.reserve(sample_count);
            for (uint16_t val : encoded) {
                temp_frequencies.push_back(static_cast<uint32_t>(val));
            }
            break;
        }
        
        case FrequencyDecompressor::SPARSE: {
            uint16_t count;
            file_stream.read(reinterpret_cast<char*>(&count), sizeof(uint16_t));
            if (file_stream.fail()) return false;
            
            FrequencyDecompressor::SparseEncoding sparse;
            sparse.indices.resize(count);
            sparse.values.resize(count);
            
            file_stream.read(reinterpret_cast<char*>(sparse.indices.data()), 
                           count * sizeof(uint16_t));
            file_stream.read(reinterpret_cast<char*>(sparse.values.data()), 
                           count * sizeof(uint32_t));
            
            if (file_stream.fail()) return false;
            temp_frequencies = FrequencyDecompressor::decode_sparse(sparse, sample_count);
            break;
        }
        
        case FrequencyDecompressor::FULL_32BIT:
        default:
            temp_frequencies.resize(sample_count);
            file_stream.read(reinterpret_cast<char*>(temp_frequencies.data()),
                           sizeof(uint32_t) * sample_count);
            if (file_stream.fail()) return false;
            break;
    }
    
    kmer_data.frequencies.resize(temp_frequencies.size());
    for (size_t i = 0; i < temp_frequencies.size(); ++i) {
        kmer_data.frequencies[i] = static_cast<double>(temp_frequencies[i]);
    }
    
    return true;
}

bool CompressedBinaryKmerReader::read_next_kmer(BimbamKmerData& kmer_data) {
    if (!is_valid || !has_more_data()) {
        return false;
    }
    
    kmer_data = BimbamKmerData(kmer_size, sample_count);
    
    bool success;
    if (is_compressed) {
        success = read_compressed_kmer(kmer_data);
    } else {
        success = read_uncompressed_kmer(kmer_data);
    }
    
    if (!success) {
        std::cerr << "Error: Failed to read k-mer data at position " << current_position << "\n";
        return false;
    }
    
    kmer_data.kmer_sequence = kmer_data.binary_to_string();
    
    ++current_position;
    return true;
}

void CompressedBinaryKmerReader::reset_to_beginning() {
    if (is_valid) {
        file_stream.clear();
        size_t header_size = sizeof(uint32_t) * 2 + sizeof(uint64_t);
        if (is_compressed) {
            header_size += sizeof(uint8_t);
        }
        file_stream.seekg(header_size, std::ios::beg);
        current_position = 0;
    }
}

// ============================================================================
// BinaryKmerReader implementation
// ============================================================================
/*
BinaryKmerReader::BinaryKmerReader(const std::string& filename) 
    : current_position(0), is_valid(false) {
    
    file_stream.open(filename, std::ios::binary);
    if (!file_stream.is_open()) {
        std::cerr << "Error: Cannot open binary file " << filename << "\n";
        return;
    }
    
    file_stream.read(reinterpret_cast<char*>(&kmer_size), sizeof(uint32_t));
    file_stream.read(reinterpret_cast<char*>(&sample_count), sizeof(uint32_t));
    file_stream.read(reinterpret_cast<char*>(&total_kmers), sizeof(uint64_t));
    
    if (file_stream.fail()) {
        std::cerr << "Error: Failed to read header from " << filename << "\n";
        return;
    }
    
    words_per_kmer = (kmer_size + 31) >> 5;
    is_valid = true;
}

BinaryKmerReader::~BinaryKmerReader() {
    if (file_stream.is_open()) {
        file_stream.close();
    }
}
*/


bool BinaryKmerReader::read_next_kmer(BimbamKmerData& kmer_data) {
    if (!is_valid || !has_more_data()) {
        return false;
    }
    
    kmer_data = BimbamKmerData(kmer_size, sample_count);
    
    file_stream.read(reinterpret_cast<char*>(kmer_data.kmer_binary.data()),
                    sizeof(uint64_t) * words_per_kmer);
    
    std::vector<uint32_t> temp_frequencies(sample_count);
    file_stream.read(reinterpret_cast<char*>(temp_frequencies.data()),
                    sizeof(uint32_t) * sample_count);
    
    if (file_stream.fail()) {
        std::cerr << "Error: Failed to read k-mer data at position " << current_position << "\n";
        return false;
    }
    
    for (size_t i = 0; i < sample_count; ++i) {
        kmer_data.frequencies[i] = static_cast<double>(temp_frequencies[i]);
    }
    
    kmer_data.kmer_sequence = kmer_data.binary_to_string();
    
    ++current_position;
    return true;
}

void BinaryKmerReader::reset_to_beginning() {
    if (is_valid) {
        file_stream.clear();
        file_stream.seekg(sizeof(uint32_t) * 2 + sizeof(uint64_t), std::ios::beg);
        current_position = 0;
    }
}

// ============================================================================
// BimbamUtils implementation
// ============================================================================

std::vector<double> BimbamUtils::normalize_to_range(const std::vector<double>& numbers,
                                                   double min_range, double max_range) {
    if (numbers.empty()) {
        return {};
    }
    
    auto minmax = std::minmax_element(numbers.begin(), numbers.end());
    double data_min = *minmax.first;
    double data_max = *minmax.second;
    double data_range = data_max - data_min;
    
    if (data_range == 0.0) {
        return std::vector<double>(numbers.size(), max_range);
    }
    
    std::vector<double> normalized;
    normalized.reserve(numbers.size());
    double target_range = max_range - min_range;
    
    for (double number : numbers) {
        double normalized_val = min_range + (number - data_min) / data_range * target_range;
        normalized.push_back(normalized_val);
    }
    
    return normalized;
}

double BimbamUtils::calculate_mean(const std::vector<double>& numbers) {
    if (numbers.empty()) return 0.0;
    double sum = std::accumulate(numbers.begin(), numbers.end(), 0.0);
    return sum / numbers.size();
}

double BimbamUtils::calculate_std_dev(const std::vector<double>& numbers, double mean) {
    if (numbers.empty()) return 0.0;
    double sq_sum = std::accumulate(numbers.begin(), numbers.end(), 0.0,
                                   [mean](double acc, double val) {
                                       double diff = val - mean;
                                       return acc + diff * diff;
                                   });
    return std::sqrt(sq_sum / numbers.size());
}

double BimbamUtils::calculate_quantile(const std::vector<double>& sorted_data, double quantile) {
    if (sorted_data.empty()) return 0.0;
    if (quantile <= 0.0) return sorted_data.front();
    if (quantile >= 1.0) return sorted_data.back();
    
    // 使用线性插值计算分位数
    double index = quantile * (sorted_data.size() - 1);
    size_t lower_index = static_cast<size_t>(std::floor(index));
    size_t upper_index = static_cast<size_t>(std::ceil(index));
    
    if (lower_index == upper_index) {
        return sorted_data[lower_index];
    }
    
    double weight = index - lower_index;
    return sorted_data[lower_index] * (1.0 - weight) + sorted_data[upper_index] * weight;
}

std::vector<double> BimbamUtils::quantile_normalize(const std::vector<double>& numbers,
                                                    double lower_quantile,
                                                    double upper_quantile,
                                                    double min_range,
                                                    double max_range) {
    if (numbers.empty()) {
        return {};
    }
    
    std::vector<double> sorted_data = numbers;
    std::sort(sorted_data.begin(), sorted_data.end());
    
    double q_lower = calculate_quantile(sorted_data, lower_quantile);
    double q_upper = calculate_quantile(sorted_data, upper_quantile);
    
    if (q_lower == q_upper) {
        return std::vector<double>(numbers.size(), max_range);
    }
    
    std::vector<double> normalized;
    normalized.reserve(numbers.size());
    
    double quantile_range = q_upper - q_lower;
    double target_range = max_range - min_range;
    
    for (double number : numbers) {
        double normalized_val;
        if (number <= q_lower) {
            normalized_val = min_range; 
        } else if (number >= q_upper) {
            normalized_val = max_range;
        } else {
            
            normalized_val = min_range + (number - q_lower) / quantile_range * target_range;
        }
        normalized.push_back(normalized_val);
    }
    
    return normalized;
}




std::pair<double, double> BimbamUtils::get_min_max(const std::vector<double>& numbers) {
    if (numbers.empty()) return {0.0, 0.0};
    auto minmax = std::minmax_element(numbers.begin(), numbers.end());
    return {*minmax.first, *minmax.second};
}


// TextBimbamProcessor implementation (BGZF)
bool TextBimbamProcessor::process_text_file(const std::string& input_file,
                                           const std::string& output_file,
                                           const BimbamConfig& config,
                                           uint64_t& processed_count) {
    std::ifstream input_stream(input_file);
    if (!input_stream) {
        std::cerr << "Error: Cannot open input file " << input_file << "\n";
        return false;
    }
    
    processed_count = 0;
    
    if (config.compress_output) {
        BgzfOutputStream bgzf_output(output_file, config.compression_level, config.bgzf_threads);
        if (!bgzf_output.is_valid()) {
            std::cerr << "Error: Cannot create BGZF output file " << output_file << "\n";
            return false;
        }
        
        bgzf_output.set_buffer_size(config.write_buffer_size);
        
        std::string line;
        std::vector<double> frequencies;
        frequencies.reserve(1000);
        
        while (std::getline(input_stream, line)) {
            if (line.empty() || line[0] == '#') {
                continue;
            }
            
            std::stringstream ss(line);
            std::string kmer;
            std::getline(ss, kmer, '\t');
            
            frequencies.clear();
            double frequency;
            while (ss >> frequency) {
                frequencies.push_back(frequency);
                if (ss.peek() == '\t') ss.ignore();
            }
            
            if (frequencies.empty()) {
                continue;
            }

            if (config.use_quantile_normalization) {
                frequencies = BimbamUtils::quantile_normalize(frequencies,
                                                             config.lower_quantile,
                                                             config.upper_quantile,
                                                             config.min_range, 
                                                             config.max_range);            

            } else if (config.normalize_values) {
                frequencies = BimbamUtils::normalize_to_range(frequencies, 
                                                             config.min_range, 
                                                             config.max_range);
            }
            
            std::ostringstream line_out;
            line_out << kmer << ", " << config.allele_info << ", ";
            for (size_t i = 0; i < frequencies.size(); ++i) {
                line_out << std::fixed << std::setprecision(6) << frequencies[i];
                if (i < frequencies.size() - 1) {
                    line_out << ", ";
                }
            }
            line_out << "\n";
            
            bgzf_output << line_out.str();
            ++processed_count;
        }
        
    } else {
        std::ofstream output_stream(output_file);
        if (!output_stream) {
            std::cerr << "Error: Cannot create output file " << output_file << "\n";
            return false;
        }
        
        std::string line;
        while (std::getline(input_stream, line)) {
            if (line.empty() || line[0] == '#') {
                continue;
            }
            
            std::stringstream ss(line);
            std::string kmer;
            std::getline(ss, kmer, '\t');
            
            std::vector<double> frequencies;
            double frequency;
            while (ss >> frequency) {
                frequencies.push_back(frequency);
                if (ss.peek() == '\t') ss.ignore();
            }
            
            if (frequencies.empty()) {
                continue;
            }

            if (config.use_quantile_normalization) {
                frequencies = BimbamUtils::quantile_normalize(frequencies,
                                                             config.lower_quantile,
                                                             config.upper_quantile,
                                                             config.min_range, 
                                                             config.max_range);
            } else if (config.normalize_values) {
                frequencies = BimbamUtils::normalize_to_range(frequencies, 
                                                             config.min_range, 
                                                             config.max_range);
            }
            
            output_stream << kmer << ", " << config.allele_info << ", ";
            for (size_t i = 0; i < frequencies.size(); ++i) {
                output_stream << std::fixed << std::setprecision(6) << frequencies[i];
                if (i < frequencies.size() - 1) {
                    output_stream << ", ";
                }
            }
            output_stream << "\n";
            
            ++processed_count;
        }
    }
    
    return true;
}

// BinaryBimbamProcessor implementation

bool BinaryBimbamProcessor::process_binary_file(const std::string& input_file,
                                               const std::string& output_file,
                                               const BimbamConfig& config,
                                               uint64_t& processed_count) {
    CompressedBinaryKmerReader reader(input_file);
    if (!reader.is_open()) {
        std::cerr << "Error: Failed to open binary file " << input_file << "\n";
        return false;
    }
    
    if (config.verbose && reader.compressed()) {
        std::cout << "Reading compressed binary file: " << input_file << "\n";
    }
    
    BimbamKmerData kmer_data;
    processed_count = 0;
    
    if (config.compress_output) {
        BgzfOutputStream bgzf_output(output_file, config.compression_level, config.bgzf_threads);
        if (!bgzf_output.is_valid()) {
            std::cerr << "Error: Cannot create BGZF output file " << output_file << "\n";
            return false;
        }
        
        bgzf_output.set_buffer_size(config.write_buffer_size);
        
        while (reader.read_next_kmer(kmer_data)) {

            if (config.use_quantile_normalization) {
                kmer_data.quantile_normalize_frequencies(config.lower_quantile,
                                                        config.upper_quantile,
                                                        config.min_range,
                                                        config.max_range);
            } else if (config.normalize_values) {
                kmer_data.frequencies = BimbamUtils::normalize_to_range(
                    kmer_data.frequencies, 
                    config.min_range, 
                    config.max_range);
            }
            
            std::ostringstream line_out;
            line_out << kmer_data.kmer_sequence << ", " << config.allele_info << ", ";
            for (size_t i = 0; i < kmer_data.frequencies.size(); ++i) {
                line_out << std::fixed << std::setprecision(6) << kmer_data.frequencies[i];
                if (i < kmer_data.frequencies.size() - 1) {
                    line_out << ", ";
                }
            }
            line_out << "\n";
            
            bgzf_output << line_out.str();
            ++processed_count;
        }
        
    } else {
        std::ofstream output_stream(output_file);
        if (!output_stream) {
            std::cerr << "Error: Cannot create output file " << output_file << "\n";
            return false;
        }
        
        while (reader.read_next_kmer(kmer_data)) {

            if (config.use_quantile_normalization) {
                kmer_data.quantile_normalize_frequencies(config.lower_quantile,
                                                        config.upper_quantile,
                                                        config.min_range,
                                                        config.max_range);
            }else if (config.normalize_values) {
                kmer_data.frequencies = BimbamUtils::normalize_to_range(
                    kmer_data.frequencies, 
                    config.min_range, 
                    config.max_range);
            }
            
            output_stream << kmer_data.kmer_sequence << ", " << config.allele_info << ", ";
            for (size_t i = 0; i < kmer_data.frequencies.size(); ++i) {
                output_stream << std::fixed << std::setprecision(6) << kmer_data.frequencies[i];
                if (i < kmer_data.frequencies.size() - 1) {
                    output_stream << ", ";
                }
            }
            output_stream << "\n";
            
            ++processed_count;
        }
    }
    
    return true;
}


EnhancedBimbamConverter::EnhancedBimbamConverter(const BimbamConfig& cfg) : config(cfg) {}

std::vector<std::string> EnhancedBimbamConverter::discover_input_files() {
    std::vector<std::string> extensions;
    
    if (config.input_format == "auto" || config.input_format == "text") {
        extensions.push_back(".txt");
        extensions.push_back(".tsv");
    }
    if (config.input_format == "auto" || config.input_format == "binary" || 
        config.input_format == "compressed") {
        extensions.push_back(".bin");
    }
    
    return get_files_with_prefix(config.input_dir, config.input_prefix, extensions);
}

std::string EnhancedBimbamConverter::determine_file_format(const std::string& filename) {
    if (config.input_format != "auto") {
        return config.input_format;
    }
    
    if (filename.find(".bin") != std::string::npos) {
        return "binary";
    } else {
        return "text";
    }
}

std::string EnhancedBimbamConverter::generate_output_filename(const std::string& input_file) {
    std::string base_name = input_file.substr(input_file.find_last_of('/') + 1);
    
    size_t last_dot = base_name.find_last_of('.');
    if (last_dot != std::string::npos) {
        base_name = base_name.substr(0, last_dot);
    }
    
    if (config.compress_output) {
        return config.output_dir + "/" + base_name + ".bimbam.gz";
    } else {
        return config.output_dir + "/" + base_name + ".bimbam";
    }
}

void EnhancedBimbamConverter::process_single_file(const std::string& input_file, 
                                                 const std::string& output_file) {
    std::string format = determine_file_format(input_file);
    uint64_t processed_count = 0;
    bool success = false;
    
    try {
        if (format == "binary" || format == "compressed") {
            success = BinaryBimbamProcessor::process_binary_file(input_file, output_file, 
                                                               config, processed_count);
        } else {
            success = TextBimbamProcessor::process_text_file(input_file, output_file, 
                                                           config, processed_count);
        }
        
        update_progress(input_file, processed_count, success);
        
    } catch (const std::exception& e) {
        std::cerr << "Error processing " << input_file << ": " << e.what() << "\n";
        update_progress(input_file, 0, false);
    }
}

void EnhancedBimbamConverter::update_progress(const std::string& filename, 
                                             uint64_t processed, bool success) {
    std::lock_guard<std::mutex> lock(progress_mutex);
    
    processed_files.fetch_add(1);
    total_kmers_processed.fetch_add(processed);
    
    if (!success) {
        failed_files.fetch_add(1);
    }
    
    if (config.verbose) {
        std::string base_name = filename.substr(filename.find_last_of('/') + 1);
        if (success) {
            std::cout << "✓ Completed " << base_name << ": " << processed << " k-mers converted\n";
        } else {
            std::cout << "✗ Failed to process " << base_name << "\n";
        }
    }
}

bool EnhancedBimbamConverter::execute_conversion() {
    if (!config.validate()) {
        return false;
    }
    
    if (!create_output_directory(config.output_dir)) {
        std::cerr << "Error: Failed to create output directory\n";
        return false;
    }
    
    auto input_files = discover_input_files();
    if (input_files.empty()) {
        std::cerr << "Error: No input files found with prefix '" << config.input_prefix 
                  << "' in " << config.input_dir << "\n";
        return false;
    }
    
    if (config.verbose) {
        std::cout << "Found " << input_files.size() << " input files\n";
        std::cout << "Processing with " << config.max_threads << " file threads\n";
        if (config.compress_output) {
            std::cout << "BGZF compression: Level " << config.compression_level 
                      << " with " << config.bgzf_threads << " threads per file\n";
            std::cout << "Write buffer: " << (config.write_buffer_size / 1024) << " KB\n";
        } else {
            std::cout << "Compression: Disabled\n";
        }
        std::cout << "\n";
    }
    
    auto start_time = std::chrono::high_resolution_clock::now();
    
    std::vector<std::future<void>> futures;
    
    for (const auto& input_file : input_files) {
        std::string output_file = generate_output_filename(input_file);
        
        if (futures.size() >= static_cast<size_t>(config.max_threads)) {
            futures[0].wait();
            futures.erase(futures.begin());
        }
        
        futures.emplace_back(std::async(std::launch::async,
                                       &EnhancedBimbamConverter::process_single_file,
                                       this, input_file, output_file));
    }
    
    for (auto& future : futures) {
        future.wait();
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    
    std::cout << "\n=== BIMBAM Conversion Complete ===\n";
    std::cout << "Processing time: " << std::fixed << std::setprecision(2) 
              << (duration.count() / 1000.0) << " seconds\n";
    std::cout << "Files processed: " << processed_files.load() << "\n";
    std::cout << "Files failed: " << failed_files.load() << "\n";
    std::cout << "Total k-mers converted: " << total_kmers_processed.load() << "\n";
    
    if (total_kmers_processed.load() > 0 && duration.count() > 0) {
        double kmers_per_sec = (total_kmers_processed.load() * 1000.0) / duration.count();
        std::cout << "Average speed: " 
                  << std::fixed << std::setprecision(0) << kmers_per_sec << " k-mers/second\n";
    }
    
    std::cout << "==================================\n";
    
    return failed_files.load() == 0;
}

void EnhancedBimbamConverter::display_statistics() const {
    std::cout << "\n=== Detailed Statistics ===\n";
    std::cout << "Total files: " << processed_files.load() << "\n";
    std::cout << "Successful: " << (processed_files.load() - failed_files.load()) << "\n";
    std::cout << "Failed: " << failed_files.load() << "\n";
    std::cout << "Total k-mers: " << total_kmers_processed.load() << "\n";
    
    if (processed_files.load() > 0) {
        double success_rate = 100.0 * (processed_files.load() - failed_files.load()) / processed_files.load();
        std::cout << "Success rate: " << std::fixed << std::setprecision(2)
                  << success_rate << "%\n";
    }
    std::cout << "===========================\n";
}


void display_bimbam_help_message(const char* program_name) {
    std::cout << R"(
#--------------------------------------------------------------------------------------------------------------#
#
# @Prog:              K-mer Matrix to BIMBAM Converter
# @Version:           v3.0.0
#
# Usage: 
#             )" << program_name << R"( [options]
#
# Required options:
#             --in  STR       Directory of k-mer matrices
#             --out STR       Output directory for BIMBAM files
#
# Optional parameters:
#             --threads INT       Number of file processing threads (default: 8)
#             --prefix STR        Input file prefix filter (default: "filtered_")
#
# Format control:
#             --input-format STR  Input format: auto|text|binary|compressed (default: auto)
#             --allele-info STR   Allele information string (default: "X, Y")
#
# Normalization options:
#             --no-normalize      Disable value normalization (default: enabled with linear)
#             --min-range FLOAT   Minimum output value (default: 0.0)
#             --max-range FLOAT   Maximum output value (default: 2.0)
#
#             --quantile-norm     Enable quantile-based normalization (handles outliers)
#                                 This method converts extreme values to boundaries:
#                                 - Values ≤ lower quantile → min-range (0.0)
#                                 - Values ≥ upper quantile → max-range (2.0)
#                                 - Values in between → linear transformation [0, 2]
# 
#             --lower-quantile FLOAT  Lower quantile threshold (default: 0.05, i.e., 5th percentile)
#                                      Values below this are set to min-range
# 
#             --upper-quantile FLOAT  Upper quantile threshold (default: 0.95, i.e., 95th percentile)
#                                      Values above this are set to max-range
# 
#
#
# BGZF Compression options:
#             --no-compress       Disable compression (default: BGZF enabled)
#             --level INT         Compression level 0-9 (default: 6)
#                                 0=no compression, 1=fastest, 9=best compression
#             --bgzf-threads INT  BGZF compression threads per file (default: 4)
#             --buffer-size INT   Write buffer size in KB (default: 128)
#                                 Recommended: 64-256 KB for best performance
#
# Other options:
#             --verbose           Verbose output
#             --stats             Output detailed statistics
#             --help              Show this help message
#
# Examples:
#             # Basic conversion with linear normalization (default)
#             )" << program_name << R"( --in filtered_matrices/ --out bimbam_output/
#
#             # Quantile-based normalization (robust to outliers)
#             )" << program_name << R"( --in matrices/ --out bimbam/ --quantile-norm
#
#             # Custom quantile thresholds (1st and 99th percentiles)
#             )" << program_name << R"( --in matrices/ --out bimbam/ --quantile-norm \
#                   --lower-quantile 0.01 --upper-quantile 0.99
#
#             # Quantile normalization with custom output range
#             )" << program_name << R"( --in matrices/ --out bimbam/ --quantile-norm \
#                   --min-range -1.0 --max-range 1.0
#
#             # High-throughput with quantile normalization
#             )" << program_name << R"( --in large_data/ --out output/ --quantile-norm \
#                   --threads 16 --bgzf-threads 8 --level 3 --verbose
#
#             # No normalization (preserve original values)
#             )" << program_name << R"( --in matrices/ --out bimbam/ --no-normalize
#
#             # Maximum compression with quantile normalization
#             )" << program_name << R"( --in matrices/ --out bimbam/ --quantile-norm --level 9
#
#             # No compression (maximum speed) with quantile normalization
#             )" << program_name << R"( --in matrices/ --out bimbam/ --quantile-norm --no-compressFormat:
#
# Output Format:
#             - With BGZF (default): .bimbam.gz files (BGZF format)
#             - Without compression: .bimbam files (plain text)
#
#--------------------------------------------------------------------------------------------------------------#  
)" << std::endl;
}

BimbamConfig parse_bimbam_arguments(int argc, char** argv) {
    BimbamConfig config;
    
    static struct option long_options[] = {
        {"in", required_argument, 0, 1001},
        {"out", required_argument, 0, 1002},
        {"threads", required_argument, 0, 1003},
        {"prefix", required_argument, 0, 1004},
        {"input-format", required_argument, 0, 1005},
        {"no-normalize", no_argument, 0, 1006},
        {"min-range", required_argument, 0, 1007},
        {"max-range", required_argument, 0, 1008},
        {"allele-info", required_argument, 0, 1009},
        {"verbose", no_argument, 0, 1010},
        {"stats", no_argument, 0, 1011},
        {"help", no_argument, 0, 1012},
        {"no-compress", no_argument, 0, 1013},
        {"level", required_argument, 0, 1014},
        {"bgzf-threads", required_argument, 0, 1015},
        {"buffer-size", required_argument, 0, 1016},
        {"quantile-norm", no_argument, 0, 1017},
        {"lower-quantile", required_argument, 0, 1018},
        {"upper-quantile", required_argument, 0, 1019},
        {0, 0, 0, 0}
    };
    
    int opt;
    int option_index = 0;
    
    while ((opt = getopt_long(argc, argv, "", long_options, &option_index)) != -1) {
        switch (opt) {
            case 1001:
                config.input_dir = optarg;
                break;
            case 1002:
                config.output_dir = optarg;
                break;
            case 1003:
                config.max_threads = std::stoi(optarg);
                break;
            case 1004:
                config.input_prefix = optarg;
                break;
            case 1005:
                config.input_format = optarg;
                if (config.input_format != "auto" && config.input_format != "text" && 
                    config.input_format != "binary" && config.input_format != "compressed") {
                    std::cerr << "Error: Invalid input format\n";
                    exit(EXIT_FAILURE);
                }
                break;
            case 1006:
                config.normalize_values = false;
                break;
            case 1007:
                config.min_range = std::stod(optarg);
                break;
            case 1008:
                config.max_range = std::stod(optarg);
                break;
            case 1009:
                config.allele_info = optarg;
                break;
            case 1010:
                config.verbose = true;
                break;
            case 1011:
                config.output_stats = true;
                break;
            case 1012:
                display_bimbam_help_message(argv[0]);
                exit(EXIT_SUCCESS);
            case 1013:
                config.compress_output = false;
                break;
            case 1014:
                config.compression_level = std::stoi(optarg);
                break;
            case 1015:
                config.bgzf_threads = std::stoi(optarg);
                break;
            case 1016: {
                int kb = std::stoi(optarg);
                config.write_buffer_size = kb * 1024;
                break;
            }
            case 1017:
                config.use_quantile_normalization = true;
                config.normalize_values = true;  
                break;
            case 1018:
                config.lower_quantile = std::stod(optarg);
                break;
            case 1019:
                config.upper_quantile = std::stod(optarg);
                break;
            default:
                display_bimbam_help_message(argv[0]);
                exit(EXIT_FAILURE);
        }
    }
    
    return config;
}

bool create_output_directory(const std::string& dir_path) {
    struct stat st;
    if (stat(dir_path.c_str(), &st) == 0) {
        return S_ISDIR(st.st_mode);
    }
    
    return mkdir(dir_path.c_str(), 0755) == 0;
}

std::vector<std::string> get_files_with_prefix(const std::string& dir, 
                                              const std::string& prefix,
                                              const std::vector<std::string>& extensions) {
    std::vector<std::string> files;
    
    struct stat buffer;
    if (stat(dir.c_str(), &buffer) != 0 || !S_ISDIR(buffer.st_mode)) {
        std::cerr << "Warning: Input directory does not exist: " << dir << "\n";
        return files;
    }
    
    DIR* dp = opendir(dir.c_str());
    if (dp == nullptr) {
        std::cerr << "Warning: Cannot open directory " << dir << "\n";
        return files;
    }
    
    struct dirent* entry;
    while ((entry = readdir(dp)) != nullptr) {
        std::string filename = entry->d_name;
        std::string full_path = dir + "/" + filename;
        
        struct stat file_stat;
        if (lstat(full_path.c_str(), &file_stat) == 0) {
            if (S_ISREG(file_stat.st_mode) || S_ISLNK(file_stat.st_mode)) {
                if (filename.find(prefix) == 0) {
                    if (extensions.empty()) {
                        files.push_back(full_path);
                    } else {
                        for (const auto& ext : extensions) {
                            if (filename.find(ext) != std::string::npos) {
                                files.push_back(full_path);
                                break;
                            }
                        }
                    }
                }
            }
        }
    }
    closedir(dp);
    
    std::sort(files.begin(), files.end());
    return files;
}

std::string detect_file_format(const std::string& filename) {
    std::ifstream file(filename, std::ios::binary);
    if (!file.is_open()) {
        return "unknown";
    }
    
    uint32_t kmer_size, sample_count;
    uint64_t total_kmers;
    
    file.read(reinterpret_cast<char*>(&kmer_size), sizeof(uint32_t));
    file.read(reinterpret_cast<char*>(&sample_count), sizeof(uint32_t));
    file.read(reinterpret_cast<char*>(&total_kmers), sizeof(uint64_t));
    
    file.close();
    
    if (kmer_size > 0 && kmer_size <= 255 && 
        sample_count > 0 && sample_count <= 10000 && 
        total_kmers > 0 && total_kmers <= 1000000000ULL) {
        return "binary";
    }
    
    return "text";
}

bool validate_bimbam_output(const std::string& filename) {
    bool is_compressed = (filename.find(".gz") != std::string::npos);
    
    if (is_compressed) {
        BGZF* bgzf_file = bgzf_open(filename.c_str(), "r");
        if (bgzf_file == nullptr) {
            return false;
        }
        
        //char buffer[1024];
        //int line_count = 0;
        kstring_t line = {0, 0, NULL};
        int line_count = 0;
/*
        while (line_count < 10) {
            ssize_t len = bgzf_getline(bgzf_file, '\n', buffer, sizeof(buffer));
            if (len <= 0) break;
            
            std::string line(buffer, len);
            if (line.empty()) continue;
            
            size_t comma_count = std::count(line.begin(), line.end(), ',');
            if (comma_count < 2) {
                bgzf_close(bgzf_file);
                return false;
            }
            ++line_count;
        }
*/
        while (line_count < 10) {
            int ret = bgzf_getline(bgzf_file, '\n', &line);
            if (ret < 0) break; 
            
            if (line.l == 0) continue;
            
            std::string line_str(line.s, line.l);
            size_t comma_count = std::count(line_str.begin(), line_str.end(), ',');
            if (comma_count < 2) {
                ks_free(&line);
                bgzf_close(bgzf_file);
                return false;
            }
            ++line_count;
        }

        ks_free(&line);
        bgzf_close(bgzf_file);
        return line_count > 0;
        
    } else {
        std::ifstream file(filename);
        if (!file.is_open()) {
            return false;
        }
        
        std::string line;
        int line_count = 0;
        while (std::getline(file, line) && line_count < 10) {
            if (line.empty()) continue;
            
            size_t comma_count = std::count(line.begin(), line.end(), ',');
            if (comma_count < 2) {
                return false;
            }
            ++line_count;
        }
        
        return line_count > 0;
    }
}

void FileProcessingStats::display() const {
    std::cout << "\n=== Detailed Processing Statistics ===\n";
    std::cout << "Total files: " << total_files << "\n";
    std::cout << "Successful files: " << successful_files << "\n";
    std::cout << "Failed files: " << failed_files << "\n";
    std::cout << "Total k-mers processed: " << total_kmers << "\n";
    std::cout << "Processing time: " << std::fixed << std::setprecision(2) 
              << processing_time_seconds << " seconds\n";
    
    if (successful_files > 0) {
        std::cout << "Average k-mers per file: " 
                  << (total_kmers / successful_files) << "\n";
    }
    if (processing_time_seconds > 0) {
        std::cout << "Processing rate: " 
                  << std::fixed << std::setprecision(0)
                  << (total_kmers / processing_time_seconds) << " k-mers/second\n";
    }
    if (total_files > 0) {
        double success_rate = 100.0 * successful_files / total_files;
        std::cout << "Success rate: " << std::fixed << std::setprecision(2)
                  << success_rate << "%\n";
    }
    std::cout << "=====================================\n";
}


extern "C" int kmtob(int argc, char* argv[]) {
    try {
        if (argc < 2) {
            display_bimbam_help_message(argv[0]);
            return EXIT_FAILURE;
        }
        
        BimbamConfig config = parse_bimbam_arguments(argc, argv);
        
        if (config.input_dir.empty() || config.output_dir.empty()) {
            std::cerr << "Error: Both input and output directories must be specified\n";
            std::cerr << "Use --help for usage information\n";
            return EXIT_FAILURE;
        }
        
        if (config.verbose) {
            config.display_config();
        }
        
        EnhancedBimbamConverter converter(config);
        
        if (!converter.execute_conversion()) {
            std::cerr << "Error: Conversion process failed\n";
            return EXIT_FAILURE;
        }
        
        if (config.output_stats) {
            converter.display_statistics();
        }
        
        std::cout << "\n✓ All BIMBAM conversion tasks completed successfully.\n";
        return EXIT_SUCCESS;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return EXIT_FAILURE;
    } catch (...) {
        std::cerr << "Error: Unknown exception occurred\n";
        return EXIT_FAILURE;
    }
}

//#ifndef NO_MAIN
//int main(int argc, char* argv[]) {
//    return kmtob(argc, argv);
//}
//#endif
