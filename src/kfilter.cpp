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

#include "kfilter.h"


void FilterConfig::display_config() const {
    std::cout << "=== Filter Configuration ===\n";
    std::cout << "Input directory: " << input_dir << "\n";
    std::cout << "Output directory: " << output_dir << "\n";
    std::cout << "Depth file: " << depth_file << "\n";
    std::cout << "Max processes: " << max_processes << "\n";
    std::cout << "Max coverage: " << max_coverage << "\n";
    std::cout << "Missing rate threshold: " << missing_rate << "\n";
    std::cout << "Genome ploidy: " << ploidy << "\n";
    std::cout << "Input format: " << input_format << "\n";
    std::cout << "Output format: " << output_format << "\n";
    std::cout << "Verbose logging: " << (verbose ? "Enabled" : "Disabled") << "\n";
    std::cout << "============================\n\n";
}

bool FilterConfig::validate() const {
    if (input_dir.empty()) {
        std::cerr << "Error: Input directory not specified\n";
        return false;
    }
    if (output_dir.empty()) {
        std::cerr << "Error: Output directory not specified\n";
        return false;
    }
    if (depth_file.empty()) {
        std::cerr << "Error: Depth file not specified\n";
        return false;
    }
    if (max_coverage <= 0) {
        std::cerr << "Error: Max coverage must be positive\n";
        return false;
    }
    if (missing_rate < 0.0 || missing_rate > 1.0) {
        std::cerr << "Error: Missing rate must be between 0 and 1\n";
        return false;
    }
    if (ploidy <= 0) {
        std::cerr << "Error: Ploidy must be positive\n";
        return false;
    }
    return true;
}

// ============================================================================
// KmerData implementation
// ============================================================================

/*
std::string KmerData::binary_to_string() const {
    std::string result;
    result.reserve(kmer_length);
    
    for (uint32_t i = 0; i < kmer_length; ++i) {
        uint32_t word_idx = i / 32;
        uint32_t bit_offset = (i % 32) * 2;
        uint32_t symbol = (kmer_binary[word_idx] >> bit_offset) & 3;
        
        char nucleotide;
        switch (symbol) {
            case 0: nucleotide = 'A'; break;
            case 1: nucleotide = 'C'; break;
            case 2: nucleotide = 'G'; break;
            case 3: nucleotide = 'T'; break;
            default: nucleotide = 'N'; break;
        }
        result += nucleotide;
    }
    return result;
}
*/

static const char LOOKUP_TABLE[4] = {'A', 'C', 'G', 'T'};
std::string KmerData::binary_to_string() const {

    std::string result(kmer_length, 'N');
    char* data = &result[0]; 
    for (int32_t i = kmer_length - 1; i >= 0; --i) {
        uint32_t word_idx = i >> 5;            
        uint32_t bit_offset = (i & 31) << 1; 
        uint32_t symbol = (kmer_binary[word_idx] >> bit_offset) & 0x03;
        data[i] = LOOKUP_TABLE[symbol];    
    }
    
    return result;
}


static const uint8_t NUCLEOTIDE_TO_CODE[256] = {
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,1,0,0,0,2,0,0,0,0,0,0,0,0,  // 64-79:  A=65, C=67, G=71
    0,0,0,0,3,0,0,0,0,0,0,0,0,0,0,0,  // 80-95:  T=84  
    0,0,0,1,0,0,0,2,0,0,0,0,0,0,0,0,  // 96-111:  a=97, c=99, g=103
    0,0,0,0,3,0,0,0,0,0,0,0,0,0,0,0,  // 112-127: t=116     
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
 };

void KmerData::string_to_binary(const std::string& seq) {
   
    kmer_length = seq.length();

    uint32_t words = (kmer_length + 31) >> 5;
    kmer_binary.assign(words, 0);
    
    
    const char* seq_data = seq.data();
//    uint32_t* binary_data = kmer_binary.data();
      
    for (uint32_t i = 0; i < kmer_length; ++i) {
        
        uint32_t symbol = NUCLEOTIDE_TO_CODE[static_cast<uint8_t>(seq_data[i])];
        
        uint32_t word_idx = i >> 5;           // i / 32
        uint32_t bit_offset = (i & 31) << 1;  // (i % 32) * 2
        
        kmer_binary[word_idx] |= (symbol << bit_offset);
    }
    
    kmer_sequence = seq;
}


/*
void KmerData::string_to_binary(const std::string& seq) {
    kmer_length = seq.length();
    uint32_t words = (kmer_length + 31) / 32;
    kmer_binary.assign(words, 0);
    
    for (uint32_t i = 0; i < kmer_length; ++i) {
        uint32_t symbol;
        switch (seq[i]) {
            case 'A': case 'a': symbol = 0; break;
            case 'C': case 'c': symbol = 1; break;
            case 'G': case 'g': symbol = 2; break;
            case 'T': case 't': symbol = 3; break;
            default: symbol = 0; break;
        }
        
        uint32_t word_idx = i / 32;
        uint32_t bit_offset = (i % 32) * 2;
        kmer_binary[word_idx] |= (static_cast<uint64_t>(symbol) << bit_offset);
    }
    kmer_sequence = seq;
}
*/

bool KmerData::apply_filtering(const std::vector<double>& depths, 
                              const std::vector<int>& ploidies,
                              const FilterConfig& config, 
                              double min_haploid_depth) {
    // Depth correction
    for (size_t i = 0; i < frequencies.size() && i < depths.size(); ++i) {
        double haploid_depth = depths[i] / ploidies[i];
        double corrected_count = static_cast<double>(frequencies[i]) / haploid_depth;
        frequencies[i] = static_cast<uint32_t>(std::round(corrected_count));
    }
    
    int zero_count = std::count_if(frequencies.begin(), frequencies.end(), 
                                  [](uint32_t x) { return x == 0; });
    
    uint32_t max_value = *std::max_element(frequencies.begin(), frequencies.end());
    double threshold_depth = config.max_coverage / min_haploid_depth;
    double zero_percentage = static_cast<double>(zero_count) / frequencies.size();
    
    if (zero_percentage > config.missing_rate) {
        return false;
    }
    if (max_value > threshold_depth) {
        return false;
    }
    
    return true;
}

// ============================================================================
// FrequencyCompressor implementation
// ============================================================================
/*
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
*/
std::vector<uint32_t> FrequencyCompressor::decode_8bit(const std::vector<uint8_t>& encoded) {
    std::vector<uint32_t> decoded;
    decoded.reserve(encoded.size());
    for (uint8_t val : encoded) {
        decoded.push_back(static_cast<uint32_t>(val));
    }
    return decoded;
}

std::vector<uint32_t> FrequencyCompressor::decode_16bit(const std::vector<uint16_t>& encoded) {
    std::vector<uint32_t> decoded;
    decoded.reserve(encoded.size());
    for (uint16_t val : encoded) {
        decoded.push_back(static_cast<uint32_t>(val));
    }
    return decoded;
}

/*
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
    
    size_t min_size = size_32bit;
    EncodingType best_type = FULL_32BIT;
    
    if (max_val <= 255 && size_8bit < min_size) {
        min_size = size_8bit;
        best_type = FULL_8BIT;
    }
    
    if (max_val <= 65535 && size_16bit < min_size) {
        min_size = size_16bit;
        best_type = FULL_16BIT;
    }
    
    if (non_zero_count < sample_count * 0.3 && size_sparse < min_size) {
        min_size = size_sparse;
        best_type = SPARSE;
    }
    
    estimated_size = min_size;
    return best_type;
}
*/



/*
CompressedBinaryKmerReader::CompressedBinaryKmerReader(const std::string& filename) 
    : current_position(0), is_valid(false), is_compressed(false) {
    
    file_stream.open(filename, std::ios::binary);
    if (!file_stream.is_open()) {
        std::cerr << "Error: Cannot open binary file " << filename << "\n";
        return;
    }
    
    // Read header
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
        // No compression flag, this is an old format file
        file_stream.clear();
        file_stream.seekg(sizeof(uint32_t) * 2 + sizeof(uint64_t), std::ios::beg);
        is_compressed = false;
    } else {
        is_compressed = (compression_flag != 0);
    }
    
    //words_per_kmer = (kmer_size + 31) / 32;
    words_per_kmer = (kmer_size + 31) >> 5;
    is_valid = true;
}

CompressedBinaryKmerReader::~CompressedBinaryKmerReader() {
    if (file_stream.is_open()) {
        file_stream.close();
    }
}
*/

bool CompressedBinaryKmerReader::read_uncompressed_kmer(KmerData& kmer_data) {
    file_stream.read(reinterpret_cast<char*>(kmer_data.kmer_binary.data()),
                    sizeof(uint64_t) * words_per_kmer);
    
    file_stream.read(reinterpret_cast<char*>(kmer_data.frequencies.data()),
                    sizeof(uint32_t) * sample_count);
    
    return !file_stream.fail();
}

bool CompressedBinaryKmerReader::read_compressed_kmer(KmerData& kmer_data) {
    file_stream.read(reinterpret_cast<char*>(kmer_data.kmer_binary.data()),
                    sizeof(uint64_t) * words_per_kmer);
    
    if (file_stream.fail()) {
        return false;
    }
    
    // Read encoding type
    uint8_t encoding_type;
    file_stream.read(reinterpret_cast<char*>(&encoding_type), sizeof(uint8_t));
    
    if (file_stream.fail()) {
        return false;
    }
    
    // Decode based on encoding type
    switch (static_cast<FrequencyCompressor::EncodingType>(encoding_type)) {
        case FrequencyCompressor::FULL_8BIT: {
            std::vector<uint8_t> encoded(sample_count);
            file_stream.read(reinterpret_cast<char*>(encoded.data()), sample_count);
            if (file_stream.fail()) return false;
            kmer_data.frequencies = FrequencyCompressor::decode_8bit(encoded);
            break;
        }
        
        case FrequencyCompressor::FULL_16BIT: {
            std::vector<uint16_t> encoded(sample_count);
            file_stream.read(reinterpret_cast<char*>(encoded.data()), 
                           sample_count * sizeof(uint16_t));
            if (file_stream.fail()) return false;
            kmer_data.frequencies = FrequencyCompressor::decode_16bit(encoded);
            break;
        }
        
        case FrequencyCompressor::SPARSE: {
            uint16_t count;
            file_stream.read(reinterpret_cast<char*>(&count), sizeof(uint16_t));
            if (file_stream.fail()) return false;
            
            FrequencyCompressor::SparseEncoding sparse;
            sparse.indices.resize(count);
            sparse.values.resize(count);
            
            file_stream.read(reinterpret_cast<char*>(sparse.indices.data()), 
                           count * sizeof(uint16_t));
            file_stream.read(reinterpret_cast<char*>(sparse.values.data()), 
                           count * sizeof(uint32_t));
            
            if (file_stream.fail()) return false;
            kmer_data.frequencies = FrequencyCompressor::decode_sparse(sparse, sample_count);
            break;
        }
        
        case FrequencyCompressor::FULL_32BIT:
        default:
            file_stream.read(reinterpret_cast<char*>(kmer_data.frequencies.data()),
                           sizeof(uint32_t) * sample_count);
            if (file_stream.fail()) return false;
            break;
    }
    
    return true;
}

bool CompressedBinaryKmerReader::read_next_kmer(KmerData& kmer_data) {
    if (!is_valid || !has_more_data()) {
        return false;
    }
    
    // Initialize kmer_data
    kmer_data = KmerData(kmer_size, sample_count);
    
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

/*
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
*/

// ============================================================================
// CompressedBinaryKmerWriter implementation
// ============================================================================

CompressedBinaryKmerWriter::CompressedBinaryKmerWriter(const std::string& filename, 
                                                      uint32_t k_size, 
                                                      uint32_t num_samples,
                                                      bool compress)
    : kmer_size(k_size), sample_count(num_samples), written_kmers(0), 
      header_written(false), enable_compression(compress),
      total_uncompressed_bytes(0), total_compressed_bytes(0) {
    
    file_stream.open(filename, std::ios::binary);
    if (!file_stream.is_open()) {
        std::cerr << "Error: Cannot create output file " << filename << "\n";
        return;
    }
    
    //words_per_kmer = (kmer_size + 31) / 32;
    words_per_kmer = (kmer_size + 31) >> 5;
    file_stream.write(reinterpret_cast<const char*>(&kmer_size), sizeof(uint32_t));
    file_stream.write(reinterpret_cast<const char*>(&sample_count), sizeof(uint32_t));
    
    kmer_count_position = file_stream.tellp();
    uint64_t placeholder = 0;
    file_stream.write(reinterpret_cast<const char*>(&placeholder), sizeof(uint64_t));
    
    uint8_t compression_flag = enable_compression ? 1 : 0;
    file_stream.write(reinterpret_cast<const char*>(&compression_flag), sizeof(uint8_t));
    
    header_written = true;
}

CompressedBinaryKmerWriter::~CompressedBinaryKmerWriter() {
    if (file_stream.is_open() && header_written) {
        finalize();
    }
}

void CompressedBinaryKmerWriter::write_uncompressed_kmer(const KmerData& kmer_data) {
    file_stream.write(reinterpret_cast<const char*>(kmer_data.kmer_binary.data()),
                     sizeof(uint64_t) * words_per_kmer);
    file_stream.write(reinterpret_cast<const char*>(kmer_data.frequencies.data()),
                     sizeof(uint32_t) * sample_count);
    
    size_t size = sizeof(uint64_t) * words_per_kmer + sizeof(uint32_t) * sample_count;
    total_uncompressed_bytes += size;
    total_compressed_bytes += size;
}

void CompressedBinaryKmerWriter::write_compressed_kmer(const KmerData& kmer_data) {
    file_stream.write(reinterpret_cast<const char*>(kmer_data.kmer_binary.data()),
                     sizeof(uint64_t) * words_per_kmer);
    
    size_t estimated_size;
    auto encoding_type = FrequencyCompressor::choose_best_encoding(kmer_data.frequencies, 
                                                                   estimated_size);
    
    file_stream.write(reinterpret_cast<const char*>(&encoding_type), sizeof(uint8_t));
    
    size_t uncompressed_size = sizeof(uint64_t) * words_per_kmer + sizeof(uint32_t) * sample_count;
    size_t compressed_size = sizeof(uint64_t) * words_per_kmer + 1;
    
    // Encode and write based on chosen type
    switch (encoding_type) {
        case FrequencyCompressor::FULL_8BIT: {
            auto encoded = FrequencyCompressor::encode_frequencies_8bit(kmer_data.frequencies);
            file_stream.write(reinterpret_cast<const char*>(encoded.data()), encoded.size());
            compressed_size += encoded.size();
            break;
        }
        
        case FrequencyCompressor::FULL_16BIT: {
            auto encoded = FrequencyCompressor::encode_frequencies_16bit(kmer_data.frequencies);
            file_stream.write(reinterpret_cast<const char*>(encoded.data()), 
                            encoded.size() * sizeof(uint16_t));
            compressed_size += encoded.size() * sizeof(uint16_t);
            break;
        }
        
        case FrequencyCompressor::SPARSE: {
            auto sparse = FrequencyCompressor::encode_sparse(kmer_data.frequencies);
            uint16_t count = static_cast<uint16_t>(sparse.indices.size());
            file_stream.write(reinterpret_cast<const char*>(&count), sizeof(uint16_t));
            file_stream.write(reinterpret_cast<const char*>(sparse.indices.data()),
                            sparse.indices.size() * sizeof(uint16_t));
            file_stream.write(reinterpret_cast<const char*>(sparse.values.data()),
                            sparse.values.size() * sizeof(uint32_t));
            compressed_size += sizeof(uint16_t) + sparse.indices.size() * sizeof(uint16_t)
                             + sparse.values.size() * sizeof(uint32_t);
            break;
        }
        
        case FrequencyCompressor::FULL_32BIT:
        default:
            file_stream.write(reinterpret_cast<const char*>(kmer_data.frequencies.data()),
                            sizeof(uint32_t) * sample_count);
            compressed_size += sizeof(uint32_t) * sample_count;
            break;
    }
    
    total_uncompressed_bytes += uncompressed_size;
    total_compressed_bytes += compressed_size;
}

bool CompressedBinaryKmerWriter::write_kmer(const KmerData& kmer_data) {
    if (!file_stream.is_open() || !header_written) {
        return false;
    }
    
    if (enable_compression) {
        write_compressed_kmer(kmer_data);
    } else {
        write_uncompressed_kmer(kmer_data);
    }
    
    if (file_stream.fail()) {
        return false;
    }
    
    ++written_kmers;
    return true;
}

void CompressedBinaryKmerWriter::finalize() {
    if (file_stream.is_open() && header_written) {
        // Update k-mer count in header
        auto current_pos = file_stream.tellp();
        file_stream.seekp(kmer_count_position);
        file_stream.write(reinterpret_cast<const char*>(&written_kmers), sizeof(uint64_t));
        file_stream.seekp(current_pos);
        file_stream.close();
        header_written = false;
    }
}

double CompressedBinaryKmerWriter::get_compression_ratio() const {
    if (total_uncompressed_bytes == 0) return 0.0;
    double saved = static_cast<double>(total_uncompressed_bytes - total_compressed_bytes);
    return (saved / total_uncompressed_bytes) * 100.0;
}

// ============================================================================
// Legacy BinaryKmerReader implementation (for backward compatibility)
// ============================================================================

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
    
    //words_per_kmer = (kmer_size + 31) / 32;
    words_per_kmer = (kmer_size + 31) >> 5;
    is_valid = true;
}

BinaryKmerReader::~BinaryKmerReader() {
    if (file_stream.is_open()) {
        file_stream.close();
    }
}

bool BinaryKmerReader::read_next_kmer(KmerData& kmer_data) {
    if (!is_valid || !has_more_data()) {
        return false;
    }
    
    kmer_data = KmerData(kmer_size, sample_count);
    
    file_stream.read(reinterpret_cast<char*>(kmer_data.kmer_binary.data()),
                    sizeof(uint64_t) * words_per_kmer);
    file_stream.read(reinterpret_cast<char*>(kmer_data.frequencies.data()),
                    sizeof(uint32_t) * sample_count);
    
    if (file_stream.fail()) {
        std::cerr << "Error: Failed to read k-mer data at position " << current_position << "\n";
        return false;
    }
    
    kmer_data.kmer_sequence = kmer_data.binary_to_string();
    
    ++current_position;
    return true;
}

/*
void BinaryKmerReader::reset_to_beginning() {
    if (is_valid) {
        file_stream.clear();
        file_stream.seekg(sizeof(uint32_t) * 2 + sizeof(uint64_t), std::ios::beg);
        current_position = 0;
    }
}
*/


// ============================================================================
// Legacy BinaryKmerWriter implementation (for backward compatibility)
// ============================================================================

BinaryKmerWriter::BinaryKmerWriter(const std::string& filename, uint32_t k_size, uint32_t num_samples)
    : kmer_size(k_size), sample_count(num_samples), written_kmers(0), header_written(false) {
    
    file_stream.open(filename, std::ios::binary);
    if (!file_stream.is_open()) {
        std::cerr << "Error: Cannot create output file " << filename << "\n";
        return;
    }
    
    //words_per_kmer = (kmer_size + 31) / 32;
    words_per_kmer = (kmer_size + 31) >> 5;
    file_stream.write(reinterpret_cast<const char*>(&kmer_size), sizeof(uint32_t));
    file_stream.write(reinterpret_cast<const char*>(&sample_count), sizeof(uint32_t));
    
    kmer_count_position = file_stream.tellp();
    uint64_t placeholder = 0;
    file_stream.write(reinterpret_cast<const char*>(&placeholder), sizeof(uint64_t));
    
    header_written = true;
}

BinaryKmerWriter::~BinaryKmerWriter() {
    if (file_stream.is_open() && header_written) {
        finalize();
    }
}

bool BinaryKmerWriter::write_kmer(const KmerData& kmer_data) {
    if (!file_stream.is_open() || !header_written) {
        return false;
    }
    
    file_stream.write(reinterpret_cast<const char*>(kmer_data.kmer_binary.data()),
                     sizeof(uint64_t) * words_per_kmer);
    file_stream.write(reinterpret_cast<const char*>(kmer_data.frequencies.data()),
                     sizeof(uint32_t) * sample_count);
    
    if (file_stream.fail()) {
        return false;
    }
    
    ++written_kmers;
    return true;
}

void BinaryKmerWriter::finalize() {
    if (file_stream.is_open() && header_written) {
        auto current_pos = file_stream.tellp();
        file_stream.seekp(kmer_count_position);
        file_stream.write(reinterpret_cast<const char*>(&written_kmers), sizeof(uint64_t));
        file_stream.seekp(current_pos);
        file_stream.close();
        header_written = false;
    }
}

// ============================================================================
// TextKmerProcessor implementation
// ============================================================================

bool TextKmerProcessor::process_text_file(const std::string& input_file,
                                         const std::string& output_file,
                                         const std::vector<double>& depths,
                                         const std::vector<int>& ploidies,
                                         double min_haploid_depth,
                                         const FilterConfig& config) {
    std::ifstream input_stream(input_file);
    if (!input_stream) {
        std::cerr << "Error: Cannot open input file " << input_file << "\n";
        return false;
    }
    
    std::ofstream output_stream(output_file);
    if (!output_stream) {
        std::cerr << "Error: Cannot create output file " << output_file << "\n";
        return false;
    }
    
    std::string line;
    uint64_t processed = 0, kept = 0;
    
    while (std::getline(input_stream, line)) {
        if (line.find("#KMER") == 0 || line.empty()) {
            continue;
        }
        
        std::istringstream iss(line);
        std::string kmer_seq;
        iss >> kmer_seq;
        
        KmerData kmer_data(kmer_seq.length(), depths.size());
        kmer_data.kmer_sequence = kmer_seq;
        
        uint32_t count;
        size_t idx = 0;
        while (iss >> count && idx < kmer_data.frequencies.size()) {
            kmer_data.frequencies[idx++] = count;
        }
        
        ++processed;
        
        if (kmer_data.apply_filtering(depths, ploidies, config, min_haploid_depth)) {
            output_stream << kmer_data.kmer_sequence;
            for (const auto& freq : kmer_data.frequencies) {
                output_stream << "\t" << freq;
            }
            output_stream << "\n";
            ++kept;
        }
    }
    
    if (config.verbose) {
        std::cout << "Processed " << input_file << ": " << kept << "/" << processed << " k-mers kept\n";
    }
    
    return true;
}

// ============================================================================
// EnhancedKmerProcessor implementation
// ============================================================================

EnhancedKmerProcessor::EnhancedKmerProcessor(const FilterConfig& cfg) 
    : config(cfg), min_haploid_depth(0.0) {
    load_depth_information();
}

void EnhancedKmerProcessor::load_depth_information() {
    std::ifstream depth_file(config.depth_file);
    if (!depth_file) {
        throw std::runtime_error("Cannot open depth file: " + config.depth_file);
    }
    
    std::string line, sample;
    double depth;
    int ploidy;
    bool has_ploidy_column = false;
    while (std::getline(depth_file, line)) {
        std::istringstream iss(line);
        iss >> sample >> depth;
        ploidy = config.ploidy; // default
        std::string ploidy_str;
        if (iss >> ploidy_str) {
            ploidy = std::stoi(ploidy_str);
            has_ploidy_column = true;
        }
        sample_depths.push_back(depth);
        sample_ploidies.push_back(ploidy);
    }
    
    if (sample_depths.empty()) {
        throw std::runtime_error("No depth information found in " + config.depth_file);
    }
    
    min_haploid_depth = std::numeric_limits<double>::max();
    bool using_variable_ploidy = false;
    for (size_t i = 0; i < sample_depths.size(); ++i) {
        double hd = sample_depths[i] / sample_ploidies[i];
        if (hd < min_haploid_depth) {
            min_haploid_depth = hd;
        }
        if (sample_ploidies[i] != config.ploidy) {
            using_variable_ploidy = true;
        }
    }
    
    if (min_haploid_depth <= 0.0) {
        throw std::runtime_error("Invalid haploid depth calculated (division by zero or negative)");
    }
    
    if (config.verbose) {
        std::cout << "Loaded depth information for " << sample_depths.size() << " samples\n";
        std::cout << "Using " << (using_variable_ploidy || has_ploidy_column ? "variable" : "uniform") << " ploidy\n";
        std::cout << "Minimum haploid depth: " << min_haploid_depth << "\n";
    }
}

std::vector<std::string> EnhancedKmerProcessor::discover_input_files() {
    std::vector<std::string> extensions;
    
    if (config.input_format == "auto" || config.input_format == "text") {
        extensions.push_back(".txt");
        extensions.push_back(".tsv");
    }
    if (config.input_format == "auto" || config.input_format == "binary" || 
        config.input_format == "compressed") {
        extensions.push_back(".bin");
    }
    
    return get_files_with_extensions(config.input_dir, extensions);
}

std::string EnhancedKmerProcessor::determine_file_format(const std::string& filename) {
    if (config.input_format != "auto") {
        return config.input_format;
    }
    
    if (filename.find(".bin") != std::string::npos) {
        return "binary"; // Will auto-detect compression when reading
    } else {
        return "text";
    }
}

std::string EnhancedKmerProcessor::generate_output_filename(const std::string& input_file) {
    std::string base_name = input_file.substr(input_file.find_last_of('/') + 1);
    
    size_t last_dot = base_name.find_last_of('.');
    if (last_dot != std::string::npos) {
        base_name = base_name.substr(0, last_dot);
    }
    
    std::string extension;
    if (config.output_format == "compressed" || config.output_format == "binary") {
        extension = ".bin";
    } else {
        extension = ".txt";
    }
    
    return config.output_dir + "/filtered_" + base_name + extension;
}

void EnhancedKmerProcessor::process_single_file(const std::string& input_file, 
                                               const std::string& output_file) {
    std::string format = determine_file_format(input_file);
    uint64_t processed = 0, kept = 0;
    
    try {
        if (format == "binary" || format == "compressed") {
            // Use compressed reader (auto-detects compression)
            CompressedBinaryKmerReader reader(input_file);
            if (!reader.is_open()) {
                throw std::runtime_error("Failed to open binary file");
            }
            
            if (config.verbose && reader.compressed()) {
                std::cout << "Reading compressed file: " << input_file << "\n";
            }
            
            KmerData kmer_data;
            bool use_compression = (config.output_format == "compressed");
            
            if (config.output_format == "binary" || config.output_format == "compressed") {
                // Binary/compressed to binary/compressed
                CompressedBinaryKmerWriter writer(output_file, reader.get_kmer_size(), 
                                                 reader.get_sample_count(), use_compression);
                if (!writer.is_open()) {
                    throw std::runtime_error("Failed to create binary output file");
                }
                
                while (reader.read_next_kmer(kmer_data)) {
                    ++processed;
                    if (kmer_data.apply_filtering(sample_depths, sample_ploidies, config, min_haploid_depth)) {
                        writer.write_kmer(kmer_data);
                        ++kept;
                    }
                }
                writer.finalize();
                
                if (config.verbose && use_compression) {
                    double ratio = writer.get_compression_ratio();
                    std::cout << "  Compression ratio: " << std::fixed << std::setprecision(1) 
                              << ratio << "%\n";
                }
            } else {
                // Binary/compressed to text
                std::ofstream output_stream(output_file);
                if (!output_stream) {
                    throw std::runtime_error("Failed to create text output file");
                }
                
                while (reader.read_next_kmer(kmer_data)) {
                    ++processed;
                    if (kmer_data.apply_filtering(sample_depths, sample_ploidies, config, min_haploid_depth)) {
                        output_stream << kmer_data.kmer_sequence;
                        for (const auto& freq : kmer_data.frequencies) {
                            output_stream << "\t" << freq;
                        }
                        output_stream << "\n";
                        ++kept;
                    }
                }
            }
        } else {
            // Text input
            if (config.output_format == "binary" || config.output_format == "compressed") {
                // Text to binary/compressed
                std::ifstream input_stream(input_file);
                if (!input_stream) {
                    throw std::runtime_error("Failed to open text input file");
                }
                
                // Determine k-mer size
                std::string line, kmer_seq;
                uint32_t kmer_size = 0;
                while (std::getline(input_stream, line) && kmer_size == 0) {
                    if (line.find("#KMER") == 0 || line.empty()) continue;
                    std::istringstream iss(line);
                    iss >> kmer_seq;
                    kmer_size = kmer_seq.length();
                }
                
                if (kmer_size == 0) {
                    throw std::runtime_error("Could not determine k-mer size from file");
                }
                
                input_stream.clear();
                input_stream.seekg(0, std::ios::beg);
                
                bool use_compression = (config.output_format == "compressed");
                CompressedBinaryKmerWriter writer(output_file, kmer_size, 
                                                 sample_depths.size(), use_compression);
                if (!writer.is_open()) {
                    throw std::runtime_error("Failed to create binary output file");
                }
                
                while (std::getline(input_stream, line)) {
                    if (line.find("#KMER") == 0 || line.empty()) continue;
                    
                    std::istringstream iss(line);
                    iss >> kmer_seq;
                    
                    KmerData kmer_data(kmer_size, sample_depths.size());
                    kmer_data.string_to_binary(kmer_seq);
                    
                    uint32_t count;
                    size_t idx = 0;
                    while (iss >> count && idx < kmer_data.frequencies.size()) {
                        kmer_data.frequencies[idx++] = count;
                    }
                    
                    ++processed;
                    if (kmer_data.apply_filtering(sample_depths, sample_ploidies, config, min_haploid_depth)) {
                        writer.write_kmer(kmer_data);
                        ++kept;
                    }
                }
                writer.finalize();
                
                if (config.verbose && use_compression) {
                    double ratio = writer.get_compression_ratio();
                    std::cout << "  Compression ratio: " << std::fixed << std::setprecision(1) 
                              << ratio << "%\n";
                }
            } else {
                // Text to text
                TextKmerProcessor::process_text_file(input_file, output_file, 
                                                   sample_depths, sample_ploidies, min_haploid_depth, config);
                return;
            }
        }
        
        update_progress(input_file, processed, kept);
        
    } catch (const std::exception& e) {
        std::cerr << "Error processing " << input_file << ": " << e.what() << "\n";
    }
}

void EnhancedKmerProcessor::update_progress(const std::string& filename, 
                                           uint64_t processed, uint64_t kept) {
    std::lock_guard<std::mutex> lock(progress_mutex);
    
    processed_files.fetch_add(1);
    total_kmers_processed.fetch_add(processed);
    total_kmers_kept.fetch_add(kept);
    
    if (config.verbose) {
        std::cout << "Completed " << filename << ": " << kept << "/" << processed 
                  << " k-mers kept (" << std::fixed << std::setprecision(2) 
                  << (100.0 * kept / processed) << "%)\n";
    }
}

bool EnhancedKmerProcessor::execute_filtering() {
    if (!config.validate()) {
        return false;
    }
    
    if (!create_output_directory(config.output_dir)) {
        std::cerr << "Error: Failed to create output directory\n";
        return false;
    }
    
    auto input_files = discover_input_files();
    if (input_files.empty()) {
        std::cerr << "Error: No input files found in " << config.input_dir << "\n";
        return false;
    }
    
    if (config.verbose) {
        std::cout << "Found " << input_files.size() << " input files\n";
        std::cout << "Processing with " << config.max_processes << " threads\n\n";
    }
    
    auto start_time = std::chrono::high_resolution_clock::now();
    
    std::vector<std::future<void>> futures;
    
    for (const auto& input_file : input_files) {
        std::string output_file = generate_output_filename(input_file);
        
        if (futures.size() >= static_cast<size_t>(config.max_processes)) {
            futures[0].wait();
            futures.erase(futures.begin());
        }
        
        futures.emplace_back(std::async(std::launch::async,
                                       &EnhancedKmerProcessor::process_single_file,
                                       this, input_file, output_file));
    }
    
    for (auto& future : futures) {
        future.wait();
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
    
    std::cout << "\n=== Filtering Complete ===\n";
    std::cout << "Processing time: " << duration.count() << " seconds\n";
    std::cout << "Files processed: " << processed_files.load() << "\n";
    std::cout << "Total k-mers processed: " << total_kmers_processed.load() << "\n";
    std::cout << "Total k-mers kept: " << total_kmers_kept.load() << "\n";
    if (total_kmers_processed.load() > 0) {
        std::cout << "Overall retention rate: " << std::fixed << std::setprecision(2)
                  << (100.0 * total_kmers_kept.load() / total_kmers_processed.load()) << "%\n";
    }
    
    return true;
}

void EnhancedKmerProcessor::display_statistics() const {
    std::cout << "\n=== Processing Statistics ===\n";
    std::cout << "Files processed: " << processed_files.load() << "\n";
    std::cout << "Total k-mers processed: " << total_kmers_processed.load() << "\n";
    std::cout << "Total k-mers kept: " << total_kmers_kept.load() << "\n";
    if (total_kmers_processed.load() > 0) {
        double retention_rate = 100.0 * total_kmers_kept.load() / total_kmers_processed.load();
        std::cout << "Overall retention rate: " << std::fixed << std::setprecision(2)
                  << retention_rate << "%\n";
    }
    std::cout << "============================\n";
}

// ============================================================================
// Utility functions
// ============================================================================

void display_usage(const char* program_name) {
    std::cout << R"(
#--------------------------------------------------------------------------------------------------------------#
#
# @Prog:              K-mer Matrix Filter with Compression Support                                                                               
# @Version:           v3.0.0
#
# Usage: 
#             )" << program_name << R"( [options]
#
# Required options:
#             -i  STR     Directory of k-mer abundance matrices <input_dir>
#             -o  STR     Output directory <output_dir>
#             -d  STR     Sample depth file <sample_id\tsequence_depth[\tploidy]>
#
# Optional parameters:
#             -t  INT     Number of threads (default: 8)
#             -c  INT     Max abundance of k-mer (default: 1000)
#             -s  FLOAT   Missing ratio threshold (default: 0.8)
#             -p  INT     Genome ploidy (default: 4)
#
# Format control:
#             --input-format STR   Input format: auto|text|binary|compressed (default: auto)
#             --output-format STR  Output format: text|binary|compressed (default: text)
#
# Other options:
#             -v          Verbose output
#             -h          Show this help message
#
# Format Details:
#             auto        - Auto-detect input format (recommended)
#             text        - Plain text tab-delimited format
#             binary      - Uncompressed binary format
#             compressed  - Compressed binary format (30-90% space savings)
#
# Examples:
#             # Process compressed files, output compressed format (preserves compression)
#             )" << program_name << R"( -i matrices/ -o filtered/ -d depths.txt --output-format compressed -v
#
#             # Process any format, output text
#             )" << program_name << R"( -i matrices/ -o filtered/ -d depths.txt
#
#             # Convert compressed to text while filtering
#             )" << program_name << R"( -i matrices/ -o filtered/ -d depths.txt --input-format compressed
#
#             # Process with custom parameters
#             )" << program_name << R"( -i matrices/ -o filtered/ -d depths.txt -t 16 -c 500 -s 0.9 --output-format compressed
#
#--------------------------------------------------------------------------------------------------------------#  
)" << std::endl;
}

FilterConfig parse_filter_arguments(int argc, char** argv) {
    FilterConfig config;
    
    static struct option long_options[] = {
        {"input-format", required_argument, 0, 1001},
        {"output-format", required_argument, 0, 1002},
        {"help", no_argument, 0, 'h'},
        {"verbose", no_argument, 0, 'v'},
        {0, 0, 0, 0}
    };
    
    int opt;
    int option_index = 0;
    
    while ((opt = getopt_long(argc, argv, "hi:c:o:p:s:d:t:v", long_options, &option_index)) != -1) {
        switch (opt) {
            case 'h':
                display_usage(argv[0]);
                exit(EXIT_SUCCESS);
            case 'i':
                config.input_dir = optarg;
                break;
            case 'c':
                config.max_coverage = std::stoi(optarg);
                break;
            case 'o':
                config.output_dir = optarg;
                break;
            case 'p':
                config.ploidy = std::stoi(optarg);
                break;
            case 's':
                config.missing_rate = std::stod(optarg);
                break;
            case 'd':
                config.depth_file = optarg;
                break;
            case 't':
                config.max_processes = std::stoi(optarg);
                break;
            case 'v':
                config.verbose = true;
                break;
            case 1001: // --input-format
                config.input_format = optarg;
                if (config.input_format != "auto" && config.input_format != "text" && 
                    config.input_format != "binary" && config.input_format != "compressed") {
                    std::cerr << "Error: Invalid input format. Use auto, text, binary, or compressed\n";
                    exit(EXIT_FAILURE);
                }
                break;
            case 1002: // --output-format
                config.output_format = optarg;
                if (config.output_format != "text" && config.output_format != "binary" && 
                    config.output_format != "compressed") {
                    std::cerr << "Error: Invalid output format. Use text, binary, or compressed\n";
                    exit(EXIT_FAILURE);
                }
                break;
            default:
                display_usage(argv[0]);
                exit(EXIT_FAILURE);
        }
    }
    
    return config;
}

/*
bool create_output_directory(const std::string& dir_path) {
    struct stat st;
    if (stat(dir_path.c_str(), &st) == 0) {
        return S_ISDIR(st.st_mode);
    }
    
    return mkdir(dir_path.c_str(), 0755) == 0;
}
*/
std::vector<std::string> get_files_with_extensions(const std::string& dir, 
                                                  const std::vector<std::string>& extensions) {
    std::vector<std::string> files;
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
                for (const auto& ext : extensions) {
                    if (filename.find(ext) != std::string::npos) {
                        files.push_back(full_path);
                        break;
                    }
                }
            }
        }
    }
    closedir(dp);
    
    std::sort(files.begin(), files.end());
    return files;
}

/*
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
*/

extern "C" int kfilter(int argc, char* argv[]) {
//int main(int argc, char* argv[]) {
    try {
        if (argc < 2) {
            display_usage(argv[0]);
            return EXIT_FAILURE;
        }
        
        FilterConfig config = parse_filter_arguments(argc, argv);
        
        if (config.verbose) {
            config.display_config();
        }
        
        EnhancedKmerProcessor processor(config);
        
        if (!processor.execute_filtering()) {
            std::cerr << "Error: Filtering process failed\n";
            return EXIT_FAILURE;
        }
        
        std::cout << "All filtering tasks completed successfully.\n";
        return EXIT_SUCCESS;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return EXIT_FAILURE;
    } catch (...) {
        std::cerr << "Error: Unknown exception occurred\n";
        return EXIT_FAILURE;
    }
}
