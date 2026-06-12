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

#ifndef BIMBAM_VCF_H
#define BIMBAM_VCF_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cstdlib>
#include <getopt.h>
#include <cstdio>
#include <algorithm>
#include <memory>
#include <atomic>
#include <mutex>
#include <thread>
#include <future>
#include <chrono>
#include <iomanip>
#include <map>
#include <set>
#include <cmath>
#include <zlib.h>  // Add gzip support


class GzipInputStream {
private:
    gzFile gz_file;
    std::string filename;
    bool is_open;
    int64_t current_line;
    
public:
    explicit GzipInputStream(const std::string& fname);
    ~GzipInputStream();
    
    bool open(const std::string& fname);
    void close();
    bool is_valid() const { return is_open && gz_file != nullptr; }
    
    // Read a line from gzip file
    bool getline(std::string& line);
    
    // Reset to beginning of file
    void reset();
    
    // Check if end of file reached
    bool eof() const;
    
    int64_t get_line_number() const { return current_line; }
};


struct VcfConfig {
    std::string input_file;
    std::string sample_file;
    std::string output_file = "out.vcf";
    std::string chromosome = "1";
    int64_t start_position = 1;
    int64_t position_increment = 1;
    std::string reference_allele = "A";
    std::string alternative_allele = "T";
    double het_threshold_low = 0.33;     // Lower threshold for heterozygous
    double het_threshold_high = 1.33;    // Upper threshold for heterozygous
    double hom_alt_threshold = 1.33;     // Threshold for homozygous alternative
    bool include_header = true;
    bool verbose = false;
    bool validate_input = true;
    bool output_stats = false;
    std::string filter_value = "PASS";
    std::string info_field = "PR";
    std::string format_field = "GT";
    bool include_genotype_quality = false;
    bool include_allele_depth = false;
    int num_threads = 1;  // Single-threaded by default for VCF output
    bool auto_detect_gzip = true;  // NEW: Auto-detect gzip format
    
    void display_config() const;
    bool validate() const;
};

// Structure to hold variant data
struct VariantData {
    std::string variant_id;
    std::string chromosome;
    int64_t position;
    std::string ref_allele;
    std::string alt_allele;
    std::string quality;
    std::string filter;
    std::string info;
    std::string format;
    std::vector<double> raw_values;
    std::vector<std::string> genotypes;
    
    VariantData() : position(0), quality("."), filter("PASS"), info("PR"), format("GT") {}
    
    void convert_to_genotypes(const VcfConfig& config);
    
    struct VariantStats {
        size_t total_samples;
        size_t hom_ref_count;
        size_t het_count;
        size_t hom_alt_count;
        size_t missing_count;
        double maf;  // Minor allele frequency
        double het_rate;
        double missing_rate;
    };
    VariantStats get_statistics() const;
    
    bool is_valid() const;
};

class BimbamReader {
private:
    std::ifstream file_stream;
    std::unique_ptr<GzipInputStream> gz_stream;
    std::string filename;
    uint64_t current_line;
    bool is_valid;
    bool is_gzipped;

    bool detect_gzip_format(const std::string& file_path);

public:
    explicit BimbamReader(const std::string& file_path);
    ~BimbamReader();
    
    bool is_open() const { return is_valid; }
    bool is_compressed() const { return is_gzipped; }
    bool read_next_variant(VariantData& variant);
    void reset_to_beginning();
    uint64_t get_current_line() const { return current_line; }
};

class VcfWriter {
private:
    std::ofstream file_stream;
    std::string filename;
    bool header_written;
    uint64_t variants_written;
    std::vector<std::string> sample_names;

    void write_vcf_header(const VcfConfig& config);
    void write_meta_information(const VcfConfig& config);

public:
    VcfWriter(const std::string& file_path, const std::vector<std::string>& samples);
    ~VcfWriter();
    
    bool is_open() const { return file_stream.is_open(); }
    bool write_variant(const VariantData& variant, const VcfConfig& config);
    bool write_header(const VcfConfig& config);
    uint64_t get_variants_written() const { return variants_written; }
    void finalize();
};

// Sample information manager
class SampleManager {
private:
    std::vector<std::string> sample_names;
    std::map<std::string, size_t> sample_indices;
    std::string sample_file;

public:
    explicit SampleManager(const std::string& sample_file_path);
    
    bool load_samples();
    const std::vector<std::string>& get_sample_names() const { return sample_names; }
    size_t get_sample_count() const { return sample_names.size(); }
    bool is_valid_sample_index(size_t index) const { return index < sample_names.size(); }
    
    // Get sample statistics
    void display_sample_info() const;
};

class EnhancedBimbamVcfConverter {
private:
    VcfConfig config;
    std::unique_ptr<SampleManager> sample_manager;
    std::atomic<uint64_t> variants_processed{0};
    std::atomic<uint64_t> variants_written{0};
    std::atomic<uint64_t> invalid_variants{0};
    std::mutex progress_mutex;

    bool validate_configuration();
    void update_progress(bool variant_valid);
    void display_processing_progress(uint64_t current, uint64_t total = 0);

public:
    explicit EnhancedBimbamVcfConverter(const VcfConfig& cfg);
    
    bool execute_conversion();
    void display_statistics() const;
};

class GenotypeUtils {
public:
    // Convert numeric value to genotype string
    static std::string value_to_genotype(double value, const VcfConfig& config);
    
    // Convert genotype string to numeric value (reverse operation)
    static double genotype_to_value(const std::string& genotype);
    
    static bool is_valid_genotype(const std::string& genotype);
    
    static int calculate_genotype_quality(double value, const VcfConfig& config);
    
    static std::pair<int, int> calculate_allele_depths(double value, int total_depth = 30);
};

struct ConversionStats {
    uint64_t total_variants;
    uint64_t valid_variants;
    uint64_t invalid_variants;
    uint64_t total_genotypes;
    uint64_t hom_ref_genotypes;
    uint64_t het_genotypes;
    uint64_t hom_alt_genotypes;
    uint64_t missing_genotypes;
    double processing_time_seconds;
    
    void display() const;
    double get_conversion_rate() const;
    double get_missing_rate() const;
};

void display_vcf_help_message(const char* program_name);
VcfConfig parse_vcf_arguments(int argc, char** argv);
bool validate_bimbam_format(const std::string& filename);
bool validate_vcf_output(const std::string& filename);

template<typename Container>
std::string join_strings(const Container& container, const std::string& delimiter);

bool file_exists(const std::string& filename);
size_t count_lines(const std::string& filename);
std::string get_file_extension(const std::string& filename);
bool is_gzipped_file(const std::string& filename);

extern "C" int btog(int argc, char* argv[]);

#endif
