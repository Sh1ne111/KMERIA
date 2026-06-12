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

#include "kbtog.h"


GzipInputStream::GzipInputStream(const std::string& fname) 
    : gz_file(nullptr), filename(fname), is_open(false), current_line(0) {
    open(fname);
}

GzipInputStream::~GzipInputStream() {
    close();
}

bool GzipInputStream::open(const std::string& fname) {
    if (is_open) {
        close();
    }
    
    filename = fname;
    gz_file = gzopen(filename.c_str(), "rb");
    
    if (gz_file == nullptr) {
        std::cerr << "Error: Cannot open gzip file " << filename << "\n";
        is_open = false;
        return false;
    }
    
    is_open = true;
    current_line = 0;
    return true;
}

void GzipInputStream::close() {
    if (gz_file != nullptr) {
        gzclose(gz_file);
        gz_file = nullptr;
    }
    is_open = false;
}

bool GzipInputStream::getline(std::string& line) {
    if (!is_valid()) {
        return false;
    }
    
    line.clear();
    char buffer[8192];
    
    // Read line using gzgets
    if (gzgets(gz_file, buffer, sizeof(buffer)) == nullptr) {
        return false; // End of file or error
    }
    
    line = buffer;
    
    // Remove trailing newline if present
    if (!line.empty() && line.back() == '\n') {
        line.pop_back();
    }
    if (!line.empty() && line.back() == '\r') {
        line.pop_back();
    }
    
    ++current_line;
    return true;
}

void GzipInputStream::reset() {
    if (is_open && gz_file != nullptr) {
        gzrewind(gz_file);
        current_line = 0;
    }
}

bool GzipInputStream::eof() const {
    if (!is_open || gz_file == nullptr) {
        return true;
    }
    return gzeof(gz_file) != 0;
}


void VcfConfig::display_config() const {
    std::cout << "=== BIMBAM to VCF Converter Configuration ===\n";
    std::cout << "Input BIMBAM file: " << input_file << "\n";
    std::cout << "Sample list file: " << sample_file << "\n";
    std::cout << "Output VCF file: " << output_file << "\n";
    std::cout << "Chromosome: " << chromosome << "\n";
    std::cout << "Start position: " << start_position << "\n";
    std::cout << "Position increment: " << position_increment << "\n";
    std::cout << "Reference allele: " << reference_allele << "\n";
    std::cout << "Alternative allele: " << alternative_allele << "\n";
    std::cout << "Heterozygous thresholds: [" << het_threshold_low << ", " << het_threshold_high << "]\n";
    std::cout << "Homozygous alt threshold: " << hom_alt_threshold << "\n";
    std::cout << "Include header: " << (include_header ? "Yes" : "No") << "\n";
    std::cout << "Validate input: " << (validate_input ? "Yes" : "No") << "\n";
    std::cout << "Auto-detect gzip: " << (auto_detect_gzip ? "Yes" : "No") << "\n";
    std::cout << "Verbose logging: " << (verbose ? "Enabled" : "Disabled") << "\n";
    std::cout << "============================================\n\n";
}

bool VcfConfig::validate() const {
    if (input_file.empty()) {
        std::cerr << "Error: Input BIMBAM file not specified\n";
        return false;
    }
    if (sample_file.empty()) {
        std::cerr << "Error: Sample list file not specified\n";
        return false;
    }
    if (output_file.empty()) {
        std::cerr << "Error: Output VCF file not specified\n";
        return false;
    }
    if (het_threshold_low < 0.0 || het_threshold_low > 2.0) {
        std::cerr << "Error: Invalid heterozygous low threshold (should be 0.0-2.0)\n";
        return false;
    }
    if (het_threshold_high < het_threshold_low) {
        std::cerr << "Error: Heterozygous high threshold must be >= heterozygous low threshold\n";
        return false;
    }
    if (hom_alt_threshold < het_threshold_high) {
        std::cerr << "Error: Homozygous alternative threshold should be >= heterozygous high threshold\n";
        std::cerr << "Current values: het_low=" << het_threshold_low 
                  << ", het_high=" << het_threshold_high 
                  << ", hom_alt=" << hom_alt_threshold << "\n";
        std::cerr << "Suggestion: Use default thresholds or specify: --het-low 0.33 --het-high 1.33 --hom-alt 1.34\n";
        return false;
    }
    if (start_position < 1) {
        std::cerr << "Error: Start position must be positive\n";
        return false;
    }
    return true;
}


void VariantData::convert_to_genotypes(const VcfConfig& config) {
    genotypes.clear();
    genotypes.reserve(raw_values.size());
    
    for (double value : raw_values) {
        genotypes.push_back(GenotypeUtils::value_to_genotype(value, config));
    }
}

VariantData::VariantStats VariantData::get_statistics() const {
    VariantStats stats;
    stats.total_samples = genotypes.size();
    stats.hom_ref_count = 0;
    stats.het_count = 0;
    stats.hom_alt_count = 0;
    stats.missing_count = 0;
    
    for (const std::string& gt : genotypes) {
        if (gt == "0/0") {
            stats.hom_ref_count++;
        } else if (gt == "0/1" || gt == "1/0") {
            stats.het_count++;
        } else if (gt == "1/1") {
            stats.hom_alt_count++;
        } else {
            stats.missing_count++;
        }
    }
    
    // Calculate frequencies
    size_t total_valid = stats.total_samples - stats.missing_count;
    if (total_valid > 0) {
        size_t total_alleles = total_valid * 2;
        size_t alt_alleles = stats.het_count + stats.hom_alt_count * 2;
        stats.maf = static_cast<double>(std::min(alt_alleles, total_alleles - alt_alleles)) / total_alleles;
        stats.het_rate = static_cast<double>(stats.het_count) / total_valid;
    } else {
        stats.maf = 0.0;
        stats.het_rate = 0.0;
    }
    
    stats.missing_rate = static_cast<double>(stats.missing_count) / stats.total_samples;
    
    return stats;
}

bool VariantData::is_valid() const {
    if (variant_id.empty()) return false;
    if (raw_values.empty()) return false;
    if (position < 1) return false;
    
    // Check if all raw values are valid (not NaN, not infinite)
    for (double value : raw_values) {
        if (std::isnan(value) || std::isinf(value)) {
            return false;
        }
    }
    
    return true;
}


bool BimbamReader::detect_gzip_format(const std::string& file_path) {
    if (file_path.size() > 3) {
        std::string ext = file_path.substr(file_path.size() - 3);
        if (ext == ".gz") {
            return true;
        }
    }
    
    // Check magic number (gzip files start with 0x1f 0x8b)
    std::ifstream file(file_path, std::ios::binary);
    if (!file.is_open()) {
        return false;
    }
    
    unsigned char magic[2];
    file.read(reinterpret_cast<char*>(magic), 2);
    file.close();
    
    return (magic[0] == 0x1f && magic[1] == 0x8b);
}

BimbamReader::BimbamReader(const std::string& file_path) 
    : filename(file_path), current_line(0), is_valid(false), is_gzipped(false) {
    
    is_gzipped = detect_gzip_format(file_path);
    
    if (is_gzipped) {
        gz_stream = std::make_unique<GzipInputStream>(filename);
        if (!gz_stream->is_valid()) {
            std::cerr << "Error: Cannot open gzipped BIMBAM file " << filename << "\n";
            return;
        }
    } else {
        file_stream.open(filename);
        if (!file_stream.is_open()) {
            std::cerr << "Error: Cannot open BIMBAM file " << filename << "\n";
            return;
        }
    }
    
    is_valid = true;
}

BimbamReader::~BimbamReader() {
    if (is_gzipped) {
        if (gz_stream) {
            gz_stream->close();
        }
    } else {
        if (file_stream.is_open()) {
            file_stream.close();
        }
    }
}

bool BimbamReader::read_next_variant(VariantData& variant) {
    if (!is_valid) return false;
    
    std::string line;
    bool read_success = false;
    
    // Read line from appropriate stream
    if (is_gzipped) {
        read_success = gz_stream->getline(line);
    } else {
        read_success = static_cast<bool>(std::getline(file_stream, line));
    }
    
    if (!read_success) {
        return false; // End of file or error
    }
    
    ++current_line;
    
    if (line.empty()) {
        return read_next_variant(variant); // Skip empty lines
    }
    
    // Parse BIMBAM format: ID, allele1, allele2, value1, value2, ...
    std::stringstream ss(line);
    std::string token;
    std::vector<std::string> tokens;
    
    while (std::getline(ss, token, ',')) {
        // Trim whitespace
        token.erase(0, token.find_first_not_of(" \t"));
        token.erase(token.find_last_not_of(" \t") + 1);
        tokens.push_back(token);
    }
    
    if (tokens.size() < 3) {
        std::cerr << "Warning: Invalid BIMBAM format at line " << current_line << "\n";
        return read_next_variant(variant); // Skip invalid lines
    }
    
    variant = VariantData();
    variant.variant_id = tokens[0];
    variant.ref_allele = tokens[1];
    variant.alt_allele = tokens[2];
    
    for (size_t i = 3; i < tokens.size(); ++i) {
        try {
            double value = std::stod(tokens[i]);
            variant.raw_values.push_back(value);
        } catch (const std::exception& e) {
            std::cerr << "Warning: Invalid numeric value '" << tokens[i] 
                      << "' at line " << current_line << "\n";
            variant.raw_values.push_back(0.0); // Use default value
        }
    }
    
    return true;
}

void BimbamReader::reset_to_beginning() {
    if (is_valid) {
        if (is_gzipped) {
            gz_stream->reset();
        } else {
            file_stream.clear();
            file_stream.seekg(0, std::ios::beg);
        }
        current_line = 0;
    }
}

// ============================================================================
// VcfWriter implementation
// ============================================================================

VcfWriter::VcfWriter(const std::string& file_path, const std::vector<std::string>& samples)
    : filename(file_path), header_written(false), variants_written(0), sample_names(samples) {
    
    file_stream.open(filename);
    if (!file_stream.is_open()) {
        std::cerr << "Error: Cannot create VCF file " << filename << "\n";
    }
}

VcfWriter::~VcfWriter() {
    finalize();
}

bool VcfWriter::write_header(const VcfConfig& config) {
    if (!file_stream.is_open() || header_written) {
        return false;
    }
    
    write_meta_information(config);
    write_vcf_header(config);
    header_written = true;
    return true;
}

void VcfWriter::write_meta_information(const VcfConfig& config) {
    file_stream << "##fileformat=VCFv4.2\n";
    
    auto now = std::chrono::system_clock::now();
    auto time_t_now = std::chrono::system_clock::to_time_t(now);
    file_stream << "##fileDate=" << std::put_time(std::gmtime(&time_t_now), "%Y%m%d") << "\n";
    
    file_stream << "##source=Enhanced_BIMBAM_to_VCF_Converter_v2.1_with_gzip\n";
    
    file_stream << "##reference=Unknown\n";
    
    file_stream << "##contig=<ID=" << config.chromosome << ">\n";
    
    file_stream << "##INFO=<ID=" << config.info_field << ",Number=0,Type=Flag,Description=\"Presence/absence variant\">\n";
    
    file_stream << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
    
    if (config.include_genotype_quality) {
        file_stream << "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n";
    }
    
    if (config.include_allele_depth) {
        file_stream << "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths for the ref and alt alleles\">\n";
        file_stream << "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth\">\n";
    }
    
    // Filter fields
    file_stream << "##FILTER=<ID=PASS,Description=\"All filters passed\">\n";
}

void VcfWriter::write_vcf_header(const VcfConfig& config) {
    file_stream << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    for (const std::string& sample : sample_names) {
        file_stream << "\t" << sample;
    }
    file_stream << "\n";
}

bool VcfWriter::write_variant(const VariantData& variant, const VcfConfig& config) {
    if (!file_stream.is_open() || !header_written) {
        return false;
    }
    
    // Write fixed fields
    file_stream << variant.chromosome << "\t"
                << variant.position << "\t"
                << variant.variant_id << "\t"
                << variant.ref_allele << "\t"
                << variant.alt_allele << "\t"
                << variant.quality << "\t"
                << variant.filter << "\t"
                << variant.info << "\t";
    
    // Write format field
    std::string format = config.format_field;
    if (config.include_genotype_quality) {
        format += ":GQ";
    }
    if (config.include_allele_depth) {
        format += ":AD:DP";
    }
    file_stream << format;
    
    // Write genotype data
    for (size_t i = 0; i < variant.genotypes.size(); ++i) {
        file_stream << "\t" << variant.genotypes[i];
        
        if (config.include_genotype_quality) {
            int gq = GenotypeUtils::calculate_genotype_quality(variant.raw_values[i], config);
            file_stream << ":" << gq;
        }
        
        if (config.include_allele_depth) {
            auto ad = GenotypeUtils::calculate_allele_depths(variant.raw_values[i]);
            int dp = ad.first + ad.second;
            file_stream << ":" << ad.first << "," << ad.second << ":" << dp;
        }
    }
    file_stream << "\n";
    
    ++variants_written;
    return true;
}

void VcfWriter::finalize() {
    if (file_stream.is_open()) {
        file_stream.close();
    }
}

// ============================================================================
// SampleManager implementation
// ============================================================================

SampleManager::SampleManager(const std::string& sample_file_path) : sample_file(sample_file_path) {}

bool SampleManager::load_samples() {
    std::ifstream file(sample_file);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open sample file " << sample_file << "\n";
        return false;
    }
    
    sample_names.clear();
    sample_indices.clear();
    
    std::string line;
    size_t index = 0;
    while (std::getline(file, line)) {
        // Trim whitespace
        line.erase(0, line.find_first_not_of(" \t\r\n"));
        line.erase(line.find_last_not_of(" \t\r\n") + 1);
        
        if (!line.empty()) {
            sample_names.push_back(line);
            sample_indices[line] = index++;
        }
    }
    
    return !sample_names.empty();
}

void SampleManager::display_sample_info() const {
    std::cout << "=== Sample Information ===\n";
    std::cout << "Number of samples: " << sample_names.size() << "\n";
    std::cout << "Sample names: ";
    if (sample_names.size() <= 10) {
        for (size_t i = 0; i < sample_names.size(); ++i) {
            std::cout << sample_names[i];
            if (i < sample_names.size() - 1) std::cout << ", ";
        }
    } else {
        for (size_t i = 0; i < 5; ++i) {
            std::cout << sample_names[i] << ", ";
        }
        std::cout << "... (and " << (sample_names.size() - 5) << " more)";
    }
    std::cout << "\n==========================\n\n";
}

// ============================================================================
// EnhancedBimbamVcfConverter implementation
// ============================================================================

EnhancedBimbamVcfConverter::EnhancedBimbamVcfConverter(const VcfConfig& cfg) : config(cfg) {
    sample_manager = std::make_unique<SampleManager>(config.sample_file);
}

bool EnhancedBimbamVcfConverter::validate_configuration() {
    if (!config.validate()) {
        return false;
    }
    
    if (!file_exists(config.input_file)) {
        std::cerr << "Error: Input file does not exist: " << config.input_file << "\n";
        return false;
    }
    
    if (!file_exists(config.sample_file)) {
        std::cerr << "Error: Sample file does not exist: " << config.sample_file << "\n";
        return false;
    }
    
    if (!sample_manager->load_samples()) {
        std::cerr << "Error: Failed to load sample information\n";
        return false;
    }
    
    if (config.verbose) {
        // Check if input is gzipped
        bool is_gz = is_gzipped_file(config.input_file);
        std::cout << "Input format: " << (is_gz ? "gzip compressed" : "plain text") << "\n\n";
        sample_manager->display_sample_info();
    }
    
    return true;
}

bool EnhancedBimbamVcfConverter::execute_conversion() {
    if (!validate_configuration()) {
        return false;
    }
    
    if (config.verbose) {
        config.display_config();
    }
    
    auto start_time = std::chrono::high_resolution_clock::now();
    
    // Initialize reader and writer
    BimbamReader reader(config.input_file);
    if (!reader.is_open()) {
        std::cerr << "Error: Failed to open input file\n";
        return false;
    }
    
    if (config.verbose && reader.is_compressed()) {
        std::cout << "Reading gzip compressed BIMBAM file...\n";
    }
    
    VcfWriter writer(config.output_file, sample_manager->get_sample_names());
    if (!writer.is_open()) {
        std::cerr << "Error: Failed to create output file\n";
        return false;
    }
    
    // Write VCF header
    if (config.include_header) {
        writer.write_header(config);
    }
    
    // Process variants
    VariantData variant;
    int64_t current_position = config.start_position;
    
    if (config.verbose) {
        std::cout << "Starting conversion...\n";
    }
    
    while (reader.read_next_variant(variant)) {
        variants_processed.fetch_add(1);
        
        // Validate variant
        bool valid = variant.is_valid();
        if (config.validate_input && !valid) {
            invalid_variants.fetch_add(1);
            if (config.verbose) {
                std::cerr << "Warning: Skipping invalid variant at line " << reader.get_current_line() << "\n";
            }
            continue;
        }
        
        // Set position and chromosome information
        variant.chromosome = config.chromosome;
        variant.position = current_position;
        current_position += config.position_increment;
        
        // Set alleles if not specified in BIMBAM file
        if (variant.ref_allele.empty() || variant.ref_allele == "X") {
            variant.ref_allele = config.reference_allele;
        }
        if (variant.alt_allele.empty() || variant.alt_allele == "Y") {
            variant.alt_allele = config.alternative_allele;
        }
        
        // Convert values to genotypes
        variant.convert_to_genotypes(config);
        
        // Write variant to VCF
        if (writer.write_variant(variant, config)) {
            variants_written.fetch_add(1);
        }
        
        // Progress reporting
        if (config.verbose && variants_processed.load() % 10000 == 0) {
            display_processing_progress(variants_processed.load());
        }
        
        update_progress(valid);
    }
    
    writer.finalize();
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
    
    std::cout << "\n=== Conversion Complete ===\n";
    std::cout << "Processing time: " << duration.count() << " seconds\n";
    std::cout << "Variants processed: " << variants_processed.load() << "\n";
    std::cout << "Variants written: " << variants_written.load() << "\n";
    std::cout << "Invalid variants: " << invalid_variants.load() << "\n";
    
    if (variants_processed.load() > 0) {
        double success_rate = 100.0 * variants_written.load() / variants_processed.load();
        std::cout << "Success rate: " << std::fixed << std::setprecision(2) << success_rate << "%\n";
    }
    
    if (duration.count() > 0) {
        std::cout << "Processing rate: " << (variants_processed.load() / duration.count()) << " variants/second\n";
    }
    
    return variants_written.load() > 0;
}

void EnhancedBimbamVcfConverter::update_progress(bool variant_valid) {
    if (!variant_valid) {
        invalid_variants.fetch_add(1);
    }
}

void EnhancedBimbamVcfConverter::display_processing_progress(uint64_t current, uint64_t total) {
    std::lock_guard<std::mutex> lock(progress_mutex);
    if (total > 0) {
        double percentage = 100.0 * current / total;
        std::cout << "Progress: " << current << "/" << total 
                  << " (" << std::fixed << std::setprecision(1) << percentage << "%)\n";
    } else {
        std::cout << "Processed: " << current << " variants\n";
    }
}

void EnhancedBimbamVcfConverter::display_statistics() const {
    std::cout << "\n=== Conversion Statistics ===\n";
    std::cout << "Total variants processed: " << variants_processed.load() << "\n";
    std::cout << "Variants written to VCF: " << variants_written.load() << "\n";
    std::cout << "Invalid variants skipped: " << invalid_variants.load() << "\n";
    
    if (variants_processed.load() > 0) {
        double success_rate = 100.0 * variants_written.load() / variants_processed.load();
        double error_rate = 100.0 * invalid_variants.load() / variants_processed.load();
        std::cout << "Success rate: " << std::fixed << std::setprecision(2) << success_rate << "%\n";
        std::cout << "Error rate: " << std::fixed << std::setprecision(2) << error_rate << "%\n";
    }
    std::cout << "=============================\n";
}

// ============================================================================
// GenotypeUtils implementation
// ============================================================================

std::string GenotypeUtils::value_to_genotype(double value, const VcfConfig& config) {
    if (std::isnan(value) || value < 0) {
        return "./.";  // Missing genotype
    }
    
    // Standard BIMBAM genotype calling (0-2 scale)
    // 0 = homozygous reference (0/0)
    // 1 = heterozygous (0/1)
    // 2 = homozygous alternative (1/1)
    
    if (value < config.het_threshold_low) {
        return "0/0";  // Homozygous reference
    } else if (value < config.het_threshold_high) {
        return "0/1";  // Heterozygous
    } else {
        return "1/1";  // Homozygous alternative
    }
}

double GenotypeUtils::genotype_to_value(const std::string& genotype) {
    if (genotype == "0/0" || genotype == "0|0") {
        return 0.0;
    } else if (genotype == "0/1" || genotype == "1/0" || genotype == "0|1" || genotype == "1|0") {
        return 1.0;
    } else if (genotype == "1/1" || genotype == "1|1") {
        return 2.0;
    } else {
        return std::nan(""); // Missing or invalid
    }
}

bool GenotypeUtils::is_valid_genotype(const std::string& genotype) {
    std::set<std::string> valid_genotypes = {
        "0/0", "0/1", "1/0", "1/1", 
        "0|0", "0|1", "1|0", "1|1", 
        "./.", ".|."
    };
    return valid_genotypes.count(genotype) > 0;
}

int GenotypeUtils::calculate_genotype_quality(double value, const VcfConfig& config) {
    // Simple quality calculation based on distance from thresholds
    double distance_to_nearest_threshold = std::numeric_limits<double>::max();
    
    // Distance to homozygous reference (0)
    distance_to_nearest_threshold = std::min(distance_to_nearest_threshold, 
                                           std::abs(value - config.het_threshold_low));
    
    // Distance to heterozygous (middle)
    double het_middle = (config.het_threshold_low + config.het_threshold_high) / 2.0;
    distance_to_nearest_threshold = std::min(distance_to_nearest_threshold, 
                                           std::abs(value - het_middle));
    
    // Distance to homozygous alternative
    distance_to_nearest_threshold = std::min(distance_to_nearest_threshold, 
                                           std::abs(value - config.hom_alt_threshold));
    
    // Convert distance to quality score (0-99)
    int quality = static_cast<int>(std::min(99.0, distance_to_nearest_threshold * 100));
    return std::max(1, quality); // Minimum quality of 1
}

std::pair<int, int> GenotypeUtils::calculate_allele_depths(double value, int total_depth) {
    // Mock implementation for allele depths based on BIMBAM value (0-2 scale)
    if (value < 0.67) {
        // Homozygous reference: most reads support reference
        return {static_cast<int>(total_depth * 0.9), static_cast<int>(total_depth * 0.1)};
    } else if (value < 1.33) {
        // Heterozygous: roughly equal depths
        return {total_depth / 2, total_depth / 2};
    } else {
        // Homozygous alternative: most reads support alternative
        return {static_cast<int>(total_depth * 0.1), static_cast<int>(total_depth * 0.9)};
    }
}

// ============================================================================
// ConversionStats implementation
// ============================================================================

void ConversionStats::display() const {
    std::cout << "\n=== Detailed Conversion Statistics ===\n";
    std::cout << "Total variants: " << total_variants << "\n";
    std::cout << "Valid variants: " << valid_variants << "\n";
    std::cout << "Invalid variants: " << invalid_variants << "\n";
    std::cout << "Total genotypes: " << total_genotypes << "\n";
    std::cout << "Homozygous reference (0/0): " << hom_ref_genotypes << "\n";
    std::cout << "Heterozygous (0/1): " << het_genotypes << "\n";
    std::cout << "Homozygous alternative (1/1): " << hom_alt_genotypes << "\n";
    std::cout << "Missing genotypes: " << missing_genotypes << "\n";
    std::cout << "Processing time: " << std::fixed << std::setprecision(2) 
              << processing_time_seconds << " seconds\n";
    
    if (total_genotypes > 0) {
        std::cout << "Genotype distribution:\n";
        std::cout << "  0/0: " << std::setprecision(1) << (100.0 * hom_ref_genotypes / total_genotypes) << "%\n";
        std::cout << "  0/1: " << std::setprecision(1) << (100.0 * het_genotypes / total_genotypes) << "%\n";
        std::cout << "  1/1: " << std::setprecision(1) << (100.0 * hom_alt_genotypes / total_genotypes) << "%\n";
        std::cout << "  ./.: " << std::setprecision(1) << (100.0 * missing_genotypes / total_genotypes) << "%\n";
    }
    
    std::cout << "Conversion rate: " << std::setprecision(2) << get_conversion_rate() << " variants/second\n";
    std::cout << "Missing data rate: " << std::setprecision(2) << get_missing_rate() << "%\n";
    std::cout << "=====================================\n";
}

double ConversionStats::get_conversion_rate() const {
    if (processing_time_seconds > 0) {
        return total_variants / processing_time_seconds;
    }
    return 0.0;
}

double ConversionStats::get_missing_rate() const {
    if (total_genotypes > 0) {
        return 100.0 * missing_genotypes / total_genotypes;
    }
    return 0.0;
}

// ============================================================================
// Utility functions implementation
// ============================================================================

void display_vcf_help_message(const char* program_name) {
    std::cout << R"(
#--------------------------------------------------------------------------------------------------------------#
#
# @Prog:              BIMBAM to VCF Converter                                                                                
# @Version:           v2.1.0
#
# Usage: 
#             " << program_name << R"  [options]
#
# Required options:
#             -i, --input FILE      Input BIMBAM format file (.bimbam or .bimbam.gz)
#             -s, --samples FILE    Sample list file (one sample name per line)
#             -o, --output FILE     Output VCF file (default: out.vcf)
#
# Genotype calling options:
#             --het-low FLOAT       Lower threshold for heterozygous calls (default: 0.33)
#             --het-high FLOAT      Upper threshold for heterozygous calls (default: 1.33)
#             --hom-alt FLOAT       Threshold for homozygous alternative calls (default: 1.34)
#
#             Note: BIMBAM values are on 0-2 scale:
#                   0.0 = homozygous reference (0/0)
#                   1.0 = heterozygous (0/1)
#                   2.0 = homozygous alternative (1/1)
#
# VCF format options:
#             --chr STRING          Chromosome identifier (default: "1")
#             --start-pos INT       Starting position (default: 1)
#             --pos-inc INT         Position increment (default: 1)
#             --ref-allele STRING   Reference allele (default: "A")
#             --alt-allele STRING   Alternative allele (default: "T")
#
# Output options:
#             --no-header           Skip VCF header
#             --include-gq          Include genotype quality scores
#             --include-ad          Include allele depths
#             --filter STRING       Filter field value (default: "PASS")
#
# Quality control:
#             --no-validate         Skip input validation
#             --verbose             Verbose output
#             --stats               Output detailed statistics
#
# Other options:
#             -h, --help            Show this help message
#
# Input Format Support (NEW):
#             - Plain text BIMBAM files (.bimbam, .txt)
#             - Gzip compressed BIMBAM files (.bimbam.gz, .gz)
#             - Auto-detection of compression format
#
# Examples:
#             # Basic conversion from plain text
#             )" << program_name << R"( -i data.bimbam -s samples.txt -o output.vcf
#
#             # Convert from gzip compressed file
#             )" << program_name << R"( -i data.bimbam.gz -s samples.txt -o output.vcf
#
#             # Custom genotype thresholds (strict calling)
#             )" << program_name << R"( -i data.bimbam.gz -s samples.txt -o output.vcf \
#                       --het-low 0.25 --het-high 1.75 --hom-alt 1.75
#
#             # Include quality information with verbose output
#             )" << program_name << R"( -i data.bimbam.gz -s samples.txt -o output.vcf \
#                       --include-gq --include-ad --verbose --stats
#
#             # Custom chromosome and position
#             )" << program_name << R"( -i data.bimbam.gz -s samples.txt -o output.vcf \
#                       --chr "chr1" --start-pos 1000000 --pos-inc 100
#
#--------------------------------------------------------------------------------------------------------------#  
)" << std::endl;
}

VcfConfig parse_vcf_arguments(int argc, char** argv) {
    VcfConfig config;
    
    static struct option long_options[] = {
        {"input", required_argument, 0, 'i'},
        {"samples", required_argument, 0, 's'},
        {"output", required_argument, 0, 'o'},
        {"het-low", required_argument, 0, 1001},
        {"het-high", required_argument, 0, 1002},
        {"hom-alt", required_argument, 0, 1003},
        {"chr", required_argument, 0, 1004},
        {"start-pos", required_argument, 0, 1005},
        {"pos-inc", required_argument, 0, 1006},
        {"ref-allele", required_argument, 0, 1007},
        {"alt-allele", required_argument, 0, 1008},
        {"no-header", no_argument, 0, 1009},
        {"include-gq", no_argument, 0, 1010},
        {"include-ad", no_argument, 0, 1011},
        {"filter", required_argument, 0, 1012},
        {"no-validate", no_argument, 0, 1013},
        {"verbose", no_argument, 0, 'v'},
        {"stats", no_argument, 0, 1014},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };
    
    int opt;
    int option_index = 0;
    
    while ((opt = getopt_long(argc, argv, "i:s:o:vh", long_options, &option_index)) != -1) {
        switch (opt) {
            case 'i':
                config.input_file = optarg;
                break;
            case 's':
                config.sample_file = optarg;
                break;
            case 'o':
                config.output_file = optarg;
                break;
            case 'v':
                config.verbose = true;
                break;
            case 'h':
                display_vcf_help_message(argv[0]);
                exit(EXIT_SUCCESS);
            case 1001: // --het-low
                config.het_threshold_low = std::stod(optarg);
                break;
            case 1002: // --het-high
                config.het_threshold_high = std::stod(optarg);
                break;
            case 1003: // --hom-alt
                config.hom_alt_threshold = std::stod(optarg);
                break;
            case 1004: // --chr
                config.chromosome = optarg;
                break;
            case 1005: // --start-pos
                config.start_position = std::stoll(optarg);
                break;
            case 1006: // --pos-inc
                config.position_increment = std::stoll(optarg);
                break;
            case 1007: // --ref-allele
                config.reference_allele = optarg;
                break;
            case 1008: // --alt-allele
                config.alternative_allele = optarg;
                break;
            case 1009: // --no-header
                config.include_header = false;
                break;
            case 1010: // --include-gq
                config.include_genotype_quality = true;
                break;
            case 1011: // --include-ad
                config.include_allele_depth = true;
                break;
            case 1012: // --filter
                config.filter_value = optarg;
                break;
            case 1013: // --no-validate
                config.validate_input = false;
                break;
            case 1014: // --stats
                config.output_stats = true;
                break;
            default:
                display_vcf_help_message(argv[0]);
                exit(EXIT_FAILURE);
        }
    }
    
    return config;
}

bool validate_bimbam_format(const std::string& filename) {
    // Check if file is gzipped
    bool is_gz = is_gzipped_file(filename);
    
    if (is_gz) {
        // Validate gzipped file
        gzFile gz = gzopen(filename.c_str(), "rb");
        if (gz == nullptr) {
            return false;
        }
        
        char buffer[8192];
        int line_count = 0;
        while (line_count < 10 && gzgets(gz, buffer, sizeof(buffer)) != nullptr) {
            std::string line(buffer);
            if (line.empty()) continue;
            
            // Count commas in BIMBAM format
            size_t comma_count = std::count(line.begin(), line.end(), ',');
            if (comma_count < 2) {
                gzclose(gz);
                return false;
            }
            ++line_count;
        }
        gzclose(gz);
        return line_count > 0;
    } else {
        // Validate plain text file
        std::ifstream file(filename);
        if (!file.is_open()) {
            return false;
        }
        
        std::string line;
        int line_count = 0;
        while (std::getline(file, line) && line_count < 10) {
            if (line.empty()) continue;
            
            // Count commas in BIMBAM format
            size_t comma_count = std::count(line.begin(), line.end(), ',');
            if (comma_count < 2) {
                return false;
            }
            ++line_count;
        }
        return line_count > 0;
    }
}

bool validate_vcf_output(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        return false;
    }
    
    std::string line;
    bool found_header = false;
    while (std::getline(file, line)) {
        if (line.find("#CHROM") == 0) {
            found_header = true;
            break;
        }
    }
    
    return found_header;
}

// Template function implementation
template<typename Container>
std::string join_strings(const Container& container, const std::string& delimiter) {
    std::ostringstream result;
    for (auto it = container.begin(); it != container.end(); ++it) {
        result << *it;
        if (std::next(it) != container.end()) {
            result << delimiter;
        }
    }
    return result.str();
}

// Explicit template instantiation for common types
template std::string join_strings<std::vector<std::string>>(const std::vector<std::string>&, const std::string&);

// File utility functions
bool file_exists(const std::string& filename) {
    std::ifstream file(filename);
    return file.good();
}

size_t count_lines(const std::string& filename) {
    bool is_gz = is_gzipped_file(filename);
    
    if (is_gz) {
        gzFile gz = gzopen(filename.c_str(), "rb");
        if (gz == nullptr) return 0;
        
        size_t count = 0;
        char buffer[8192];
        while (gzgets(gz, buffer, sizeof(buffer)) != nullptr) {
            ++count;
        }
        gzclose(gz);
        return count;
    } else {
        std::ifstream file(filename);
        size_t count = 0;
        std::string line;
        while (std::getline(file, line)) {
            ++count;
        }
        return count;
    }
}

std::string get_file_extension(const std::string& filename) {
    size_t last_dot = filename.find_last_of('.');
    if (last_dot == std::string::npos) {
        return "";
    }
    return filename.substr(last_dot + 1);
}

bool is_gzipped_file(const std::string& filename) {
    if (filename.size() > 3) {
        std::string ext = filename.substr(filename.size() - 3);
        if (ext == ".gz") {
            return true;
        }
    }
    
    std::ifstream file(filename, std::ios::binary);
    if (!file.is_open()) {
        return false;
    }
    
    unsigned char magic[2];
    file.read(reinterpret_cast<char*>(magic), 2);
    file.close();
    
    return (magic[0] == 0x1f && magic[1] == 0x8b);
}


//int main(int argc, char* argv[]) {
extern "C" int btog(int argc, char* argv[]) {
    try {
        if (argc < 2) {
            display_vcf_help_message(argv[0]);
            return EXIT_FAILURE;
        }
        
        VcfConfig config = parse_vcf_arguments(argc, argv);
        
        if (config.input_file.empty() || config.sample_file.empty()) {
            std::cerr << "Error: Both input BIMBAM file and sample list must be specified\n";
            std::cerr << "Use --help for usage information\n";
            return EXIT_FAILURE;
        }
        
        EnhancedBimbamVcfConverter converter(config);
        
        if (!converter.execute_conversion()) {
            std::cerr << "Error: Conversion process failed\n";
            return EXIT_FAILURE;
        }
        
        if (config.output_stats) {
            converter.display_statistics();
        }
        
        std::cout << "BIMBAM to VCF conversion completed successfully.\n";
        std::cout << "Output written to: " << config.output_file << "\n";
        
        return EXIT_SUCCESS;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return EXIT_FAILURE;
    } catch (...) {
        std::cerr << "Error: Unknown exception occurred\n";
        return EXIT_FAILURE;
    }
}
