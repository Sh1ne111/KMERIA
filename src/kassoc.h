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


#ifndef ASSOCIATION_H
#define ASSOCIATION_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <thread>
#include <mutex>
#include <atomic>
#include <future>
#include <chrono>
#include <dirent.h>
#include <sys/stat.h>
#include <cstdlib>
#include <unistd.h>
#include <algorithm>
#include <map>
#include <set>
#include <memory>
#include <iomanip>
#include <regex>
#include <getopt.h>

enum class AssociationTool {
    GEMMA,
    BIMBAMASSO
};

struct AssocConfig {
    std::string input_dir;
    std::string covariate_file;
    std::string kinship_file;
    std::string phenotype_file;
    std::string sample_file;  // New: for bimbam format
    std::string output_dir;
    int max_threads = 8;
    int phenotype_column = 1;
    
    AssociationTool tool = AssociationTool::GEMMA;
    
    std::string analysis_method = "lmm";
    bool use_relatedness_matrix = true;
    bool no_snp_format = true;
    double minor_allele_freq = 0.01;
    double missing_threshold = 0.05;
    bool verbose = false;
    bool dry_run = false;
    
    bool bimbam_gzip = false;
    int kinship_precision = 10;
    bool kinship_use_ibs = false;  // false=BN, true=IBS
    bool kinship_random_fill = false;
    bool kinship_hetero_div = false;
    int kinship_method = 3;  // 1=IBS+mean, 2=IBS+random, 3=BN
    bool write_eigen_files = false;
    bool disable_gls = false;
    int output_precision = 5;
    int start_marker = 0;
    int end_marker = -1;  // -1 means all markers
    
    bool validate_inputs = true;
    bool check_dependencies = true;
    bool cleanup_temp_files = true;
    bool generate_kinship = false;  // Generate kinship before association
    
    bool compress_output = false;
    bool generate_plots = false;
    std::string log_level = "INFO";
    
    void display_config() const;
    bool validate() const;
    std::string get_tool_name() const;
};

class BimbamFormatConverter {
public:
    // Convert phenotype file to bimbam format (FAMID INDID PHENO)
    static bool convert_phenotype_file(
        const std::string& input_file,
        const std::string& output_file,
        int pheno_column
    );
    
    // Generate sample file from bimbam genotype file (FAMID INDID)
    static bool generate_sample_file(
        const std::string& bimbam_file,
        const std::string& output_file,
        bool is_gzipped = false
    );
    
    static std::vector<std::string> extract_sample_names(
        const std::string& bimbam_file,
        bool is_gzipped = false
    );
    
    static bool validate_bimbam_format(const std::string& filename, bool is_gzipped);
};

class DependencyChecker {
private:
    std::map<std::string, std::string> required_tools;
    std::map<std::string, bool> tool_status;

public:
    DependencyChecker();
    
    bool check_all_dependencies(AssociationTool tool);
    bool check_tool(const std::string& tool_name, const std::string& version_cmd);
    void display_dependency_status() const;
    bool are_all_dependencies_met() const;
    
    std::string get_tool_version(const std::string& tool_name) const;
    void add_dependency(const std::string& tool_name, const std::string& check_command);
};

class FileValidator {
public:
    struct ValidationResult {
        bool is_valid;
        std::string error_message;
        size_t line_count;
        size_t column_count;
        
        ValidationResult() : is_valid(false), line_count(0), column_count(0) {}
    };
    
    static ValidationResult validate_bimbam_file(const std::string& filename, bool is_gzipped = false);
    static ValidationResult validate_phenotype_file(const std::string& filename, int expected_columns);
    static ValidationResult validate_covariate_file(const std::string& filename);
    static ValidationResult validate_kinship_file(const std::string& filename);
    static ValidationResult validate_sample_file(const std::string& filename);
    
    static bool check_sample_consistency(const std::vector<std::string>& file_list);
    
private:
    static std::vector<std::string> extract_sample_names(const std::string& filename, const std::string& format);
};

struct AssociationJob {
    std::string input_file;
    std::string output_prefix;
    std::string full_command;
    std::string kinship_file;  // For bimbam: generated or provided
    std::string sample_file;   // For bimbam: generated sample file
    std::chrono::system_clock::time_point start_time;
    std::chrono::system_clock::time_point end_time;
    int exit_code;
    bool completed;
    bool kinship_generated;
    std::string error_message;
    
    AssociationJob() : exit_code(-1), completed(false), kinship_generated(false) {}
    
    double get_execution_time() const;
    bool was_successful() const { return completed && exit_code == 0; }
};

class JobManager {
private:
    std::vector<AssociationJob> job_queue;
    std::vector<std::future<void>> active_jobs;
    std::atomic<size_t> completed_jobs{0};
    std::atomic<size_t> failed_jobs{0};
    std::mutex progress_mutex;
    std::mutex queue_mutex;
    
    const AssocConfig& config;

public:
    explicit JobManager(const AssocConfig& cfg);
    
    void add_job(const AssociationJob& job);
    void execute_all_jobs();
    void display_progress();
    void display_final_statistics() const;
    
    size_t get_total_jobs() const { return job_queue.size(); }
    size_t get_completed_jobs() const { return completed_jobs.load(); }
    size_t get_failed_jobs() const { return failed_jobs.load(); }
    
    // Command construction for different tools
    std::string construct_gemma_command(const AssociationJob& job) const;
    std::string construct_bimbam_kinship_command(const AssociationJob& job) const;
    std::string construct_bimbam_assoc_command(const AssociationJob& job) const;
    
private:
    void execute_single_job(AssociationJob& job);
    void update_progress(const AssociationJob& job);
    bool generate_kinship_for_job(AssociationJob& job);
};

class ResultProcessor {
private:
    std::string output_directory;
    std::vector<std::string> result_files;
    AssociationTool tool;

public:
    explicit ResultProcessor(const std::string& output_dir, AssociationTool t = AssociationTool::GEMMA);
    
    bool collect_results();
    bool generate_summary_statistics();
    bool create_manhattan_plot_data();
    bool compress_results();
    
    void display_result_summary() const;
    
private:
    struct AssociationResult {
        std::string marker_id;
        std::string chromosome;
        int position;
        double p_value;
        double beta;
        double se;
        
        bool is_significant(double threshold = 5e-8) const {
            return p_value < threshold;
        }
    };
    
    std::vector<AssociationResult> parse_result_file(const std::string& filename) const;
    std::vector<AssociationResult> parse_gemma_result(const std::string& filename) const;
    std::vector<AssociationResult> parse_bimbam_result(const std::string& filename) const;
    void write_summary_file(const std::vector<AssociationResult>& results);
};

class EnhancedAssociationAnalyzer {
private:
    AssocConfig config;
    std::unique_ptr<DependencyChecker> dependency_checker;
    std::unique_ptr<JobManager> job_manager;
    std::unique_ptr<ResultProcessor> result_processor;
    
    std::vector<std::string> input_files;
    std::atomic<bool> analysis_running{false};

public:
    explicit EnhancedAssociationAnalyzer(const AssocConfig& cfg);
    ~EnhancedAssociationAnalyzer();
    
    bool initialize();
    bool execute_analysis();
    void cleanup();
    
    void display_progress() const;
    bool is_analysis_running() const { return analysis_running.load(); }
    
private:
    bool prepare_analysis();
    bool validate_all_inputs();
    bool create_output_directory();
    bool prepare_bimbam_files();
    bool generate_global_kinship();
    std::vector<std::string> discover_input_files();
    std::string generate_output_prefix(const std::string& input_file) const;
};

void display_association_help_message(const char* program_name);
AssocConfig parse_association_arguments(int argc, char** argv);
bool create_directory(const std::string& dir_path);
std::string get_base_filename(const std::string& file_path, bool remove_extension = true);
std::vector<std::string> get_files_with_extension(const std::string& dir, const std::string& extension);
bool file_exists_and_readable(const std::string& filename);
std::string get_timestamp_string();
void setup_logging(const std::string& log_level, const std::string& log_file);

class StatUtils {
public:
    static double bonferroni_correction(double p_value, size_t num_tests);
    static double fdr_correction(const std::vector<double>& p_values, size_t index);
    static bool is_genome_wide_significant(double p_value, double threshold = 5e-8);
    static bool is_suggestive(double p_value, double threshold = 1e-5);
    
    static double calculate_r_squared(double beta, double se, size_t n);
    static std::string interpret_effect_size(double r_squared);
};

extern "C" int assoc(int argc, char* argv[]);

#endif
