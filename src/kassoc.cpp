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

#include "kassoc.h"


std::string AssocConfig::get_tool_name() const {
    return (tool == AssociationTool::GEMMA) ? "GEMMA" : "bimbamAsso";
}

void AssocConfig::display_config() const {
    std::cout << "=== Association Analysis Configuration ===\n";
    std::cout << "Tool: " << get_tool_name() << "\n";
    std::cout << "Input directory: " << input_dir << "\n";
    std::cout << "Phenotype file: " << phenotype_file << "\n";
    std::cout << "Covariate file: " << covariate_file << "\n";
    std::cout << "Kinship file: " << kinship_file << "\n";
    if (tool == AssociationTool::BIMBAMASSO) {
        std::cout << "Sample file: " << sample_file << "\n";
        std::cout << "Bimbam gzipped: " << (bimbam_gzip ? "Yes" : "No") << "\n";
        std::cout << "Kinship method: " << kinship_method << "\n";
        std::cout << "Output precision: " << output_precision << "\n";
    }
    std::cout << "Output directory: " << output_dir << "\n";
    std::cout << "Max threads: " << max_threads << "\n";
    std::cout << "Phenotype column: " << phenotype_column << "\n";
    if (tool == AssociationTool::GEMMA) {
        std::cout << "Analysis method: " << analysis_method << "\n";
        std::cout << "MAF threshold: " << minor_allele_freq << "\n";
        std::cout << "Missing threshold: " << missing_threshold << "\n";
    }
    std::cout << "Use relatedness matrix: " << (use_relatedness_matrix ? "Yes" : "No") << "\n";
    std::cout << "Generate kinship: " << (generate_kinship ? "Yes" : "No") << "\n";
    std::cout << "Verbose logging: " << (verbose ? "Enabled" : "Disabled") << "\n";
    std::cout << "==========================================\n\n";
}

bool AssocConfig::validate() const {
    if (input_dir.empty()) {
        std::cerr << "Error: Input directory not specified\n";
        return false;
    }
    if (phenotype_file.empty()) {
        std::cerr << "Error: Phenotype file not specified\n";
        return false;
    }
    if (output_dir.empty()) {
        std::cerr << "Error: Output directory not specified\n";
        return false;
    }
    if (max_threads < 1 || max_threads > 64) {
        std::cerr << "Error: Invalid number of threads (1-64)\n";
        return false;
    }
    if (phenotype_column < 1) {
        std::cerr << "Error: Phenotype column must be positive\n";
        return false;
    }
    
    if (tool == AssociationTool::GEMMA) {
        if (minor_allele_freq < 0.0 || minor_allele_freq > 0.5) {
            std::cerr << "Error: Invalid MAF threshold (0.0-0.5)\n";
            return false;
        }
    } else if (tool == AssociationTool::BIMBAMASSO) {
        if (sample_file.empty() && !generate_kinship) {
            std::cerr << "Warning: Sample file not specified. Will auto-generate.\n";
        }
        if (kinship_method < 1 || kinship_method > 3) {
            std::cerr << "Error: Invalid kinship method (1-3)\n";
            return false;
        }
    }
    
    return true;
}


bool BimbamFormatConverter::convert_phenotype_file(
    const std::string& input_file,
    const std::string& output_file,
    int pheno_column
) {
    std::ifstream infile(input_file);
    if (!infile.is_open()) {
        std::cerr << "Error: Cannot open phenotype file: " << input_file << "\n";
        return false;
    }
    
    std::ofstream outfile(output_file);
    if (!outfile.is_open()) {
        std::cerr << "Error: Cannot create output file: " << output_file << "\n";
        return false;
    }
    
    std::string line;
    int line_num = 0;
    
    while (std::getline(infile, line)) {
        if (line.empty()) continue;
        
        std::istringstream iss(line);
        std::vector<std::string> tokens;
        std::string token;
        
        while (iss >> token) {
            tokens.push_back(token);
        }
        
        if (tokens.size() < static_cast<size_t>(pheno_column + 1)) {
            std::cerr << "Warning: Line " << line_num << " has insufficient columns\n";
            continue;
        }
        
        outfile << tokens[0] << " " << tokens[0] << " " << tokens[pheno_column] << "\n";
        line_num++;
    }
    
    std::cout << "Converted " << line_num << " samples to bimbam phenotype format\n";
    return true;
}

bool BimbamFormatConverter::generate_sample_file(
    const std::string& bimbam_file,
    const std::string& output_file,
    bool is_gzipped
) {
    auto sample_names = extract_sample_names(bimbam_file, is_gzipped);
    
    if (sample_names.empty()) {
        std::cerr << "Error: Could not extract sample names from " << bimbam_file << "\n";
        return false;
    }
    
    std::ofstream outfile(output_file);
    if (!outfile.is_open()) {
        std::cerr << "Error: Cannot create sample file: " << output_file << "\n";
        return false;
    }
    
    for (size_t i = 0; i < sample_names.size(); ++i) {
        outfile << "sample_" << i << " sample_" << i << "\n";
    }
    
    std::cout << "Generated sample file with " << sample_names.size() << " samples\n";
    return true;
}

std::vector<std::string> BimbamFormatConverter::extract_sample_names(
    const std::string& bimbam_file,
    bool is_gzipped
) {
    std::vector<std::string> samples;
    
    FILE* fp = nullptr;
    if (is_gzipped) {
        std::string cmd = "zcat " + bimbam_file + " 2>/dev/null | head -n 1";
        fp = popen(cmd.c_str(), "r");
    } else {
        fp = fopen(bimbam_file.c_str(), "r");
    }
    
    if (!fp) {
        return samples;
    }
    
    char buffer[65536];
    if (is_gzipped) {
        if (fgets(buffer, sizeof(buffer), fp)) {
            std::string line(buffer);
            size_t comma_count = std::count(line.begin(), line.end(), ',');
            if (comma_count >= 3) {
                size_t n_samples = comma_count - 2;
                for (size_t i = 0; i < n_samples; ++i) {
                    samples.push_back("sample_" + std::to_string(i));
                }
            }
        }
        pclose(fp);
    } else {
        if (fgets(buffer, sizeof(buffer), fp)) {
            std::string line(buffer);
            size_t comma_count = std::count(line.begin(), line.end(), ',');
            if (comma_count >= 3) {
                size_t n_samples = comma_count - 2;
                for (size_t i = 0; i < n_samples; ++i) {
                    samples.push_back("sample_" + std::to_string(i));
                }
            }
        }
        fclose(fp);
    }
    
    return samples;
}

bool BimbamFormatConverter::validate_bimbam_format(const std::string& filename, bool is_gzipped) {
    FILE* fp = nullptr;
    if (is_gzipped) {
        std::string cmd = "zcat " + filename + " 2>/dev/null | head -n 5";
        fp = popen(cmd.c_str(), "r");
    } else {
        fp = fopen(filename.c_str(), "r");
    }
    
    if (!fp) {
        return false;
    }
    
    char buffer[65536];
    int lines_checked = 0;
    size_t expected_commas = 0;
    
    while (lines_checked < 5 && fgets(buffer, sizeof(buffer), fp)) {
        std::string line(buffer);
        if (line.empty()) continue;
        
        size_t comma_count = std::count(line.begin(), line.end(), ',');
        
        if (lines_checked == 0) {
            expected_commas = comma_count;
            if (expected_commas < 3) {
                if (is_gzipped) pclose(fp); else fclose(fp);
                return false;
            }
        } else if (comma_count != expected_commas) {
            if (is_gzipped) pclose(fp); else fclose(fp);
            return false;
        }
        
        lines_checked++;
    }
    
    if (is_gzipped) pclose(fp); else fclose(fp);
    return lines_checked > 0;
}


DependencyChecker::DependencyChecker() {}

bool DependencyChecker::check_all_dependencies(AssociationTool tool) {
    std::cout << "\n============================================\n";
    std::cout << "Checking software dependencies...\n";
    std::cout << "============================================\n";
    
    required_tools.clear();
    
    if (tool == AssociationTool::GEMMA) {
        required_tools["gemma"] = "gemma -h 2>&1";
        required_tools["R"] = "R --version 2>&1";
    } else if (tool == AssociationTool::BIMBAMASSO) {
        required_tools["bimbamKin"] = "bimbamKin 2>&1 | head -1";
        required_tools["bimbamAsso"] = "bimbamAsso 2>&1 | head -1";
        required_tools["R"] = "R --version 2>&1";
    }
    
    bool all_met = true;
    for (const auto& tool_entry : required_tools) {
        bool status = check_tool(tool_entry.first, tool_entry.second);
        tool_status[tool_entry.first] = status;
        
        if (!status) {
            if (tool == AssociationTool::GEMMA && tool_entry.first == "gemma") {
                all_met = false;
            } else if (tool == AssociationTool::BIMBAMASSO && 
                      (tool_entry.first == "bimbamKin" || tool_entry.first == "bimbamAsso")) {
                all_met = false;
            }
        }
    }
    
    display_dependency_status();
    std::cout << "============================================\n\n";
    
    return all_met;
}

bool DependencyChecker::check_tool(const std::string& tool_name, const std::string& version_cmd) {
    FILE* pipe = popen(version_cmd.c_str(), "r");
    if (!pipe) {
        std::cout << tool_name << ":\t\tFAILED (not found)\n";
        return false;
    }
    
    std::string output;
    char buffer[256];
    while (fgets(buffer, sizeof(buffer), pipe) != nullptr) {
        output += buffer;
    }
    int status = pclose(pipe);
    
    bool tool_found = false;
    std::string version_info = "Unknown";
    
    if (tool_name == "gemma") {
        if (output.find("GEMMA") != std::string::npos || status == 0) {
            tool_found = true;
            std::regex version_regex(R"(GEMMA\s+([0-9]+\.[0-9]+(?:\.[0-9]+)?))");
            std::smatch match;
            if (std::regex_search(output, match, version_regex)) {
                version_info = match[1].str();
            }
        }
    } else if (tool_name == "bimbamKin" || tool_name == "bimbamAsso") {
        if (output.find("Usage:") != std::string::npos || output.find(tool_name) != std::string::npos) {
            tool_found = true;
            version_info = "Available";
        }
    } else if (tool_name == "R") {
        if (output.find("R version") != std::string::npos || status == 0) {
            tool_found = true;
            std::regex version_regex(R"(R version ([0-9]+\.[0-9]+\.[0-9]+))");
            std::smatch match;
            if (std::regex_search(output, match, version_regex)) {
                version_info = match[1].str();
            }
        }
    }
    
    if (tool_found) {
        std::cout << tool_name << ":\t\tOK (v" << version_info << ")\n";
    } else {
        std::cout << tool_name << ":\t\tFAILED\n";
    }
    
    return tool_found;
}

void DependencyChecker::display_dependency_status() const {
    std::cout << "\nDependency Summary:\n";
    for (const auto& status : tool_status) {
        std::cout << "  " << status.first << ": " 
                  << (status.second ? "Available" : "Missing") << "\n";
    }
}

bool DependencyChecker::are_all_dependencies_met() const {
    for (const auto& status : tool_status) {
        if (status.first != "R" && !status.second) {
            return false;
        }
    }
    return true;
}

std::string DependencyChecker::get_tool_version(const std::string& tool_name) const {
    return "Unknown";
}

void DependencyChecker::add_dependency(const std::string& tool_name, const std::string& check_command) {
    required_tools[tool_name] = check_command;
}


FileValidator::ValidationResult FileValidator::validate_bimbam_file(const std::string& filename, bool is_gzipped) {
    ValidationResult result;
    
    FILE* fp = nullptr;
    if (is_gzipped) {
        std::string cmd = "zcat " + filename + " 2>/dev/null | head -n 100";
        fp = popen(cmd.c_str(), "r");
    } else {
        fp = fopen(filename.c_str(), "r");
    }
    
    if (!fp) {
        result.error_message = "Cannot open file: " + filename;
        return result;
    }
    
    char buffer[65536];
    size_t line_count = 0;
    size_t expected_columns = 0;
    
    while (line_count < 100 && fgets(buffer, sizeof(buffer), fp)) {
        std::string line(buffer);
        if (line.empty()) continue;
        
        size_t comma_count = std::count(line.begin(), line.end(), ',');
        size_t current_columns = comma_count + 1;
        
        if (line_count == 0) {
            expected_columns = current_columns;
            if (expected_columns < 3) {
                result.error_message = "Invalid BIMBAM format: too few columns";
                if (is_gzipped) pclose(fp); else fclose(fp);
                return result;
            }
        } else if (current_columns != expected_columns) {
            result.error_message = "Inconsistent column count at line " + std::to_string(line_count + 1);
            if (is_gzipped) pclose(fp); else fclose(fp);
            return result;
        }
        
        line_count++;
    }
    
    if (is_gzipped) pclose(fp); else fclose(fp);
    
    result.is_valid = true;
    result.line_count = line_count;
    result.column_count = expected_columns;
    
    return result;
}

FileValidator::ValidationResult FileValidator::validate_phenotype_file(const std::string& filename, int expected_columns) {
    ValidationResult result;
    
    std::ifstream file(filename);
    if (!file.is_open()) {
        result.error_message = "Cannot open phenotype file: " + filename;
        return result;
    }
    
    std::string line;
    size_t line_count = 0;
    
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        
        std::istringstream iss(line);
        std::string token;
        size_t column_count = 0;
        
        while (iss >> token) {
            column_count++;
        }
        
        if (line_count == 0) {
            result.column_count = column_count;
            if (column_count < static_cast<size_t>(expected_columns + 1)) {
                result.error_message = "Insufficient columns in phenotype file";
                return result;
            }
        }
        
        line_count++;
    }
    
    result.is_valid = true;
    result.line_count = line_count;
    
    return result;
}

FileValidator::ValidationResult FileValidator::validate_covariate_file(const std::string& filename) {
    ValidationResult result;
    
    if (filename.empty()) {
        result.is_valid = true;
        return result;
    }
    
    std::ifstream file(filename);
    if (!file.is_open()) {
        result.error_message = "Cannot open covariate file: " + filename;
        return result;
    }
    
    std::string line;
    size_t line_count = 0;
    size_t expected_columns = 0;
    
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        
        std::istringstream iss(line);
        std::string token;
        size_t column_count = 0;
        
        while (iss >> token) {
            column_count++;
        }
        
        if (line_count == 0) {
            expected_columns = column_count;
            if (expected_columns < 2) {
                result.error_message = "Invalid covariate file format";
                return result;
            }
        } else if (column_count != expected_columns) {
            result.error_message = "Inconsistent column count in covariate file";
            return result;
        }
        
        line_count++;
    }
    
    result.is_valid = true;
    result.line_count = line_count;
    result.column_count = expected_columns;
    
    return result;
}

FileValidator::ValidationResult FileValidator::validate_kinship_file(const std::string& filename) {
    ValidationResult result;
    
    if (filename.empty()) {
        result.is_valid = true;
        return result;
    }
    
    std::ifstream file(filename);
    if (!file.is_open()) {
        result.error_message = "Cannot open kinship file: " + filename;
        return result;
    }
    
    std::string line;
    size_t line_count = 0;
    size_t first_row_columns = 0;
    
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        
        std::istringstream iss(line);
        std::string token;
        size_t column_count = 0;
        
        while (iss >> token) {
            column_count++;
        }
        
        if (line_count == 0) {
            first_row_columns = column_count;
        } else if (column_count != first_row_columns) {
            result.error_message = "Non-square kinship matrix";
            return result;
        }
        
        line_count++;
    }
    
    if (line_count != first_row_columns) {
        result.error_message = "Kinship matrix is not square";
        return result;
    }
    
    result.is_valid = true;
    result.line_count = line_count;
    result.column_count = first_row_columns;
    
    return result;
}

FileValidator::ValidationResult FileValidator::validate_sample_file(const std::string& filename) {
    ValidationResult result;
    
    if (filename.empty()) {
        result.is_valid = true;
        return result;
    }
    
    std::ifstream file(filename);
    if (!file.is_open()) {
        result.error_message = "Cannot open sample file: " + filename;
        return result;
    }
    
    std::string line;
    size_t line_count = 0;
    
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        
        std::istringstream iss(line);
        std::string famid, indid;
        
        if (!(iss >> famid >> indid)) {
            result.error_message = "Invalid sample file format at line " + std::to_string(line_count + 1);
            return result;
        }
        
        line_count++;
    }
    
    if (line_count == 0) {
        result.error_message = "Empty sample file";
        return result;
    }
    
    result.is_valid = true;
    result.line_count = line_count;
    result.column_count = 2;
    
    return result;
}

bool FileValidator::check_sample_consistency(const std::vector<std::string>& file_list) {
    return true;
}

std::vector<std::string> FileValidator::extract_sample_names(const std::string& filename, const std::string& format) {
    return std::vector<std::string>();
}


double AssociationJob::get_execution_time() const {
    if (!completed) return 0.0;
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    return duration.count() / 1000.0;
}


JobManager::JobManager(const AssocConfig& cfg) : config(cfg) {}

void JobManager::add_job(const AssociationJob& job) {
    std::lock_guard<std::mutex> lock(queue_mutex);
    job_queue.push_back(job);
}

void JobManager::execute_all_jobs() {
    if (config.verbose) {
        std::cout << "Starting execution of " << job_queue.size() << " association jobs...\n";
    }
    
    for (size_t i = 0; i < job_queue.size(); ++i) {
//        if (active_jobs.size() >= static_cast<size_t>(config.max_threads)) {
        while (active_jobs.size() >= static_cast<size_t>(config.max_threads)) {
            for (auto it = active_jobs.begin(); it != active_jobs.end();) {
                if (it->wait_for(std::chrono::milliseconds(100)) == std::future_status::ready) {
                    it = active_jobs.erase(it);
                } else {
                    ++it;
                }
            }
        }
        
        active_jobs.emplace_back(
            std::async(std::launch::async, &JobManager::execute_single_job, this, std::ref(job_queue[i]))
        );
        
        if (config.verbose && (i + 1) % 10 == 0) {
            display_progress();
        }
    }
    
    for (auto& job : active_jobs) {
        job.wait();
    }
    
    if (config.verbose) {
        display_final_statistics();
    }
}

void JobManager::execute_single_job(AssociationJob& job) {
    job.start_time = std::chrono::system_clock::now();
    
    if (config.tool == AssociationTool::BIMBAMASSO && config.generate_kinship) {
        if (!generate_kinship_for_job(job)) {
            job.exit_code = 1;
            job.completed = true;
            job.error_message = "Failed to generate kinship matrix";
            job.end_time = std::chrono::system_clock::now();
            update_progress(job);
            return;
        }
    }
    
    if (config.tool == AssociationTool::BIMBAMASSO) {
        job.full_command = construct_bimbam_assoc_command(job);
    } else {
        job.full_command = construct_gemma_command(job);
    }
    
    if (config.dry_run) {
        std::cout << "DRY RUN: " << job.full_command << "\n";
        job.exit_code = 0;
        job.completed = true;
    } else {
        job.exit_code = system(job.full_command.c_str());
        job.completed = true;
    }
    
    job.end_time = std::chrono::system_clock::now();
    update_progress(job);
}

void JobManager::update_progress(const AssociationJob& job) {
    std::lock_guard<std::mutex> lock(progress_mutex);
    
    completed_jobs.fetch_add(1);
    
    if (!job.was_successful()) {
        failed_jobs.fetch_add(1);
        if (config.verbose) {
            std::cerr << "Job failed: " << job.input_file 
                      << " (exit code: " << job.exit_code << ")\n";
        }
    } else if (config.verbose) {
        std::cout << "Completed: " << job.output_prefix 
                  << " (" << std::fixed << std::setprecision(2) 
                  << job.get_execution_time() << "s)\n";
    }
}

void JobManager::display_progress() {
    size_t completed = completed_jobs.load();
    size_t total = job_queue.size();
    size_t failed = failed_jobs.load();
    
    double percentage = (total > 0) ? (100.0 * completed / total) : 0.0;
    
    std::cout << "Progress: " << completed << "/" << total 
              << " (" << std::fixed << std::setprecision(1) << percentage << "%)";
    
    if (failed > 0) {
        std::cout << ", Failed: " << failed;
    }
    
    std::cout << "\n";
}

void JobManager::display_final_statistics() const {
    size_t total = job_queue.size();
    size_t completed = completed_jobs.load();
    size_t failed = failed_jobs.load();
    size_t successful = completed - failed;
    
    std::cout << "\n=== Job Execution Summary ===\n";
    std::cout << "Total jobs: " << total << "\n";
    std::cout << "Completed: " << completed << "\n";
    std::cout << "Successful: " << successful << "\n";
    std::cout << "Failed: " << failed << "\n";
    
    if (total > 0) {
        double success_rate = 100.0 * successful / total;
        std::cout << "Success rate: " << std::fixed << std::setprecision(2) 
                  << success_rate << "%\n";
    }
    
    double total_time = 0.0;
    size_t time_count = 0;
    for (const auto& job : job_queue) {
        if (job.completed) {
            total_time += job.get_execution_time();
            time_count++;
        }
    }
    
    if (time_count > 0) {
        double avg_time = total_time / time_count;
        std::cout << "Average execution time: " << std::fixed << std::setprecision(2) 
                  << avg_time << " seconds\n";
    }
    
    std::cout << "=============================\n\n";
}

std::string JobManager::construct_bimbam_kinship_command(const AssociationJob& job) const {
    std::ostringstream cmd;
    
    cmd << "bimbamKin";
    cmd << " -b";
    
    if (config.bimbam_gzip) {
        cmd << " -g";
    }
    
    cmd << " -d " << config.kinship_precision;
    
    if (config.kinship_use_ibs) {
        cmd << " -s";
    }
    
    if (config.kinship_random_fill) {
        cmd << " -r";
    }
    
    if (config.kinship_hetero_div) {
        cmd << " -h";
    }
    
    if (config.verbose) {
        cmd << " -v";
    }
    
    cmd << " " << job.input_file;
    cmd << " " << job.kinship_file;
    
    return cmd.str();
}

std::string JobManager::construct_bimbam_assoc_command(const AssociationJob& job) const {
    std::ostringstream cmd;
    
    cmd << "bimbamAsso";
    cmd << " -b";
    cmd << " -t " << job.input_file;
    cmd << " -p " << config.phenotype_file;
    cmd << " -s " << job.sample_file;
    cmd << " -o " << job.output_prefix;
    
    if (config.bimbam_gzip) {
        cmd << " -g";
    }
    
    if (!job.kinship_file.empty() && config.use_relatedness_matrix) {
        cmd << " -k " << job.kinship_file;
    } else if (config.generate_kinship && config.kinship_method > 0) {
        cmd << " -K " << config.kinship_method;
    }
    
    if (!config.covariate_file.empty()) {
        cmd << " -c " << config.covariate_file;
    }
    
    cmd << " -d " << config.output_precision;
    
    if (config.start_marker > 0) {
        cmd << " -S " << config.start_marker;
    }
    if (config.end_marker > 0) {
        cmd << " -E " << config.end_marker;
    }
    
    if (config.write_eigen_files) {
        cmd << " -w";
    }
    
    if (config.disable_gls) {
        cmd << " -N";
    }
    
    if (config.verbose) {
        cmd << " -v";
    }
    
    return cmd.str();
}

std::string JobManager::construct_gemma_command(const AssociationJob& job) const {
    std::ostringstream cmd;
    
    cmd << "gemma";
    cmd << " -g " << job.input_file;
    cmd << " -p " << config.phenotype_file;
    
    if (!config.covariate_file.empty()) {
        cmd << " -c " << config.covariate_file;
    }
    
    if (!config.kinship_file.empty() && config.use_relatedness_matrix) {
        cmd << " -k " << config.kinship_file;
    }
    
    cmd << " -n " << (config.phenotype_column + 1);
    
    if (config.analysis_method == "lmm") {
        cmd << " -lmm";
    } else if (config.analysis_method == "bslmm") {
        cmd << " -bslmm";
    } else if (config.analysis_method == "loco") {
        cmd << " -loco";
    }
    
    if (config.no_snp_format) {
        cmd << " -notsnp";
    }
    
    cmd << " -maf " << config.minor_allele_freq;
    cmd << " -miss " << config.missing_threshold;
    cmd << " -outdir " << config.output_dir;
    cmd << " -o " << job.output_prefix;
    
    return cmd.str();
}

bool JobManager::generate_kinship_for_job(AssociationJob& job) {
    if (!config.generate_kinship || config.tool != AssociationTool::BIMBAMASSO) {
        return true;
    }
    
    if (!job.kinship_file.empty() && file_exists_and_readable(job.kinship_file)) {
        return true;
    }
    
    job.kinship_file = config.output_dir + "/" + job.output_prefix + ".kinship.txt";
    
    std::string kin_cmd = construct_bimbam_kinship_command(job);
    
    if (config.verbose) {
        std::cout << "Generating kinship matrix: " << kin_cmd << "\n";
    }
    
    if (config.dry_run) {
        std::cout << "DRY RUN (kinship): " << kin_cmd << "\n";
        job.kinship_generated = true;
        return true;
    }
    
    int exit_code = system(kin_cmd.c_str());
    job.kinship_generated = (exit_code == 0);
    
    if (!job.kinship_generated) {
        std::cerr << "Error: Failed to generate kinship matrix for " << job.input_file << "\n";
        return false;
    }
    
    if (config.verbose) {
        std::cout << "Kinship matrix generated: " << job.kinship_file << "\n";
    }
    
    return true;
}


EnhancedAssociationAnalyzer::EnhancedAssociationAnalyzer(const AssocConfig& cfg) : config(cfg) {
    dependency_checker = std::make_unique<DependencyChecker>();
    job_manager = std::make_unique<JobManager>(config);
    result_processor = std::make_unique<ResultProcessor>(config.output_dir, config.tool);
}

EnhancedAssociationAnalyzer::~EnhancedAssociationAnalyzer() = default;

bool EnhancedAssociationAnalyzer::initialize() {
    if (config.verbose) {
        config.display_config();
    }
    
    if (!config.validate()) {
        return false;
    }
    
    if (config.check_dependencies && !dependency_checker->check_all_dependencies(config.tool)) {
        std::cerr << "Error: Required dependencies not met\n";
        return false;
    }
    
    if (!create_output_directory()) {
        std::cerr << "Error: Failed to create output directory\n";
        return false;
    }
    
    return prepare_analysis();
}

bool EnhancedAssociationAnalyzer::prepare_analysis() {
    input_files = discover_input_files();
    
    if (input_files.empty()) {
        std::cerr << "Error: No BIMBAM files found in " << config.input_dir << "\n";
        return false;
    }
    
    if (config.verbose) {
        std::cout << "Found " << input_files.size() << " BIMBAM files for analysis\n";
    }
    
    if (config.tool == AssociationTool::BIMBAMASSO) {
        if (!prepare_bimbam_files()) {
            return false;
        }
    }
    
    if (config.generate_kinship && config.kinship_file.empty() && 
        config.tool == AssociationTool::BIMBAMASSO && input_files.size() == 1) {
        if (!generate_global_kinship()) {
            std::cerr << "Warning: Failed to generate global kinship matrix\n";
        }
    }
    
    for (const auto& file : input_files) {
        AssociationJob job;
        job.input_file = file;
        job.output_prefix = generate_output_prefix(file);
        
        if (!config.kinship_file.empty()) {
            job.kinship_file = config.kinship_file;
        }
        
        if (config.tool == AssociationTool::BIMBAMASSO) {
            if (!config.sample_file.empty()) {
                job.sample_file = config.sample_file;
            } else {
                job.sample_file = config.output_dir + "/" + job.output_prefix + ".samples.txt";
                if (!BimbamFormatConverter::generate_sample_file(file, job.sample_file, config.bimbam_gzip)) {
                    std::cerr << "Error: Failed to generate sample file for " << file << "\n";
                    return false;
                }
            }
        }
        
        job_manager->add_job(job);
    }
    
    return true;
}

bool EnhancedAssociationAnalyzer::prepare_bimbam_files() {
    if (!config.phenotype_file.empty()) {
        std::string bimbam_pheno = config.output_dir + "/phenotype.bimbam.txt";
        
        if (config.verbose) {
            std::cout << "Converting phenotype file to bimbam format...\n";
        }
        
        if (!BimbamFormatConverter::convert_phenotype_file(
                config.phenotype_file, 
                bimbam_pheno, 
                config.phenotype_column)) {
            std::cerr << "Error: Failed to convert phenotype file\n";
            return false;
        }
        
        config.phenotype_file = bimbam_pheno;
    }
    
    if (config.sample_file.empty() && input_files.size() == 1) {
        config.sample_file = config.output_dir + "/samples.txt";
        
        if (config.verbose) {
            std::cout << "Generating sample file from genotype data...\n";
        }
        
        if (!BimbamFormatConverter::generate_sample_file(
                input_files[0], 
                config.sample_file, 
                config.bimbam_gzip)) {
            std::cerr << "Warning: Could not generate sample file automatically\n";
            config.sample_file.clear();
        }
    }
    
    return true;
}

bool EnhancedAssociationAnalyzer::generate_global_kinship() {
    if (input_files.empty()) return false;
    
    std::string kinship_output = config.output_dir + "/global.kinship.txt";
    
    if (config.verbose) {
        std::cout << "Generating global kinship matrix...\n";
    }
    
    AssociationJob kin_job;
    kin_job.input_file = input_files[0];
    kin_job.kinship_file = kinship_output;
    
    std::string kin_cmd = job_manager->construct_bimbam_kinship_command(kin_job);
    
    if (config.verbose) {
        std::cout << "Kinship command: " << kin_cmd << "\n";
    }
    
    if (config.dry_run) {
        std::cout << "DRY RUN (global kinship): " << kin_cmd << "\n";
        config.kinship_file = kinship_output;
        return true;
    }
    
    int exit_code = system(kin_cmd.c_str());
    
    if (exit_code == 0) {
        config.kinship_file = kinship_output;
        if (config.verbose) {
            std::cout << "Global kinship matrix generated successfully\n";
        }
        return true;
    }
    
    return false;
}

bool EnhancedAssociationAnalyzer::validate_all_inputs() {
    if (!config.validate_inputs) {
        return true;
    }
    
    std::cout << "Validating input files...\n";
    
    auto pheno_result = FileValidator::validate_phenotype_file(config.phenotype_file, config.phenotype_column);
    if (!pheno_result.is_valid) {
        std::cerr << "Phenotype file validation failed: " << pheno_result.error_message << "\n";
        return false;
    }
    
    if (!config.covariate_file.empty()) {
        auto covar_result = FileValidator::validate_covariate_file(config.covariate_file);
        if (!covar_result.is_valid) {
            std::cerr << "Covariate file validation failed: " << covar_result.error_message << "\n";
            return false;
        }
    }
    
    if (!config.kinship_file.empty()) {
        auto kinship_result = FileValidator::validate_kinship_file(config.kinship_file);
        if (!kinship_result.is_valid) {
            std::cerr << "Kinship file validation failed: " << kinship_result.error_message << "\n";
            return false;
        }
    }
    
    if (config.tool == AssociationTool::BIMBAMASSO && !config.sample_file.empty()) {
        auto sample_result = FileValidator::validate_sample_file(config.sample_file);
        if (!sample_result.is_valid) {
            std::cerr << "Sample file validation failed: " << sample_result.error_message << "\n";
            return false;
        }
    }
    
    size_t files_to_check = std::min(static_cast<size_t>(5), input_files.size());
    for (size_t i = 0; i < files_to_check; ++i) {
        auto bimbam_result = FileValidator::validate_bimbam_file(input_files[i], config.bimbam_gzip);
        if (!bimbam_result.is_valid) {
            std::cerr << "BIMBAM file validation failed (" << input_files[i] << "): " 
                      << bimbam_result.error_message << "\n";
            return false;
        }
    }
    
    std::cout << "Input validation completed successfully\n";
    return true;
}

bool EnhancedAssociationAnalyzer::execute_analysis() {
    if (!validate_all_inputs()) {
        return false;
    }
    
    analysis_running.store(true);
    
    auto start_time = std::chrono::high_resolution_clock::now();
    
    job_manager->execute_all_jobs();
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
    
    analysis_running.store(false);
    
    std::cout << "\n=== Analysis Complete ===\n";
    std::cout << "Total execution time: " << duration.count() << " seconds\n";
/*    
    if (result_processor->collect_results()) {
        result_processor->display_result_summary();
        
        if (config.generate_plots) {
            result_processor->create_manhattan_plot_data();
        }
        
        if (config.compress_output) {
            result_processor->compress_results();
        }
    }
    
    if (config.cleanup_temp_files) {
        cleanup();
    }
*/
  
    return job_manager->get_failed_jobs() == 0;
}

void EnhancedAssociationAnalyzer::cleanup() {
    if (config.verbose) {
        std::cout << "Cleaning up temporary files...\n";
    }
}

void EnhancedAssociationAnalyzer::display_progress() const {
    if (job_manager) {
        job_manager->display_progress();
    }
}

bool EnhancedAssociationAnalyzer::create_output_directory() {
    return create_directory(config.output_dir);
}

std::vector<std::string> EnhancedAssociationAnalyzer::discover_input_files() {
    return get_files_with_extension(config.input_dir, ".bimbam");
}

std::string EnhancedAssociationAnalyzer::generate_output_prefix(const std::string& input_file) const {
    return get_base_filename(input_file, true);
}



ResultProcessor::ResultProcessor(const std::string& output_dir, AssociationTool t) 
    : output_directory(output_dir), tool(t) {}

bool ResultProcessor::collect_results() {
    if (tool == AssociationTool::GEMMA) {
        result_files = get_files_with_extension(output_directory, ".assoc.txt");
    } else {
        result_files = get_files_with_extension(output_directory, ".ps");
    }
    
    if (result_files.empty()) {
        std::cout << "Warning: No result files found in " << output_directory << "\n";
        return false;
    }
    
    return true;
}

bool ResultProcessor::generate_summary_statistics() {
    if (result_files.empty()) {
        return false;
    }
    
    std::vector<AssociationResult> all_results;
    
    for (const auto& file : result_files) {
        auto file_results = parse_result_file(file);
        all_results.insert(all_results.end(), file_results.begin(), file_results.end());
    }
    
    if (!all_results.empty()) {
        write_summary_file(all_results);
        return true;
    }
    
    return false;
}

void ResultProcessor::display_result_summary() const {
    std::cout << "\n=== Result Summary ===\n";
    std::cout << "Result files found: " << result_files.size() << "\n";
    
    size_t total_tests = 0;
    size_t genome_wide_significant = 0;
    size_t suggestive = 0;
    
    for (const auto& file : result_files) {
        auto results = parse_result_file(file);
        total_tests += results.size();
        
        for (const auto& result : results) {
            if (StatUtils::is_genome_wide_significant(result.p_value)) {
                genome_wide_significant++;
            } else if (StatUtils::is_suggestive(result.p_value)) {
                suggestive++;
            }
        }
    }
    
    std::cout << "Total associations tested: " << total_tests << "\n";
    std::cout << "Genome-wide significant (p < 5e-8): " << genome_wide_significant << "\n";
    std::cout << "Suggestive (p < 1e-5): " << suggestive << "\n";
    
    if (total_tests > 0) {
        double bonferroni_threshold = 0.05 / total_tests;
        std::cout << "Bonferroni threshold: " << std::scientific << bonferroni_threshold << "\n";
    }
    
    std::cout << "======================\n\n";
}

std::vector<ResultProcessor::AssociationResult> ResultProcessor::parse_result_file(const std::string& filename) const {
    if (tool == AssociationTool::GEMMA) {
        return parse_gemma_result(filename);
    } else {
        return parse_bimbam_result(filename);
    }
}

std::vector<ResultProcessor::AssociationResult> ResultProcessor::parse_gemma_result(const std::string& filename) const {
    std::vector<AssociationResult> results;
    
    std::ifstream file(filename);
    if (!file.is_open()) {
        return results;
    }
    
    std::string line;
    bool header_skipped = false;
    
    while (std::getline(file, line)) {
        if (!header_skipped) {
            header_skipped = true;
            continue;
        }
        
        std::istringstream iss(line);
        AssociationResult result;
        
        iss >> result.chromosome >> result.marker_id >> result.position 
            >> result.beta >> result.se >> result.p_value;
        
        if (!iss.fail()) {
            results.push_back(result);
        }
    }
    
    return results;
}

std::vector<ResultProcessor::AssociationResult> ResultProcessor::parse_bimbam_result(const std::string& filename) const {
    std::vector<AssociationResult> results;
    
    std::ifstream file(filename);
    if (!file.is_open()) {
        return results;
    }
    
    std::string line;
    
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        
        std::istringstream iss(line);
        AssociationResult result;
        
        iss >> result.marker_id >> result.beta >> result.p_value;
        
        if (!iss.fail()) {
            result.chromosome = "NA";
            result.position = 0;
            result.se = 0.0;
            results.push_back(result);
        }
    }
    
    return results;
}

void ResultProcessor::write_summary_file(const std::vector<AssociationResult>& results) {
    std::string summary_file = output_directory + "/association_summary.txt";
    
    std::ofstream out(summary_file);
    if (!out.is_open()) {
        std::cerr << "Warning: Cannot create summary file " << summary_file << "\n";
        return;
    }
    
    out << "chromosome\tmarker_id\tposition\tbeta\tse\tp_value\tsignificance\n";
    
    for (const auto& result : results) {
        std::string significance = "NS";
        if (StatUtils::is_genome_wide_significant(result.p_value)) {
            significance = "GWS";
        } else if (StatUtils::is_suggestive(result.p_value)) {
            significance = "SUG";
        }
        
        out << result.chromosome << "\t" << result.marker_id << "\t" 
            << result.position << "\t" << result.beta << "\t" 
            << result.se << "\t" << std::scientific << result.p_value 
            << "\t" << significance << "\n";
    }
}

bool ResultProcessor::create_manhattan_plot_data() {
    std::string script_file = output_directory + "/manhattan_plot.R";
    
    std::ofstream script(script_file);
    if (!script.is_open()) {
        return false;
    }
    
    script << R"R(
library(ggplot2)
library(dplyr)

data <- read.table("association_summary.txt", header=TRUE, stringsAsFactors=FALSE)

manhattan_plot <- ggplot(data, aes(x=position, y=-log10(p_value), color=chromosome)) +
    geom_point(alpha=0.6, size=0.8) +
    geom_hline(yintercept=-log10(5e-8), color="red", linetype="dashed", alpha=0.7) +
    geom_hline(yintercept=-log10(1e-5), color="blue", linetype="dashed", alpha=0.7) +
    labs(title="Manhattan Plot", x="Position", y="-log10(P-value)") +
    theme_minimal() +
    theme(legend.position="none")

ggsave("manhattan_plot.png", manhattan_plot, width=12, height=6, dpi=300)
)R";
    
    script.close();
    
    std::string r_command = "cd " + output_directory + " && R --vanilla < manhattan_plot.R > manhattan_plot.log 2>&1";
    int result = system(r_command.c_str());
    
    return result == 0;
}

bool ResultProcessor::compress_results() {
    std::string tar_command = "cd " + output_directory + " && tar -czf association_results.tar.gz *.ps *.assoc.txt *.log 2>/dev/null";
    return system(tar_command.c_str()) == 0;
}


double StatUtils::bonferroni_correction(double p_value, size_t num_tests) {
    return std::min(1.0, p_value * num_tests);
}

double StatUtils::fdr_correction(const std::vector<double>& p_values, size_t index) {
    return 0.0;
}

bool StatUtils::is_genome_wide_significant(double p_value, double threshold) {
    return p_value < threshold;
}

bool StatUtils::is_suggestive(double p_value, double threshold) {
    return p_value < threshold;
}

double StatUtils::calculate_r_squared(double beta, double se, size_t n) {
    double t_stat = beta / se;
    return (t_stat * t_stat) / (t_stat * t_stat + n - 2);
}

std::string StatUtils::interpret_effect_size(double r_squared) {
    if (r_squared < 0.01) return "negligible";
    if (r_squared < 0.09) return "small";
    if (r_squared < 0.25) return "medium";
    return "large";
}


void display_association_help_message(const char* program_name) {
    std::cout << R"(
#########################################################################################################################
#
# @Prog:              KmersGWAS
# @Version:           v2.1.0
#
# Usage: 
#             )" << program_name << R"( [options]
#
# Required options:
#             -i, --input DIR          Directory containing BIMBAM format files
#             -p, --pheno FILE         Phenotype file
#             -o, --output DIR         Output directory for results
#
# Tool selection:
#             --tool STR               Tool to use: bimbamAsso | gemma
#
# Optional files:
#             -c, --covar FILE         Covariate file
#             -k, --kinship FILE       Kinship/relatedness matrix
#             -s, --sample FILE        Sample file (for bimbamAsso tool)
#             --auto-sample            Auto-generate sample file from genotype data
#
# Analysis options:
#             -n, --ncol INT           Phenotype column number (default: 1)
#             -m, --method STR         Analysis method for GEMMA: lmm|bslmm|loco (default: lmm)
#             --maf FLOAT              Minor allele frequency threshold (default: 0.01)
#             --miss FLOAT             Missing data threshold (default: 0.05)
#             --no-kinship             Don't use kinship matrix even if provided
#
# Bimbam-specific options:
#             --bimbam-gzip            Input bimbam files are gzipped
#             --gen-kinship            Generate kinship matrix before association
#             --kin-method INT         Kinship method: 1=IBS+mean, 2=IBS+random, 3=BN (default: 3)
#             --kin-precision INT      Kinship precision digits (default: 10)
#             --out-precision INT      Output precision digits (default: 5)
#             --start-marker INT       Start marker index (default: 0)
#             --end-marker INT         End marker index (default: all)
#             --write-eigen            Write eigenvalue/eigenvector files
#             --disable-gls            Disable GLS, use OLS instead
#
# Performance options:
#             -t, --threads INT        Number of parallel threads (default: 8)
#             --dry-run                Show commands without executing
#
# Quality control:
#             --no-validate            Skip input file validation
#             --no-check-deps          Skip dependency checking
#
# Output options:
#             --verbose                Verbose output and progress reporting
#             --compress               Compress output files
#             --no-cleanup             Keep temporary files
#
# Other options:
#             -h, --help               Show this help message
#
# Examples:
#             
#             # Using kmeria association tools (bimbamAsso) with pre-computed kinship
#             )" << program_name << R"( --tool bimbamAsso -i bimbam_files/ -p pheno.txt \
#                       -s samples.txt -k kinship.txt --bimbam-gzip -t num_threads -o results/
#
#########################################################################################################################
)" << std::endl;
}

AssocConfig parse_association_arguments(int argc, char** argv) {
    AssocConfig config;
    
    static struct option long_options[] = {
        {"input", required_argument, 0, 'i'},
        {"pheno", required_argument, 0, 'p'},
        {"covar", required_argument, 0, 'c'},
        {"kinship", required_argument, 0, 'k'},
        {"sample", required_argument, 0, 's'},
        {"output", required_argument, 0, 'o'},
        {"threads", required_argument, 0, 't'},
        {"ncol", required_argument, 0, 'n'},
        {"method", required_argument, 0, 'm'},
        {"tool", required_argument, 0, 2001},
        {"maf", required_argument, 0, 1001},
        {"miss", required_argument, 0, 1002},
        {"no-kinship", no_argument, 0, 1003},
        {"dry-run", no_argument, 0, 1004},
        {"no-validate", no_argument, 0, 1005},
        {"no-check-deps", no_argument, 0, 1006},
        {"verbose", no_argument, 0, 1007},
        {"generate-plots", no_argument, 0, 1008},
        {"compress", no_argument, 0, 1009},
        {"no-cleanup", no_argument, 0, 1010},
        {"bimbam-gzip", no_argument, 0, 2002},
        {"gen-kinship", no_argument, 0, 2003},
        {"auto-sample", no_argument, 0, 2004},
        {"kin-method", required_argument, 0, 2005},
        {"kin-precision", required_argument, 0, 2006},
        {"out-precision", required_argument, 0, 2007},
        {"start-marker", required_argument, 0, 2008},
        {"end-marker", required_argument, 0, 2009},
        {"write-eigen", no_argument, 0, 2010},
        {"disable-gls", no_argument, 0, 2011},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };
    
    int opt;
    int option_index = 0;
    
    while ((opt = getopt_long(argc, argv, "i:p:c:k:s:o:t:n:m:h", long_options, &option_index)) != -1) {
        switch (opt) {
            case 'i': config.input_dir = optarg; break;
            case 'p': config.phenotype_file = optarg; break;
            case 'c': config.covariate_file = optarg; break;
            case 'k': config.kinship_file = optarg; break;
            case 's': config.sample_file = optarg; break;
            case 'o': config.output_dir = optarg; break;
            case 't': config.max_threads = std::stoi(optarg); break;
            case 'n': config.phenotype_column = std::stoi(optarg); break;
            case 'm': config.analysis_method = optarg; break;
            case 2001:
                if (std::string(optarg) == "gemma") {
                    config.tool = AssociationTool::GEMMA;
                } else if (std::string(optarg) == "bimbamAsso") {
                    config.tool = AssociationTool::BIMBAMASSO;
                } else {
                    std::cerr << "Error: Unknown tool '" << optarg << "'\n";
                    exit(EXIT_FAILURE);
                }
                break;
            case 1001: config.minor_allele_freq = std::stod(optarg); break;
            case 1002: config.missing_threshold = std::stod(optarg); break;
            case 1003: config.use_relatedness_matrix = false; break;
            case 1004: config.dry_run = true; break;
            case 1005: config.validate_inputs = false; break;
            case 1006: config.check_dependencies = false; break;
            case 1007: config.verbose = true; break;
            case 1008: config.generate_plots = true; break;
            case 1009: config.compress_output = true; break;
            case 1010: config.cleanup_temp_files = false; break;
            case 2002: config.bimbam_gzip = true; break;
            case 2003: config.generate_kinship = true; break;
            case 2004: break; // --auto-sample handled in prepare_bimbam_files
            case 2005: config.kinship_method = std::stoi(optarg); break;
            case 2006: config.kinship_precision = std::stoi(optarg); break;
            case 2007: config.output_precision = std::stoi(optarg); break;
            case 2008: config.start_marker = std::stoi(optarg); break;
            case 2009: config.end_marker = std::stoi(optarg); break;
            case 2010: config.write_eigen_files = true; break;
            case 2011: config.disable_gls = true; break;
            case 'h':
                display_association_help_message(argv[0]);
                exit(EXIT_SUCCESS);
            default:
                display_association_help_message(argv[0]);
                exit(EXIT_FAILURE);
        }
    }
    
    return config;
}

bool create_directory(const std::string& dir_path) {
    struct stat st;
    if (stat(dir_path.c_str(), &st) == 0) {
        return S_ISDIR(st.st_mode);
    }
    return mkdir(dir_path.c_str(), 0755) == 0;
}

std::string get_base_filename(const std::string& file_path, bool remove_extension) {
    size_t last_slash = file_path.find_last_of("/\\");
    std::string filename = (last_slash != std::string::npos) ? 
                          file_path.substr(last_slash + 1) : file_path;
    
    if (remove_extension) {
        size_t last_dot = filename.find_last_of('.');
        if (last_dot != std::string::npos) {
            filename = filename.substr(0, last_dot);
        }
    }
    
    return filename;
}

std::vector<std::string> get_files_with_extension(const std::string& dir, const std::string& extension) {
    std::vector<std::string> files;
    
    DIR* directory = opendir(dir.c_str());
    if (!directory) {
        // If dir is a file, check if it matches extension
        struct stat st;
        if (stat(dir.c_str(), &st) == 0 && S_ISREG(st.st_mode)) {
            if (dir.find(extension) != std::string::npos) {
                files.push_back(dir);
            }
        }
        return files;
    }
    
    struct dirent* entry;
    while ((entry = readdir(directory)) != nullptr) {
        std::string filename = entry->d_name;
        std::string full_path = dir + "/" + filename;
        
        struct stat file_stat;
        if (lstat(full_path.c_str(), &file_stat) == 0) {
            if ((S_ISREG(file_stat.st_mode) || S_ISLNK(file_stat.st_mode)) &&
                filename.find(extension) != std::string::npos) {
                files.push_back(full_path);
            }
        }
    }
    
    closedir(directory);
    std::sort(files.begin(), files.end());
    
    return files;
}

bool file_exists_and_readable(const std::string& filename) {
    std::ifstream file(filename);
    return file.good();
}

std::string get_timestamp_string() {
    auto now = std::chrono::system_clock::now();
    auto time_t_now = std::chrono::system_clock::to_time_t(now);
    std::ostringstream oss;
    oss << std::put_time(std::localtime(&time_t_now), "%Y%m%d_%H%M%S");
    return oss.str();
}

void setup_logging(const std::string& log_level, const std::string& log_file) {
    // Placeholder for logging setup
}

//int main(int argc, char* argv[]) {
extern "C" int assoc(int argc, char* argv[]) {
    try {
        if (argc < 2) {
            display_association_help_message(argv[0]);
            return EXIT_FAILURE;
        }
        
        AssocConfig config = parse_association_arguments(argc, argv);
        
        EnhancedAssociationAnalyzer analyzer(config);
        
        if (!analyzer.initialize()) {
            std::cerr << "Error: Failed to initialize analysis\n";
            return EXIT_FAILURE;
        }
        
        if (!analyzer.execute_analysis()) {
            std::cerr << "Error: Analysis execution failed\n";
            return EXIT_FAILURE;
        }
        
        std::cout << "Association analysis completed successfully!\n";
        std::cout << "Results written to: " << config.output_dir << "\n";
        std::cout << "Tool used: " << config.get_tool_name() << "\n";
        
        return EXIT_SUCCESS;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return EXIT_FAILURE;
    } catch (...) {
        std::cerr << "Error: Unknown exception occurred\n";
        return EXIT_FAILURE;
    }
}
