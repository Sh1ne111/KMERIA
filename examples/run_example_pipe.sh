#!/bin/bash

# This script is an example for downloading sweet potato data and simulating phenotypes for KMERIA analyses.
# Date: $(date)

set -euo pipefail  # Exit on error, undefined variables, and pipe failures

# =============================================================================
# Configuration and Default Parameters
# =============================================================================

# Default parameters - modify these according to your system and requirements
THREADS=8
MEMORY="32G"
KMER_SIZE=31
MIN_ABUNDANCE=5
MAX_ABUNDANCE=1000
BATCH_SIZE=1
PLOIDY=6
PHENO_COL=2  # Column number for phenotype values
KMC_MEMORY=32

# Software paths - adjust these to your system
KMERIA_WRAPPER="$HOME/software/KMERIA/scripts/kmeria_wrapper.pl"
SRA_TOOLKIT_PATH="/usr/local/bin"  # Adjust to your SRA toolkit installation

# Directory structure
WORK_DIR=$(pwd)
SRA_DIR="$WORK_DIR/01_sra_data"
FASTQ_DIR="$WORK_DIR/02_fastq_data"
KMERIA_RESULTS_DIR="$WORK_DIR/03_kmeria_results"
LOG_DIR="$WORK_DIR/logs"
mkdir -p $LOG_DIR

# Input files
SRA_LIST_FILE=""
PHENOTYPE_FILE=""
DEPTH_FILE=""

# =============================================================================
# Function Definitions
# =============================================================================

# Display usage information
show_usage() {
    cat << EOF
This script is an example for downloading sweet potato SRA data and simulating phenotypes for KMERIA analyses

Usage: $0 -s <sra_list_file> -p <phenotype_file> -d <depth_file> [options]

Required Arguments:
    -s, --sra-list      File containing SRA accession numbers (one per line)
                        Example format:
                        SRR28578193
                        SRR28578210
                        ...

    -p, --phenotype     File containing phenotype data (sample_name phenotype_value)
                        Example format:
                        SRR28578193 10
                        SRR28578210 20
                        ...

    -d, --depth         File containing sequencing depth information (sample_name depth)
                        Example format:
                        SRR28578193 25
                        SRR28578210 30
                        SRR28578407 28
                        ...

Optional Arguments:
    -t, --threads       Number of threads to use (default: $THREADS)
    -m, --memory        Memory allocation (default: $MEMORY)
    -k, --kmer          K-mer size (default: $KMER_SIZE)
    --batch-size        Batch size for processing (default: $BATCH_SIZE)
    --ploidy            Genome ploidy (default: $PLOIDY)
    --pheno-col         Phenotype column number (default: $PHENO_COL)
    -h, --help          Show this help message

Example Usage:
    $0 -s sra_list.txt -p phenotypes.txt -d sample_depths.txt -t 16 -m 16G

Steps:
    1. Download SRA files from NCBI
    2. Convert SRA files to FASTQ format
    3. Prepare input files for KMERIA analysis
    4. Generate KMERIA job scripts

Note: This pipeline generates job scripts that need to be submitted manually to your cluster system.

EOF
}

# Logging function with timestamps
log_message() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$LOG_DIR/pipeline.log"
}

# Error handling function
handle_error() {
    log_message "ERROR: $1"
    exit 1
}

# Check if required software is available
check_dependencies() {
    log_message "Checking software dependencies..."
    
    # Check SRA toolkit components
    local sra_tools=("prefetch" "fasterq-dump")
    for tool in "${sra_tools[@]}"; do
        if ! command -v "$SRA_TOOLKIT_PATH/$tool" &> /dev/null; then
            handle_error "SRA toolkit ($tool) not found at $SRA_TOOLKIT_PATH. Please install SRA toolkit and update SRA_TOOLKIT_PATH."
        fi
    done
    
    # Check KMERIA wrapper
    if [[ ! -f "$KMERIA_WRAPPER" ]]; then
        handle_error "KMERIA wrapper not found at $KMERIA_WRAPPER. Please check the path."
    fi
    
    # Check essential tools
    local essential_tools=("pigz" "perl")
    for tool in "${essential_tools[@]}"; do
        if ! command -v "$tool" &> /dev/null; then
            handle_error "$tool not found. Please install $tool."
        fi
    done
    
    log_message "All dependencies check passed!"
}

# Validate input files format and content
validate_input_files() {
    log_message "Validating input files..."
    
    # Validate SRA list file
    if [[ ! -s "$SRA_LIST_FILE" ]]; then
        handle_error "SRA list file is empty or does not exist: $SRA_LIST_FILE"
    fi
    
    # Check SRA ID format (basic validation)
    while IFS= read -r sra_id || [[ -n "$sra_id" ]]; do
        [[ -z "$sra_id" || "$sra_id" =~ ^[[:space:]]*# ]] && continue
        if [[ ! "$sra_id" =~ ^[SED]RR[0-9]+$ ]]; then
            log_message "WARNING: SRA ID '$sra_id' may not be in correct format (expected: SRR/ERR/DRR followed by numbers)"
        fi
    done < "$SRA_LIST_FILE"
    
    # Validate phenotype file
    if [[ ! -s "$PHENOTYPE_FILE" ]]; then
        handle_error "Phenotype file is empty or does not exist: $PHENOTYPE_FILE"
    fi
    
    local pheno_samples=()
    while IFS= read -r line || [[ -n "$line" ]]; do
        [[ -z "$line" || "$line" =~ ^[[:space:]]*# ]] && continue
        read -r sample_name pheno_value <<< "$line"
        if [[ -z "$sample_name" || -z "$pheno_value" ]]; then
            handle_error "Invalid line in phenotype file: '$line'. Expected format: sample_name phenotype_value"
        fi
        pheno_samples+=("$sample_name")
    done < "$PHENOTYPE_FILE"
    
    # Validate depth file
    if [[ ! -s "$DEPTH_FILE" ]]; then
        handle_error "Depth file is empty or does not exist: $DEPTH_FILE"
    fi
    
    local depth_samples=()
    while IFS= read -r line || [[ -n "$line" ]]; do
        [[ -z "$line" || "$line" =~ ^[[:space:]]*# ]] && continue
        read -r sample_name depth_value <<< "$line"
        if [[ -z "$sample_name" || -z "$depth_value" ]]; then
            handle_error "Invalid line in depth file: '$line'. Expected format: sample_name depth_value"
        fi
        if ! [[ "$depth_value" =~ ^[0-9]+([.][0-9]+)?$ ]]; then
            handle_error "Invalid depth value '$depth_value' for sample '$sample_name'. Depth must be a number."
        fi
        depth_samples+=("$sample_name")
    done < "$DEPTH_FILE"
    
    # Cross-validate that samples exist in both phenotype and depth files
    log_message "Cross-validating sample consistency between files..."
    local missing_in_depth=0
    for pheno_sample in "${pheno_samples[@]}"; do
        if [[ ! " ${depth_samples[*]} " =~ " ${pheno_sample} " ]]; then
            log_message "WARNING: Sample '$pheno_sample' found in phenotype file but missing in depth file"
            missing_in_depth=$((missing_in_depth + 1))
        fi
    done
    
    if [[ $missing_in_depth -gt 0 ]]; then
        log_message "WARNING: $missing_in_depth sample(s) have phenotype data but no depth information"
    fi
    
    log_message "Input file validation completed!"
}

# Create directory structure
setup_directories() {
    log_message "Setting up directory structure..."
    
    # Create all necessary directories
    local directories=("$SRA_DIR" "$FASTQ_DIR" "$KMERIA_RESULTS_DIR" "$LOG_DIR")
    for dir in "${directories[@]}"; do
        mkdir -p "$dir"
        log_message "Created directory: $dir"
    done
    
    log_message "Directory structure created successfully!"
}

# Download SRA data with improved error handling and progress tracking
download_sra_data() {
    local sra_file="$1"
    log_message "Starting SRA data download from NCBI..."
    
    local total_samples=0
    local downloaded_samples=0
    
    # Count total samples first
    while IFS= read -r sra_id || [[ -n "$sra_id" ]]; do
        [[ -z "$sra_id" || "$sra_id" =~ ^[[:space:]]*# ]] && continue
        total_samples=$((total_samples + 1))
    done < "$sra_file"
    
    log_message "Total samples to download: $total_samples"
    
    # Download each SRA file
    while IFS= read -r sra_id || [[ -n "$sra_id" ]]; do
        # Skip empty lines and comments
        [[ -z "$sra_id" || "$sra_id" =~ ^[[:space:]]*# ]] && continue
        
        downloaded_samples=$((downloaded_samples + 1))
        log_message "[$downloaded_samples/$total_samples] Downloading SRA data for: $sra_id"
        
        # Check if SRA file already exists
        if [[ -f "$SRA_DIR/$sra_id/$sra_id.sra" ]]; then
            log_message "SRA file for $sra_id already exists, skipping download"
            continue
        fi
        
        # Download SRA file using prefetch with retry mechanism
        local retry_count=0
        local max_retries=3
        local success=false
        
        while [[ $retry_count -lt $max_retries ]]; do
            if "$SRA_TOOLKIT_PATH/prefetch" --output-directory "$SRA_DIR" --progress "$sra_id" 2>/dev/null; then
                success=true
                break
            else
                retry_count=$((retry_count + 1))
                log_message "Download attempt $retry_count failed for $sra_id, retrying..."
                sleep 5  # Wait 5 seconds before retry
            fi
        done
        
        if [[ "$success" != true ]]; then
            handle_error "Failed to download SRA data for $sra_id after $max_retries attempts"
        fi
        
        log_message "Successfully downloaded: $sra_id"
    done < "$sra_file"
    
    log_message "SRA data download completed successfully!"
}

# Convert SRA to FASTQ format with enhanced processing
convert_sra_to_fastq() {
    log_message "Converting SRA files to FASTQ format..."
    
    local total_conversions=0
    local completed_conversions=0
    
    # Count total conversions needed
    for sra_path in "$SRA_DIR"/*/*.sra; do
        [[ ! -f "$sra_path" ]] && continue
        total_conversions=$((total_conversions + 1))
    done
    
    if [[ $total_conversions -eq 0 ]]; then
        handle_error "No SRA files found in $SRA_DIR. Please check if download step completed successfully."
    fi
    
    log_message "Total SRA files to convert: $total_conversions"
    
    # Process each SRA file
    for sra_path in "$SRA_DIR"/*/*.sra; do
        [[ ! -f "$sra_path" ]] && continue
        
        # Extract SRA ID from path
        local sra_id=$(basename "$(dirname "$sra_path")")
        completed_conversions=$((completed_conversions + 1))
        log_message "[$completed_conversions/$total_conversions] Converting $sra_id to FASTQ format..."
        
        # Check if FASTQ files already exist
        if ls "$FASTQ_DIR"/${sra_id}*.fastq.gz >/dev/null 2>&1; then
            log_message "FASTQ files for $sra_id already exist, skipping conversion"
            continue
        fi
        
        # Convert SRA to FASTQ using fasterq-dump
        # The --split-files option handles both single-end and paired-end data appropriately
        if ! "$SRA_TOOLKIT_PATH/fasterq-dump" \
            --outdir "$FASTQ_DIR" \
            --threads "$THREADS" \
            --split-files \
            --progress \
            --skip-technical \
            "$sra_path"; then
            handle_error "Failed to convert $sra_id to FASTQ format"
        fi
        
        # Compress FASTQ files to save disk space
        log_message "Compressing FASTQ files for $sra_id..."
        for fastq_file in "$FASTQ_DIR"/${sra_id}*.fastq; do
            if [[ -f "$fastq_file" ]]; then
                pigz -p "$THREADS" "$fastq_file"
                log_message "Compressed: $(basename "$fastq_file")"
            fi
        done
        
        log_message "Successfully processed: $sra_id"
    done
    
    log_message "SRA to FASTQ conversion completed successfully!"
}

# Generate required files for KMERIA analysis
prepare_kmeria_files() {
    log_message "Preparing input files for KMERIA analysis..."
    
    # Generate sample list file from FASTQ files
    log_message "Generating sample list from FASTQ files..."
    find "$FASTQ_DIR" -name "*.fastq.gz" -type f | \
        sed 's/.*\///' | \
        sed 's/_[12]\.fastq\.gz$//' | \
        sed 's/\.fastq\.gz$//' | \
        sort -u > "$KMERIA_RESULTS_DIR/sample.list"
    
    local sample_count=$(wc -l < "$KMERIA_RESULTS_DIR/sample.list")
    log_message "Generated sample list with $sample_count samples"
    
    # Copy and format the depth file for KMERIA
    log_message "Processing depth file for KMERIA..."
    
    # Create header for depth file if it doesn't exist
    if ! head -n1 "$DEPTH_FILE" | grep -q "sample.*depth" 2>/dev/null; then
        # Add header if the depth file doesn't have one
        echo -e "sample\tdepth" > "$KMERIA_RESULTS_DIR/sample.depth"
        # Convert space-separated to tab-separated format
        awk '{print $1 "\t" $2}' "$DEPTH_FILE" >> "$KMERIA_RESULTS_DIR/sample.depth"
    else
        # File already has header, just copy it
        cp "$DEPTH_FILE" "$KMERIA_RESULTS_DIR/sample.depth"
    fi
    
    # Process phenotype file for KMERIA format
    log_message "Processing phenotype file for KMERIA..."
    
    # KMERIA expects phenotype file without header
    # Convert the input format to the required format
    awk -v col="$PHENO_COL" '{print $1, $col}' "$PHENOTYPE_FILE" > "$KMERIA_RESULTS_DIR/pheno_imputed_noheader.txt"
    
    local pheno_count=$(wc -l < "$KMERIA_RESULTS_DIR/pheno_imputed_noheader.txt")
    log_message "Processed phenotype data for $pheno_count samples"
    
    # Validate that we have matching samples across all files
    log_message "Validating sample consistency across input files..."
    
    local samples_with_fastq=$(wc -l < "$KMERIA_RESULTS_DIR/sample.list")
    local samples_with_pheno=$(wc -l < "$KMERIA_RESULTS_DIR/pheno_imputed_noheader.txt")
    local samples_with_depth=$(tail -n +2 "$KMERIA_RESULTS_DIR/sample.depth" | wc -l)
    
    log_message "Sample counts - FASTQ: $samples_with_fastq, Phenotype: $samples_with_pheno, Depth: $samples_with_depth"
    
    if [[ $samples_with_fastq -eq 0 ]]; then
        handle_error "No samples found in FASTQ directory. Check if SRA conversion completed successfully."
    fi
    
    log_message "KMERIA input files prepared successfully!"
}

# Run KMERIA analysis pipeline
run_kmeria_analysis() {
    log_message "Generating KMERIA k-mer analysis job scripts..."
    
    # Change to the KMERIA results directory for output organization
    cd "$KMERIA_RESULTS_DIR"
    
    # Prepare the KMERIA command with all optimized parameters
    local kmeria_cmd="perl \"$KMERIA_WRAPPER\" --step all \
        --input \"$FASTQ_DIR\" \
        --output . \
        --threads $THREADS \
        --memory \"$MEMORY\" \
        --kmer $KMER_SIZE \
        --use-kmc \
        --kmc-memory $KMC_MEMORY \
        --use-kctm2 \
        --min-abund $MIN_ABUNDANCE \
        --max-abund $MAX_ABUNDANCE \
        --batch_size $BATCH_SIZE \
        --ploidy $PLOIDY \
        --depth-file sample.depth \
        --pheno pheno_imputed_noheader.txt \
        --pheno-col 2 \
        --samples sample.list"
    
    # Execute KMERIA wrapper to generate job scripts
    log_message "Executing KMERIA wrapper to generate job scripts..."
    
    if eval "$kmeria_cmd" 2>&1 | tee -a "$LOG_DIR/kmeria_analysis.log"; then
        log_message "KMERIA job scripts generated successfully!"
    else
        handle_error "Failed to generate KMERIA job scripts. Check the log for details."
    fi
    
    # Return to original directory
    cd "$WORK_DIR"
    
    # List generated job scripts for user reference
    log_message "Generated job scripts:"
    find "$KMERIA_RESULTS_DIR" -name "*.sh" -type f | while read -r script; do
        log_message "  - $(realpath "$script")"
    done
    
    log_message "KMERIA job script generation completed!"
    log_message "Please submit the generated job scripts to your cluster system manually."
}

# Generate comprehensive summary report
generate_summary() {
    log_message "Generating comprehensive analysis summary..."
    
    local summary_file="$WORK_DIR/analysis_summary.txt"
    local fastq_files_count=$(find "$FASTQ_DIR" -name "*.fastq.gz" | wc -l)
    local sample_count=$(wc -l < "$KMERIA_RESULTS_DIR/sample.list" 2>/dev/null || echo "0")
    local job_scripts_count=$(find "$KMERIA_RESULTS_DIR" -name "*.sh" -type f 2>/dev/null | wc -l)
    
    cat > "$summary_file" << EOF
=================================================================
KMERIA SRA Analysis Pipeline - Comprehensive Summary Report
=================================================================
Generated on: $(date)
Working Directory: $WORK_DIR

INPUT PARAMETERS
================
SRA List File: $SRA_LIST_FILE
Phenotype File: $PHENOTYPE_FILE  
Depth File: $DEPTH_FILE
Processing Threads: $THREADS
Memory Allocation: $MEMORY
K-mer Size: $KMER_SIZE
Minimum K-mer Abundance: $MIN_ABUNDANCE  
Maximum K-mer Abundance: $MAX_ABUNDANCE
Batch Size: $BATCH_SIZE
Genome Ploidy: $PLOIDY
Phenotype Column: $PHENO_COL
KMC Memory: $KMC_MEMORY GB

DIRECTORY STRUCTURE
===================
Working Directory: $WORK_DIR
├── 01_sra_data/          # Downloaded SRA files
├── 02_fastq_data/        # Converted FASTQ files  
├── 03_kmeria_results/    # KMERIA analysis files and job scripts
└── logs/                 # Pipeline execution logs

PROCESSING RESULTS
==================
Total SRA Files Downloaded: $(find "$SRA_DIR" -name "*.sra" 2>/dev/null | wc -l)
Total FASTQ Files Generated: $fastq_files_count
Unique Samples Processed: $sample_count
KMERIA Job Scripts Generated: $job_scripts_count

GENERATED FILES FOR ANALYSIS
=============================
- sample.list: Complete list of processed samples
- sample.depth: Sequencing depth information (from input file)  
- pheno_imputed_noheader.txt: Phenotype data in KMERIA format

NEXT STEPS
==========
1. Review the generated job scripts in: $KMERIA_RESULTS_DIR
2. Submit job scripts to your cluster system in the following order:
   a) K-mer counting jobs (count_batch_*.sh)
   b) K-mer matrix building (kctm_job.sh)  
   c) Matrix filtering (filter_job.sh)
   d) Format conversion (convert_job.sh)
   e) Association analysis (asso_job.sh)

3. Monitor job progress through your cluster's job management system
4. Results will be organized in subdirectories within: $KMERIA_RESULTS_DIR

TROUBLESHOOTING
===============
- Pipeline logs: $LOG_DIR/pipeline.log
- KMERIA analysis logs: $LOG_DIR/kmeria_analysis.log  
- Individual job logs will be created when jobs are executed

For questions about KMERIA parameters or results interpretation,
refer to the KMERIA documentation and user manual.

=================================================================
Pipeline completed successfully! 
Review this summary and submit job scripts to continue analysis.
=================================================================
EOF

    log_message "Comprehensive summary report generated: $summary_file"
    
    # Display key information to user
    echo
    echo "=========================================="
    echo "KMERIA Pipeline Preparation Complete!"
    echo "=========================================="
    echo "Samples processed: $sample_count"
    echo "Job scripts generated: $job_scripts_count"
    echo "Summary report: $summary_file"
    echo "=========================================="
}

# =============================================================================
# Main Pipeline Execution
# =============================================================================

# Parse command line arguments with enhanced validation
parse_arguments() {
    while [[ $# -gt 0 ]]; do
        case $1 in
            -s|--sra-list)
                SRA_LIST_FILE="$2"
                shift 2
                ;;
            -p|--phenotype)
                PHENOTYPE_FILE="$2"  
                shift 2
                ;;
            -d|--depth)
                DEPTH_FILE="$2"
                shift 2
                ;;
            -t|--threads)
                if [[ "$2" =~ ^[1-9][0-9]*$ ]]; then
                    THREADS="$2"
                else
                    handle_error "Invalid thread count: $2. Must be a positive integer."
                fi
                shift 2
                ;;
            -m|--memory)
                MEMORY="$2"
                shift 2
                ;;
            -k|--kmer)
                if [[ "$2" =~ ^[1-9][0-9]*$ ]] && [[ "$2" -ge 15 ]] && [[ "$2" -le 127 ]]; then
                    KMER_SIZE="$2"
                else
                    handle_error "Invalid k-mer size: $2. Must be between 15 and 127."
                fi
                shift 2
                ;;
            --min-abund)
                if [[ "$2" =~ ^[0-9]+$ ]]; then
                    MIN_ABUNDANCE="$2"
                else
                    handle_error "Invalid minimum abundance: $2. Must be a non-negative integer."
                fi
                shift 2
                ;;
            --max-abund)
                if [[ "$2" =~ ^[1-9][0-9]*$ ]]; then
                    MAX_ABUNDANCE="$2"
                else
                    handle_error "Invalid maximum abundance: $2. Must be a positive integer."
                fi
                shift 2
                ;;
            --batch-size)
                if [[ "$2" =~ ^[1-9][0-9]*$ ]]; then
                    BATCH_SIZE="$2"
                else
                    handle_error "Invalid batch size: $2. Must be a positive integer."
                fi
                shift 2
                ;;
            --ploidy)
                if [[ "$2" =~ ^[1-9][0-9]*$ ]]; then
                    PLOIDY="$2"
                else
                    handle_error "Invalid ploidy: $2. Must be a positive integer."
                fi
                shift 2
                ;;
            --pheno-col)
                if [[ "$2" =~ ^[1-9][0-9]*$ ]]; then
                    PHENO_COL="$2"
                else
                    handle_error "Invalid phenotype column: $2. Must be a positive integer."
                fi
                shift 2
                ;;
            --kmc-memory)
                if [[ "$2" =~ ^[1-9][0-9]*$ ]]; then
                    KMC_MEMORY="$2"
                else
                    handle_error "Invalid KMC memory: $2. Must be a positive integer."
                fi
                shift 2
                ;;
            -h|--help)
                show_usage
                exit 0
                ;;
            *)
                echo "Unknown option: $1"
                show_usage
                exit 1
                ;;
        esac
    done
}

# Main function orchestrating the entire pipeline
main() {
    # Parse command line arguments first
    parse_arguments "$@"
    
    # Validate that all required arguments are provided
    if [[ -z "$SRA_LIST_FILE" || -z "$PHENOTYPE_FILE" || -z "$DEPTH_FILE" ]]; then
        echo "Error: SRA list file (-s), phenotype file (-p), and depth file (-d) are all required!"
        echo
        show_usage
        exit 1
    fi
    
    # Check if input files exist and convert to absolute paths
    for file_var in "SRA_LIST_FILE" "PHENOTYPE_FILE" "DEPTH_FILE"; do
        local file_path="${!file_var}"
        if [[ ! -f "$file_path" ]]; then
            handle_error "${file_var//_/ } not found: $file_path"
        fi
        # Convert to absolute path for consistency
        declare -g "$file_var=$(realpath "$file_path")"
    done
    
    # Display pipeline startup information
    echo "=============================================================="
    echo "KMERIA SRA Analysis Pipeline (Simplified Version)"  
    echo "=============================================================="
    echo "SRA List File: $SRA_LIST_FILE"
    echo "Phenotype File: $PHENOTYPE_FILE"
    echo "Depth File: $DEPTH_FILE"
    echo "Working Directory: $WORK_DIR"
    echo "Processing Threads: $THREADS"
    echo "Memory Allocation: $MEMORY"
    echo "K-mer Size: $KMER_SIZE"
    echo "Genome Ploidy: $PLOIDY"
    echo "=============================================================="
    echo "Pipeline Steps:"
    echo "  1. Setup directories and validate inputs"
    echo "  2. Download SRA data from NCBI"
    echo "  3. Convert SRA files to FASTQ format"  
    echo "  4. Prepare KMERIA input files"
    echo "  5. Generate KMERIA analysis job scripts"
    echo "=============================================================="
    echo
    
    # Execute the complete pipeline
    setup_directories
    check_dependencies
    validate_input_files
    download_sra_data "$SRA_LIST_FILE"
    convert_sra_to_fastq
    prepare_kmeria_files
    run_kmeria_analysis
    generate_summary
    
    # Final success message
    log_message "Complete pipeline execution finished successfully!"
    echo
    echo "=============================================================="
    echo "Pipeline Preparation Complete!"
    echo "=============================================================="
    echo "✓ SRA data downloaded and converted to FASTQ"
    echo "✓ KMERIA input files prepared"  
    echo "✓ Job scripts generated and ready for submission"
    echo
    echo "Next step: Submit the generated job scripts to your cluster"
    echo "Check the summary report for detailed information:"
    echo "$WORK_DIR/analysis_summary.txt"
    echo "=============================================================="
}

# Execute the main function with all provided arguments
main "$@"
