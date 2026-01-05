#!/usr/bin/env perl
#
# Author: Chen Shuai (chensss1209@gmail.com)
#
# DEBUG (2025-11-11): Fixing the covariate format for the bimbamAsso analysis .
# DEBUG (2025-11-09): Removed the temporary BIMBAM file after the PCA and kinship calculations were completed.
#

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use File::Path qw(make_path);
use Cwd qw(abs_path);
use threads;
use Pod::Usage;

our $VERSION = "2.0.1";
our $PROGRAM = "kmeria_wrapper.pl";


my $threads = 16;         
my $memory = "32G";       # Default memory per job
my $kmer_size = 31;       
my $min_abundance = 4;    
my $max_abundance = 1000; # Default maximum k-mer abundance
my $JOB_SCHEDULER = "local"; # Default scheduler (local, slurm, sge, pbs)
my $queue = "share";     # Default queue for job submission
my $WALLTIME = "720:00:00"; # Default time limit
my $batch_size = 4;      
my $input_dir = ""; 
my $output_dir = "";      # Default output directory
my $help = 0;
my $step = "";            # Step to execute
my $sample_list = "";     # List of samples
my $kctm_missing = 0.6;   
my $ploidy = 4;           # Default ploidy
my $depth_file = "";      # Sample depth file
my $covar_file = "";      # Covariate file for association testing
my $kinship_file = "";    # Kinship matrix file
my $pheno_file = "";      # Phenotype file
my $pheno_col = 1;        # Phenotype column
my $kmc_memory = "16";    
my $kctm_batch = 10000;   # Batch size for kctm (new parameter)
my $bgzf_threads = 16;    # Threads for bgzf compression in m2b
my $sketch_size = 8000000; 
my $use_bimbam_tools = 0; # Use bimbamKin and bimbamAsso instead of kmeria asso
my $kinship_precision = 10; # Precision for kinship matrix (default: 10)
my $output_precision = 5;   # Precision for association output (default: 5)


my $count_dir;
my $kctm_dir;
my $filter_dir;
my $bimbam_dir;
my $asso_dir;

my $use_kmc;
my $count_separate_strands;
my $text_output; 
my $partition_bits;
my $compress_homopolymers;


GetOptions(
    "threads|t=i"       => \$threads,
    "memory|m=s"        => \$memory,
    "kmer|k=i"          => \$kmer_size,
    "min-abund|c=i"     => \$min_abundance,
    "max-abund|C=i"     => \$max_abundance,
    "scheduler|s=s"     => \$JOB_SCHEDULER,
    "queue|q=s"         => \$queue,
    "time=s"            => \$WALLTIME,
    "batch-size|b=i"    => \$batch_size,
    "output|o=s"        => \$output_dir,
    "help|h"            => \$help,
    "step=s"            => \$step,
    "samples=s"         => \$sample_list,
    "missing=f"         => \$kctm_missing,
    "ploidy|p=i"        => \$ploidy,
    "depth-file|d=s"    => \$depth_file,
    "input|i=s"         => \$input_dir,
    "covar=s"           => \$covar_file,
    "kinship=s"         => \$kinship_file,
    "pheno=s"           => \$pheno_file,
    "pheno-col|n=i"     => \$pheno_col,
    "kmc-memory=i"      => \$kmc_memory,
    "kctm-batch=i"      => \$kctm_batch,
    "bgzf-threads=i"    => \$bgzf_threads,
    "sketch-size=i"     => \$sketch_size,
    "use-bimbam-tools"  => \$use_bimbam_tools,
    "kinship-precision=i" => \$kinship_precision,
    "output-precision=i"  => \$output_precision,
    "use-kmc"           => \$use_kmc,
    "count-separate-strands" => \$count_separate_strands,
    "text-output|T"     => \$text_output,
    "partition-bits=i"  => \$partition_bits,
    "compress-homopolymers" => \$compress_homopolymers,
) or pod2usage(2);


pod2usage(-verbose => 2) if $help;


unless ($step) {
    print "Error: --step is required.\n";
    pod2usage(1);
}


if ($output_dir) {
    make_path($output_dir) unless -d $output_dir;
    $output_dir = abs_path($output_dir);
    
    # Initialize standard directory paths based on output directory
    initialize_directories($output_dir);
}


sub initialize_directories {
    my ($base_dir) = @_;
    
    $count_dir = "$base_dir/01_kmer_counts";
    $kctm_dir = "$base_dir/02_kmer_matrices";
    $filter_dir = "$base_dir/03_filtered_matrices";
    $bimbam_dir = "$base_dir/04_bimbam";
    $asso_dir = "$base_dir/05_association";
    
    return ($count_dir, $kctm_dir, $filter_dir, $bimbam_dir, $asso_dir);
}

# ===== NEW FUNCTIONS FOR SAMPLE ORDER MANAGEMENT =====

sub read_sample_list {
    my ($sample_file) = @_;
    my @samples;
    
    if ($sample_file && -f $sample_file) {
        open(my $fh, '<', $sample_file) or die "Cannot open sample list $sample_file: $!";
        while (my $line = <$fh>) {
            chomp $line;
            $line =~ s/^\s+|\s+$//g;  # Remove leading/trailing whitespace
            next if $line =~ /^\s*$/;  # Skip empty lines
            next if $line =~ /^#/;     # Skip comment lines
            push @samples, $line;
        }
        close($fh);
    }
    
    return @samples;
}

sub save_sample_order {
    my ($output_dir, @samples) = @_;
    
    my $order_file = "$output_dir/sample_order.txt";
    open(my $fh, '>', $order_file) or die "Cannot create sample order file $order_file: $!";
    foreach my $sample (@samples) {
        print $fh "$sample\n";
    }
    close($fh);
    
    print "Sample order saved to: $order_file\n";
    print "Total samples: " . scalar(@samples) . "\n";
    
    return $order_file;
}

sub verify_sample_order {
    my ($dir1, $file1, $dir2, $file2) = @_;
    
    my @samples1 = read_sample_list("$dir1/$file1");
    my @samples2 = read_sample_list("$dir2/$file2");
    
    if (scalar(@samples1) != scalar(@samples2)) {
        warn "WARNING: Sample count mismatch between $file1 and $file2\n";
        return 0;
    }
    
    for (my $i = 0; $i < scalar(@samples1); $i++) {
        if ($samples1[$i] ne $samples2[$i]) {
            warn "WARNING: Sample order mismatch at position $i: $samples1[$i] vs $samples2[$i]\n";
            return 0;
        }
    }
    
    return 1;
}


if ($step eq "count") {
    run_count_step();
} elsif ($step eq "kctm") {
    run_kctm_step();
} elsif ($step eq "filter" || $step eq "flt") {
    run_filter_step();
} elsif ($step eq "m2b") {
    run_convert_step();
} elsif ($step eq "asso") {
    run_asso_step();
} elsif ($step eq "all") {
    run_full_pipeline();
} else {
    print "Error: Unknown step '$step'. Valid steps are: count, kctm, filter, m2b, asso, all\n";
    exit 1;
}


sub run_full_pipeline {
    
    if (!$input_dir || !$output_dir) {
        die "Error: Input directory (-i) and output directory (-o) are required for pipeline command\n";
    }
    
    if (!$depth_file || !-f $depth_file) {
        die "Error: Depth file (-d) is required and must exist for the complete pipeline\n";
    }
    
  
    make_path($count_dir, $kctm_dir, $filter_dir, $bimbam_dir, $asso_dir);
    
  
    print "Generating scripts for KMERIA pipeline...\n";
    
    # Step 1: Count k-mers
    print "Step 1: Generating scripts for counting k-mers...\n";
    run_count_step($input_dir, $count_dir);
    
    # Step 2: Build k-mer matrix
    print "Step 2: Generating script for building k-mer matrix...\n";
    run_kctm_step($count_dir, $kctm_dir);
    
    # Step 3: Filter k-mer matrix
    print "Step 3: Generating script for filtering k-mer matrix...\n";
    run_filter_step($kctm_dir, $filter_dir);
    
    # Step 4: Convert to BIMBAM format
    print "Step 4: Generating script for converting to BIMBAM format...\n";
    run_convert_step($filter_dir, $bimbam_dir);
    
    # Step 5: Association analysis (if phenotype data is available)
    if ($pheno_file) {
        print "Step 5: Generating script for association analysis...\n";
        run_asso_step($bimbam_dir, $asso_dir);
    } else {
        print "Step 5: Skipping association analysis script generation (no phenotype file provided)\n";
    }
    
    print "\n=== Pipeline Script Generation Complete ===\n";
    print "All scripts have been generated successfully!\n";
    print "Please submit the job scripts manually in the following order:\n";
    print "  1. $count_dir/count_batch_*.sh\n";
    print "  2. $kctm_dir/kctm_job.sh\n";
    print "  3. $filter_dir/filter_job.sh\n";
    print "  4. $bimbam_dir/convert_job.sh\n";
    if ($pheno_file) {
        print "  5. $asso_dir/asso_job.sh\n";
    }
    print "===========================================\n";
}


sub create_job_header {
    my ($fh, $job_name, $log_prefix) = @_;
    
    if ($JOB_SCHEDULER eq "slurm") {
        print $fh "#!/bin/bash\n";
        print $fh "#SBATCH --job-name=$job_name\n";
        print $fh "#SBATCH --output=${log_prefix}.log\n";
        print $fh "#SBATCH --error=${log_prefix}.err\n";
        print $fh "#SBATCH --ntasks=1\n";
        print $fh "#SBATCH --cpus-per-task=$threads\n";
        #print $fh "#SBATCH --mem=$memory\n";
        print $fh "#SBATCH --time=$WALLTIME\n";
        print $fh "#SBATCH --partition=$queue\n\n";
    } elsif ($JOB_SCHEDULER eq "sge") {
        print $fh "#!/bin/bash\n";
        print $fh "#\$ -N $job_name\n";
        print $fh "#\$ -o ${log_prefix}.log\n";
        print $fh "#\$ -e ${log_prefix}.err\n";
        print $fh "#\$ -pe smp $threads\n";
        print $fh "#\$ -l h_vmem=$memory\n";
        print $fh "#\$ -l h_rt=$WALLTIME\n";
        print $fh "#\$ -q $queue\n\n";
    } elsif ($JOB_SCHEDULER eq "pbs") {
        print $fh "#!/bin/bash\n";
        print $fh "#PBS -N $job_name\n";
        print $fh "#PBS -o ${log_prefix}.log\n";
        print $fh "#PBS -e ${log_prefix}.err\n";
        print $fh "#PBS -l select=1:ncpus=$threads:mem=$memory\n";
        print $fh "#PBS -l walltime=$WALLTIME\n";
        print $fh "#PBS -q $queue\n\n";
    } else {
        # Local execution, just add a basic header
        print $fh "#!/bin/bash\n";
        print $fh "# Local execution job: $job_name\n";
        print $fh "# Log files: ${log_prefix}.log, ${log_prefix}.err\n\n";
    }
    
    print $fh "set -e  # Exit on error\n";
    print $fh "set -u  # Exit on undefined variable\n";
    print $fh "set -o pipefail  # Exit on pipe failure\n\n";
}


sub run_count_step {
    my ($input, $output) = @_;
    
    
    $input = $input_dir if !$input;
    
    if (!$output) {
        $output = $count_dir;
    }
    
    if (!$input || !$output) {
        die "Error: Input directory (-i) and output directory (-o) are required for count command\n";
    }
    
    # Create output directory if it doesn't exist
    make_path($output) unless -d $output;
    $output = abs_path($output);
    
    
    my @samples;
    if ($sample_list && -f $sample_list) {
        @samples = read_sample_list($sample_list);
    } elsif (-d $input) {
        $input = abs_path($input);
        # If input is a directory, find all FASTQ files
        opendir(DIR, $input) or die "Cannot open directory $input: $!";
        my %sample_files;
        while (my $file = readdir(DIR)) {
            if ($file =~ /^(.*?)(?:_R?[12])?\.(?:fastq|fq|fasta|fa)(?:\.gz)?$/i) {
                my $sample = $1;
                $sample_files{$sample} = 1;
            }
        }
        closedir(DIR);
        @samples = sort keys %sample_files;
    } else {
        die "Error: Input ($input) is neither a directory nor a file\n";
    }
    
    if (scalar(@samples) == 0) {
        die "Error: No samples found to process\n";
    }
    
    # ===== SAVE SAMPLE ORDER =====
    save_sample_order($output, @samples);
    
    # Process samples in batches
    my $total_samples = scalar(@samples);
    my $num_batches = int(($total_samples + $batch_size - 1) / $batch_size);
    
    my $count_method = $use_kmc ? "KMC" : "kmeria count";
    print "Found $total_samples samples, generating $num_batches batch scripts using $count_method\n";
    
    # Create a job script for each batch
    for (my $batch = 0; $batch < $num_batches; $batch++) {
        my $start_idx = $batch * $batch_size;
        my $end_idx = ($start_idx + $batch_size - 1) > $#samples ? $#samples : ($start_idx + $batch_size - 1);
        
        # Get samples for this batch
        my @batch_samples = @samples[$start_idx..$end_idx];
        
        # Create job script
        my $job_script = "$output/count_batch_${batch}.sh";
        open(my $fh, '>', $job_script) or die "Cannot create job script $job_script: $!";
        
        # Write job header
        my $job_name = $use_kmc ? "kmc_count_$batch" : "kmeria_count_$batch";
        create_job_header($fh, $job_name, "$output/count_batch_$batch");
        
        print $fh "echo '========================================'\n";
        print $fh "echo 'K-mer Counting Batch $batch ($count_method)'\n";
        print $fh "echo 'Samples: " . join(", ", @batch_samples) . "'\n";
        print $fh "echo '========================================'\n\n";
        
    
        foreach my $sample (@batch_samples) {
            # Find corresponding FASTQ/FASTA files
            my @input_files;
            if (-d $input) {
                $input = abs_path($input);
                opendir(DIR, $input) or die "Cannot open directory $input: $!";
                while (my $file = readdir(DIR)) {
                    if ($file =~ /^${sample}(?:_R?[12])?\.(?:fastq|fq|fasta|fa)(?:\.gz)?$/i) {
                        push @input_files, "$input/$file";
                    }
                }
                closedir(DIR);
                @input_files = sort @input_files;  # Ensure consistent ordering
            }
            
            # Skip if no input files found
            if (scalar(@input_files) == 0) {
                print $fh "echo 'WARNING: No input files found for sample $sample, skipping'\n\n";
                next;
            }
            
            if ($use_kmc) {
                # Use KMC for k-mer counting
                print $fh "# KMC k-mer counting for sample: $sample\n";
                print $fh "mkdir -p $output/$sample\n";
                print $fh "echo 'Processing sample: $sample with KMC'\n";
                
                # Create a temporary file list for KMC input
                print $fh "cat > $output/${sample}_fastq_list.txt << 'EOF'\n";
                foreach my $input_file (@input_files) {
                    print $fh "$input_file\n";
                }
                print $fh "EOF\n\n";
                
                # Run KMC with the file list
                print $fh "kmc -k$kmer_size -t$threads -m$kmc_memory -b -ci$min_abundance -cs$max_abundance \\\n";
                print $fh "    \@$output/${sample}_fastq_list.txt $output/${sample}_k${kmer_size} $output/$sample\n\n";
                
                # Sort k-mers (required for kctm)
                print $fh "# Sort k-mers for matrix construction\n";
                print $fh "kmc_tools -t$threads transform $output/${sample}_k${kmer_size} sort $output/${sample}_sort_k${kmer_size}\n\n";
                
                # Cleanup temporary files
                print $fh "# Cleanup temporary files\n";
                print $fh "rm -f $output/${sample}_fastq_list.txt\n";
                print $fh "rm -f $output/${sample}_k${kmer_size}.kmc_pre $output/${sample}_k${kmer_size}.kmc_suf\n";
                print $fh "rm -rf $output/$sample\n";
                
            } else {
                # Use kmeria count for k-mer counting
                print $fh "# kmeria count for sample: $sample\n";
                print $fh "echo 'Processing sample: $sample with kmeria count'\n";
                
                # Build kmeria count command
                my $kmeria_cmd = "kmeria count -k $kmer_size -t $threads";
                
                # Add optional flags
                $kmeria_cmd .= " -C" if $count_separate_strands;
                $kmeria_cmd .= " -T" if $text_output;
                $kmeria_cmd .= " -p $partition_bits" if defined $partition_bits;
                $kmeria_cmd .= " -c" if $compress_homopolymers;
                
                # Determine output format extension
                my $output_ext = $text_output ? "txt" : "bin";
                
                # Process each input file
                foreach my $input_file (@input_files) {
                    my $file_basename = basename($input_file);
                    $file_basename =~ s/\.(?:fastq|fq|fasta|fa)(?:\.gz)?$//i;
                    
                    print $fh "$kmeria_cmd -o $output/${sample}_${file_basename}_k${kmer_size}.$output_ext $input_file\n";
                }
                
                # If multiple files exist (paired-end), merge them if needed
                if (scalar(@input_files) > 1) {
                    print $fh "\n# Note: Multiple input files detected for $sample\n";
                    print $fh "# Consider merging count results if necessary\n\n";
                } else {
                    print $fh "\n";
                }
                
                # For compatibility with kctm step, create a sorted version marker
                # Note: kmeria count output may need to be converted to KMC format for kctm
                print $fh "# Note: For kctm compatibility, kmeria count output may need format conversion\n";
                print $fh "# Current output: $output/${sample}_*_k${kmer_size}.$output_ext\n\n";
            }
            
            print $fh "echo 'Finished processing sample: $sample'\n";
            print $fh "echo '----------------------------------------'\n\n";
        }
        
        print $fh "echo 'Batch $batch completed successfully'\n";
        
        close($fh);
        chmod 0755, $job_script;
        
        print "Generated script: $job_script\n";
    }
    
    print "\nAll k-mer counting job scripts have been generated.\n";
    
    if (!$use_kmc) {
        print "\nIMPORTANT NOTE:\n";
        print "You are using 'kmeria count' which outputs in binary/text format.\n";
        print "The kctm step currently expects KMC sorted database format.\n";
        print "You may need to:\n";
        print "  1. Convert kmeria count output to KMC format, OR\n";
        print "  2. Use a different kctm workflow compatible with kmeria count output\n\n";
    }
    
    print "Please submit them manually to your cluster system.\n\n";
}


sub run_kctm_step {
    my ($input, $output) = @_;
    
    # Handle input/output directory
    $input = $input_dir if !$input;
    
    if (!$output) {
        $output = $kctm_dir;
    }
    
   
    if (!$input || !$output) {
        die "Error: Input directory (-i) and output directory (-o) are required for kctm command\n";
    }
    
    # Create output directory if it doesn't exist
    make_path($output) unless -d $output;
    $output = abs_path($output);
    $input = abs_path($input);
    
    # ===== CHECK FOR SAMPLE ORDER FILE =====
    my $sample_order_file = "$input/sample_order.txt";
    if (!-f $sample_order_file) {
        warn "WARNING: Sample order file not found at $sample_order_file\n";
        warn "Will use alphabetical order from ls command (may cause order mismatch)\n";
    }
    
    # Create job script
    my $job_script = "$output/kctm_job.sh";
    open(my $fh, '>', $job_script) or die "Cannot create job script $job_script: $!";
    
    # Write job header
    create_job_header($fh, "kmeria_kctm", "$output/kctm_job");
    
    print $fh "echo '========================================'\n";
    print $fh "echo 'K-mer Matrix Construction'\n";
    print $fh "echo '========================================'\n\n";
    
    # ===== CRITICAL FIX: USE SAMPLE ORDER FILE =====
    print $fh "# Create list of sorted k-mer database files in CORRECT ORDER\n";
    
    if (-f $sample_order_file) {
        print $fh "# Using sample order from: $sample_order_file\n";
        print $fh "if [ ! -f $sample_order_file ]; then\n";
        print $fh "    echo 'ERROR: Sample order file not found: $sample_order_file'\n";
        print $fh "    exit 1\n";
        print $fh "fi\n\n";
        
        print $fh "# Generate sample list based on original order\n";
        print $fh "> $output/sample_sort_k${kmer_size}.list\n";
        print $fh "while IFS= read -r sample; do\n";
        print $fh "    kmer_file=\"$input/\${sample}_sort_k${kmer_size}\"\n";
        print $fh "    if [ -f \"\${kmer_file}.kmc_pre\" ] && [ -f \"\${kmer_file}.kmc_suf\" ]; then\n";
        print $fh "        echo \"\${kmer_file}\" >> $output/sample_sort_k${kmer_size}.list\n";
        print $fh "    else\n";
        print $fh "        echo \"WARNING: K-mer files not found for sample: \$sample\"\n";
        print $fh "        echo \"Expected: \${kmer_file}.kmc_pre and \${kmer_file}.kmc_suf\"\n";
        print $fh "    fi\n";
        print $fh "done < $sample_order_file\n\n";
    } else {
        print $fh "# WARNING: No sample order file found, using alphabetical order\n";
        print $fh qq{ls $input/*_sort_k${kmer_size}.kmc_pre | sed "s/\\.kmc_pre//g" | sort > $output/sample_sort_k${kmer_size}.list\n\n};
    }
    
    print $fh "# Check if sample list is not empty\n";
    print $fh "if [ ! -s $output/sample_sort_k${kmer_size}.list ]; then\n";
    print $fh "    echo 'ERROR: No sorted k-mer files found in $input'\n";
    print $fh "    exit 1\n";
    print $fh "fi\n\n";
    
    print $fh "echo \"Found \$(wc -l < $output/sample_sort_k${kmer_size}.list) samples for matrix construction\"\n\n";
    
    # ===== PRINT SAMPLE LIST FOR VERIFICATION =====
    print $fh "# Print sample list for verification\n";
    print $fh "echo \"Sample order:\"\n";
    print $fh "cat $output/sample_sort_k${kmer_size}.list\n";
    print $fh "echo \"\"\n\n";
    
    # Run kmeria kctm
    print $fh "# Build k-mer count matrix\n";
    print $fh "kmeria kctm -i $output/sample_sort_k${kmer_size}.list \\\n";
    print $fh "    -o $output/sample_k${kmer_size} \\\n";
    print $fh "    -v --no-header -b $kctm_batch\n\n";
    
    # ===== SAVE KCTM SAMPLE ORDER =====
    print $fh "# Save sample order used in kctm for verification\n";
    print $fh "cat $output/sample_sort_k${kmer_size}.list | xargs -n1 basename | sed 's/_sort_k${kmer_size}\$//' > $output/kctm_sample_order.txt\n\n";
    
    print $fh "echo 'K-mer matrix construction completed successfully'\n";
    print $fh "echo 'Sample order saved to: $output/kctm_sample_order.txt'\n";
    
    close($fh);
    chmod 0755, $job_script;
    
    print "Generated script: $job_script\n";
    print "K-mer matrix building job script has been generated.\n";
    
    if (-f $sample_order_file) {
        print "✓ Sample order file found: $sample_order_file\n";
        print "✓ K-mer matrix will use the same sample order as count step\n";
    } else {
        print "✗ WARNING: No sample order file found!\n";
        print "  The sample order may not match the original sample list\n";
    }
    
    print "Please submit it manually to your cluster system after all counting jobs complete.\n\n";
}


sub run_filter_step {
    my ($input, $output) = @_;
    
    # Handle input/output directory
    $input = $input_dir if !$input;
    
    if (!$output) {
        $output = $filter_dir;
    }
    
    
    if (!$input || !$output) {
        die "Error: Input directory (-i) and output directory (-o) are required for filter command\n";
    }
    
   
    if (!$depth_file || !-f $depth_file) {
        die "Error: Depth file (-d) is required and must exist for filter command\n";
    }
    
    # Create output directory if it doesn't exist
    make_path($output) unless -d $output;
    $output = abs_path($output);
    $input = abs_path($input);
    $depth_file = abs_path($depth_file);
    
    # Create job script
    my $job_script = "$output/filter_job.sh";
    open(my $fh, '>', $job_script) or die "Cannot create job script $job_script: $!";
    
    # Write job header
    create_job_header($fh, "kmeria_filter", "$output/filter_job");
    
    print $fh "echo '========================================'\n";
    print $fh "echo 'K-mer Matrix Filtering'\n";
    print $fh "echo '========================================'\n\n";
    
    # ===== VERIFY SAMPLE ORDER =====
    print $fh "# Verify sample order consistency\n";
    my $count_order = abs_path("$count_dir/sample_order.txt");
    my $kctm_order = abs_path("$input/kctm_sample_order.txt");
    
    print $fh "if [ -f $count_order ] && [ -f $kctm_order ]; then\n";
    print $fh "    echo 'Checking sample order consistency...'\n";
    print $fh "    if ! diff -q $count_order $kctm_order > /dev/null 2>&1; then\n";
    print $fh "        echo 'WARNING: Sample order mismatch detected!'\n";
    print $fh "        echo 'Count order: $count_order'\n";
    print $fh "        echo 'KCTM order: $kctm_order'\n";
    print $fh "        diff $count_order $kctm_order || true\n";
    print $fh "    else\n";
    print $fh "        echo '✓ Sample order is consistent'\n";
    print $fh "    fi\n";
    print $fh "else\n";
    print $fh "    echo 'WARNING: Cannot verify sample order (order files not found)'\n";
    print $fh "fi\n";
    print $fh "echo ''\n\n";
    
    # Write actual commands
    print $fh "# Filter k-mer matrix\n";
    print $fh "kmeria filter -i $input \\\n";
    print $fh "    -o $output \\\n";
    print $fh "    -t $threads \\\n";
    print $fh "    -c $max_abundance \\\n";
    print $fh "    -p $ploidy \\\n";
    print $fh "    -s $kctm_missing \\\n";
    print $fh "    --output-format compressed \\\n";
    print $fh "    -v \\\n";
    print $fh "    -d $depth_file\n\n";
    
    print $fh "echo 'K-mer matrix filtering completed successfully'\n";
    
    close($fh);
    chmod 0755, $job_script;
    
    print "Generated script: $job_script\n";
    print "K-mer matrix filtering job script has been generated.\n";
    print "Please submit it manually to your cluster system after kctm completes.\n\n";
}


sub run_convert_step {
    my ($input, $output) = @_;
    
  
    $input = $input_dir if !$input;
    
    if (!$output) {
        $output = $bimbam_dir;
    }
    
   
    if (!$input || !$output) {
        die "Error: Input directory (-i) and output directory (-o) are required for convert command\n";
    }
    
    if (!$depth_file || !-f $depth_file) {
        die "Error: Depth file (-d) is required for generating sample list in convert step\n";
    }
    
    
    make_path($output) unless -d $output;
    $output = abs_path($output);
    $input = abs_path($input);
    $depth_file = abs_path($depth_file);
    
    
    my $job_script = "$output/convert_job.sh";
    open(my $fh, '>', $job_script) or die "Cannot create job script $job_script: $!";
    
    
    create_job_header($fh, "kmeria_m2b", "$output/convert_job");
    
    print $fh "echo '========================================'\n";
    print $fh "echo 'Converting to BIMBAM Format'\n";
    print $fh "echo '========================================'\n\n";
    
    # Convert to BIMBAM format
    print $fh "# Convert k-mer matrices to BIMBAM format\n";
    print $fh "kmeria m2b --in $input \\\n";
    print $fh "    --out $output \\\n";
    print $fh "    --threads $threads \\\n";
    print $fh "    --bgzf-threads $bgzf_threads \\\n";
    print $fh "    --verbose \\\n";
    print $fh "    --stats\n\n";
    
    
    print $fh "# Create a sketch sample for PCA and kinship calculation\n";
    print $fh "kmeria sketch $output/*.bimbam.gz -n $sketch_size > $output/sampling.bimbam\n\n";
    
    # ===== USE ORIGINAL SAMPLE ORDER =====
    print $fh "# Generate sample list from original sample order\n";
    my $original_order = abs_path("$count_dir/sample_order.txt");
    
    print $fh "if [ -f $original_order ]; then\n";
    print $fh "    echo 'Using original sample order from: $original_order'\n";
    print $fh "    cp $original_order $output/tmp_sample.list\n";
    print $fh "else\n";
    print $fh "    echo 'WARNING: Original sample order not found, using depth file'\n";
    print $fh "    cut -f1 $depth_file | sort | uniq > $output/tmp_sample.list\n";
    print $fh "fi\n\n";
    
    # Convert to VCF
    print $fh "# Convert BIMBAM to VCF format\n";
    print $fh "kmeria b2g -i $output/sampling.bimbam \\\n";
    print $fh "    -s $output/tmp_sample.list \\\n";
    print $fh "    --no-validate \\\n";
    print $fh "    -o $output/sampling.vcf\n\n";
    # remove the temp bimbam file
    # print $fh "rm $output/sampling.bimbam\n\n";
    
    # Convert to PLINK format
    print $fh "# Convert VCF to PLINK binary format\n";
    print $fh "plink --vcf $output/sampling.vcf \\\n";
    print $fh "    --make-bed \\\n";
    print $fh "    --out $output/sampling\n\n";
    
    print $fh "echo 'Conversion to BIMBAM format completed successfully'\n";
    
    close($fh);
    chmod 0755, $job_script;
    
    print "Generated script: $job_script\n";
    print "Conversion to BIMBAM format job script has been generated.\n";
    print "Please submit it manually to your cluster system after filtering completes.\n\n";
}


sub run_asso_step {
    my ($input, $output) = @_;
    
    # Handle input/output directory
    $input = $input_dir if !$input;
    
    if (!$output) {
        $output = $asso_dir;
    }
    
 
    if (!$input || !$output || !$pheno_file) {
        die "Error: Input directory (-i), output directory (-o), and phenotype file (--pheno) are required for asso command\n";
    }
    
   
    make_path($output) unless -d $output;
    $output = abs_path($output);
    $input = abs_path($input);
    $pheno_file = abs_path($pheno_file);
    
    # Create PCA and kinship directories
    make_path("$output/pca", "$output/kinship");
    
    # Create job script
    my $job_script = "$output/asso_job.sh";
    open(my $fh, '>', $job_script) or die "Cannot create job script $job_script: $!";
    
    # Write job header
    create_job_header($fh, "kmeria_asso", "$output/asso_job");
    
    print $fh "echo '========================================'\n";
    print $fh "echo 'Association Analysis'\n";
    print $fh "echo '========================================'\n\n";
    
    my $tool = $use_bimbam_tools ? "bimbamAsso" : "gemma";
    
    # Generate PCA if covariate file not provided
    if (!$covar_file || !-f $covar_file) {
        print $fh "# Generate PCA for covariates\n";
        print $fh "echo 'Calculating PCA...'\n";
        print $fh "plink --bfile $input/sampling \\\n";
        print $fh "    --allow-extra-chr \\\n";
        print $fh "    --allow-no-sex \\\n";
        print $fh "    --pca 10 \\\n";
        print $fh "    --out $output/pca/samplPCA\n\n";
        
        print $fh "# Format PCA output as covariate file\n";
        if ($use_bimbam_tools) {
            print $fh "cat $output/pca/samplPCA.eigenvec | awk '{print \$1,\$2,\"1\",\$3,\$4,\$5}' > $output/pca/PCA.txt\n\n";
        } else {
            print $fh "cat $output/pca/samplPCA.eigenvec | awk '{print \"1\",\$3,\$4,\$5}' > $output/pca/PCA.txt\n\n";
        }
        
        $covar_file = "$output/pca/PCA.txt";
    } else {
        $covar_file = abs_path($covar_file);
    }
    
    # Generate kinship matrix if not provided
    if (!$kinship_file || !-f $kinship_file) {
        if ($use_bimbam_tools) {
            print $fh "# Calculate kinship matrix using bimbamKin\n";
            print $fh "echo 'Calculating kinship matrix...'\n";
            print $fh "bimbamKin $input/sampling.bimbam \\\n";
            print $fh "    $output/kinship/sampling.BN.kinMat \\\n";
            print $fh "    -b -d $kinship_precision\n\n";
            $kinship_file = "$output/kinship/sampling.BN.kinMat";
        } else {
            print $fh "# Calculate kinship matrix using GEMMA\n";
            print $fh "echo 'Calculating kinship matrix...'\n";
            print $fh "gemma -bfile $input/sampling \\\n";
            print $fh "    -gk \\\n";
            print $fh "    -p $pheno_file \\\n";
            print $fh "    -outdir $output/kinship \\\n";
            print $fh "    -o kinship\n\n";
            $kinship_file = "$output/kinship/kinship.cXX.txt";
        }
    } else {
        $kinship_file = abs_path($kinship_file);
    }
    
    # Run association analysis with bimbamAsso
    if ($use_bimbam_tools) {
        print $fh "# Generate sample list file for bimbamAsso\n";
        print $fh "echo 'Generating sample list...'\n";
        
        # ===== USE ORIGINAL SAMPLE ORDER =====
        my $original_order = abs_path("$count_dir/sample_order.txt");
        
        print $fh "if [ -f $original_order ]; then\n";
        print $fh "    echo 'Using original sample order from: $original_order'\n";
        print $fh "    cat $original_order | awk '{print \$1,\$1}' > $output/sample.list\n";
        print $fh "else\n";
        if ($depth_file) {
            print $fh "    echo 'WARNING: Using depth file for sample order'\n";
            print $fh "    cut -f1 $depth_file | awk '{print \$1,\$1}' > $output/sample.list\n";
        }
        print $fh "fi\n\n";
    }
    
    print $fh "# Run association analysis\n";
    print $fh "echo 'Running association analysis...'\n";
    print $fh "kmeria asso --tool $tool \\\n";
    print $fh "    -i $input \\\n";
    print $fh "    -p $pheno_file \\\n";
    print $fh "    -n $pheno_col \\\n";
    print $fh "    -o $output \\\n";
    print $fh "    -t $threads \\\n";
    print $fh "    --verbose \\\n";
    print $fh "    --bimbam-gzip \\\n";
    if ($use_bimbam_tools) {
#        print $fh "    --bimbam-gzip \\\n";
        print $fh "    --out-precision $output_precision \\\n";
        print $fh "    -s $output/sample.list \\\n";
    }
    if ($covar_file) {
        print $fh "    -c $covar_file \\\n";
    }
    if ($kinship_file) {
        print $fh "    -k $kinship_file \\\n";
    }
    print $fh "    2>&1 | tee $output/asso.log\n\n";
    
    print $fh "echo 'Association analysis completed successfully'\n";
    
    close($fh);
    chmod 0755, $job_script;
    
    print "Generated script: $job_script\n";
    print "Association job script has been generated.\n";
    print "Please submit it manually to your cluster system after conversion completes.\n\n";
}

__END__

=head1 NAME

kmeria_wrapper.pl - A parallel wrapper for KMERIA pipeline (v2.0.1)

=head1 SYNOPSIS

kmeria_wrapper.pl --step <step> [options]

 Steps:
    count            Count k-mers from FASTQ files using KMC
    kctm             Build population-level k-mer matrices
    filter           Filter raw k-mer matrices
    m2b              Convert k-mer matrices to BIMBAM dosage format
    asso             Conduct k-mer association study
    all              Run the complete KMERIA pipeline

 Global Options:
   --help|-h                    Show detailed help message
   --threads|-t      [INT]      Number of threads per job [default: 16]
   --step            [STR]      Pipeline step to run (required)
   --scheduler|-s    [STR]      Job scheduler (local, slurm, sge, pbs) [default: local]
   --samples         [FILE]     File containing list of samples (one per line)
   --queue|-q        [STR]      Queue/Partition for job submission [default: share]
   --memory|-m       [STR]      Memory per job [default: 32G]
   --time            [STR]      Time limit for job [default: 720:00:00]

 Command-specific options:

 For 'count' (K-mer counting with KMC or kmeria count):
   --input|-i        [DIR]      Directory with FASTQ/FASTA files
   --output|-o       [DIR]      Output directory [default: 01_kmer_counts]
   --kmer|-k         [INT]      K-mer size [default: 31]
   --min-abund       [INT]      Minimum k-mer abundance (KMC only) [default: 4]
   --max-abund       [INT]      Maximum k-mer abundance (KMC only) [default: 1000]
   --batch-size|-b   [INT]      Number of samples per batch script [default: 4]
   
   KMC-specific options:
   --use-kmc                    Use KMC instead of kmeria count
   --kmc-memory      [INT]      Memory allocation for KMC in GB [default: 16]
   
   kmeria count-specific options (default):
   --count-separate-strands     Count strands separately (no canonical)
   -T|--text-output             Text output instead of binary
   --partition-bits  [INT]      Partitioning bits [default: 16]
   --compress-homopolymers      Compress homopolymers (experimental)

 For 'kctm' (K-mer matrix construction):
   --input|-i        [DIR]      Directory with sorted k-mer databases
   --output|-o       [DIR]      Output directory [default: 02_kmer_matrices]
   --kctm-batch      [INT]      Batch size for kctm processing [default: 10000]

 For 'filter' (Matrix filtering):
   --input|-i        [DIR]      Directory with k-mer matrices
   --output|-o       [DIR]      Output directory [default: 03_filtered_matrices]
   --max-abund       [INT]      Maximum k-mer abundance [default: 1000]
   --missing         [FLOAT]    Missing ratio threshold [default: 0.6]
   --ploidy|-p       [INT]      Genome ploidy [default: 4]
   --depth-file|-d   [FILE]     Sample depth file (REQUIRED)

 For 'm2b' (Matrix to BIMBAM conversion):
   --input|-i        [DIR]      Directory with filtered matrices
   --output|-o       [DIR]      Output directory [default: 04_bimbam]
   --bgzf-threads    [INT]      Threads for BGZF compression [default: 16]
   --sketch-size     [INT]      Number of k-mers for sampling [default: 8000000]
   --depth-file|-d   [FILE]     Sample depth file for sample list generation

 For 'asso' (Association analysis):
   --input|-i        [DIR]      Directory with BIMBAM files
   --output|-o       [DIR]      Output directory [default: 05_association]
   --pheno           [FILE]     Phenotype file (REQUIRED)
   --pheno-col|-n    [INT]      Phenotype column [default: 1]
   --covar           [FILE]     Covariate file (if not provided, PCA will be calculated)
   --kinship         [FILE]     Kinship matrix file (if not provided, will be calculated)
   --use-bimbam-tools           Use bimbamAsso mode in kassoc (default: gemma)
   --kinship-precision [INT]    Precision for kinship matrix [default: 10]
   --output-precision  [INT]    Precision for association output [default: 5]

=head1 DESCRIPTION

kmeria_wrapper.pl is a wrapper script for generating job scripts for the KMERIA v2.0 pipeline.
It generates bash scripts that can be manually submitted to cluster job schedulers (SLURM, SGE, PBS)
or executed locally.

Key changes in v2.0.1:
- **CRITICAL FIX**: Maintains consistent sample order across all pipeline steps
- Sample order is saved from count step and reused in all downstream steps
- Added verification checks to detect sample order mismatches
- Prevents genotype-phenotype association errors due to sample order inconsistency

=head1 WORKFLOW

The complete KMERIA pipeline consists of 5 steps:

1. count  - Count k-mers from FASTQ/FASTA files (kmeria count by default, or KMC with --use-kmc)
2. kctm   - Build k-mer count matrices from count results
3. filter - Filter matrices by abundance, ploidy, and missing rate
4. m2b    - Convert matrices to BIMBAM format and generate VCF/PLINK files
5. asso   - Perform association analysis with kinship correction using kassoc

=head1 EXAMPLES

=head2 Generate scripts for complete pipeline

perl kmeria_wrapper.pl --step all \
  --input /data/fastq_files \
  --output /results/kmeria_analysis \
  --samples sample.list \
  --depth-file sample_depth.tsv \
  --pheno phenotypes.txt \
  --threads 32 \
  --memory 32G \
  --kmer 31 \
  --min-abund 4 \
  --max-abund 1000 \
  --batch-size 4 \
  --ploidy 4 \
  --missing 0.6 \
  --scheduler slurm \
  --queue normal

=head2 Count k-mers using KMC (recommended)

perl kmeria_wrapper.pl --step count \
  --input /data/fastq_files \
  --output /results/01_kmer_counts \
  --samples sample.list \
  --threads 32 \
  --kmer 31 \
  --min-abund 4 \
  --batch-size 4 \
  --use-kmc \
  --kmc-memory 16 \
  --scheduler slurm

=head1 IMPORTANT NOTES ON SAMPLE ORDER

This version (2.0.1) addresses a **CRITICAL** issue where sample order could become
inconsistent between steps, leading to incorrect phenotype-genotype associations.

The fix works by:
1. Saving the sample order from the count step to sample_order.txt
2. Using this file in kctm to maintain the exact same order
3. Verifying order consistency in subsequent steps
4. Using the original order for all downstream analyses

**IMPORTANT**: If you're using pre-existing count results, ensure the samples are
processed in the same order as your phenotype and depth files!

=head2 Verifying Sample Order

After running the pipeline, you can verify sample order consistency:

  # Check if sample order is consistent
  diff 01_kmer_counts/sample_order.txt 02_kmer_matrices/kctm_sample_order.txt
  
  # Should show no differences if order is consistent

=head2 Sample Order Files

- **01_kmer_counts/sample_order.txt**: Original sample order from count step
- **02_kmer_matrices/kctm_sample_order.txt**: Sample order used in kctm
- **04_bimbam/tmp_sample.list**: Sample order used for VCF/PLINK conversion
- **05_association/sample.list**: Sample order used for association analysis

All these files should contain samples in the SAME order!

=head1 FILE FORMATS

=head2 Sample list file (--samples)
Plain text file with one sample name per line (no header):
  sample1
  sample2
  sample3

**IMPORTANT**: The order in this file will be preserved throughout the entire pipeline!

=head2 Depth file (--depth-file)
Tab-separated file with sample names and depth values:
  sample1    45.2
  sample2    52.8
  sample3    38.9

**IMPORTANT**: Sample order should match the sample list file!

=head2 Phenotype file (--pheno)
Tab or space-separated file with format:
  sample1  1.5    0
  sample2  2.3    1
  sample3  1.8    1

**IMPORTANT**: Sample order should match the sample list file!

=head1 TROUBLESHOOTING

=head2 Sample Order Mismatch

If you see warnings about sample order mismatch:

1. Check the sample_order.txt files in each step directory
2. Verify your input files (sample list, depth file, phenotype file) have consistent order
3. If using pre-existing count results, regenerate them with the correct sample order
4. Use the verification commands:
   
   diff 01_kmer_counts/sample_order.txt 02_kmer_matrices/kctm_sample_order.txt

=head2 Missing sample_order.txt

If sample_order.txt is missing from the count step:

1. The kctm step will fall back to alphabetical order (ls command)
2. This may cause sample order mismatch
3. **Solution**: Re-run the count step to generate sample_order.txt

=head1 VERSION

Version 2.0.1 - Sample order consistency fix (2025)

=head1 CHANGELOG

v2.0.1 (2025):
- **CRITICAL FIX**: Ensures consistent sample order across all pipeline steps
- Added sample_order.txt generation in count step
- Modified kctm step to use sample_order.txt instead of ls output
- Added sample order verification in filter and subsequent steps
- Prevents mismatch between sample order in matrices and phenotype/depth files
- All downstream steps now use the original sample order from count step

v2.0.0 (2025):
- Added support for all kmeria count parameters (-C, -T, -p, -c)
- Updated filter step to use new compressed output format
- Enhanced m2b step with BGZF compression and statistics
- Updated association step to use kassoc tool

=head1 AUTHOR

Updated for KMERIA v2.0.1 - Sample order consistency fix
Original concept based on KMERIA by Chen Shuai (chensss1209@gmail.com)

=head1 SEE ALSO

KMERIA documentation: https://github.com/Sh1ne111/KMERIA

=cut
