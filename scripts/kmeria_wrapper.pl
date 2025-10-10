#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use File::Path qw(make_path);
use Cwd qw(abs_path);
use threads;
use Pod::Usage;

our $VERSION = "1.0.2";
our $PROGRAM = "kmeria_wrapper.pl";

# Default parameters
my $threads = 16;         # Default number of threads
my $memory = "32G";       # Default memory per job
my $kmer_size = 31;       # Default k-mer size：31
my $min_abundance = 5;    # Default minimum k-mer abundance
my $max_abundance = 1000; # Default maximum k-mer abundance
my $JOB_SCHEDULER = "local"; # Default scheduler (local, slurm, sge, pbs)
my $queue = "share";     # Default queue for job submission
my $WALLTIME = "720:00:00"; # Default time limit
my $batch_size = 4;      # Number of samples to process in each batch (default: 4)
my $input_dir = ""; 
my $output_dir = "";      # Default output directory
my $help = 0;
my $step = "";            # Step to execute
my $sample_list = "";     # List of samples
my $min_recurrence = 0;   #
my $kctm_missing = 0.8;   # Default missing ratio for 'kmeria kctm'
my $ploidy = 4;           # Default ploidy
my $depth_file = "";      # Sample depth file
my $covar_file = "";      # Covariate file for association testing
my $kinship_file = "";    # Kinship matrix file
my $pheno_file = "";      # Phenotype file
my $pheno_col = 1;        # Phenotype column
my $use_kmc = 0;          # Whether to use 'KMC' instead of 'kmeria count'
my $kmc_memory = "32";    # Default memory for KMC (in GB)
my $use_kctm2 = 0;

# Define standard directory structure 
my $count_dir;
my $kctm_dir;
my $filter_dir;
my $bimbam_dir;
my $asso_dir;

# Parse command-line options
GetOptions(
    "threads|t=i"       => \$threads,
    "memory|m=s"        => \$memory,
    "kmer|k=i"          => \$kmer_size,
    "min-abund|c=i"     => \$min_abundance,
    "max-abund|C=i"     => \$max_abundance,
    "scheduler|s=s"     => \$JOB_SCHEDULER,
    "queue|q=s"         => \$queue,
    "time=s"            => \$WALLTIME,
    "batch_size|b=i"    => \$batch_size,
    "output|o=s"        => \$output_dir,
    "help|h"            => \$help,
    "step=s"            => \$step,
    "samples=s"         => \$sample_list,
    "min-recur=i"       => \$min_recurrence,
    "missing=f"         => \$kctm_missing,
    "ploidy|p=i"        => \$ploidy,
    "depth-file|d=s"    => \$depth_file,
    "input|i=s"         => \$input_dir,
    "covar|c=s"         => \$covar_file,
    "kinship|k=s"       => \$kinship_file,
    "pheno=s"           => \$pheno_file,
    "pheno-col|n=i"     => \$pheno_col,
    "use-kmc"           => \$use_kmc,          # Option for using 'KMC' to count k-mer
    "kmc-memory=i"      => \$kmc_memory,       # Memory for KMC (in GB)
    "use-kctm2"         => \$use_kctm2, 
) or pod2usage(2);

# Display help message if requested
pod2usage(-verbose => 2) if $help;

# Check for required arguments
unless ($step) {
    print "Error: --step is required.\n";
    pod2usage(1);
}

# Create output directory if it doesn't exist
if ($output_dir) {
    make_path($output_dir) unless -d $output_dir;
    $output_dir = abs_path($output_dir);
    
    # Initialize standard directory paths based on output directory
    initialize_directories($output_dir);
}

# Function to initialize standard directory paths
sub initialize_directories {
    my ($base_dir) = @_;
    
    $count_dir = "$base_dir/01_kmer_counts";
    $kctm_dir = "$base_dir/02_kmer_matrices";
    $filter_dir = "$base_dir/03_filtered_matrices";
    $bimbam_dir = "$base_dir/04_bimbam";
    $asso_dir = "$base_dir/05_association";
    
    return ($count_dir, $kctm_dir, $filter_dir, $bimbam_dir, $asso_dir);
}

# Main workflow controller
if ($step eq "count") {
    run_count_step();
} elsif ($step eq "kctm") {
    run_kctm_step();
} elsif ($step eq "flt") {
    run_filter_step();
} elsif ($step eq "m2b") {
    run_convert_step();
} elsif ($step eq "asso") {
    run_asso_step();
} elsif ($step eq "all") {
    run_full_pipeline();
} else {
    print "Error: Unknown step '$step'. Valid steps are: count, kctm, flt, m2b, asso, all\n";
    exit 1;
}

# Function to run the complete pipeline
sub run_full_pipeline {
    # Check required parameters
    if (!$input_dir || !$output_dir) {
        die "Error: Input directory (-i) and output directory (-o) are required for pipeline command\n";
    }
    
    # Create output directory structure (already initialized in the main script)
    make_path($count_dir, $kctm_dir, $filter_dir, $bimbam_dir, $asso_dir);
    
    # Run each step sequentially
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
    
    print "All scripts have been generated successfully!\n";
    print "Please submit the job scripts manually to your cluster system.\n";
}

# Function to create job script header
sub create_job_header {
    my ($fh, $job_name, $log_prefix) = @_;
    
    if ($JOB_SCHEDULER eq "slurm") {
        print $fh "#!/bin/bash\n";
        print $fh "#SBATCH --job-name=$job_name\n";
        print $fh "#SBATCH --output=${log_prefix}.log\n";
        print $fh "#SBATCH --error=${log_prefix}.err\n";
        print $fh "#SBATCH --ntasks=1\n";
        print $fh "#SBATCH --cpus-per-task=$threads\n";
        print $fh "#SBATCH --mem=$memory\n";
        print $fh "#SBATCH --time=$WALLTIME\n";
        print $fh "#SBATCH --partition=$queue\n";
    } elsif ($JOB_SCHEDULER eq "sge") {
        print $fh "#!/bin/bash\n";
        print $fh "#\$ -N $job_name\n";
        print $fh "#\$ -o ${log_prefix}.log\n";
        print $fh "#\$ -e ${log_prefix}.err\n";
        print $fh "#\$ -pe smp $threads\n";
        print $fh "#\$ -l h_vmem=$memory\n";
        print $fh "#\$ -l h_rt=$WALLTIME\n";
        print $fh "#\$ -q $queue\n";
    } elsif ($JOB_SCHEDULER eq "pbs") {
        print $fh "#!/bin/bash\n";
        print $fh "#PBS -N $job_name\n";
        print $fh "#PBS -o ${log_prefix}.log\n";
        print $fh "#PBS -e ${log_prefix}.err\n";
        print $fh "#PBS -l select=1:ncpus=$threads:mem=$memory\n";
        print $fh "#PBS -l walltime=$WALLTIME\n";
        print $fh "#PBS -q $queue\n";
    } else {
        # Local execution, just add a basic header
        print $fh "#!/bin/bash\n";
        print $fh "# Local execution job: $job_name\n";
        print $fh "# Log files: ${log_prefix}.log, ${log_prefix}.err\n\n";
    }
}

# Function to run k-mer counting
sub run_count_step {
    my ($input, $output) = @_;
    
    # Handle input/output directory
    $input = $input_dir if !$input;
    
    if (!$output) {
        $output = $count_dir;
    }
    
    # Check required parameters
    if (!$input || !$output) {
        die "Error: Input directory (-i) and output directory (-o) are required for count command\n";
    }
    
    # Create output directory if it doesn't exist
    make_path($output) unless -d $output;
    $output = abs_path($output);
    
    # Get list of samples to process
    my @samples;
    if (-d $input) {
        $input = abs_path($input);
        # If input is a directory, find all FASTQ files
        opendir(DIR, $input) or die "Cannot open directory $input: $!";
        my %sample_files;
        while (my $file = readdir(DIR)) {
            if ($file =~ /^(.*?)(?:_R?[12])?\.(?:fastq|fq)(?:\.gz)?$/i) {
                my $sample = $1;
                $sample_files{$sample} = 1;
            }
        }
        closedir(DIR);
        @samples = sort keys %sample_files;
    } elsif ($sample_list && -f $sample_list) {
        # If a sample list file is provided, read samples from it
        open(my $fh, '<', $sample_list) or die "Cannot open file $sample_list: $!";
        while (my $line = <$fh>) {
            chomp $line;
            push @samples, $line if $line;
        }
        close($fh);
    } elsif (-f $input) {
        # If input is a file, assume it's a list of samples
        open(my $fh, '<', $input) or die "Cannot open file $input: $!";
        while (my $line = <$fh>) {
            chomp $line;
            push @samples, $line if $line;
        }
        close($fh);
    } else {
        die "Error: Input ($input) is neither a directory nor a file\n";
    }
    
    if (scalar(@samples) == 0) {
        die "Error: No samples found to process\n";
    }
    
    # Process samples in batches
    my $total_samples = scalar(@samples);
    my $num_batches = int(($total_samples + $batch_size - 1) / $batch_size);
    
    print "Found $total_samples samples, generating $num_batches batch scripts\n";
    
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
        create_job_header($fh, "kmeria_count_$batch", "$output/count_batch_$batch");
        
        # Write actual commands for each sample in the batch
        foreach my $sample (@batch_samples) {
            # Find corresponding FASTQ files
            my @fastq_files;
            if (-d $input) {
                $input = abs_path($input);
                opendir(DIR, $input) or die "Cannot open directory $input: $!";
                while (my $file = readdir(DIR)) {
                    if ($file =~ /^${sample}(?:_R?[12])?\.(?:fastq|fq)(?:\.gz)?$/i) {
                        push @fastq_files, "$input/$file";
                    }
                }
                closedir(DIR);
            } else {
                push @fastq_files, $sample;
            }
            
            # Skip if no FASTQ files found
            if (scalar(@fastq_files) == 0) {
                print $fh "echo 'No FASTQ files found for sample $sample, skipping'\n";
                next;
            }
 
            if ($use_kmc) {
                # Create temporary directory for KMC and a file list for input files
                print $fh "mkdir -p $output/$sample\n";
                print $fh "echo 'Processing sample $sample...'\n";
                
                # Create a temporary file list for KMC input
                print $fh "# Create a temporary file listing all FASTQ files for this sample\n";
                print $fh "cat > $output/${sample}_fastq_list.txt << EOF\n";
                foreach my $fastq_file (@fastq_files) {
                    print $fh "$fastq_file\n";
                }
                print $fh "EOF\n\n";
                
                # Run KMC using the file list
                print $fh "# Run KMC with the file list\n";
				# Turn off the cononical mode and fixed the error of max abundance parameter 
                print $fh "kmc -k$kmer_size -t$threads -m$kmc_memory -b -ci$min_abundance -cs$max_abundance \\\n";
                print $fh "    \@$output/${sample}_fastq_list.txt $output/${sample}_k${kmer_size} $output/$sample\n";
                
                # Dump k-mers
                print $fh "# Dump k-mers to TSV format\n";
                print $fh "kmc_tools transform $output/${sample}_k${kmer_size} dump -s $output/${sample}_k${kmer_size}.kc.tsv\n";
                
                # Compress the tsv file
                print $fh "# Compress the TSV file\n";
                print $fh "pigz -p $threads $output/${sample}_k${kmer_size}.kc.tsv\n";
                
                # Cleanup temporary files
                print $fh "# Cleanup temporary files\n";
                print $fh "rm -f $output/${sample}_fastq_list.txt\n";
                print $fh "rm -rf $output/$sample\n";
            } else {
                print $fh "echo 'Processing sample $sample...'\n";
                print $fh "kmeria count @fastq_files -t $threads -k $kmer_size -o $output/${sample}_k${kmer_size}.kc\n";
                
                # Dump k-mers
                print $fh "kmeria dump $output/${sample}_k${kmer_size}.kc -c $min_abundance -C $max_abundance -o $output/${sample}_k${kmer_size}.kc.tsv\n";
                
                # Compress the tsv file
                print $fh "pigz -p $threads $output/${sample}_k${kmer_size}.kc.tsv\n";
            }
            
            print $fh "echo 'Finished processing sample $sample'\n\n";
        }
        
        close($fh);
        chmod 0755, $job_script;
        
        print "Generated script: $job_script\n";
    }
    
    print "All k-mer counting job scripts have been generated. Please submit them manually to your cluster system.\n";
}

# Function to build k-mer matrices
sub run_kctm_step {
    my ($input, $output) = @_;
    
    # Handle input/output directory
    $input = $input_dir if !$input;
    
    if (!$output) {
        $output = $kctm_dir;
    }
    
    # Check required parameters
    if (!$input || !$output) {
        die "Error: Input directory (-i) and output directory (-o) are required for kctm command\n";
    }
    
    # Create output directory if it doesn't exist
    make_path($output) unless -d $output;
    $output = abs_path($output);
    
    # Create job script
    my $job_script = "$output/kctm_job.sh";
    open(my $fh, '>', $job_script) or die "Cannot create job script $job_script: $!";
    
    # Write job header
    create_job_header($fh, "kmeria_kctm", "$output/kctm_job");
    
    # Write actual commands
    print $fh "echo 'Building k-mer matrix...'\n";
    
    if ($use_kctm2) {
        # Create a sample list file
        print $fh "ls $input/*_k${kmer_size}.kc.tsv.gz > $output/sample.list\n";
        print $fh "# Using kctm2 - faster method\n";
        print $fh "kmeria kctm2 -i $output/sample.list -o $output/kmer_matrices\n";
    } else {
        # Use kctm
        print $fh "kmeria kctm -m $min_recurrence -c $min_abundance -x $max_abundance $input/*_k${kmer_size}.kc.tsv.gz -o $output/kmer_matrices\n";
    }
    
    print $fh "echo 'Finished building k-mer matrix'\n";
    close($fh);
    chmod 0755, $job_script;
    
    print "Generated script: $job_script\n";
    print "K-mer matrix building job script has been generated. Please submit it manually to your cluster system.\n";
}

# Function to filter k-mer matrix
sub run_filter_step {
    my ($input, $output) = @_;
    
    # Handle input/output directory
    $input = $input_dir if !$input;
    
    if (!$output) {
        $output = $filter_dir;
    }
    
    # Check required parameters
    if (!$input || !$output) {
        die "Error: Input directory (-i) and output directory (-o) are required for filter command\n";
    }
    
    # Create output directory if it doesn't exist
    make_path($output) unless -d $output;
    $output = abs_path($output);
    
    # Verify depth file
    if (!$depth_file) {
        die "Error: Depth file (-d) is required and must exist for filter command\n";
    }
    
    # Create job script
    my $job_script = "$output/filter_job.sh";
    open(my $fh, '>', $job_script) or die "Cannot create job script $job_script: $!";
    
    # Write job header
    create_job_header($fh, "kmeria_filter", "$output/filter_job");
    
    # Write actual commands
    print $fh "echo 'Filtering k-mer matrix...'\n";
    $depth_file = abs_path($depth_file);
    print $fh "kmeria flt -i $input -o $output -t $threads -c $max_abundance \\\n";
    print $fh "            -s $kctm_missing -p $ploidy -d $depth_file\n";
    print $fh "echo 'Finished filtering k-mer matrix'\n";
    
    close($fh);
    chmod 0755, $job_script;
    
    print "Generated script: $job_script\n";
    print "K-mer matrix filtering job script has been generated. Please submit it manually to your cluster system.\n";
}

# Function to convert k-mer matrix to BIMBAM format
sub run_convert_step {
    my ($input, $output) = @_;
    
    # Handle input/output directory
    $input = $input_dir if !$input;
    
    if (!$output) {
        $output = $bimbam_dir;
    }
    
    # Check required parameters
    if (!$input || !$output) {
        die "Error: Input directory (-i) and output directory (-o) are required for convert command\n";
    }
    
    # Create output directory if it doesn't exist
    make_path($output) unless -d $output;
    $output = abs_path($output);
    
    # Create job script
    my $job_script = "$output/convert_job.sh";
    open(my $fh, '>', $job_script) or die "Cannot create job script $job_script: $!";
    
    # Write job header
    create_job_header($fh, "kmeria_convert", "$output/convert_job");
    
    # Write actual commands
    print $fh "echo 'Converting to BIMBAM format...'\n";
    print $fh "kmeria m2b --in $input --out $output --threads $threads\n";
    print $fh "kmeria sketch $output/*.bimbam -n 10000000 > $output/sampling.bimbam\n";
    
    # Check if sample list is provided
    if ($sample_list && -f $sample_list) {
       $sample_list = abs_path($sample_list); 
       print $fh "kmeria b2g -i $output/sampling.bimbam -s $sample_list -o $output/sampling.vcf\n";
    } else {
        print $fh "# Generating sample list from input directory\n";
        $depth_file = abs_path($depth_file);
        print $fh "cut -f1 $depth_file| sort | uniq > $output/tmp_sample.list\n";
        print $fh "kmeria b2g -i $output/sampling.bimbam -s $output/tmp_sample.list -o $output/sampling.vcf\n";
    }
    
    print $fh "plink --vcf $output/sampling.vcf --make-bed --out $output/sampling\n";
    print $fh "echo 'Finished converting to BIMBAM format'\n";
    
    close($fh);
    chmod 0755, $job_script;
    
    print "Generated script: $job_script\n";
    print "Conversion to BIMBAM format job script has been generated. Please submit it manually to your cluster system.\n";
}

# Function to run association analysis
sub run_asso_step {
    my ($input, $output) = @_;
    
    # Handle input/output directory
    $input = $input_dir if !$input;
    
    if (!$output) {
        $output = $asso_dir;
    }
    
    # Check required parameters
    if (!$input || !$output || !$pheno_file) {
        die "Error: Input directory (-i), output directory (-o), and phenotype file (--pheno) are required for asso command\n";
    }
    
    # Create output directory if it doesn't exist
    make_path($output) unless -d $output;
    $output = abs_path($output);
    
    # Create PCA and kinship directories
    make_path("$output/pca");
    
    # Create job script
    my $job_script = "$output/asso_job.sh";
    open(my $fh, '>', $job_script) or die "Cannot create job script $job_script: $!";
    
    # Write job header
    create_job_header($fh, "kmeria_asso", "$output/asso_job");
    
    # Write actual commands
    print $fh "echo 'Association study Starting...'\n";
    
    # Generate PCA if needed
    if (!$covar_file) {
        print $fh "# Generating PCA for covariates\n";
        print $fh "plink --bfile $input/sampling --allow-extra-chr --allow-no-sex --pca 10 --out $output/pca/samplPCA\n";
        print $fh "cat $output/pca/samplPCA.eigenvec | awk '{print \"1\",\$3,\$4,\$5}' > $output/pca/PCA.txt\n";
        $covar_file = "$output/pca/PCA.txt";
    }
    
    # Generate kinship matrix if needed
    if (!$kinship_file) {
        $kinship_file = abs_path($kinship_file);
	$pheno_file   = abs_path($pheno_file);
        print $fh "# Generating kinship matrix\n";
        print $fh "gemma -bfile $input/sampling -gk -p $pheno_file -outdir $output/kinship -o kinship\n";
        $kinship_file = "$output/kinship/kinship.cXX.txt";
    }
    
    # Run association analysis
    print $fh "echo 'Running association analysis...'\n";
    my $asso_output = "$output/asso_output";
    print $fh "kmeria asso --input $input --output $asso_output --threads $threads \\\n";
    print $fh "           --covar $covar_file --kinship $kinship_file \\\n";
    print $fh "           --pheno $pheno_file --pheno_cols $pheno_col\n";
    print $fh "echo 'Finished association analysis.'\n";
    
    close($fh);
    chmod 0755, $job_script;
    
    print "Generated script: $job_script\n";
    print "Association job script has been generated. Please submit it manually to your cluster system.\n";
}

__END__

=head1 NAME

kmeria_wrapper.pl - A parallel wrapper for KMERIA pipeline.

=head1 SYNOPSIS

kmeria_wrapper.pl --step <step> [options]

 Steps:
    count            Count k-mers from FASTA/FASTQ files
    kctm             Build a population-level k-mer matries
    flt              Filter the raw k-mer matries
    m2b              Convert k-mer count matries to dosage (BIMBAM) format
    asso             Conduct a k-mer association study
    all              Run the complete KMERIA pipeline

 Global Options:
   --help|-h                    Show detailed help message
   --threads|-t      [INT]      Number of threads per job [default: 16]
   --step            [STR]      Pipeline step to run (required)
   --scheduler|-s    [STR]      Job scheduler (local, slurm, sge, pbs) [default: local]
   --samples         [FILE]     File containing list of samples (one per line)
   --queue|-q        [STR]      Queue/Partition for job submission [default: normal]
   --memory|-m       [STR]      Memory per job [default: 32G]
   --time            [STR]      Time limit for job [default: 720:00:00]

Command-specific options:

 For 'count':   
   --input|-i        [DIR]       Directory with FASTQ files or a list file of samples
   --output|-o       [DIR]       Output directory [default: current directory]
   --kmer|-k         [INT]       K-mer size [default: 31]
   --min-abund|-c    [INT]       Minimum k-mer abundance [default: 5]
   --max-abund|-C    [INT]       Maximum k-mer abundance [default: 1000]
   --batch_size|-b   [INT]       Number of samples to process in each batch (default: 50)
   --use-kmc                     Use KMC instead of KMERIA for k-mer counting
   --kmc-memory      [INT]       Memory allocation for KMC in GB [default: 32]
 
 For 'kctm' or 'kctm2':
   --input|-i        [DIR]       Directory with k-mer count files
   --output|-o       [DIR]       Output directory for k-mer matrices
   --min-recur       [INT]       Minimum recurrence number of k-mer [1], common recurrence ratio: 20-50%
   --use-kctm2                   Same function as kctm but faster

 For 'flt':
   --input|-i        [DIR]       Directory with k-mer abundance matrices
   --output|-o       [DIR]       Output directory for filtered matrices
   --max-abund|-C    [INT]       Maximum k-mer abundance [default: 1000]
   --missing         [FLOAT]     Missing ratio for kctm [default: 0.8]
   --ploidy|-p       [INT]       Genome ploidy [default: 4]
   --depth-file|-d   [FILE]      Sample depth file

 For 'm2b'
   --input|-i        [DIR]       Directory with filtered k-mer matrices
   --output|-o       [DIR]       Output directory

 For 'asso' 
   --input|-i        [DIR]       Directory with BIMBAM files
   --output|-o       [DIR]       Output directory for association results
   --covar|-c        [FILE]      Covariate file for association testing
   --kinship|-k      [STR]       Kinship matrix file
   --pheno|-s        [FILE]      Phenotype file
   --pheno-col|-n    [INT]       Phenotype column [default: 1]


=head1 DESCRIPTION

kmeria_wrapper.pl is a wrapper script for generating job scripts for the KMERIA pipeline,
with support for various job schedulers (SLURM, SGE, PBS). These scripts need to be manually
submitted to the cluster system.

=head1 EXAMPLES

=head2 Generate scripts for the complete pipeline using KMERIA

perl /path/to/kmeria_wrapper.pl --step all \
  --input /path/to/fastq_files \
  --output /path/to/kmeria_results \ 
  --samples sample.list \
  --threads 32 \
  --memory 32G \
  --kmer 31 \
  --min-abund 5 \
  --max-abund 1000 \
  --batch_size 2 \
  --use-kmc \
  --use-kctm2 \   # Optional
  --kmc-memory 32 \  #  Optional
  --ploidy 4 \
  --depth-file /path/to/sample_depths.txt \
  --pheno /path/to/phenotypes.txt \
  --pheno-col 1 \
  --scheduler slurm \
  --queue hebhcnormal01

=head2 Count k-mers for multiple samples in parallel using KMERIA

perl /path/to/kmeria_wrapper.pl --step count \
  --input /path/to/fastq_files \
  --output /path/to/output/kmer_counts \
  --threads 32 \
  --memory 32G \
  --kmer 31 \
  --min-abund 5 \
  --max-abund 1000 \
  --batch_size 2 \
  --use-kmc \  Optional
  --kmc-memory 32 \
  --scheduler slurm \
  --queue hebhcnormal01   # queue name


=head1 REQUIREMENTS

- KMERIA (https://github.com/Sh1ne111/KMERIA)
- KMC (https://github.com/refresh-bio/KMC) - optional, for KMC counting mode
- Perl modules: Getopt::Long, File::Basename, File::Path, Cwd, Pod::Usage
- Job schedulers (optional): SLURM, SGE, or PBS
- pigz (for parallel compression)

=head1 AUTHOR

Based on KMERIA by Chen Shuai (chensss1209@gmail.com)
