
## This script is an example for downloading sweet potato SRA data and simulating phenotypes for KMERIA analyses

```bash

Usage: run_example_pipe.sh -s <sra_list_file> -p <phenotype_file> -d <depth_file> [options]

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
    -t, --threads       Number of threads to use (default: 8)
    -m, --memory        Memory allocation (default: 32G)
    -k, --kmer          K-mer size (default: 31)
    --batch-size        Batch size for processing (default: 1)
    --ploidy            Genome ploidy (default: 6)
    --pheno-col         Phenotype column number (default: 2)
    -h, --help          Show this help message

Example Usage:
    run_example_pipeline.sh -s sra_list.txt -p phenotypes.txt -d sample_depths.txt -t 16 -m 16G

Steps:
    1. Download SRA files from NCBI
    2. Convert SRA files to FASTQ format
    3. Prepare input files for KMERIA analysis
    4. Generate KMERIA job scripts

Note: This pipeline generates job scripts that need to be submitted manually to your cluster system.
```
