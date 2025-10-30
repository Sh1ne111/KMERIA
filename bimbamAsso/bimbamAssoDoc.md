## bimbamKin usage
```
Usage: bimbamKin [input_file] [output_file]
Required parameters
        [input_file] : bimbam dosage file (support .gz)
        [output_file]: output file name (required for bimbam)
        -b           : input is bimbam dosage format
Optional parameters
        -g           : input bimbam file is gzip compressed
        -d [# digits]: precision of the kinship values (default : 10)
        -s           : compute IBS matrix instead of BN
        -r           : randomly fill missing genotypes
        -h           : use hetero division for missing data
        -v           : turn on verbose mode
```

## bimbamAsso usage
```
Usage: bimbamKin [input_file] [output_file]
Required parameters
        [input_file] : bimbam dosage file (support .gz)
        [output_file]: output file name (required for bimbam)
        -b           : input is bimbam dosage format
Optional parameters
        -g           : input bimbam file is gzip compressed
        -d [# digits]: precision of the kinship values (default : 10)
        -s           : compute IBS matrix instead of BN
        -r           : randomly fill missing genotypes
        -h           : use hetero division for missing data
        -v           : turn on verbose mode
(kmeria_env) [agis_chenshuai@login02 software]$ bimbamAsso
Usage: bimbamAsso [options]

Required parameters:
  -p [phenotype_file] : Phenotype file with columns (FAMID INDID PHENO)
  -o [out_prefix]     : Output file prefix
  -t [genotype_file]  : Genotype file (TPED or Bimbam format)
                        Use -b flag for Bimbam format

Bimbam format options:
  -b                  : Input is in Bimbam dosage format
  -s [sample_file]    : Sample information file (FAMID INDID)
                        Required when using -b flag
  -g                  : Input genotype file is gzip compressed
                        Works with both TPED and Bimbam formats

Bimbam file format:
  Column 1: Marker ID (e.g., k-mer sequence)
  Column 2: Chromosome (can be placeholder like 'X')
  Column 3: Position (can be placeholder like 'Y')
  Column 4+: Dosage values for each sample (0-2 or dosages)
  Delimiter: comma (,)
  Example: AACCGAAA,X,Y,0.27,0.22,0.29,0.33,...

Sample file format (for -s):
  Column 1: Family ID
  Column 2: Individual ID
  Delimiter: space or tab

Kinship options:
  -k [kinship_file]   : Pre-computed kinship matrix (n×n)
  -K [method]         : Generate kinship matrix
                        1 = IBS with mean fill-in
                        2 = IBS with random fill-in
                        3 = Balding-Nichols

Optional parameters:
  -c [covar_file]     : Covariate file (FAMID INDID COV1 COV2 ...)
  -i [in_prefix]      : Input prefix for pre-computed eigenvectors
  -d [digits]         : Output precision (default: 5)
  -S [start_index]    : Start marker index (default: 0)
  -E [end_index]      : End marker index (default: all markers)
  -v                  : Verbose mode
  -w                  : Write eigenvalue/eigenvector files
  -N                  : Disable GLS (use OLS instead)

Output files:
  [prefix].ps         : Association results
                        Format: MARKER_ID BETA P_VALUE
  [prefix].reml       : REML estimates
  [prefix].log        : Log file
  [prefix].kinf       : Kinship matrix (if -K or -w used)
  [prefix].eLvals     : Eigenvalues (if -w used)
  [prefix].eLvecs     : Eigenvectors (if -w used)

Example usage:
  # Bimbam format with gzip
  bimbamAsso -b -g -t kmers.bimbam.gz -s samples.txt -p pheno.txt -k kinship.txt -o results
```
