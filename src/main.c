/**
    The MIT License

   Copyright (c) Chen Shuai, 2024-

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

 * @details
 * KMERIA is a comprehensive toolkit for k-mer based GWAS, specifically
 * designed for polyploid organisms.
 * 
 * @author Chen Shuai (chensss1209@gmail.com)
 * 
 * @version 2.0.4
 * @date 2026-06-12
 * 
 * @see https://github.com/Sh1ne111/KMERIA
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

#ifndef VERSION
#define VERSION "v2.0.4"
#define DATE "2026-06-12"
#endif

#define PROGRAM_NAME "KMERIA"
#define AUTHOR "Chen Shuai"
#define EMAIL "chensss1209@gmail.com"
#define GITHUB_URL "https://github.com/Sh1ne111/KMERIA"

/* ============================================================================
 * Documentation
 * ============================================================================ */

/**
 * @brief Count k-mers from FASTA/FASTQ files
 *
 * @details
 * This function processes sequencing data files and generates k-mer count tables.
 * Supports both FASTA and FASTQ formats with automatic format detection.
 */
int kcount(int argc, char *argv[]);

/**
 * @brief Dump binary k-mers into text format
 *
 * @details
 * Converts binary k-mer database to Plain text format for inspection
 * and downstream analysis.
 */
int kdump(int argc, char *argv[]);

/**
 * @brief Build a population-level k-mer count matrix
 * 
 * @details
 * Constructs a comprehensive k-mer count matrix across all samples in the
 * population, essential for association analysis.
 */
extern int kctm(int argc, char *argv[]);

/**
 * @brief Filter raw k-mer matrix based on quality criteria
 * 
 * @details
 * Applies various filtering criteria (e.g., frequency, coverage, quality) to
 * remove low-quality or uninformative k-mers.
 */
extern int kfilter(int argc, char *argv[]);

/**
 * @brief Convert k-mer matrix to BIMBAM dosage format
 * 
 * @details
 * Transforms k-mer count matrix into BIMBAM format, compatible with various
 * association testing tools like GEMMA.
 */
extern int kmtob(int argc, char *argv[]);

/**
 * @brief Convert BIMBAM format to genotype format
 * 
 * @details
 * Converts dosage data in BIMBAM format to discrete genotype calls for
 * alternative analysis methods.
 */
extern int btog(int argc, char *argv[]);

/**
 * @brief Conduct association study using k-mers
 * 
 * @details
 * Performs genome-wide association study using k-mer markers, supporting
 * multiple statistical models and covariate correction.
 */
extern int assoc(int argc, char *argv[]);

/**
 * @brief Fetch reads containing specific k-mers from FASTQ files
 * 
 * @details
 * Extracts sequencing reads that contain specified k-mers for validation
 * and further analysis.
 */
int fetch_reads(int argc, char *argv[]);
int fetch_reads_tgs(int argc, char *argv[]);
/**
 * @brief Random sampling of k-mers for PCA and kinship calculation
 * 
 * @details
 * Performs random sampling of k-mers to construct a representative subset
 * for population structure analysis (PCA) and kinship matrix calculation.
 */
int ksketch(int argc, char *argv[]);

/**
 * @brief Extract reads containing k-mers from BAM files
 * 
 * @details
 * Searches aligned reads in BAM format for specific k-mers, useful for
 * validating k-mer associations in aligned context.
 */
int kbam(int argc, char* argv[]);

/**
 * @brief Add p-values to BAM file for significant k-mers
 * 
 * Annotates BAM files with p-values from association analysis, enabling
 * visualization of significant associations in genome browsers.
 */
int kaddp(int argc, char* argv[]);


static int usage() {
    fprintf(stderr, "\n");
    fprintf(stderr, "#===============================================================================#\n");
    fprintf(stderr, "#                                                                               #\n");
    fprintf(stderr, 
        "#                 _  ____  __ ______ _____  _____                               #\n"
        "#                | |/ /  \\/  |  ____|  __ \\|_   _|   /\\                         #\n"
        "#                | ' /| \\  / | |__  | |__) | | |    /  \\                        #\n"
        "#                |  < | |\\/| |  __| |  _  /  | |   / /\\ \\                       #\n"
        "#                | . \\| |  | | |____| | \\ \\ _| |_ / ____ \\                      #\n"
        "#                |_|\\_|_|  |_|______|_|  \\_\\_____/_/    \\_\\                     #\n"
        "#                                                                               #\n");
    fprintf(stderr, "#===============================================================================#\n\n");
    
    fprintf(stderr, "Program:  \e[1;31m%s\e[0m - A \e[1;31mKMER\e[0m-based genome-w\e[1;31mI\e[0mde \e[1;31mA\e[0;0mssociation testing approach\n", PROGRAM_NAME);
    fprintf(stderr, "          for polyploids\n\n");
    fprintf(stderr, "Version:  %s (%s)\n", VERSION, DATE);
    fprintf(stderr, "Author:   %s <%s>\n", AUTHOR, EMAIL);
    fprintf(stderr, "GitHub:   %s\n\n", GITHUB_URL);
    
    fprintf(stderr, "\e[1;33mUsage:\e[0m    \e[1;31mkmeria\e[0m <command> [options]\n\n");
    fprintf(stderr, "\e[1;33mCommands:\e[0m\n\n");
    
    fprintf(stderr, "  \e[1;32mData Processing:\e[0m\n");
    fprintf(stderr, "    count      Count k-mers from FASTA/FASTQ files\n");
    fprintf(stderr, "    dump       Convert binary k-mer file to plain text\n");
    fprintf(stderr, "    kctm       Build population k-mer counting matrix\n");
    fprintf(stderr, "    filter     Filter k-mer matrix by frequency and quality\n\n");
    
    fprintf(stderr, "  \e[1;32mFormat Conversion:\e[0m\n");
    fprintf(stderr, "    m2b        Convert k-mer matrix to BIMBAM dosage format\n");
    fprintf(stderr, "    b2g        Convert BIMBAM format to genotype format\n\n");
    
    fprintf(stderr, "  \e[1;32mAnalysis:\e[0m\n");
    fprintf(stderr, "    sketch     Random sampling for PCA and kinship calculation\n");
    fprintf(stderr, "    asso       Conduct k-mer genome-wide association study\n\n");
    
    fprintf(stderr, "  \e[1;32mUtilities:\e[0m\n");
    fprintf(stderr, "    fkr        Fetch reads associated k-mers from FASTQ files\n");
    fprintf(stderr, "    fkrtgs     Fetch reads associated k-mers from TGS FASTQ files\n");
    fprintf(stderr, "    kbam       Extract reads associated k-mers from BAM files\n");
    fprintf(stderr, "    addp       Annotate BAM with association p-values\n\n");
    
    fprintf(stderr, "\e[1;33mAdditional Help:\e[0m\n");
    fprintf(stderr, "    kmeria <command> -h     Show detailed help for specific command\n");
    fprintf(stderr, "    Visit %s for documentation\n\n", GITHUB_URL);
    
    fprintf(stderr, "#===========================================================================#\n");
    fprintf(stderr, "#  Citation: If you use KMERIA, please cite our paper at                    #\n");
    fprintf(stderr, "#         [https://doi.org/10.1038/s41588-026-02641-8]                      #\n");
    fprintf(stderr, "#===========================================================================#\n\n");
    
    return 1;
}


int main(int argc, char *argv[]) {
    if (argc < 2) {
        return usage();
    }
    
    const char *command = argv[1];
    
    if (strcmp(command, "-h") == 0 || strcmp(command, "--help") == 0) {
        return usage();
    }
    
    if (strcmp(command, "-v") == 0 || strcmp(command, "--version") == 0) {
        fprintf(stdout, "%s %s\n", PROGRAM_NAME, VERSION);
        return 0;
    }
    
    int ret = 0;
    
    if (strcmp(command, "count") == 0) {
        ret = kcount(argc - 1, argv + 1);
    } 
    else if (strcmp(command, "dump") == 0) {
        ret = kdump(argc - 1, argv + 1);
    } 
    else if (strcmp(command, "kctm") == 0) {
        ret = kctm(argc - 1, argv + 1);
    } 
    else if (strcmp(command, "filter") == 0) {
        ret = kfilter(argc - 1, argv + 1);
    } 
    else if (strcmp(command, "m2b") == 0) {
        ret = kmtob(argc - 1, argv + 1);
    } 
    else if (strcmp(command, "b2g") == 0) {
        ret = btog(argc - 1, argv + 1);
    } 
    else if (strcmp(command, "sketch") == 0) {
        ret = ksketch(argc - 1, argv + 1);
    } 
    else if (strcmp(command, "asso") == 0) {
        ret = assoc(argc - 1, argv + 1);
    } 
    else if (strcmp(command, "fkr") == 0) {
        ret = fetch_reads(argc - 1, argv + 1);
    }
    else if (strcmp(command, "fkrtgs") == 0) {
        ret = fetch_reads_tgs(argc - 1, argv + 1);
    } 
    else if (strcmp(command, "kbam") == 0) {
        ret = kbam(argc - 1, argv + 1);
    } 
    else if (strcmp(command, "addp") == 0) {
        ret = kaddp(argc - 1, argv + 1);
    } 
    else {
        fprintf(stderr, "\n\e[1;31m[ERROR]\e[0m Unrecognized command: '%s'\n", command);
        fprintf(stderr, "    Run 'kmeria -h' for available commands\n\n");
        return 1;
    }
    
    if (ret != 0) {
        fprintf(stderr, "\n\e[1;31m[ERROR]\e[0m Command '%s' failed with exit code %d\n\n", command, ret);
    }
    
    return ret;
}

