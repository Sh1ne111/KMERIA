#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <pthread.h>
#include "khash.h"
#include "sam.h"

// The version of the program
#define KATCHER_VERSION "0.1"

// The initial size of the array of strings created to hold k-mers by read_kmer()
#define N_KMERS_INIT 1000

// The length of the k-mers
#define KMER_LENGTH 31

// The maximum number of k-mers that the program can process
#define MAX_KMERS 100000

// Initializing a hash table type named "kmer_hash" using klib's khash functions
KHASH_SET_INIT_STR(kmer_hash);

// Declaring the lock used to prevent multiple threads from writing to
// the output file simultaneously
extern pthread_mutex_t lock;

// A struct that holds the data required for a single analysis thread to execute
//   kmer_table: a pointer to the hash table containing the k-mers
//   hash_end:   the end of the hash table ; will not change after the table has been computed
//   output:     a pointer to the sam output file
//   header:     a pointer to the header
//   bambuf:     a pointer to the buffered records being processed
//   n_records:  the number of records that will be processed by this thread
//   seq:        pointer to the sequence being queried
typedef struct thread_data {
	khash_t(kmer_hash) *kmer_table;
	khint_t hash_end;
	samFile *output;
	sam_hdr_t *header;
	bam1_t **bambuf;
	int n_records;
	char **seq;
} thread_data;

// A struct that holds the data for the decompression thread reading the bam file
//   bamfile:  a pointer to the input bam file
//   header:   a pointer to the header of the input bam file
//   tmpbuf:   a pointer to the temporary buffer into which records are read
//   max_read: the maximum number of records to read from the file (equal to buffer size)
//   n_read:   the number of records read
typedef struct buffer_data {
	samFile *bamfile;
	sam_hdr_t *header;
	bam1_t **tmpbuf;
	int max_read;
	int n_read;
} buffer_data;

// A struct that holds the data for copying the records from the temporary buffer to the main buffer
//   tmpbuf:     a pointer to the temporary buffer into which records are read
//   bambuf:     a pointer to the main buffer that is being processed by analysis threads
//   seq:        a pointer to the character sequences of the reads associated with bambuf
//   n_records:  the number of records that this thread must process
//   max_length: the maximum length of reads; used to check that it does not exceed seq_alloc
typedef struct copy_data {
	bam1_t **tmpbuf;
	bam1_t **bambuf;
	char **seq;
	int n_records;
	int max_length;
} copy_data;

// A struct that holds the parameters as read from the command line
//   input_file:  name of the input bam file
//   kmer_list:   name of the file containing the list of k-mers
//   output_file: name of the output sam file
//   n_threads:   number of analysis threads used
//   bufsize:     number of records read at once into the buffer
//   seq_alloc:   number of bytes to allocate for the read sequence (read length + 1 for ending NULL)
typedef struct params {
	char* input_file;
	char* kmer_list;
	char* output_file;
	int n_threads;
	int bufsize;
	int seq_alloc;
} params;

// A function that initializes the buffers used by the program
// bambuf: a pointer to the buffer that will store the alignments being processd
// tmpbuf: a pointer to the temporary buffer of alignments into which the records are read from file
// seq: a pointer to an array of strings that will hold the read sequences
// bufsize: the maximum number of records or sequences that can be in the buffer at any given time
// max_length: the maximum length that a read can be (determines the memory allocated to each sequence)
void init_buffers(bam1_t **bambuf, bam1_t **tmpbuf, char **seq, int bufsize, int max_length);

// A function that takes the read sequence as encoded by 4-bit nucleotides
// in the BAM format and converts it to a string
// bamseq: a pointer to the start of the sequence in the BAM file
// seq: a pointer to an array of characters large enough to hold the read sequence
// seqlength: the length of the read sequence (in nucleotides)
void bamseq_to_char(uint8_t *bamseq, char *seq, int seqlength);

// A function that generates the reverse-complemented version of a k-mer
// forward: a pointer to an array of characters holding the forward k-mer
// reverse: a pointer to an array of characters that will hold the reverse k-mer (must be sufficiently large)
// kmer_length: the length the the k-mers (must be the same for all k-mers)
// n_kmers: the adress of an int variable where the number of k-mers read will be stored
void revcomp(char *forward, char *reverse, int kmer_length);

// A function that allows matching read sequence and k-mers with multithreading
void *match_kmers(void *input);

// A function that fills a temporary buffer with bam records
void *fillbuf(void* input);

// A function that copies alignments from the temporary buffer to the main one
void *copy_buf(void *input);

// A function meant for multithreading that generates a character string sequence from the bam representation
void *prepare_seq(void *input);

// A function to parse the command line arguments
int parse_args(int argc, char* argv[], params *args);

// A function to read a list of k-mers from a text file
// kmer_list: a pointer to a file open for reading which contains the list of k-mers
// n_init: the size to initialize the array of kmers with
// kmer_length: the length of the k-mers (will throw an error if it does not correspond)
char **read_kmers(FILE *kmer_list, int n_init, int kmer_length, int *n_kmers);

#endif
