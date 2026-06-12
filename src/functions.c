#include "../include/functions.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <pthread.h>
#include <getopt.h>

#include "../include/sam.h"
#include "../include/khash.h"

// Global static variables used in this file
// Creating a global table used to translate 4-bit integers into chars
static char nt_table[16] = "*AC*G***T******N";

// Defining the lock declared in the header
// Used to prevent threads from writing to the output file simultaneously
pthread_mutex_t lock;

// bamseq_to_char
void bamseq_to_char(uint8_t *bamseq, char *seq, int seqlength) {

	int i;

	for(i = 0; i < seqlength; i++) {
		seq[i] = nt_table[bam_seqi(bamseq, i)];
	}

	// Appending a NULL character so that this is understood as a string
	seq[i] = '\0';

	return;
}

// revcomp
void revcomp(char *forward, char *reverse, int kmer_length) {

	for(int i = 0; i < kmer_length; i++) {
		switch (forward[i]) {
			case 'A':
				reverse[kmer_length - 1 - i] = 'T';
				break;

			case 'T':
				reverse[kmer_length - 1 - i] = 'A';
				break;

			case 'G':
				reverse[kmer_length - 1 - i] = 'C';
				break;

			case 'C':
				reverse[kmer_length - 1 - i] = 'G';
				break;

			default:
				fprintf(stderr, "ERROR: Unknown nucleotide found in k-mer sequence %s\n", forward);
				exit(1);
		}
	}

	// Ending the array with a NULL character to indicate the end of the string
	reverse[kmer_length] = '\0';
}

// A function that checks whether the k-mers in a given read are found in the k-mer hash table
// This function is set for multithreading
// CAUTION: this function modifies the seq char** by setting the last character to a \0 character until the whole
// read was processed. This is OK in the current implementation because each record is accessed only by one thread
// However, this should be taken into account if the implementation were to change
void *match_kmers(void *input) {
	thread_data* chunk = (thread_data*) input;
	int write_op, seqlength, match, tag_size;
	char *j, *seqstart, *km_tag;

	// Iterating over all the records
	for(int i = 0; i < chunk->n_records; i++) {

		// Resetting the match flag
		match = 0;

		// Initializing the KM tag to an empty string
		tag_size = (KMER_LENGTH + 2) * 10;
		km_tag = (char*) malloc(tag_size * sizeof(char));
		km_tag[0] = '\0';

		// Getting the length of the read and a pointer to its start
		seqlength = strlen(chunk->seq[i]);
		seqstart = chunk->seq[i];

		// Iterating over all the k-mers in the read
		// j is the last character of the read
		for(j = seqstart + seqlength - 1; j > seqstart + KMER_LENGTH - 2; j--) {

			// Checking if the k-mer until the end of the read is found in the k-mer hash table
			if(kh_get(kmer_hash, chunk->kmer_table, j - KMER_LENGTH + 1) != chunk->hash_end) {
				// Then we first mark this record as having to be written
				match = 1;

				// Then we add the k-mer being queried to the aux KM tag that will be written
				if(strlen(km_tag) + KMER_LENGTH + 2 >= tag_size) {
					tag_size *= 2;
					km_tag = (char*) realloc(km_tag, tag_size * sizeof(char));
				}

				if(strlen(km_tag)) {
					strcat(km_tag, ",");
				}

				if(strcat(km_tag, j - KMER_LENGTH + 1) == NULL) {
					fprintf(stderr, "Error concatenating k-mer %s to auxiliary tag.\n", j - KMER_LENGTH + 1);
					exit(1);
				}
			}

			// Setting the last character in the string to \0 so we can process the next k-mer
			*j = '\0';
		}

		// At the end of this record we check if it should be written to output
		if(match == 1) {
			// If it does we append the KM tag to the record and then write it to file
			if(bam_aux_append(chunk->bambuf[i], "KM", 'Z', strlen(km_tag) + 1, km_tag) != 0) {
				fprintf(stderr, "Error appending KM tag %s to BAM record. Aborting.\n", km_tag);
				exit(1);
			}

			pthread_mutex_lock(&lock);
			write_op = sam_write1(chunk->output, chunk->header, chunk->bambuf[i]);
			pthread_mutex_unlock(&lock);
		}

		free(km_tag);
	}
}

// A function that fills a buffer with bam records and counts the number of records that were read
// This function is set such that it can be run as an independent thread that will decompress
// the input while the records that have already being read are being processed
void *fillbuf(void* input) {
	buffer_data *data = (buffer_data*) input;
        int read_op;

	// Resetting the number of records read
        data->n_read = 0;

	// Reading the records
        while(data->n_read < data->max_read && (read_op = sam_read1(data->bamfile, data->header, data->tmpbuf[data->n_read])) > 0) {
                data->n_read++;
        }
}

// A function that copies alignment records from the temporary buffer to the main one
void *copy_buf(void *input) {
	copy_data *data = (copy_data*) input;
	bam1_t *result;

	for(int i = 0; i < data->n_records; i++) {
		result = bam_copy1(data->bambuf[i], data->tmpbuf[i]);
	}
}

// A function meant for multithreading that generates a character string sequence from the bam representation
void *prepare_seq(void *input) {
	copy_data *data = (copy_data*) input;
	uint8_t *bamseq = NULL;
	int32_t seqlength;

	for(int i = 0; i < data->n_records; i++) {
		// Also preparing the seq array by converting the native bam format to a character string
		bamseq = bam_get_seq(data->bambuf[i]); // Pointer to the start of the sequence
		seqlength = data->bambuf[i]->core.l_qseq; // Length of the sequence
		if(seqlength > data->max_length) {
			fprintf(stderr, "Read length %d greater than maximum permissible value (%d)\n", seqlength, data->max_length);
			exit(1);
		}

		// Filling the seq memory segment with the character sequence
		bamseq_to_char(bamseq, data->seq[i], seqlength);
	}
}

int parse_args(int argc, char* argv[], params *args) {

	int longindex;
	char param;

	char *help = "\nThis program was copied from the kather.\n" 
		"Usage: kmeria kbam -i <aligned.bam> -k <signalkmers.txt>\n"
		"\t-i, --input:     Input aligned bam files (required)\n"
		"\t-k, --kmers:     List of k-mers to look for in the reads (required)\n"
		"\t-o, --output:    Name of the output file in bam format (default: stdout)\n"
		"\t-t, --threads:   Number of analysis threads used in addition to decompression thread (default: 1)\n"
		"\t-b, --bufsize:   Number of records stored in the buffers (default: 100,000)\n"
		"\t-m, --maxlength: Maximum read length allowed, used for memory allocation (default: 300)\n"
		"\t-h, --help:      Print this help message\n";

	// Using getopt_long for parsing the arguments
	struct option long_options[] = {
		{"input",     required_argument, 0,  'i'},
		{"kmers",     required_argument, 0,  'k'},
		{"output",    required_argument, 0,  'o'},
		{"threads",   required_argument, 0,  't'},
		{"bufsize",   required_argument, 0,  'b'},
		{"maxlength", required_argument, 0,  'm'},
		{"help",      no_argument,       0,  'h'},
		{0,           0,                 0,   0 }
	};

	while((param = getopt_long(argc, argv, "i:k:o:t:b:m:h", long_options, &longindex)) != -1) {
		switch(param) {
			case 'i':
				args->input_file = strdup(optarg);
				break;
			case 'k':
				args->kmer_list = strdup(optarg);
				break;
			case 'o':
				args->output_file = strdup(optarg);
				break;
			case 't':
				args->n_threads = atoi(optarg);
				break;
			case 'b':
				args->bufsize = atoi(optarg);
				break;
			case 'm':
				args->seq_alloc = atoi(optarg);
				break;
			case 'h':
			case '?':
				fprintf(stderr, help, argv[0]);
				return 1;
				break;
			default:
				fprintf(stderr, "Unknown option, check usage by running %s --help", argv[0]);
		}
	}

	// Checking that all the required parameters are present
	if(args->input_file == NULL) {
		fprintf(stderr, "ERROR: Input .bam file must be specified as argument -i or --input\n\n");
		fprintf(stderr, help, argv[0]);
		return 1;
	}

	if(args->kmer_list == NULL) {
		fprintf(stderr, "ERROR: The list of k-mers must be specified as argument -k or --kmers\n");
		fprintf(stderr, help, argv[0]);
		return 1;
	}

	// Checking that the bufsize is a multiple of the number of threads
	if(args->bufsize % args->n_threads != 0) {
		fprintf(stderr, "ERROR: Buffer size (-b) must be a multiple of the number of threads (-t)\n");
		fprintf(stderr, help, argv[0]);
		return 1;
	}

	return 0;
}

void init_buffers(bam1_t **bambuf, bam1_t **tmpbuf, char **seq, int bufsize, int max_length) {

        for(int i = 0; i < bufsize; i++) {
                bambuf[i] = bam_init1();
                tmpbuf[i] = bam_init1();
                seq[i] = (char*) malloc((max_length + 1) * sizeof(char));

		if(seq[i] == NULL) {
			fprintf(stderr, "Unable to allocate memory for read sequences. Aborting.\n");
			exit(1);
		}
        }
}

// A function to read a list of k-mers from a text file
char **read_kmers(FILE *kmer_list, int n_init, int kmer_length, int *n_kmers) {
	// Creating an array of strings to store the k-mers as they are read
	char **kmers = (char**) malloc(n_init * sizeof(char*));
	int array_size = n_init;

	for(int i = 0; i < array_size; i++) kmers[i] = (char*) malloc((kmer_length + 1) * sizeof(char));

	// The number of records read from file
	*n_kmers = 0;

	while(fgets(kmers[*n_kmers], kmer_length + 2, kmer_list) != NULL) {
		// Removing the newline character from the string	
		kmers[*n_kmers][strcspn(kmers[*n_kmers], "\n")] = '\0';

		// Checking that the length of the string is kmer_length
		if(strlen(kmers[*n_kmers]) != kmer_length) {
			fprintf(stderr, "The k-mer length parameter is %d but the length of k-mer %s is %ld. Aborting.\n",
					kmer_length, kmers[*n_kmers], strlen(kmers[*n_kmers]));
			exit(1);
		}
		
		// Computing the reverse complement of the k-mer and putting it in the next slot in the array
		revcomp(kmers[*n_kmers], kmers[*n_kmers + 1], kmer_length);

		*n_kmers += 2;

		// Checking that the current array can hold at least the next two k-mers (forward and reverse)
		if((*n_kmers + 1) >= array_size) {
			array_size *= 2;
			kmers = (char**) realloc(kmers, array_size * sizeof(char*));
			for(int i = *n_kmers; i < array_size; i++) kmers[i] = (char*) malloc((kmer_length + 1) * sizeof(char));
		}
	}

	return kmers;
}



