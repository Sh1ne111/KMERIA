#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>

#include <pthread.h>

#include "../include/sam.h"
#include "../include/khash.h"

#include "../include/functions.h"

// This program takes a bam file and a list of k-mers as input, and outputs all
// reads containing these kmers or their reverse complement in BAM format
// The k-mer list must be a text file with one k-mer per line

int kbam(int argc, char* argv[]) {

	// Parsing the command-line arguments
	params args = {NULL, NULL, "-", 1, 100000, 300};
	if(parse_args(argc, argv, &args) != 0) return 1;

	// Declaring simple variables
	int n_read, records_per_thread, executing_threads, khash_return, n_kmers;
	long long int n_processed = 0;
	char **all_kmers, **seq, *cl_tag;
	FILE *kmer_list = NULL;

	// Declaring structs and struct arrays containing data passed to the threads for multithreading
	thread_data *thread_chunks; // Data for analysis threads
	buffer_data bufdata; // Data for decompression thread
	copy_data *copy_chunks; // Data for threads that copy data between buffers

	// Declaring the array of threads and a thread for decompressing the file
	pthread_t *threads = NULL, bufthread;

	// Declaring and initializing the hash table that will store the hashed values of the kmers, and its end iterator hash_end
	khash_t(kmer_hash) *kmer_table;
	kmer_table = kh_init(kmer_hash);
	khint_t hash_end;

	// Declaring and initializing htslib-related files and variables
        bam1_t **bambuf, **tmpbuf;

	samFile *input = sam_open(args.input_file, "r");
	if(input == NULL) {
		fprintf(stderr, "Error opening file %s for reading. Aborting.\n", args.input_file);
		exit(1);
	}

	samFile *output = hts_open(args.output_file, "wb");
	if(output == NULL) {
		fprintf(stderr, "Error opening file %s for writing. Aborting.\n", args.output_file);
		exit(1);
	}

	sam_hdr_t *header = sam_hdr_init();
	header = sam_hdr_read(input);
	if(header == NULL) {
		fprintf(stderr, "Error reading header from file %s. Aborting.\n", args.input_file);
		exit(1);
	}

	// Adding a PG line to the header and writing it to the output file
	cl_tag = stringify_argv(argc, argv);
	if(cl_tag == NULL) {
		fprintf(stderr, "Unable to generate string representation of command line. Aborting.\n");
		exit(1);
	}

	if(sam_hdr_add_pg(header, "katcher", "VN", KATCHER_VERSION, "CL", cl_tag, NULL) != 0) {
		fprintf(stderr, "Unable to add PG line to header. Aborting.\n");
		exit(1);
	}

	if(sam_hdr_write(output, header) != 0) {
		fprintf(stderr, "Error writing header to output file. Aborting.\n");
		exit(1);
	}

        // Initializing buffers for alignment records and sequence
        bambuf = (bam1_t**) malloc(args.bufsize * sizeof(bam1_t*));
	assert(bambuf != NULL);
        tmpbuf = (bam1_t**) malloc(args.bufsize * sizeof(bam1_t*));
	assert(tmpbuf != NULL);
	seq = (char**) malloc(args.bufsize * sizeof(char*));
	assert(seq != NULL);

	init_buffers(bambuf, tmpbuf, seq, args.bufsize, args.seq_alloc);

	// Opening the file with the k-mer list
	kmer_list = fopen(args.kmer_list, "r");

	// Reading the k-mers from file
	all_kmers = read_kmers(kmer_list, N_KMERS_INIT, KMER_LENGTH, &n_kmers);
	fprintf(stderr, "Read %d k-mers from file %s.\n", n_kmers, args.kmer_list);

	// Adding all k-mers to the hash table
	for(int i = 0; i < n_kmers; i++) kh_put(kmer_hash, kmer_table, all_kmers[i], &khash_return);

	// Getting the end position of the hash table
	hash_end = kh_end(kmer_table);

	// Creating an array of thread_data structs with as many elements as there are threads
	// The copy chunks which are meant for copying the data from tmpbuf to bambuf wiht multithreading are also processed here
	thread_chunks = (thread_data*) malloc(args.n_threads * sizeof(thread_data));
	copy_chunks = (copy_data*) malloc(args.n_threads * sizeof(copy_data));
	records_per_thread = args.bufsize / args.n_threads;

	// Initializing the chunk elements
	for(int i = 0; i < args.n_threads; i++) {
		thread_chunks[i].kmer_table = kmer_table;
		thread_chunks[i].hash_end = hash_end;
		thread_chunks[i].output = output;
		thread_chunks[i].header = header;
		thread_chunks[i].bambuf = bambuf + i * records_per_thread;
		thread_chunks[i].seq = seq + i * records_per_thread;
		thread_chunks[i].n_records = records_per_thread; // this value might change after filling the buffer

		copy_chunks[i].tmpbuf = tmpbuf + i * records_per_thread;
		copy_chunks[i].bambuf = bambuf + i * records_per_thread;
		copy_chunks[i].seq = seq + i * records_per_thread;
		copy_chunks[i].n_records = records_per_thread; // we won't care about adjusting this at the end of file for now
		copy_chunks[i].max_length = args.seq_alloc;
	}

	// Creating the array of analysis threads
	threads = (pthread_t*) malloc(args.n_threads * sizeof(pthread_t));

	// Initializing the lock
	if(pthread_mutex_init(&lock, NULL) != 0) {
		fprintf(stderr, "mutex init has failed\n");
		exit(1);
	}

	// Filling the temporary buffer by preparing the data and running the associated thread
	bufdata.bamfile = input;
	bufdata.header = header;
	bufdata.tmpbuf = tmpbuf;
	bufdata.max_read = args.bufsize;

	pthread_create(&bufthread, NULL, fillbuf, &bufdata);
	pthread_join(bufthread, NULL);

	n_read = bufdata.n_read;

	// Loop over the alignments in the bam file
        while(n_read > 0) {

		// Copying data from the temporary buffer to the one used for analysis
		// With multithreading
		for(int i = 0; i < args.n_threads; i++) pthread_create(&threads[i], NULL, copy_buf, &copy_chunks[i]);
		for(int i = 0; i < args.n_threads; i++) pthread_join(threads[i], NULL);

		// We no longer need the data in tmpbuf so we can launch the decompression thread
		pthread_create(&bufthread, NULL, fillbuf, &bufdata);

		// Generate the character string representation of the sequences in the bam file
		for(int i = 0; i < args.n_threads; i++) pthread_create(&threads[i], NULL, prepare_seq, &copy_chunks[i]);
		for(int i = 0; i < args.n_threads; i++) pthread_join(threads[i], NULL);

		// If we filled the buffer then the number of executing threads is args.n_threads
		if(n_read == args.bufsize) {
			executing_threads = args.n_threads;
		} else {
			// Otherwise there will be a certain number of threads that process
			// the usual number of records, and the last thread will process
			// the remaining records
			if(n_read % records_per_thread == 0) {
				executing_threads = n_read / records_per_thread;
			} else {
				executing_threads = n_read / records_per_thread + 1;
			}

			fprintf(stderr, "Buffer did not fill completely.\n"
					"Remaining records processed by %d threads processing "
					"%d records and one thread processing %d records.\n",
					executing_threads - 1,
					records_per_thread, 
					n_read - (executing_threads - 1) * records_per_thread);

			thread_chunks[executing_threads - 1].n_records = n_read - (executing_threads - 1) * records_per_thread;
		}

		// Launching the analysis threads
		for(int i = 0; i < executing_threads; i++) pthread_create(&threads[i], NULL, match_kmers, &thread_chunks[i]);
		for(int i = 0; i < executing_threads; i++) pthread_join(threads[i], NULL);

		n_processed += n_read;
		if(n_processed % 100000 == 0) fprintf(stderr, "%lld reads processed using %d threads\n", n_processed, executing_threads);

		// We wait for the buffer to be filled before starting the loop again
		pthread_join(bufthread, NULL);
		n_read = bufdata.n_read;
	}

	// Freeing the file handles and memory allocated
	sam_close(input);
	sam_close(output);
	fclose(kmer_list);

	sam_hdr_destroy(header);

	for(int i = 0; i < args.bufsize; i++) {
		bam_destroy1(bambuf[i]);
		bam_destroy1(tmpbuf[i]);
		free(seq[i]);
	}

	for(int i = 0; i < n_kmers; i++) free(all_kmers[i]);

	free(cl_tag);
	free(bambuf);
	free(tmpbuf);
	free(seq);
	free(all_kmers);
	free(thread_chunks);

	kh_destroy(kmer_hash, kmer_table);

	return 0;
}

