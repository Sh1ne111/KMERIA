#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../include/khash.h"
#include "../include/kstring.h"
#include "../include/sam.h"

#include "../include/functions.h"

#define BUFSIZE 5000
#define N_FIELDS_INIT 15
#define PV_STRING_INIT 100

// Initializing a hash table that links k-mers to their p-value
KHASH_MAP_INIT_STR(ktable, char*);

void read_pvalues(char *input, khash_t(ktable) *hash_table);

// A program that takes a .bam file with KM tags and adds a PV tag with the p-values associated
// with the k-mers. The program also needs a file with the k-mers and their associated p-values as input

// Usage: ./kmeria addp <input.bam> <pvalues.txt> > <output.bam>

int kaddp(int argc, char *argv[]) {

	// Initializing simple variables
	char *cl_tag, *pfile = argv[2], *token, *pvalue;

	// Initializing klib-related variables
	kstring_t km_string = KS_INITIALIZE;
	khash_t(ktable) *hash_table;
	hash_table = kh_init(ktable);
	khiter_t k_iterator;

	// Initializing htslib-related variables
	samFile   *input  = hts_open(argv[1], "r");
	samFile   *output = hts_open("-", "wb");
	bam1_t    *alignment = bam_init1();
	sam_hdr_t *header = sam_hdr_init();

	// Reading the header and outputting it to file
	header = sam_hdr_read(input);
	if(header == NULL) {
		fprintf(stderr, "Error reading header from %s. Aborting.\n", argv[1]);
		exit(1);
	}

	cl_tag = stringify_argv(argc, argv);
	if(sam_hdr_add_pg(header, "add_pvalues", "CL", cl_tag, NULL) != 0) {
		fprintf(stderr, "Error adding @PG line to header. Aborting.\n");
	}

	if(sam_hdr_write(output, header) != 0) {
		fprintf(stderr, "Error writing header to file. Aborting.\n");
	}
	
	// Reading the p-values and k-mers from the file so as to generate the hash table
	read_pvalues(argv[2], hash_table);

	/* DEBUG : printing all the keys and values in the hash table to see if it worked
	for(k_iterator = kh_begin(); k_iterator != kh_end(hash_table); k_iterator++) {
		fprintf(stderr, "%s: %s\n", kh_key(hash_table, k_iterator), kh_value(hash_table, k_iterator));
	}
	*/

	// Looping over the alignments
	while(sam_read1(input, header, alignment) >= 0) {
		// Initializing the p-value string to an empty string with allocated size PV_STRING_INIT
		size_t pv_size = PV_STRING_INIT;
		char *pv_string = (char*) malloc(pv_size * sizeof(char));
		pv_string[0] = '\0';

		// Extracting the KM tag data
		bam_aux_get_str(alignment, "KM", &km_string);

		// Putting the NULL character at the end of string in km_string and stripping the first few chars (KM:Z:)
		km_string.s[ks_len(&km_string)] = '\0';
		km_string.s = km_string.s + 5;

		// Looping over the elements of the KM string to find their p-value
		token = strtok(km_string.s, ",");

		while(token != NULL) {
			k_iterator = kh_get(ktable, hash_table, token);

			// The key should always be found in the hash table
			if(k_iterator == kh_end(hash_table)) {
				fprintf(stderr, "ERROR: k-mer %s not found in hash table. Aborting.\n", token);
			}

			pvalue = kh_value(hash_table, k_iterator);

			// Now adding that p-value to the PV string
			if(strlen(pv_string) + strlen(pvalue) + 2 >= pv_size) {
				pv_size *= 2;
				pv_string = (char*) realloc(pv_string, pv_size * sizeof(char));
			}

			if(strlen(pv_string)) strcat(pv_string, ",");

			strcat(pv_string, pvalue);
			token = strtok(NULL, ",");
		}

		// Adding the tag to the record
		if(bam_aux_append(alignment, "PV", 'Z', strlen(pv_string) + 1, pv_string) != 0) {
			fprintf(stderr, "ERROR: unable to add PV tag to BAM record. Aborting.");
		}

		// Writing the record
		if(sam_write1(output, header, alignment) < 0) {
			fprintf(stderr, "ERROR: Unable to write to output file. Aborting.\n");
			exit(1);
		}

		// Reset the km_string and the p-value string for the next iteration
		ks_release(&km_string);
		free(pv_string);
	}

	// Closing files
	sam_close(input);
	sam_close(output);

	// Freeing memory
	for (k_iterator = kh_begin(hash_table); k_iterator != kh_end(hash_table); k_iterator++)
		if (kh_exist(hash_table, k_iterator)) free((char*)kh_key(hash_table, k_iterator));
	kh_destroy(ktable, hash_table);

	return 0;
}

// A function that reads the file containing the p-values and updates the
// hash table that links the k-mers to their p-values (double precision)
void read_pvalues(char *input, khash_t(ktable) *hash_table) {
	// Opening the file
	FILE *pfile = fopen(input, "r");

	// Creating simple variables
	char buffer[BUFSIZE], **fields, *token, *rev_kmer = (char*) malloc((KMER_LENGTH + 1) * sizeof(char));
	int kmer_field, pvalue_field, n_fields, kh_return_val;
	int n_col = 0, k_match = 0, p_match = 0, row = 1, field_size = N_FIELDS_INIT;

	// Initializing an iterator for the hash table
	khiter_t k_iterator;

	// Initializing an array of strings for the fields
	fields = (char**) malloc(field_size * sizeof(char*));

	// Reading the first line of the file
	if(fgets(buffer, BUFSIZE, pfile) == NULL) {
		fprintf(stderr, "Unable to read input file %s. Aborting\n", input);
		exit(1);
	}

	// Removing the newline character from the input string
	buffer[strcspn(buffer, "\n")] = '\0';

	// Looping over the fields in the first column to find the column that contains the p-values
	token = strtok(buffer, "\t");

	while(token != NULL) {

		// Adding the field to the fields array
		if(n_col + 1 >= field_size) {
			field_size *= 2;
			fields = (char**) realloc(fields, field_size * sizeof(char*));
		}

		if(strcmp(token, "rs") == 0) {
			k_match++;
			kmer_field = n_col;
		}

		if(strcmp(token, "p_wald") == 0) {
			p_match++;
			pvalue_field = n_col;
		}

		token = strtok(NULL, "\t");

		n_col++;
	}

	// Checking that one and only one column matched for k-mers
	if(k_match == 0) {
		fprintf(stderr, "ERROR: k-mer column name rs not found in %s.\n", input);
		exit(1);
	}

	if(k_match > 1) {
		fprintf(stderr, "ERROR: k-mer column name rs matches more than one column name in %s.\n", input);
		exit(1);
	}

	// Checking that one and only one column matched for p-values
	if(p_match == 0) {
		fprintf(stderr, "ERROR: p-value column name p_wald not found in %s.\n", input);
		exit(1);
	}

	if(p_match > 1) {
		fprintf(stderr, "ERROR: p-value column name p_wald matches more than one column name in %s.\n", input);
		exit(1);
	}

	// Looping over the rows containing the k-mers and their p-values
	while(fgets(buffer, BUFSIZE, pfile) != NULL) {

		// Resetting the counter for the number of fields read in the row
		n_fields = 0;

		// Removing the newline character from the input string
		buffer[strcspn(buffer, "\n")] = '\0';

		// Looping over the fields in each row to add the data to the hash table
		token = strtok(buffer, "\t");

		while(token != NULL) {

			// Checking that we are not exceeding the number of columns read in the first row
			if(n_fields == n_col) {
				fprintf(stderr, "ERROR: Row number %d contains more fields than expected (%d).\n", row, n_col);
			}

			// Adding the k-mer and its p-value to the array
			fields[n_fields] = strdup(token);
			token = strtok(NULL, "\t");
			n_fields++;
		}

		// Checking that we stopped at the right number of fields
		if(n_fields != n_col) {
			fprintf(stderr, "ERROR: Row number %d contains %d fields whereas the expected number is %d.\n", row, n_fields, n_col);
		}

		// Removing the part after (and including) the underscore in the k-mer
		fields[kmer_field][strcspn(fields[kmer_field], "_")] = '\0';

		k_iterator = kh_put(ktable, hash_table, strdup(fields[kmer_field]), &kh_return_val);
		kh_value(hash_table, k_iterator) = strdup(fields[pvalue_field]);

		// Adding also the reverse complement with the same p-value
		revcomp(fields[kmer_field], rev_kmer, KMER_LENGTH);

		k_iterator = kh_put(ktable, hash_table, strdup(rev_kmer), &kh_return_val);
		kh_value(hash_table, k_iterator) = strdup(fields[pvalue_field]);

		row++;
	}

	// Closing files and freeing memory
	fclose(pfile);

	for(int i = 0; i < n_fields; i++) {
		free(fields[i]);
	}

	free(fields);
}

