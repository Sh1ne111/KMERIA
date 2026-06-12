 /*
  The MIT License

   Copyright (c) Chen Shuai, 2024- <chensss1209@gmail.com>

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
*/

#ifndef KMER_CONVERTER_H
#define KMER_CONVERTER_H

#include <stdio.h>
#include <stdint.h>
#include "kvec.h" 

typedef kvec_t(int) int_vector_t;

#define FORMAT_KMER_ONLY 0   // Output only the k-mer sequence.
#define FORMAT_ALL_COUNTS 1  // Output k-mer followed by all counts, comma-separated.
#define FORMAT_SUM 2         // Output k-mer and the sum of all its counts.
#define FORMAT_DIFFERENCE 3  // Output k-mer and the result of subtracting subsequent counts from the first.
#define FORMAT_MIN 4         // Output k-mer and the minimum count among all samples.
#define FORMAT_MAX 5         // Output k-mer and the maximum count among all samples.
#define FORMAT_TYPE_BOUNDARY FORMAT_MAX // Boundary for predefined format types.

/*
 * Converts a binary k-mer file to a Plain text format based on specified filters.
 *
 * @param instream          Input file stream (binary k-mer format).
 * @param outstream         Output file stream (text format).
 * @param min_sample_pass   The minimum number of samples that must meet the count criteria.
 * @param max_sample_pass   The maximum number of samples allowed to meet the count criteria (-1 for no limit).
 * @param min_count_filters A vector of minimum count thresholds for each sample.
 * @param max_count_filters A vector of maximum count thresholds for each sample (-1 for no limit).
 * @param output_format     The desired output format, using the FORMAT_* constants or a column index.
 *
 */
int convert_kmer_binary_to_text(FILE *instream, FILE *outstream, int min_sample_pass, int max_sample_pass,
                                int_vector_t *min_count_filters, int_vector_t *max_count_filters, int output_format);

int kdump_main(int argc, char *argv[]);

#endif
