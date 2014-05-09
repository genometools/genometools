/*
  Copyright (c) 2014 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
  Copyright (c) 2014 Center for Bioinformatics, University of Hamburg

  Permission to use, copy, modify, and distribute this software for any
  purpose with or without fee is hereby granted, provided that the above
  copyright notice and this permission notice appear in all copies.

  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/

#ifndef N_R_ENCSEQ_H
#define N_R_ENCSEQ_H

#include "core/encseq_api.h"
#include "core/error_api.h"
#include "core/logger_api.h"
#include "core/range_api.h"
#include "core/unused_api.h"
#include "match/xdrop.h"

#define GT_NRENCSEQ_MIN_KMER_POS 5
#define GT_NRENCSEQ_FILE_SUFFIX ".nre"

/* The <GtNREncseq> class stores redundand */
typedef struct GtNREncseq GtNREncseq;

/* Return new <GtNREncseq> object filed with data read from file with name
   <basename_nre>.
   Returns <NULL> on error, or fails if read/file access failes. */
GtNREncseq* gt_n_r_encseq_new_from_file(const char *basename_nre,
                                        const GtEncseq *orig_es,
                                        GtError *err);

/* Write <nre> to File <fp>, fails hard on io-error, returns error value and
   sets <err> accordingly on data errors. */
int         gt_n_r_encseq_write_nre(GtNREncseq *nre, FILE* fp, GtError *err);

/*prints statistical info of <n_r_encseq>*/
void        gt_n_r_encseq_print_info(const GtNREncseq *n_r_encseq);

/*get length of original sequence collection*/
GtUword     gt_n_r_encseq_get_orig_length(GtNREncseq *n_r_encseq);

/*get length of unique db*/
GtUword     gt_n_r_encseq_get_unique_length(GtNREncseq *n_r_encseq);

/* free space for <n_r_encseq> */
void        gt_n_r_encseq_delete(GtNREncseq *n_r_encseq);

/* The <GtNREncseqCompressor> class is used to compress an encoded sequence
   further by compressing similar sequences. */
typedef struct GtNREncseqCompressor GtNREncseqCompressor;

/* Returns a new <GtNREncseqCompressor> object */
GtNREncseqCompressor* gt_n_r_encseq_compressor_new(
                                                 GtUword initsize,
                                                 GtUword minalignlength,
                                                 GtUword max_kmer_poss,
                                                 GtWord xdropscore,
                                                 GtXdropArbitraryscores *scores,
                                                 unsigned int kmersize,
                                                 unsigned int windowsize,
                                                 GtLogger *logger);

/* use non optimized seed extension. every seed is used for xdrop without
   filtering */
void                  gt_n_r_encseq_compressor_disable_opt(
                                        GtNREncseqCompressor *n_r_e_compressor);

/* Free memory of <n_r_e_compressor> */
void                  gt_n_r_encseq_compressor_delete(
                                        GtNREncseqCompressor *n_r_e_compressor);

/* Start the compression of <encseq> and store the results on disk to files with
   <basename> as basename. */
int                   gt_n_r_encseq_compressor_compress(
                                         GtNREncseqCompressor *n_r_e_compressor,
                                         GtStr *basename,
                                         GtEncseq *encseq,
                                         GtError *err);

typedef struct GtNREncseqDecompressor GtNREncseqDecompressor;

GtNREncseqDecompressor *gt_n_r_encseq_decompressor_new(GtNREncseq *nre);
void                    gt_n_r_encseq_decompressor_delete(
                                    GtNREncseqDecompressor *n_r_e_decompressor);

/*Given a <range> in the original sequence database, this function writes the
  sequence fragments in the requested <range> to the file pointed to by <fp>,
  either as raw sequence or, if <fasta> is TRUE, as valid fasta entries*/
int                     gt_n_r_encseq_decompressor_extract_originrange(
                                                   FILE* fp,
                                                   GtNREncseqDecompressor *nred,
                                                   GtRange *range,
                                                   bool fasta,
                                                   GtError *err);

/*calls gt_n_r_encseq_decompressor_extract_originrange with range of whole
  original sequence*/
int                     gt_n_r_encseq_decompressor_extract_origin_complete(
                                                   FILE* fp,
                                                   GtNREncseqDecompressor *nred,
                                                   bool fasta,
                                                   GtError *err);
/*adds the index of the unique entry <uentry_id> to sequences
  that are to be extracted*/
void                    gt_n_r_encseq_decompressor_add_unique_idx_to_extract(
                                                   GtNREncseqDecompressor *nred,
                                                   GtUword uentry_id);

/*starts extraction of all sequences added by
  gt_n_r_decompressor_add_unique_range_to_extract and returns number of
  extracted symbols excluding header descriptions*/
GtUword                 gt_n_r_encseq_decompressor_start_unique_extraction(
                                                   FILE *fp,
                                                   GtNREncseqDecompressor *nred,
                                                   GtError *err);

int gt_n_r_encseq_unit_test(GT_UNUSED GtError *err);

#endif
