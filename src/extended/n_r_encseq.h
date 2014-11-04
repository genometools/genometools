/*
  Copyright (c) 2014 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
  Copyright (c) 2014 Florian Markowsky <1markows@informatik.uni-hamburg.de>
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
#include "core/file_api.h"

#define GT_NRENCSEQ_MIN_KMER_POS 5
#define GT_NRENCSEQ_FILE_SUFFIX ".nre"

/* The <GtNREncseq> class efficiently stores Sequences, either DNA or protein by
   finding redundancies and storing only references with editscripts of those.
   */
typedef struct GtNREncseq GtNREncseq;

/* Return new <GtNREncseq> object filled with data read from file with name
   <basename_nre>.
   Returns <NULL> on error, or fails if read/file access failes. */
GtNREncseq* gt_n_r_encseq_new_from_file(const char *basename_nre,
                                        GtLogger *logger, GtError *err);

/* Write <nre> to File <fp>, fails hard on io-error, returns error value and
   sets <err> accordingly on data errors. */
int         gt_n_r_encseq_write(GtNREncseq *nre, FILE* fp, GtError *err);

/* Prints statistical info of <n_r_encseq> to stdout.
   TODO DW refine this, and maybe split it up to get the different values
   seperately without printing. */
void        gt_n_r_encseq_print_info(const GtNREncseq *n_r_encseq);

/* Return the total length of the original sequence collection, including
   seperators (like <GtEncseq>). */
GtUword     gt_n_r_encseq_get_orig_length(GtNREncseq *n_r_encseq);

/* Return the length of unique db, that is the length of non-redundand sequences
   stored concatenated with seperators, representing the reference for all
   redundand sequences. */
GtUword     gt_n_r_encseq_get_unique_length(GtNREncseq *n_r_encseq);

/* Free space for <n_r_encseq> */
void        gt_n_r_encseq_delete(GtNREncseq *n_r_encseq);

/* The <GtNREncseqCompressor> class is used to compress an encoded sequence
   further by compressing similar sequences. It is used to create <GtNREncseq>
   objects and store them to disk. */
typedef struct GtNREncseqCompressor GtNREncseqCompressor;

/* Returns a new <GtNREncseqCompressor> object. */
GtNREncseqCompressor* gt_n_r_encseq_compressor_new(
                                                 GtUword initsize,
                                                 GtUword minalignlength,
                                                 GtUword max_kmer_poss,
                                                 GtWord xdropscore,
                                                 GtXdropArbitraryscores *scores,
                                                 unsigned int kmersize,
                                                 unsigned int windowsize,
                                                 GtLogger *logger);

/* Free memory of <n_r_e_compressor>. */
void                  gt_n_r_encseq_compressor_delete(
                                        GtNREncseqCompressor *n_r_e_compressor);

/* Analyze and compress <encseq>, stores resulting <GtNREncseq> to disk, using
   <basename> and <GT_NRENCSEQ_FILE_SUFFIX> as filename.
   Provide <logger> for verbose output. */
/* Due to change soon!
   TODO DW don't create encseq directly, call this iteratively, add finalize FKT
   */
int                   gt_n_r_encseq_compressor_compress(
                                         GtNREncseqCompressor *n_r_e_compressor,
                                         GtStr *basename,
                                         GtEncseq *encseq,
                                         GtLogger *logger,
                                         GtError *err);

/* This option turns of optimized seed extension. Every seed is used for xdrop
   without filtering.
   Consider this to be a development option for benchmarking purpose. Option
   incompatible with diagonals option! */
void                  gt_n_r_encseq_compressor_disable_opt(
                                        GtNREncseqCompressor *n_r_e_compressor);

/* Enable fast filtering for kmer hit pairs using diagonals and last hits on
   said diagonals. Option incompatible with disable optimization option. */
void                  gt_n_r_encseq_compressor_enable_diagonal_filter(
                                        GtNREncseqCompressor *n_r_e_compressor);

/* The <GtNREncseqDecompressor> class is used to decompress parts of or whole
   <GtEncseq> objects.
   TODO DW integrate this directly into <GtNREncseq>. */
typedef struct GtNREncseqDecompressor GtNREncseqDecompressor;

/* Return new <GtNREncseqCompressor> object ready to decompress <nre>. */
GtNREncseqDecompressor *gt_n_r_encseq_decompressor_new(GtNREncseq *nre);

/* Free memory of <n_r_e_decompressor>. */
void                    gt_n_r_encseq_decompressor_delete(
                                    GtNREncseqDecompressor *n_r_e_decompressor);

/* Write positions defined by <range> (positions correspond to uncompressed
   coordinates) to <fp>, either raw sequence (seperated by '|') or, if <fasta>
   is TRUE, as valid fasta entries.
   TODO DW split to two functions. */
int                     gt_n_r_encseq_decompressor_extract_originrange(
                                                   GtFile* fp,
                                                   GtNREncseqDecompressor *nred,
                                                   GtRange *range,
                                                   bool fasta,
                                                   GtError *err);

/* Write complete uncompressed sequence to <fp>. Writes valid fasta if <fasta>
   is true, raw sequences seperated by '|' otherwise.
   TODO DW split to two functions. */
int                     gt_n_r_encseq_decompressor_extract_origin_complete(
                                                   GtFile* fp,
                                                   GtNREncseqDecompressor *nred,
                                                   bool fasta,
                                                   GtError *err);

/* Adds the index of the unique entry <uentry_id> to sequences that are to be
   extracted. */
void                    gt_n_r_encseq_decompressor_add_unique_idx_to_extract(
                                                   GtNREncseqDecompressor *nred,
                                                   GtUword uentry_id);

/* Starts extraction of all sequences added by
   <gt_n_r_decompressor_add_unique_range_to_extract()> and stores number of
   extracted symbols excluding header descriptions in chars_printed.
   returns != 0 on error. */
int                     gt_n_r_encseq_decompressor_start_unique_extraction(
                                                   GtFile *fp,
                                                   GtNREncseqDecompressor *nred,
                                                   GtUword *chars_printed,
                                                   GtError *err);

int gt_n_r_encseq_unit_test(GT_UNUSED GtError *err);

#endif
