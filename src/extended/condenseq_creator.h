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

#ifndef CONDENSEQ_CREATOR_H
#define CONDENSEQ_CREATOR_H

#include "core/types_api.h"
#include "match/xdrop.h"
#include "core/logger_api.h"
#include "extended/condenseq.h"

/* The <GtCondenseqCreator> class is used to create an encoded sequence
   further by compressing similar sequences. It is used to create <GtCondenseq>
   objects and store them to disk. */
typedef struct GtCondenseqCreator GtCondenseqCreator;

/* Returns a new <GtCondenseqCreator> object. Checks if minalignlength is to
   large to be stored. Returns NULL in that case and sets <err> accordingly. */
GtCondenseqCreator* gt_condenseq_creator_new(GtUword initsize,
                                             GtUword minalignlength,
                                             GtWord xdropscore,
                                             GtXdropArbitraryscores *scores,
                                             unsigned int kmersize,
                                             unsigned int windowsize,
                                             GtLogger *logger,
                                             GtError *err);
/* Free memory of <condenseq_creator>. */
void                gt_condenseq_creator_delete(
                                         GtCondenseqCreator *condenseq_creator);
/* Analyze and compress <encseq>, stores resulting <GtCondenseq> to disk, using
   <basename> and <GT_CONDENSEQ_FILE_SUFFIX as filename.
   Provide <logger> for verbose output. */
/* Due to change soon!
   TODO DW don't create encseq directly, call this iteratively, add finalize FKT
   */
int                 gt_condenseq_creator_create(
                                          GtCondenseqCreator *condenseq_creator,
                                          GtStr *basename,
                                          GtEncseq *encseq,
                                          GtLogger *logger,
                                          GtLogger *kdb_logger,
                                          GtError *err);
/* This option turns of optimized seed extension. Every seed is used for xdrop
   without filtering. And no cutoff is used for number of considered kmer
   positions.
   Consider this to be a development option for benchmarking purpose.
   Disables diagonals and filtered extension. */
void                gt_condenseq_creator_enable_brute_force(
                                         GtCondenseqCreator *condenseq_creator);
/* This option turns on optimized seed extension. Every seed is used for xdrop
   without filtering.
   Consider this to be a development option for benchmarking purpose.
   Disables diagonals. */
void                gt_condenseq_creator_enable_opt(
                                         GtCondenseqCreator *condenseq_creator);
/* Enable sparse diagonals. If both sparse and full diagonals are active they
   will be used in parallel and checked against each other. */
void                gt_condenseq_creator_enable_diagonals(
                                         GtCondenseqCreator *condenseq_creator);
/* Disable sparse diagonals, if full diagonals is inactive, too, it has the same
   effect as gt_condenseq_creator_enable_opt(). */
void                gt_condenseq_creator_disable_diagonals(
                                         GtCondenseqCreator *condenseq_creator);
/* Enable full diagonals. If both sparse and full diagonals are active they
   will be used in parallel and checked against each other. */
void                gt_condenseq_creator_enable_full_diagonals(
                                         GtCondenseqCreator *condenseq_creator);
/* Disable full diagonals, if sparse diagonals is inactive, too, it has the same
   effect as gt_condenseq_creator_enable_opt(). */
void                gt_condenseq_creator_disable_full_diagonals(
                                         GtCondenseqCreator *condenseq_creator);
/* If a kmer appears more often than <cutoff_value> in the unique data it won't
   be used to find new alignments, because it likely has a low chance to find
   new alignments. */
void                gt_condenseq_creator_set_cutoff(
                                          GtCondenseqCreator *condenseq_creator,
                                          GtUword cutoff_value);
/* If this option is used every kmer in the unique data will be used to find new
   alignments. */
void                gt_condenseq_creator_disable_cutoff(
                                         GtCondenseqCreator *condenseq_creator);
/* If this option is set the <cutoff_value> will be deduced by the current mean
   value of kmers in the unique data. */
void                gt_condenseq_creator_use_mean_cutoff(
                                         GtCondenseqCreator *condenseq_creator);
/* If this option is set only <cutoff_value> many kmers will be saved to find
   alignments. Only works when a cutoff is set. */
void                gt_condenseq_creator_disable_prune(
                                         GtCondenseqCreator *condenseq_creator);
/* This option specifies which fraction of the mean value of each kmer in the
   unique data is used to calculate a current <cutoff_value> (mean/<fraction>.*/
void                gt_condenseq_creator_set_mean_fraction(
                                          GtCondenseqCreator *condenseq_creator,
                                          GtUword fraction);
/* Percentage of sparse diagonals that is allowed to be outside of used ranges
   and marked for deletion. 0 <= <percent> < 100. */
void gt_condenseq_creator_set_diags_clean_limit(
                                          GtCondenseqCreator *condenseq_creator,
                                          unsigned int percent);
#endif
