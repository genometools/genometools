/*
  Copyright (c) 2015 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2015 Center for Bioinformatics, University of Hamburg

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

#ifndef SEED_EXTEND_ITER_H
#define SEED_EXTEND_ITER_H

#include "core/str_api.h"
#include "core/types_api.h"
#include "core/error_api.h"

/* This is the class name for the iterator on matches in the output format
   used by gt seed_extend and gt repfind. */

typedef struct GtSeedextendMatchIterator GtSeedextendMatchIterator;

/* The constructor takes are match file, tries to open it
   and delivers the iterator which allows to iterate over
   the matches. In case of an error, a NULL pointer is returned */

GtSeedextendMatchIterator *gt_seedextend_match_iterator_new(
                                            const GtStr *matchfilename,
                                            GtError *err);

/* The destructor removes an iterator object */

void gt_seedextend_match_iterator_delete(GtSeedextendMatchIterator *semi);

/* This function reads the next match and returns a <GtQuerymatch>-object.
   If there is no match left, then a NULL-ptr is returned. */

GtQuerymatch *gt_seedextend_match_iterator_next(
                             GtSeedextendMatchIterator *semi);

/* The following function reads all matches into an arrays and sorts the,. If
   <ascending is true, then all matches are sorted in ascending order of
   the query position they occur at. Otherwise, all matches are sorted in
   descending order of the query position they occur at. */

GtUword gt_seedextend_match_iterator_all_sorted(
                                         GtSeedextendMatchIterator *semi,
                                         bool ascending);

/* If the previous function has been called, the matches are stored in a table
   (in sorted order) and the following function allows to obtain the
   <querymatch> stored at a given index in this table */

GtQuerymatch *gt_seedextend_match_iterator_get(
                            const GtSeedextendMatchIterator *semi,
                            GtUword idx);

/* The following function sets some option related to the output of
   the matches and the corresponding alignments. If only the edit operation
   list is required, then set <generatealignment> to <true>,
   <alignmentwidth> to 0 and the other three boolean parameters to <false>. */

int gt_seedextend_match_iterator_querymatchoutoptions_set(
                    GtSeedextendMatchIterator *semi,
                    bool always_polished_ends,
                    GtExtendCharAccess a_extend_char_access,
                    GtExtendCharAccess b_extend_char_access,
                    const GtSeedExtendDisplayFlag *out_display_flag,
                    GtError *err);

/* The following function return different components of the iterator
   object. */

const GtEncseq *gt_seedextend_match_iterator_aencseq(
                        const GtSeedextendMatchIterator *semi);

const GtEncseq *gt_seedextend_match_iterator_bencseq(
                        const GtSeedextendMatchIterator *semi);

GtUword gt_seedextend_match_iterator_history_size(
                        const GtSeedextendMatchIterator *semi);

GtUword gt_seedextend_match_iterator_errorpercentage(
                        const GtSeedextendMatchIterator *semi);

bool gt_seedextend_match_iterator_bias_parameters(
                        const GtSeedextendMatchIterator *semi);

bool gt_seedextend_match_iterator_has_seed(
                        const GtSeedextendMatchIterator *semi);

bool gt_seedextend_match_iterator_has_cigar(
                        const GtSeedextendMatchIterator *semi);

bool gt_seedextend_match_iterator_has_subjectid(
                        const GtSeedextendMatchIterator *semi);

bool gt_seedextend_match_iterator_has_queryid(
                        const GtSeedextendMatchIterator *semi);

bool gt_seedextend_match_iterator_has_seqnums(
                        const GtSeedextendMatchIterator *semi);

GtUword gt_seedextend_match_iterator_trace_delta(
                        const GtSeedextendMatchIterator *semi);

bool gt_seedextend_match_iterator_dtrace(const GtSeedextendMatchIterator *semi);

double gt_seedextend_match_iterator_evalue(const GtSeedextendMatchIterator
                                             *semi);

double gt_seedextend_match_iterator_bitscore(
                      const GtSeedextendMatchIterator *semi);

const char *gt_seedextend_match_iterator_Options_line(
                    const GtSeedextendMatchIterator *semi);

void gt_seedextend_match_iterator_verify_alignment_set(
                                  GtSeedextendMatchIterator *semi);

#endif
