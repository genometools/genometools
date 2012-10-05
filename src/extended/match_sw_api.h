/*
  Copyright (c) 2012 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2012 Center for Bioinformatics, University of Hamburg

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

#ifndef MATCH_SW_API_H
#define MATCH_SW_API_H

 /* The <GtMatchSW> class, implementing the <GtMatch> interface, is meant to
    store results from Smith-Waterman matching (using the <swalign> module). */
typedef struct GtMatchSW GtMatchSW;

#include "extended/match_api.h"

/* Creates a new <GtMatch> object, storing the alignment length <length>,
   the edit distance <edist> and the sequence numbers in the sequence
   collections in addition to the generic match contents <seqid1>, <seqid2>,
   <start_seq1>, <start_seq2>, <end_seq1> and <end_seq2>. */
GtMatch*      gt_match_sw_new(const char *seqid1,
                              const char *seqid2,
                              unsigned long seqno1,
                              unsigned long seqno2,
                              unsigned long length,
                              unsigned long edist,
                              unsigned long start_seq1,
                              unsigned long start_seq2,
                              unsigned long end_seq1,
                              unsigned long end_seq2,
                              GtMatchDirection dir);

/* Returns the sequence number of the match <ms> in the first sequence set. */
unsigned long gt_match_sw_get_seqno1(const GtMatchSW *ms);

/* Returns the sequence number of the match <ms> in the second sequence set. */
unsigned long gt_match_sw_get_seqno2(const GtMatchSW *ms);

/* Returns the alignment length of the match <ms>. */
unsigned long gt_match_sw_get_alignment_length(const GtMatchSW *ms);

/* Returns the edit distance of the match <ms>. */
unsigned long gt_match_sw_get_edist(const GtMatchSW *ms);

#endif
