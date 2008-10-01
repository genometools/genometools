/*
  Copyright (c) 2008 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2008 Center for Bioinformatics, University of Hamburg

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

#ifndef PBS_H
#define PBS_H

#include "core/undef.h"
#include "core/bioseq.h"
#include "core/dlist.h"
#include "core/seq.h"
#include "core/strand.h"
#include "core/score_function.h"
#include "extended/alignment.h"
#include "ltr/ltrelement.h"

typedef struct GtPBSOptions {
  unsigned int radius,
               max_edist;
  GtRange alilen,
          offsetlen,
          trnaoffsetlen;
  int ali_score_match,
      ali_score_mismatch,
      ali_score_insertion,
      ali_score_deletion;
  GtBioseq *trna_lib;
} GtPBSOptions;

/* This struct holds information about a primer binding site (PBS). */
typedef struct GtPBS_Hit {
  unsigned long start,
                end,
                edist,
                offset,
                tstart,
                alilen;
  GtStrand strand;
  double score;
  const char *trna;
} GtPBS_Hit;

typedef struct GtPBSResults {
  GtDlist *hits_fwd, *hits_rev;
  GtPBS_Hit *best_hit;
} GtPBSResults;

/* Aligns tRNA from a library to the LTR retrotransposon candidate and
   returns highest-scoring hit (newly created). */
void  gt_pbs_find(const char *seq,
                  const char *rev_seq,
                  GtLTRElement *element,
                  GtPBSResults *results,
                  GtPBSOptions *o,
                  GtError *err);

void  gt_pbs_clear_results(GtPBSResults *results);

#endif
