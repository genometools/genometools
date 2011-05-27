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

#include "core/undef_api.h"
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

typedef struct GtPBSHit GtPBSHit;
typedef struct GtPBSResults GtPBSResults;

GtPBSResults*  gt_pbs_find(const char *seq,
                           const char *rev_seq,
                           GtLTRElement *element,
                           GtPBSOptions *o,
                           GtError *err);

GtRange        gt_pbs_hit_get_coords(const GtPBSHit*);
GtStrand       gt_pbs_hit_get_strand(const GtPBSHit*);
double         gt_pbs_hit_get_score(const GtPBSHit*);
unsigned long  gt_pbs_hit_get_edist(const GtPBSHit*);
unsigned long  gt_pbs_hit_get_offset(const GtPBSHit*);
unsigned long  gt_pbs_hit_get_tstart(const GtPBSHit*);
const char*    gt_pbs_hit_get_trna(const GtPBSHit*);
unsigned long  gt_pbs_hit_get_alignment_length(const GtPBSHit*);

unsigned long  gt_pbs_results_get_number_of_hits(const GtPBSResults*);
GtPBSHit*      gt_pbs_results_get_ranked_hit(const GtPBSResults*,
                                             unsigned long);
void           gt_pbs_results_delete(GtPBSResults*);

int            gt_pbs_unit_test(GtError*);

#endif
