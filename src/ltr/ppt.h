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

#ifndef PPT_H
#define PPT_H

#include "core/alpha.h"
#include "core/range.h"
#include "core/strand.h"
#include "core/undef.h"
#include "extended/hmm.h"
#include "ltr/ltrelement.h"

typedef struct {
  GtRange ppt_len, ubox_len;
  unsigned int radius;
} GtPPTOptions;

typedef struct GtPPTHit GtPPTHit;
typedef struct GtPPTResults GtPPTResults;

/* Searches for PPTs in the given sequence. */
GtPPTResults*   gt_ppt_find(const char *seq,
                            const char *rev_seq,
                            GtLTRElement *element,
                            GtPPTOptions*);

/* A PPT hit representation */
GtRange         gt_ppt_hit_get_coords(GtPPTHit*);
GtPPTHit*       gt_ppt_hit_get_ubox(GtPPTHit*);
GtStrand        gt_ppt_hit_get_strand(GtPPTHit*);

/* A collection of PPT hits */
unsigned long   gt_ppt_results_get_number_of_hits(GtPPTResults*);
GtPPTHit*       gt_ppt_results_get_ranked_hit(GtPPTResults*, unsigned long);
void            gt_ppt_results_delete(GtPPTResults*);

int             gt_ppt_unit_test(GtError*);

#endif
