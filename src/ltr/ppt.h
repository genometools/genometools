/*
  Copyright (c) 2008-2011 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2008-2011 Center for Bioinformatics, University of Hamburg

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

#include "core/alphabet.h"
#include "core/range.h"
#include "core/strand.h"
#include "core/undef_api.h"
#include "extended/hmm.h"
#include "ltr/ltrelement.h"

#define PPT_PURINE_PROB     0.97
#define PPT_PYRIMIDINE_PROB 0.03
#define BKG_A_PROB          0.25
#define BKG_C_PROB          0.25
#define BKG_G_PROB          0.25
#define BKG_T_PROB          0.25
#define UBOX_U_PROB         0.91

typedef struct {
  GtRange ppt_len, ubox_len;
  double ppt_pyrimidine_prob,
         ppt_purine_prob,
         bkg_a_prob,
         bkg_g_prob,
         bkg_t_prob,
         bkg_c_prob,
         ubox_u_prob;
  unsigned int radius,
               max_ubox_dist;
} GtPPTOptions;

typedef struct GtPPTHit GtPPTHit;
typedef struct GtPPTResults GtPPTResults;

GtHMM*          gt_ppt_hmm_new(const GtAlphabet *alpha, GtPPTOptions *opts);

GtPPTResults*   gt_ppt_find(const char *seq,
                            const char *rev_seq,
                            GtLTRElement *element,
                            GtPPTOptions*);

GtRange         gt_ppt_hit_get_coords(const GtPPTHit*);
GtPPTHit*       gt_ppt_hit_get_ubox(const GtPPTHit*);
GtStrand        gt_ppt_hit_get_strand(const GtPPTHit*);

unsigned long   gt_ppt_results_get_number_of_hits(GtPPTResults*);
GtPPTHit*       gt_ppt_results_get_ranked_hit(GtPPTResults*, unsigned long);
void            gt_ppt_results_delete(GtPPTResults*);

int             gt_ppt_unit_test(GtError*);

#endif
