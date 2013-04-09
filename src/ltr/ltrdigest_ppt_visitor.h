/*
  Copyright (c) 2013 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2013 Center for Bioinformatics, University of Hamburg

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

#ifndef LTRDIGEST_PPT_VISITOR_H
#define LTRDIGEST_PPT_VISITOR_H

#include "core/error_api.h"
#include "core/range_api.h"
#include "extended/node_visitor_api.h"
#include "extended/region_mapping_api.h"

#define PPT_PURINE_PROB     0.97
#define PPT_PYRIMIDINE_PROB 0.03
#define BKG_A_PROB          0.25
#define BKG_C_PROB          0.25
#define BKG_G_PROB          0.25
#define BKG_T_PROB          0.25
#define UBOX_U_PROB         0.91

/* Implements the <GtNodeVisitor> interface. */
typedef struct GtLTRdigestPPTVisitor GtLTRdigestPPTVisitor;

GtNodeVisitor* gt_ltrdigest_ppt_visitor_new(GtRegionMapping *rmap,
                                            GtRange ppt_len,
                                            GtRange ubox_len,
                                            double ppt_pyrimidine_prob,
                                            double ppt_purine_prob,
                                            double bkg_a_prob,
                                            double bkg_g_prob,
                                            double bkg_t_prob,
                                            double bkg_c_prob,
                                            double ubox_u_prob,
                                            unsigned int radius,
                                            unsigned int max_ubox_dist,
                                            GtError *err);

#endif
