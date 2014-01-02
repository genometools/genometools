/*
  Copyright (c) 2011 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2011 Center for Bioinformatics, University of Hamburg

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

#ifndef SFX_SHORTREADSORT_H
#define SFX_SHORTREADSORT_H

#include "core/readmode_api.h"
#include "core/encseq_api.h"
#include "core/unused_api.h"
#include "sfx-lcpvalues.h"
#include "sfx-suffixgetset.h"
#include "spmsuftab.h"
#include "seqnumrelpos.h"

typedef struct GtShortreadsortworkinfo GtShortreadsortworkinfo;

typedef struct
{
  GtUword *suftab_bucket;
  uint16_t *lcptab_bucket;
} GtShortreadsortresult;

size_t gt_shortreadsort_size(bool firstcodes,GtUword bucketsize,
                             GtUword maxremain);

GtShortreadsortworkinfo *gt_shortreadsort_new(GtUword maxwidth,
                                              GtUword maxremain,
                                              GtReadmode readmode,
                                              bool firstcodes,
                                              bool withmediumsizelcps);

GtUword gt_shortreadsort_sumofstoredvalues(
                                      const GtShortreadsortworkinfo *srsw);

void gt_shortreadsort_delete(GtShortreadsortworkinfo *srsw);

void gt_shortreadsort_assigntableoflcpvalues(
          GtShortreadsortworkinfo *srsw,GtLcpvalues *tableoflcpvalues);

void gt_shortreadsort_sssp_sort(GtShortreadsortworkinfo *srsw,
                                const GtEncseq *encseq,
                                GtUword maxremain,
                                GtReadmode readmode,
                                GtEncseqReader *esr,
                                GtSuffixsortspace *sssp,
                                GtUword subbucketleft,
                                GtUword width,
                                GtUword depth,
                                GtUword maxdepth);

void gt_shortreadsort_firstcodes_sort(GtShortreadsortresult *srsresult,
                                      GtShortreadsortworkinfo *srsw,
                                      const GtSeqnumrelpos *snrp,
                                      const GtEncseq *encseq,
                                      const GtSpmsuftab *spmsuftab,
                                      GtUword subbucketleft,
                                      GtUword width,
                                      GtUword depth,
                                      GtUword maxdepth);

void gt_shortreadsort_sssp_add_unsorted(const GtShortreadsortworkinfo *srsw,
                                        GtUword bucketleftidx,
                                        GtUword subbucketleft,
                                        GtUword width,
                                        GtUword maxdepth,
                                        GtSuffixsortspace *sssp,
                                        GtProcessunsortedsuffixrange
                                          processunsortedsuffixrange,
                                        void *processunsortedsuffixrangeinfo);

GtUword gt_shortreadsort_maxwidth(bool firstcodes,
                                        GtUword maxremain,
                                        size_t sizeofworkspace);

#endif
