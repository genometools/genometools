/*
  Copyright (c) 2009 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2009 Center for Bioinformatics, University of Hamburg

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

#ifndef SFX_DIFFCOV_H
#define SFX_DIFFCOV_H

#include "core/encseq.h"
#include "core/unused_api.h"
#include "core/readmode.h"
#include "core/logger.h"
#include "core/timer_api.h"
#include "core/error_api.h"
#include "sfx-strategy.h"
#include "sfx-lcpvalues.h"
#include "sfx-suffixgetset.h"

typedef struct GtDifferencecover GtDifferencecover;

/* The following function is used for test purposes only */

void gt_differencecover_check(const GtEncseq *encseq,
                               GtReadmode readmode);

GtDifferencecover *gt_differencecover_new(unsigned int vparam,
                                          const GtEncseq *encseq,
                                          GtReadmode readmode,
                                          unsigned int outerprefixlength,
                                          GtLogger *logger);

unsigned long gt_differencecover_samplesize(const GtDifferencecover *dcov);

GtDifferencecover *gt_differencecover_prepare_sample(
                                        unsigned int vparam,
                                        const GtEncseq *encseq,
                                        GtReadmode readmode,
                                        unsigned int prefixlength,
                                        const Sfxstrategy *sfxstrategy,
                                        GtOutlcpinfo *outlcpinfosample,
                                        GtLogger *logger,
                                        GtTimer *sfxprogress,
                                        GtError *err);

bool gt_differencecover_is_empty(const GtDifferencecover *dcov);

void gt_differencecover_delete(GtDifferencecover *dcov);

size_t gt_differencecover_requiredspace(const GtDifferencecover *dcov);

void gt_differencecover_sortunsortedbucket(void *data,
                                           unsigned long blisbl,
                                           unsigned long width,
                                           GT_UNUSED unsigned long depth);

void gt_differencecover_completelargelcpvalues(void *data,
                                               const GtSuffixsortspace *sssp,
                                               GtLcpvalues *tableoflcpvalues,
                                               unsigned long width,
                                               unsigned long posoffset);

void gt_differencecover_set_sssp_lcp(GtDifferencecover *dcov,
                                     GtSuffixsortspace *sssp,
                                     GtOutlcpinfo *outlcpinfo);

#endif
