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
#include "sfx-suffixgetset.h"
#include "suffixptr.h"

typedef struct Differencecover Differencecover;

/* The following function is used for test purposes only */

void gt_differencecovers_check(const GtEncseq *encseq,
                               GtReadmode readmode);

Differencecover *gt_differencecover_new(unsigned int vparam,
                                        const GtEncseq *encseq,
                                        GtReadmode readmode,
                                        GtLogger *logger);

int gt_differencecover_vparamverify(const Differencecover *dcov,GtError *err);

void gt_differencecover_sortsample(Differencecover *dcov,
                                   bool cmpcharbychar,
                                   bool withcheck);

void gt_differencecover_delete(Differencecover *dcov);

void dc_sortunsortedbucket(void *data,
                           Suffixptr *subbucket,
                           unsigned long subbucketleft,
                           unsigned long width,
                           GT_UNUSED unsigned long depth);

#endif
