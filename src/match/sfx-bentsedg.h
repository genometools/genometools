/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

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

#ifndef SFX_BENTSEDG_H
#define SFX_BENTSEDG_H
#include <stdio.h>
#include "core/error_api.h"
#include "core/str.h"
#include "core/codetype.h"
#include "core/logger_api.h"
#include "core/encseq.h"
#include "bcktab.h"
#include "sfx-bltrie.h"
#include "sfx-strategy.h"
#include "sfx-copysort.h"
#include "sfx-lcpvalues.h"
#include "sfx-suffixgetset.h"

typedef void (*GtCompletelargelcpvalues) (void *,
                                         const GtSuffixsortspace *,
                                         GtLcpvalues *,
                                         unsigned long,
                                         unsigned long);

void gt_sortallbuckets(GtSuffixsortspace *suffixsortspace,
                       unsigned long numberofsuffixes,
                       GtBucketspec2 *bucketspec2,
                       const GtEncseq *encseq,
                       GtReadmode readmode,
                       GtCodetype mincode,
                       GtCodetype maxcode,
                       GtBcktab *bcktab,
                       unsigned int numofchars,
                       unsigned int prefixlength,
                       GtOutlcpinfo *outlcpinfo,
                       unsigned int sortmaxdepth,
                       const Sfxstrategy *sfxstrategy,
                       GtProcessunsortedsuffixrange processunsortedsuffixrange,
                       void *processunsortedsuffixrangeinfo,
                       unsigned long long *bucketiterstep,
                       GtLogger *logger);

void gt_sortallsuffixesfromstart(GtSuffixsortspace *suffixsortspace,
                                 unsigned long numberofsuffixes,
                                 const GtEncseq *encseq,
                                 GtReadmode readmode,
                                 GtOutlcpinfo *outlcpinfo,
                                 unsigned int sortmaxdepth,
                                 const Sfxstrategy *sfxstrategy,
                                 GtProcessunsortedsuffixrange
                                   processunsortedsuffixrange,
                                 void *processunsortedsuffixrangeinfo,
                                 GtLogger *logger);

size_t gt_size_of_sort_workspace (const Sfxstrategy *sfxstrategy);

#endif
