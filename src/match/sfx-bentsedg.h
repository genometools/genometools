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
#include "core/defined-types.h"
#include "core/codetype.h"
#include "core/encseq.h"
#include "bcktab.h"
#include "sfx-strategy.h"
#include "sfx-copysort.h"
#include "sfx-bltrie.h"
#include "sfx-suffixgetset.h"

typedef struct Outlcpinfo Outlcpinfo;

Outlcpinfo *gt_Outlcpinfo_new(const char *indexname,
                              unsigned int numofchars,
                              unsigned int prefixlength,
                              unsigned long totallength,
                              bool assideeffect,
                              GtError *err);

void gt_Outlcpinfo_delete(Outlcpinfo *outlcpinfo);

unsigned long gt_Outlcpinfo_numoflargelcpvalues(const Outlcpinfo *outlcpinfo);

unsigned long gt_Outlcpinfo_maxbranchdepth(const Outlcpinfo *outlcpinfo);

void gt_qsufsort(GtSuffixsortspace *suffixsortspace,
                 unsigned long partwidth,
                 int mmapfiledesc,
                 GtStr *mmapfilename,
                 const GtEncseq *encseq,
                 GtReadmode readmode,
                 GtCodetype mincode,
                 GtCodetype maxcode,
                 Bcktab *bcktab,
                 unsigned int numofchars,
                 unsigned int prefixlength,
                 bool hashexceptions,
                 bool absoluteinversesuftab,
                 Outlcpinfo *outlcpinfo);

void gt_sortallbuckets(GtSuffixsortspace *suffixsortspace,
                       GtBucketspec2 *bucketspec2,
                       const GtEncseq *encseq,
                       GtReadmode readmode,
                       GtCodetype mincode,
                       GtCodetype maxcode,
                       unsigned long partwidth,
                       Bcktab *bcktab,
                       unsigned int numofchars,
                       unsigned int prefixlength,
                       Outlcpinfo *outlcpinfo,
                       const Sfxstrategy *sfxstrategy,
                       unsigned long long *bucketiterstep,
                       GtLogger *logger);

void gt_sortbucketofsuffixes(GtSuffixsortspace *suffixsortspace,
                             unsigned long numberofsuffixes,
                             GtBucketspec2 *bucketspec2,
                             const GtEncseq *encseq,
                             GtReadmode readmode,
                             GtCodetype mincode,
                             GtCodetype maxcode,
                             const Bcktab *bcktab,
                             unsigned int numofchars,
                             unsigned int prefixlength,
                             unsigned int sortmaxdepth,
                             const Sfxstrategy *sfxstrategy,
                             void *voiddcov,
                             Dc_processunsortedrange dc_processunsortedrange,
                             GtLogger *logger);

void gt_sortallsuffixesfromstart(GtSuffixsortspace *suffixsortspace,
                                 unsigned long numberofsuffixes,
                                 const GtEncseq *encseq,
                                 GtReadmode readmode,
                                 unsigned int numofchars,
                                 unsigned int sortmaxdepth,
                                 const Sfxstrategy *sfxstrategy,
                                 void *voiddcov,
                                 Dc_processunsortedrange
                                   dc_processunsortedrange,
                                 GtLogger *logger);

#endif
