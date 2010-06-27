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
#include "sfx-strategy.h"
#include "sfx-copysort.h"
#include "bcktab.h"
#include "suffixptr.h"

typedef struct Outlcpinfo Outlcpinfo;

typedef struct
{
  Suffixptr *sortspace;
  unsigned long suftaboffset;
} Suftab;

Outlcpinfo *gt_newOutlcpinfo(const char *indexname,
                             unsigned int numofchars,
                             unsigned int prefixlength,
                             unsigned long totallength,
                             bool assideeffect,
                             GtError *err);

void gt_freeOutlcptab(Outlcpinfo **outlcpinfoptr);

unsigned long getnumoflargelcpvalues(const Outlcpinfo *outlcpinfo);

unsigned long getmaxbranchdepth(const Outlcpinfo *outlcpinfo);

void gt_qsufsort(Suffixptr *leftptr,
                 unsigned long partwidth,
                 int mmapfiledesc,
                 GtStr *mmapfilename,
                 unsigned long *longest,
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

void gt_sortallbuckets(Suftab *suftab,
                       Definedunsignedlong *longest,
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

void gt_sortbucketofsuffixes(Suffixptr *suffixestobesorted,
                             unsigned long numberofsuffixes,
                             GtBucketspec2 *bucketspec2,
                             const GtEncseq *encseq,
                             GtReadmode readmode,
                             GtCodetype mincode,
                             GtCodetype maxcode,
                             const Bcktab *bcktab,
                             unsigned int numofchars,
                             unsigned int prefixlength,
                             const Sfxstrategy *sfxstrategy,
                             void *voiddcov,
                             void (*dc_processunsortedrange)(void *,
                                                             Suffixptr *,
                                                             unsigned long,
                                                             unsigned long),
                             GtLogger *logger);

#endif
