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
#include "defined-types.h"
#include "encodedsequence.h"
#include "intcode-def.h"
#include "core/seqpos.h"
#include "sfx-strategy.h"
#include "sfx-copysort.h"
#include "bcktab.h"

typedef struct Outlcpinfo Outlcpinfo;

typedef struct
{
  Seqpos *sortspace,
         offset;    /* negative offset */
  DefinedSeqpos longest;
} Suftab;

Outlcpinfo *newOutlcpinfo(const GtStr *indexname,
                          unsigned int prefixlength,
                          unsigned int numofchars,
                          Seqpos totallength,
                          bool assideeffect,
                          GtError *err);

void freeOutlcptab(Outlcpinfo **outlcpinfoptr);

Seqpos getnumoflargelcpvalues(const Outlcpinfo *outlcpinfo);

Seqpos getmaxbranchdepth(const Outlcpinfo *outlcpinfo);

void qsufsort(Seqpos *sortspace,
              int mmapfiledesc,
              Seqpos *longest,
              const GtEncodedsequence *encseq,
              GtReadmode readmode,
              Codetype mincode,
              Codetype maxcode,
              Seqpos partwidth,
              Bcktab *bcktab,
              unsigned int numofchars,
              unsigned int prefixlength,
              bool hashexceptions,
              bool absoluteinversesuftab,
              Outlcpinfo *outlcpinfo);

void sortallbuckets(Suftab *suftab,
                    GtBucketspec2 *bucketspec2,
                    const GtEncodedsequence *encseq,
                    GtReadmode readmode,
                    Codetype mincode,
                    Codetype maxcode,
                    Seqpos partwidth,
                    Bcktab *bcktab,
                    unsigned int numofchars,
                    unsigned int prefixlength,
                    Outlcpinfo *outlcpinfo,
                    const Sfxstrategy *sfxstrategy,
                    unsigned long long *bucketiterstep,
                    GtLogger *logger);

void sortbucketofsuffixes(Seqpos *suffixestobesorted,
                          GtBucketspec2 *bucketspec2,
                          unsigned long numberofsuffixes,
                          const GtEncodedsequence *encseq,
                          GtReadmode readmode,
                          Codetype mincode,
                          Codetype maxcode,
                          const Bcktab *bcktab,
                          unsigned int numofchars,
                          unsigned int prefixlength,
                          const Sfxstrategy *sfxstrategy,
                          void *voiddcov,
                          void (*dc_processunsortedrange)(void *,Seqpos *,
                                                          Seqpos *,Seqpos),
                          GtLogger *logger);

#endif
