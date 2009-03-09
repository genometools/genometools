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

#ifndef SFX_LCPSUB_H
#define SFX_LCPSUB_H

#include <stdio.h>
#include "core/symboldef.h"
#include "core/arraydef.h"
#include "seqpos-def.h"

DECLAREARRAYSTRUCT(Largelcpvalue);

typedef struct
{
  void *reservoir;
  size_t sizereservoir;
  Seqpos *spaceSeqpos, /* pointer into reservoir */
         maxbranchdepth,
         numoflargelcpvalues,
         countoutputlcpvalues;
  Uchar *smalllcpvalues; /* pointer into reservoir */
  ArrayLargelcpvalue largelcpvalues;
  const Seqpos *suftabbase;
} Lcpsubtab;

void multilcpvalue(Lcpsubtab *lcpsubtab,
                   unsigned long bucketsize,
                   Seqpos posoffset,
                   FILE *fplcptab,
                   FILE *fpllvtab);

void outmany0lcpvalues(Lcpsubtab *lcpsubtab,Seqpos totallength,
                       FILE *outfplcptab);

#endif
