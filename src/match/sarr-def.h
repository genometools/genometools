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

#ifndef SARR_DEF_H
#define SARR_DEF_H

#include <stdio.h>
#include "core/defined-types.h"
#include "core/encseq.h"
#include "core/types_api.h"
#include "core/unused_api.h"
#include "core/codetype.h"
#include "match/declare-readfunc.h"

#include "lcpoverflow.h"
#include "bcktab.h"

#define SARR_ESQTAB 1U
#define SARR_SUFTAB (1U << 1)
#define SARR_LCPTAB (1U << 2)
#define SARR_BWTTAB (1U << 3)
#define SARR_DESTAB (1U << 4)
#define SARR_SDSTAB (1U << 5)
#define SARR_BCKTAB (1U << 6)
#define SARR_SSPTAB (1U << 7)

#define SARR_ALLTAB (SARR_ESQTAB |\
                     SARR_SUFTAB |\
                     SARR_LCPTAB |\
                     SARR_BWTTAB |\
                     SARR_DESTAB |\
                     SARR_SSPTAB)

DECLAREBufferedfiletype(GtUword);
DECLAREREADFUNCTION(GtUword);

#if defined (_LP64) || defined (_WIN64)
DECLAREBufferedfiletype(uint32_t);
DECLAREREADFUNCTION(uint32_t);
#endif

DECLAREBufferedfiletype(GtUchar);
DECLAREREADFUNCTION(GtUchar);

DECLAREBufferedfiletype(Largelcpvalue);
DECLAREREADFUNCTION(Largelcpvalue);

typedef GtUword ESASuffixptr;

#define ESASUFFIXPTRGET(TAB,IDX)     TAB[IDX]

typedef struct
{
  GtEncseq *encseq;
  Definedunsignedlong numoflargelcpvalues, /* only in esa-map.c */
                      maxbranchdepth;
  Defineddouble averagelcp;
  Definedunsignedlong longest; /* for BWT */
  GtReadmode readmode; /* relevant when reading the encoded sequence */
  bool mirroredencseq;
  GtUword numberofallsortedsuffixes;
  /* either with mapped input */
  const ESASuffixptr *suftab;
  const GtUchar *lcptab;
  const Largelcpvalue *llvtab;
  const GtUchar *bwttab;
  unsigned int prefixlength;
  GtBcktab *bcktab;
  /* or with streams */
#if defined (_LP64) || defined (_WIN64)
  GtBufferedfile_uint32_t suftabstream_uint32_t;
#endif
  GtBufferedfile_GtUword suftabstream_GtUword;
  GtBufferedfile_GtUchar bwttabstream,
                         lcptabstream;
  GtBufferedfile_Largelcpvalue llvtabstream;
} Suffixarray;

/*@unused@*/ static inline const Largelcpvalue *getlargelcpvalue(
                       const Suffixarray *suffixarray,
                       GtUword pos)
{
  const Largelcpvalue *leftptr, *rightptr, *midptr;
  GtUword len;

  gt_assert(suffixarray->numoflargelcpvalues.defined);

  leftptr = suffixarray->llvtab;
  rightptr = suffixarray->llvtab +
             suffixarray->numoflargelcpvalues.valueunsignedlong - 1;

  while (leftptr<=rightptr)
  {
    len = (GtUword) (rightptr-leftptr);
    midptr = leftptr + GT_DIV2(len);
    if (pos < midptr->position)
    {
      rightptr = midptr - 1;
    } else
    {
      if (pos > midptr->position)
      {
        leftptr = midptr + 1;
      } else
      {
        return midptr;
      }
    }
  }
  return NULL;
}

/*@unused@*/ static inline GtUword lcptable_get(
                       const Suffixarray *suffixarray,
                       GtUword pos)
{
  GtUchar smalllcpvalue;
  const Largelcpvalue *largelcpvalue;

  gt_assert(pos <= gt_encseq_total_length(suffixarray->encseq));
  smalllcpvalue = suffixarray->lcptab[pos];
  if (smalllcpvalue != (GtUchar) LCPOVERFLOW)
  {
    return (GtUword) smalllcpvalue;
  }
  largelcpvalue = getlargelcpvalue(suffixarray,pos);
  gt_assert(largelcpvalue != NULL);
  return largelcpvalue->value;
}

#endif
