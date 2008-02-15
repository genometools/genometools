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

#include "libgtcore/chardef.h"
#include "libgtcore/minmax.h"
#include "libgtcore/symboldef.h"
#include "libgtcore/unused.h"
#include "sarr-def.h"
#include "seqpos-def.h"
#include "divmodmul.h"

typedef struct
{
  Seqpos left,
         right;
} Simplelcpinterval;

#define SEQUENCE(ENCSEQ,POS) (((POS) == totallength) \
                             ? (Uchar) SEPARATOR\
                             : getencodedchar(ENCSEQ,POS,Forwardmode))

static Seqpos findright(const Encodedsequence *encseq,
                        const Seqpos *suftab,
                        Uchar cc,
                        unsigned long offset,
                        Seqpos l,
                        Seqpos r)
{
  Seqpos pos, mid, totallength = getencseqtotallength(encseq);
  Uchar midcc;

  while (r > l+1)
  {
    mid = DIV2(l+r);
    pos = suftab[mid] + offset;
    midcc = SEQUENCE(encseq,pos);
    if (cc < midcc)
    {
      r = mid;
    } else
    {
      l = mid;
    }
  }
  return l;
}

static bool findcharintervalbin(const Encodedsequence *encseq,
                                const Seqpos *suftab,
                                Simplelcpinterval *itv,
                                Uchar cc,
                                unsigned long lcpvalue,
                                Seqpos i,
                                Seqpos j)
{
  Uchar leftcc, rightcc;
  Seqpos pos, rightbound, leftbound = i,
         totallength = getencseqtotallength(encseq);

  pos = suftab[j] + lcpvalue;
  rightcc = SEQUENCE(encseq,pos);
  while (true)
  {
    pos = suftab[leftbound] + lcpvalue;
    leftcc = SEQUENCE(encseq,pos);
    if (leftcc == rightcc)
    {
      break;
    }
    rightbound = findright(encseq,suftab,leftcc,lcpvalue,leftbound,j);
    if (leftcc == cc)
    {
      itv->left = leftbound;
      itv->right = rightbound;
      return true;
    }
    if (leftcc > cc)
    {
      return false;
    }
    leftbound = rightbound+1;
  }
  if (leftcc == cc)
  {
    itv->left = leftbound;
    itv->right = j;
    return true;
  }
  return false;
}

unsigned long suffixarrayuniqueforward (const void *genericindex,
                                        unsigned long offset,
                                        Seqpos left,
                                        Seqpos right,
                                        UNUSED Seqpos *witnessposition,
                                        const Uchar *qstart,
                                        const Uchar *qend)
{
  Simplelcpinterval itv;
  const Uchar *qptr;
  const Suffixarray *suffixarray = (const Suffixarray *) genericindex;

  itv.left = left;
  itv.right = right;
  for (qptr = qstart; /* Nothing */; qptr++, offset++)
  {
    if (itv.left < itv.right)
    {
      if (qptr >= qend || ISSPECIAL(*qptr) ||
          !findcharintervalbin(suffixarray->encseq,
                               suffixarray->suftab,
                               &itv,
                               *qptr,
                               offset,
                               itv.left,itv.right))
      {
        break;
      }
    } else
    {
      return offset;
    }
  }
  return 0;
}

unsigned long suffixarraymstats (const void *genericindex,
                                 unsigned long offset,
                                 Seqpos left,
                                 Seqpos right,
                                 Seqpos *witnessposition,
                                 const Uchar *qstart,
                                 const Uchar *qend)
{
  Simplelcpinterval itv;
  const Uchar *qptr;
  const Suffixarray *suffixarray = (const Suffixarray *) genericindex;

  itv.left = left;
  itv.right = right;
  for (qptr = qstart; /* Nothing */; qptr++, offset++)
  {
    assert(itv.left <= itv.right);
    if (qptr >= qend || ISSPECIAL(*qptr) ||
        !findcharintervalbin(suffixarray->encseq,
                             suffixarray->suftab,
                             &itv,
                             *qptr,
                             offset,
                             itv.left,itv.right))
    {
      if (witnessposition != NULL)
      {
        *witnessposition = suffixarray->suftab[itv.left];
      }
      break;
    }
  }
  return offset;
}
