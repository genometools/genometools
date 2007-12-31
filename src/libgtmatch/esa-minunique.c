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
#include "sarr-def.h"
#include "seqpos-def.h"
#include "divmodmul.h"

typedef struct
{
  Seqpos left,
         right;
} Simplelcpinterval;

static Seqpos findright(const Encodedsequence *encseq,
                        const Seqpos *suftab,
                        Seqpos offset,
                        Uchar cc,
                        Seqpos l,
                        Seqpos r)
{
  Seqpos pos, mid;
  Uchar midcc;

  while (r > l+1)
  {
    mid = DIV2(l+r);
    pos = suftab[mid] + offset;
    midcc = getencodedchar(encseq,pos,Forwardmode);
    if (midcc > cc)
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
                                Seqpos lcpvalue,
                                Seqpos i,
                                Seqpos j)
{
  Uchar leftcc, rightcc;
  Seqpos pos, rightbound, leftbound = i;

  pos = suftab[j] + lcpvalue;
  rightcc = getencodedchar(encseq,pos,Forwardmode);
  while (true)
  {
    pos = suftab[leftbound] + lcpvalue;
    leftcc = getencodedchar(encseq,pos,Forwardmode);
    if (leftcc == rightcc)
    {
      break;
    }
    rightbound = findright(encseq,suftab,lcpvalue,leftcc,leftbound,j);
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
                                        const Uchar *qstart,
                                        const Uchar *qend)
{
  Simplelcpinterval itv;
  const Uchar *qptr;
  const Suffixarray *suffixarray = (const Suffixarray *) genericindex;

  itv.left = 0;
  itv.right = getencseqtotallength(suffixarray->encseq);
  for (qptr = qstart; /* Nothing */; qptr++)
  {
    if (itv.left < itv.right)
    {
      if (qptr >= qend)
      {
        break;
      }
      if (ISSPECIAL(*qptr))
      {
        break;
      }
      if (!findcharintervalbin(suffixarray->encseq,
                               suffixarray->suftab,
                               &itv,
                               *qptr,
                               (Seqpos) (qptr - qstart),
                               itv.left,itv.right))
      {
        break;
      }
    } else
    {
      return (unsigned long) (qptr - qstart);
    }
  }
  return 0;
}
