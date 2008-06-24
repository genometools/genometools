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
#include "libgtcore/symboldef.h"
#include "libgtcore/unused.h"
#include "sarr-def.h"
#include "seqpos-def.h"
#include "esa-splititv.h"

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
  Seqpos totallength;

  itv.left = left;
  itv.right = right;
  totallength = getencseqtotallength(suffixarray->encseq);
  for (qptr = qstart; /* Nothing */; qptr++, offset++)
  {
    if (itv.left < itv.right)
    {
      if (qptr >= qend || ISSPECIAL(*qptr) ||
          !lcpintervalfindcharchildintv(suffixarray->encseq,
                                        suffixarray->readmode,
                                        totallength,
                                        suffixarray->suftab,
                                        &itv,
                                        *qptr,
                                        (Seqpos) offset,
                                        itv.left,
                                        itv.right))
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
  Seqpos totallength;

  itv.left = left;
  itv.right = right;
  totallength = getencseqtotallength(suffixarray->encseq);
  for (qptr = qstart; /* Nothing */; qptr++, offset++)
  {
    assert(itv.left <= itv.right);
    if (qptr >= qend || ISSPECIAL(*qptr) ||
        !lcpintervalfindcharchildintv(suffixarray->encseq,
                                      suffixarray->readmode,
                                      totallength,
                                      suffixarray->suftab,
                                      &itv,
                                      *qptr,
                                      (Seqpos) offset,
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
