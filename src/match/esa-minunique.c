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

#include "core/chardef.h"
#include "core/types_api.h"
#include "core/unused_api.h"
#include "sarr-def.h"

#include "esa-splititv.h"
#include "esa-minunique.h"

unsigned long gt_suffixarrayuniqueforward (const void *genericindex,
                                       unsigned long offset,
                                       unsigned long left,
                                       unsigned long right,
                                       GT_UNUSED unsigned long
                                          *witnessposition,
                                       const GtUchar *qstart,
                                       const GtUchar *qend)
{
  Simplelcpinterval itv;
  const GtUchar *qptr;
  const Suffixarray *suffixarray = (const Suffixarray *) genericindex;
  unsigned long totallength;

  itv.left = left;
  itv.right = right;
  totallength = gt_encseq_total_length(suffixarray->encseq);
  for (qptr = qstart; /* Nothing */; qptr++, offset++)
  {
    if (itv.left < itv.right)
    {
      if (qptr >= qend || ISSPECIAL(*qptr) ||
          !gt_lcpintervalfindcharchildintv(suffixarray->encseq,
                                           suffixarray->readmode,
                                           totallength,
                                           suffixarray->suftab,
                                           &itv,
                                           *qptr,
                                           offset,
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

unsigned long gt_suffixarraymstats (const void *genericindex,
                                    unsigned long offset,
                                    unsigned long left,
                                    unsigned long right,
                                    unsigned long *witnessposition,
                                    const GtUchar *qstart,
                                    const GtUchar *qend)
{
  Simplelcpinterval itv;
  const GtUchar *qptr;
  const Suffixarray *suffixarray = (const Suffixarray *) genericindex;
  unsigned long totallength;

  itv.left = left;
  itv.right = right;
  totallength = gt_encseq_total_length(suffixarray->encseq);
  for (qptr = qstart; /* Nothing */; qptr++, offset++)
  {
    gt_assert(itv.left <= itv.right);
    if (qptr >= qend || ISSPECIAL(*qptr) ||
        !gt_lcpintervalfindcharchildintv(suffixarray->encseq,
                                      suffixarray->readmode,
                                      totallength,
                                      suffixarray->suftab,
                                      &itv,
                                      *qptr,
                                      (unsigned long) offset,
                                      itv.left,itv.right))
    {
      if (witnessposition != NULL)
      {
        *witnessposition = ESASUFFIXPTRGET(suffixarray->suftab,itv.left);
      }
      break;
    }
  }
  return offset;
}

unsigned long gt_suffixarrayfindmums (const void *genericindex,
                                      unsigned long offset,
                                      unsigned long left,
                                      unsigned long right,
                                      unsigned long *witnessposition,
                                      const GtUchar *qstart,
                                      const GtUchar *qend)
{
  const Suffixarray *suffixarray = (const Suffixarray *) genericindex;
  Simplelcpinterval itv;
  const GtUchar *qptr;
  unsigned long totallength;

  itv.left = left;
  itv.right = right;
  totallength = gt_encseq_total_length(suffixarray->encseq);
  *witnessposition = ULONG_MAX;
  for (qptr = qstart; /* Nothing */; qptr++, offset++)
  {
    gt_assert(itv.left <= itv.right);
    if (qptr >= qend || ISSPECIAL(*qptr) ||
        !gt_lcpintervalfindcharchildintv(suffixarray->encseq,
                                         suffixarray->readmode,
                                         totallength,
                                         suffixarray->suftab,
                                         &itv,
                                         *qptr,
                                         offset,
                                         itv.left,
                                         itv.right))
    {
      if (itv.left == itv.right)
      {
        *witnessposition = ESASUFFIXPTRGET(suffixarray->suftab,itv.left);
      }
      break;
    }
  }
  return offset;
}

GtRange gt_suffixarrayfindinterval (const void *genericindex,
                                    unsigned long offset,
                                    unsigned long left,
                                    unsigned long right,
                                    unsigned long *matchlength,
                                    const GtUchar *qstart,
                                    const GtUchar *qend)
{
  const Suffixarray *suffixarray = (const Suffixarray *) genericindex;
  const GtUchar *qptr;
  unsigned long totallength = gt_encseq_total_length(suffixarray->encseq);
  GtRange result;

  gt_assert(qstart < qend);
  for (qptr = qstart; /* Nothing */; qptr++, offset++)
  {
    Simplelcpinterval itv;

    if (qptr >= qend || ISSPECIAL(*qptr) ||
        !gt_lcpintervalfindcharchildintv(suffixarray->encseq,
                                         suffixarray->readmode,
                                         totallength,
                                         suffixarray->suftab,
                                         &itv,
                                         *qptr,
                                         offset,
                                         left,
                                         right))
    {
      break;
    }
    gt_assert(itv.left <= itv.right);
    left = itv.left;
    right = itv.right;
  }
  *matchlength = (unsigned long) (qptr - qstart);
  result.start = left;
  result.end = right;
  return result;
}
