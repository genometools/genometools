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

#include <math.h>
#include "core/chardef.h"
#include "core/divmodmul.h"
#include "core/unused_api.h"
#include "fmindex.h"

#include "fmi-occ.gen"
#include "fmi-locate.pr"

unsigned long gt_skfmuniqueforward (const void *genericindex,
                                 GT_UNUSED unsigned long offset,
                                 GT_UNUSED unsigned long left,
                                 GT_UNUSED unsigned long right,
                                 GT_UNUSED unsigned long *witnessposition,
                                 const GtUchar *qstart,
                                 const GtUchar *qend)
{
  GtUchar cc;
  const GtUchar *qptr;
  GtUlongBound bwtbound;
  const Fmindex *fmindex = (Fmindex *) genericindex;

  gt_assert(qstart < qend);
  qptr = qstart;
  cc = *qptr++;
  if (ISSPECIAL(cc))
  {
    return 0;
  }
  bwtbound.lbound = fmindex->tfreq[cc];
  bwtbound.ubound = fmindex->tfreq[cc+1];
  while (qptr < qend && bwtbound.lbound + 1 < bwtbound.ubound)
  {
    cc = *qptr;
    if (ISSPECIAL (cc))
    {
      return 0;
    }
    bwtbound.lbound = fmindex->tfreq[cc] +
                      fmoccurrence (fmindex, cc, bwtbound.lbound);
    bwtbound.ubound = fmindex->tfreq[cc] +
                      fmoccurrence (fmindex, cc, bwtbound.ubound);
    qptr++;
  }
  if (bwtbound.lbound + 1 == bwtbound.ubound)
  {
    return (unsigned long) (qptr - qstart);
  }
  return 0;
}

unsigned long gt_skfmmstats (const void *genericindex,
                          GT_UNUSED unsigned long offset,
                          GT_UNUSED unsigned long left,
                          GT_UNUSED unsigned long right,
                          unsigned long *witnessposition,
                          const GtUchar *qstart,
                          const GtUchar *qend)
{
  GtUchar cc;
  const GtUchar *qptr;
  unsigned long prevlbound;
  unsigned long matchlength;
  GtUlongBound bwtbound;
  const Fmindex *fmindex = (Fmindex *) genericindex;

  gt_assert(qstart < qend);
  qptr = qstart;
  cc = *qptr;
  if (ISSPECIAL(cc))
  {
    return 0;
  }
  bwtbound.lbound = fmindex->tfreq[cc];
  bwtbound.ubound = fmindex->tfreq[cc+1];
  if (bwtbound.lbound >= bwtbound.ubound)
  {
    return 0;
  }
  prevlbound = bwtbound.lbound;
  for (qptr++; qptr < qend; qptr++)
  {
    cc = *qptr;
    if (ISSPECIAL (cc))
    {
      break;
    }
    bwtbound.lbound = fmindex->tfreq[cc] +
                      fmoccurrence (fmindex, cc, bwtbound.lbound);
    bwtbound.ubound = fmindex->tfreq[cc] +
                      fmoccurrence (fmindex, cc, bwtbound.ubound);
    if (bwtbound.lbound >= bwtbound.ubound)
    {
      break;
    }
    prevlbound = bwtbound.lbound;
  }
  matchlength = (unsigned long) (qptr - qstart);
  if (witnessposition != NULL)
  {
    unsigned long startpos = gt_fmfindtextpos (fmindex,prevlbound);
    gt_assert((fmindex->bwtlength-1) >= (startpos + matchlength));
    *witnessposition = (fmindex->bwtlength-1) - (startpos + matchlength);
  }
  return matchlength;
}
