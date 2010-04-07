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

#include "core/divmodmul.h"

#include "fmindex.h"
#include "fmi-occ.gen"

static unsigned long searchsmallestgeq(const GtPairBwtidx *left,
                                const GtPairBwtidx *right,
                                unsigned long key)
{
  const GtPairBwtidx *leftptr, *midptr, *rightptr, *found = NULL;
  unsigned long len;

  gt_assert(left != NULL);
  gt_assert(right != NULL);
  leftptr = left;
  rightptr = right;
  while (leftptr<=rightptr)
  {
    len = (unsigned long) (rightptr-leftptr);
    midptr = leftptr + GT_DIV2(len);
    if (key < midptr->bwtpos)
    {
      found = midptr;
      rightptr = midptr - 1;
    } else
    {
      if (key > midptr->bwtpos)
      {
        leftptr = midptr + 1;
      } else
      {
        return midptr->suftabvalue;
      }
    }
  }
  gt_assert(found != NULL);
  return found->suftabvalue;
}

unsigned long gt_fmfindtextpos (const Fmindex *fm,unsigned long idx)
{
  unsigned long offset = 0;
  GtUchar cc;

  while ((idx & fm->markdistminus1) != 0)
  {
    if (idx == fm->longestsuffixpos || ISSPECIAL(cc = ACCESSBWTTEXT(idx)))
    {
      unsigned long smallestgeq
               = searchsmallestgeq(fm->specpos.spaceGtPairBwtidx,
                                   fm->specpos.spaceGtPairBwtidx +
                                   fm->specpos.nextfreeGtPairBwtidx - 1,
                                   idx);
      return (smallestgeq + offset) % fm->bwtlength;
    }
    idx = fm->tfreq[cc] + fmoccurrence (fm, cc, idx);
    offset++;
  }
  return (fm->markpostable[idx / fm->markdist] + offset) % fm->bwtlength;
}
