/*
  Copyright (c) 2011 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2011 Center for Bioinformatics, University of Hamburg

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
#include "firstcodes-tab.h"

unsigned long gt_firstcodes_get_leftborder(const GtFirstcodestab *fct,
                                           unsigned long idx)
{
  gt_assert(idx <= fct->differentcodes);
  return fct->countocc[idx];
}

size_t gt_firstcodes_size_to_split(const GtFirstcodestab *fct)
{
  return fct->size_to_split;
}

unsigned long gt_firstcodes_numofallcodes(const GtFirstcodestab *fct)
{
  return fct->differentcodes;
}

unsigned long gt_firstcodes_findfirstlarger(const GtFirstcodestab *fct,
                                            unsigned long suftaboffset)
{
  unsigned long left = 0, right = fct->differentcodes, mid, midval,
                found = fct->differentcodes;

  while (left+1 < right)
  {
    mid = GT_DIV2(left+right);
    midval = gt_firstcodes_get_leftborder(fct,mid);
    if (suftaboffset == midval)
    {
      return mid;
    }
    if (suftaboffset < midval)
    {
      found = mid;
      right = mid - 1;
    } else
    {
      left = mid + 1;
    }
  }
  return found;
}

unsigned long gt_firstcodes_mapped_range_size(const GtFirstcodestab *fct,
                                              unsigned long minindex,
                                              unsigned long maxindex)
{
  size_t idx;
  GtSfxmappedrange *maptab[3];
  unsigned long sumsize = 0;

  maptab[0] = fct->mappedcountocc;
  maptab[1] = fct->mappedallfirstcodes;
  maptab[2] = fct->mappedmarkprefix;;
  for (idx = 0; idx < sizeof (maptab)/sizeof (maptab[0]); idx++)
  {
    if (maptab[idx] == NULL)
    {
      return (unsigned long) fct->size_to_split;
    }
    sumsize += gt_Sfxmappedrange_size_mapped(maptab[idx],minindex,maxindex);
  }
  return sumsize;
}

unsigned long gt_firstcodes_idx2code(const GtFirstcodestab *fct,
                                     unsigned long idx)
{
  gt_assert(idx < fct->differentcodes);
  return fct->allfirstcodes[idx];
}
