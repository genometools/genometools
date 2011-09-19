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

#include <limits.h>
#include "intbits.h"
#include "error_api.h"
#include "mathsupport.h"
#include "ensure.h"
#include "assert_api.h"
#include "compactulongstore.h"

struct GtCompactUlongstore
{
  unsigned long *tab,
                numofentries,
                maskright;
  unsigned int bitsperentry,
               bitsleft;
};

GtCompactUlongstore *gt_GtCompactulongstore_new(unsigned long numofentries,
                                                unsigned int bitsperentry)
{
  GtCompactUlongstore *cus;
  unsigned long arraysize, totalbits;

  gt_assert(numofentries <= ULONG_MAX/bitsperentry);
  totalbits = numofentries * bitsperentry;
  cus = gt_malloc(sizeof (*cus));
  cus->numofentries = numofentries;
  arraysize = GT_DIVWORDSIZE(totalbits);
  if (GT_MODWORDSIZE(totalbits) > 0)
  {
    arraysize++;
  }
  cus->tab = gt_calloc((size_t) arraysize,sizeof (*cus->tab));
  cus->bitsperentry = bitsperentry;
  gt_assert(bitsperentry <= (unsigned int) GT_INTWORDSIZE);
  cus->bitsleft = (unsigned int) GT_INTWORDSIZE - cus->bitsperentry;
  cus->maskright = ~0UL >> cus->bitsleft;
  return cus;
}

size_t gt_GtCompactulongstore_size(unsigned long numofentries,
                                   unsigned int bitsperentry)
{
  unsigned long arraysize, totalbits;

  gt_assert(numofentries <= ULONG_MAX/bitsperentry);
  totalbits = numofentries * bitsperentry;
  arraysize = GT_DIVWORDSIZE(totalbits);
  if (GT_MODWORDSIZE(totalbits) > 0)
  {
    arraysize++;
  }
  return sizeof (GtCompactUlongstore) + sizeof (unsigned long) * arraysize;
}

void gt_GtCompactulongstore_delete(GtCompactUlongstore *cus)
{
  if (cus != NULL)
  {
    gt_free(cus->tab);
    gt_free(cus);
  }
}

unsigned long gt_GtCompactulongstore_get(const GtCompactUlongstore *cus,
                                         unsigned long idx)
{
  unsigned int unitoffset;
  unsigned long unitindex;

  gt_assert(idx < cus->numofentries);
  idx *= cus->bitsperentry;
  unitoffset = (unsigned int) GT_MODWORDSIZE(idx);
  unitindex = GT_DIVWORDSIZE(idx);
  if (unitoffset <= (unsigned int) cus->bitsleft)
  {
    return (unsigned long) (cus->tab[unitindex] >>
                            (cus->bitsleft - unitoffset))
           & cus->maskright;
  } else
  {
    return (unsigned long)
           ((cus->tab[unitindex] << (unitoffset + cus->bitsperentry
                                                - GT_INTWORDSIZE)) |
            (cus->tab[unitindex+1] >> (GT_INTWORDSIZE + cus->bitsleft -
                                       unitoffset)))
           & cus->maskright;
  }
}

void gt_GtCompactulongstore_update(GtCompactUlongstore *cus,
                                   unsigned long idx,unsigned long value)
{
  unsigned int unitoffset;
  unsigned long unitindex;

  gt_assert(idx < cus->numofentries);
  gt_assert(value <= cus->maskright);
  idx *= cus->bitsperentry;
  unitoffset = (unsigned int) GT_MODWORDSIZE(idx);
  unitindex = GT_DIVWORDSIZE(idx);
  if (unitoffset <= (unsigned int) cus->bitsleft)
  {
    unsigned int shiftleft = cus->bitsleft - unitoffset;

    cus->tab[unitindex]
      = (cus->tab[unitindex] & ~(cus->maskright << shiftleft)) |
        (value << shiftleft);
  } else
  {
    unsigned int shiftright = unitoffset - cus->bitsleft;

    cus->tab[unitindex]
      = (cus->tab[unitindex] & ~(cus->maskright >> shiftright)) |
        (value >> shiftright);
    cus->tab[unitindex+1]
      = (cus->tab[unitindex+1] & ~(cus->maskright <<
                                   (GT_INTWORDSIZE - shiftright))) |
        (value << (GT_INTWORDSIZE - shiftright));
  }
}

int gt_GtCompactulongstore_unit_test(GtError *err)
{
  GtCompactUlongstore *cus;
  const unsigned long constnums = 100000UL;
  unsigned int bits;
  unsigned long value, nums, idx, numforbits, *checknumbers;
  int had_err = 0;

  checknumbers = gt_malloc(sizeof (*checknumbers) * constnums);
  for (bits = 1U; bits <= (unsigned int) GT_INTWORDSIZE-1; bits++)
  {
    numforbits = 1UL << bits;
    nums = numforbits < constnums ? numforbits : constnums;
    cus = gt_GtCompactulongstore_new(nums,bits);
    for (idx = 0; idx < nums; idx++)
    {
      checknumbers[idx] = nums == constnums ? gt_rand_max(numforbits-1) : idx;
      gt_GtCompactulongstore_update(cus,idx,checknumbers[idx]);
      value = gt_GtCompactulongstore_get(cus,idx);
      ensure(had_err,checknumbers[idx] == value);
    }
    for (idx = 0; had_err == 0 && idx < nums; idx++)
    {
      value = gt_GtCompactulongstore_get(cus,idx);
      ensure(had_err,checknumbers[idx] == value);
    }
    gt_GtCompactulongstore_delete(cus);
  }
  gt_free(checknumbers);
  return had_err;
}
