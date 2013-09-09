/*
  Copyright (c) 2011 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2012 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
  Copyright (c) 2011-2012 Center for Bioinformatics, University of Hamburg

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
#include "compact_ulong_store.h"

struct GtCompactUlongStore
{
  GtUword *tab,
                numofentries,
                maskright;
  unsigned int bitsperentry,
               bitsleft;
};

GtCompactUlongStore *gt_compact_ulong_store_new(GtUword numofentries,
                                                unsigned int bitsperentry)
{
  GtCompactUlongStore *cus;
  GtUword arraysize, totalbits;

  /* if this assertion appears, then probably the 32-bit version of the
     code is used. Better use the 64-bit version by compiling with
     64bit=yes */
  gt_assert(numofentries <= ULONG_MAX/bitsperentry);
  gt_assert(bitsperentry <= (unsigned int) GT_INTWORDSIZE);
  totalbits = numofentries * bitsperentry;
  cus = gt_malloc(sizeof (*cus));
  cus->numofentries = numofentries;
  arraysize = GT_DIVWORDSIZE(totalbits);
  if (GT_MODWORDSIZE(totalbits) > 0) {
    arraysize++;
  }
  cus->tab = gt_calloc((size_t) arraysize,sizeof (*cus->tab));
  cus->bitsperentry = bitsperentry;
  cus->bitsleft = (unsigned int) GT_INTWORDSIZE - cus->bitsperentry;
  cus->maskright = ~0UL >> cus->bitsleft;
  return cus;
}

size_t gt_compact_ulong_store_size(GtUword numofentries,
                                   unsigned int bitsperentry)
{
  GtUword arraysize, totalbits;

  /* if this assertion appears, then probably the 32-bit version of the
     code is used. Better use the 64-bit version by compiling with
     64bit=yes */
  gt_assert(numofentries <= ULONG_MAX/bitsperentry);
  totalbits = numofentries * bitsperentry;
  arraysize = GT_DIVWORDSIZE(totalbits);
  if (GT_MODWORDSIZE(totalbits) > 0) {
    arraysize++;
  }
  return sizeof (GtCompactUlongStore) + sizeof (GtUword) * arraysize;
}

void gt_compact_ulong_store_delete(GtCompactUlongStore *cus)
{
  if (cus != NULL) {
    gt_free(cus->tab);
    gt_free(cus);
  }
}

GtUword gt_compact_ulong_store_get(const GtCompactUlongStore *cus,
                                         GtUword idx)
{
  unsigned int unitoffset;
  GtUword unitindex;

  gt_assert(idx < cus->numofentries);
  idx *= cus->bitsperentry;
  unitoffset = (unsigned int) GT_MODWORDSIZE(idx);
  unitindex = GT_DIVWORDSIZE(idx);
  if (unitoffset <= (unsigned int) cus->bitsleft) {
    return (GtUword) (cus->tab[unitindex] >>
                             (cus->bitsleft - unitoffset))
           & cus->maskright;
  } else {
    return (GtUword)
           ((cus->tab[unitindex] <<
                (unitoffset + cus->bitsperentry - GT_INTWORDSIZE)) |
              (cus->tab[unitindex+1] >>
                (GT_INTWORDSIZE + cus->bitsleft - unitoffset))) &
           cus->maskright;
  }
}

void gt_compact_ulong_store_update(GtCompactUlongStore *cus,
                                   GtUword idx, GtUword value)
{
  unsigned int unitoffset;
  GtUword unitindex;

  gt_assert(idx < cus->numofentries && value <= cus->maskright);
  idx *= cus->bitsperentry;
  unitoffset = (unsigned int) GT_MODWORDSIZE(idx);
  unitindex = GT_DIVWORDSIZE(idx);
  if (unitoffset <= (unsigned int) cus->bitsleft) {
    unsigned int shiftleft = cus->bitsleft - unitoffset;

    cus->tab[unitindex]
      = (GtUword) (cus->tab[unitindex] & ~(cus->maskright << shiftleft)) |
        (GtUword) (value << shiftleft);
  } else {
    unsigned int shift = unitoffset - cus->bitsleft;

    cus->tab[unitindex]
      = (GtUword) (cus->tab[unitindex] & ~(cus->maskright >> shift)) |
        (GtUword)(value >> shift);
    shift = (unsigned int) GT_INTWORDSIZE - shift;
    cus->tab[unitindex+1]
      = (GtUword) (cus->tab[unitindex+1] & ~(cus->maskright << shift)) |
        (GtUword) (value << shift);
  }
}

int gt_compact_ulong_store_unit_test(GtError *err)
{
  GtCompactUlongStore *cus;
  const GtUword constnums = 100000UL;
  unsigned int bits;
  GtUword value, nums, idx, numforbits, *checknumbers;
  int had_err = 0;

  checknumbers = gt_malloc(sizeof (*checknumbers) * constnums);
  for (bits = 1U; bits <= (unsigned int) GT_INTWORDSIZE-1; bits++) {
    numforbits = 1UL << bits;
    nums = numforbits < constnums ? numforbits : constnums;
    cus = gt_compact_ulong_store_new(nums,bits);
    for (idx = 0; idx < nums; idx++) {
      checknumbers[idx] = nums == constnums ? gt_rand_max(numforbits-1) : idx;
      gt_compact_ulong_store_update(cus,idx,checknumbers[idx]);
      value = gt_compact_ulong_store_get(cus,idx);
      gt_ensure(checknumbers[idx] == value);
    }
    for (idx = 0; had_err == 0 && idx < nums; idx++) {
      value = gt_compact_ulong_store_get(cus,idx);
      gt_ensure(checknumbers[idx] == value);
    }
    gt_compact_ulong_store_delete(cus);
  }
  gt_free(checknumbers);
  return had_err;
}
