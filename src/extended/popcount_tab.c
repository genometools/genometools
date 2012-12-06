/*
  Copyright (c) 2012 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
  Copyright (c) 2012 Center for Bioinformatics, University of Hamburg

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

#include "core/assert_api.h"
#include "core/combinatorics.h"
#include "core/compact_ulong_store.h"
#include "core/ensure.h"
#include "core/error_api.h"
#include "core/intbits.h"
#include "core/log_api.h"
#include "core/ma.h"
#include "core/unused_api.h"
#include "extended/popcount_tab.h"

struct GtPopcountTab
{
  unsigned char blocksize;
  unsigned long num_of_blocks;
  GtCompactUlongStore *blocks;
  unsigned long *offsets;
};

static void init_offset_tab(GtPopcountTab *popcount_tab)
{
  unsigned long idx, class_size,
                *offsets = popcount_tab->offsets,
                blocksize = popcount_tab->blocksize,
                num_of_blocks = popcount_tab->num_of_blocks;

  offsets[0] = 1UL;
  offsets[1] = blocksize + 1UL;
  for (idx = 2UL; idx < blocksize - 2; idx++) {
    class_size = gt_binomialCoeff_with_ln(blocksize, idx);
    offsets[idx] = offsets[idx - 1] + class_size;
  }
  offsets[blocksize - 2] = num_of_blocks - blocksize - 1;
  offsets[blocksize - 1] = num_of_blocks - 1;
}

static unsigned long next_perm(unsigned long v)
{
  unsigned long head, tail;

  head = (v | (v - 1)) + 1;
  tail = v & (((head & -head) >> 1) -1);
  if (tail != 0) {
#ifdef __SSE4_2__
    int zero_trail = __builtin_ctzl(tail);
    tail >>= zero_trail;
#else
    while (( tail & 1 ) == 0)
      tail = tail >> 1;
#endif
  }
  return head | tail;
}

static unsigned long perm_start(unsigned char bits)
{
  return (1UL << bits) - 1;
}

static unsigned long gen_blocks(unsigned char popcount, unsigned long idx,
                                unsigned long blockmask,
                                GtCompactUlongStore *blocks)
{
  unsigned long v, init;

  v = init = perm_start(popcount);
  while (v >= init) {
    gt_compact_ulong_store_update(blocks, idx++, v);
    if (popcount == 1U)
      v = (v << 1) & blockmask;
    else
      v = next_perm(v) & blockmask;
  }
  return idx;
}

static unsigned long init_blocks_tab(GtCompactUlongStore *blocks,
                                     unsigned char blocksize)
{
  unsigned long idx = 1UL,
                blockmask;
  unsigned char popcount = (unsigned char)1;
  blockmask = perm_start(blocksize);
  gt_compact_ulong_store_update(blocks, 0, 0);
  while (popcount < blocksize) {
    idx = gen_blocks(popcount++, idx, blockmask, blocks);
  }
  gt_compact_ulong_store_update(blocks, idx++, blockmask);
  return idx;
}

GtPopcountTab *gt_popcount_tab_new(unsigned char blocksize)
{
  GtPopcountTab *popcount_tab;
  GT_UNUSED unsigned long idx_check;
  gt_assert(blocksize <= (unsigned char) GT_INTWORDSIZE);
  popcount_tab = gt_malloc(sizeof (GtPopcountTab));
  popcount_tab->num_of_blocks = 1UL << blocksize;
  popcount_tab->blocksize = blocksize;
  popcount_tab->blocks = gt_compact_ulong_store_new(popcount_tab->num_of_blocks,
                                                    blocksize);
  popcount_tab->offsets = gt_malloc(sizeof (popcount_tab->offsets) * blocksize);
  init_offset_tab(popcount_tab);
  idx_check = init_blocks_tab(popcount_tab->blocks, blocksize);
  gt_assert(idx_check == popcount_tab->num_of_blocks);
  return popcount_tab;
}

void gt_popcount_tab_delete(GtPopcountTab *popcount_tab)
{
  if (popcount_tab != NULL) {
    gt_free(popcount_tab->offsets);
    gt_compact_ulong_store_delete(popcount_tab->blocks);
    gt_free(popcount_tab);
  }
}

unsigned long gt_popcount_tab_get(GtPopcountTab *popcount_tab,
                                  unsigned char popcount,
                                  unsigned long offset) {
  gt_assert(popcount <= popcount_tab->blocksize);
  if (popcount == 0) {
    gt_assert(offset == 0);
    return 0;
  }
  if (popcount < popcount_tab->blocksize)
    gt_assert(offset < popcount_tab->offsets[popcount] -
                       popcount_tab->offsets[popcount - 1]);
  else
    gt_assert(offset == 0);
  return gt_compact_ulong_store_get(popcount_tab->blocks,
                                    popcount_tab->offsets[popcount - 1] +
                                      offset);
}

size_t gt_popcount_tab_get_size(unsigned char blocksize) {
  unsigned long num_of_blocks = 1UL << blocksize;
  size_t size = gt_compact_ulong_store_size(num_of_blocks, blocksize);
  size += sizeof (GtPopcountTab);
  size += blocksize * sizeof (unsigned long);
  return size;
}

int gt_popcount_tab_unit_test(GtError *err)
{
  int had_err = 0;
  unsigned long idx, jdx;
  unsigned char popcount;
  static const unsigned char blocksize = (unsigned char) 4;
  static const unsigned long blocksize_four[] =
    {0,1UL,2UL,4UL,8UL,3UL,5UL,6UL,9UL,10UL,12UL,7UL,11UL,13UL,14UL,15UL};
  GtPopcountTab *popcount_t = gt_popcount_tab_new(blocksize);

  for (idx = 0; idx < (1UL << blocksize); idx++) {
    gt_ensure(had_err, blocksize_four[idx] ==
              gt_compact_ulong_store_get(popcount_t->blocks, idx));
  }
  for (popcount = (unsigned char) 0, idx = 0;
       !had_err && popcount <= blocksize;
       idx += gt_binomialCoeff_with_ln(blocksize, popcount), popcount++) {
    for (jdx = 0;
         !had_err && jdx < gt_binomialCoeff_with_ln(blocksize, popcount);
         jdx++) {
      gt_ensure(had_err, blocksize_four[idx + jdx] ==
                         gt_popcount_tab_get(popcount_t, popcount, jdx));
    }
  }
  gt_popcount_tab_delete(popcount_t);
  return had_err;
}
