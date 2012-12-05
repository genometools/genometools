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
#include "core/compactulongstore.h"
#include "core/error_api.h"
#include "core/intbits.h"
#include "core/log_api.h"
#include "core/ma.h"
#include "core/unused_api.h"
#include "extended/popcount_tab.h"

struct GtPopcountTab
{
  GtCompactUlongstore *blocks;
  unsigned long *offsets;
};

static void init_offset_tab(unsigned long *offsets,
                            unsigned long blocksize,
                            unsigned long num_of_blocks)
{
  unsigned long idx, class_size;

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
                                GtCompactUlongstore *blocks)
{
  unsigned long v, init;

  v = init = perm_start(popcount);
  while (v >= init) {
    gt_GtCompactulongstore_update(blocks, idx++, v);
    if (popcount == 1U)
      v = (v << 1) & blockmask;
    else
      v = next_perm(v) & blockmask;
  }
  return idx;
}

static unsigned long init_blocks_tab(GtCompactUlongstore *blocks,
                                     unsigned char blocksize)
{
  unsigned long idx = 1UL,
                blockmask;
  unsigned char popcount = (unsigned char)1;
  blockmask = perm_start(blocksize);
  gt_GtCompactulongstore_update(blocks, 0, 0);
  while (popcount < blocksize) {
    idx = gen_blocks(popcount++, idx, blockmask, blocks);
  }
  gt_GtCompactulongstore_update(blocks, idx++, blockmask);
  return idx;
}

GtPopcountTab *gt_popcount_tab_new(unsigned char blocksize)
{
  GtPopcountTab *popcount_tab;
  unsigned long num_of_blocks = 1UL << blocksize;
  GT_UNUSED unsigned long idx_check;
  gt_assert(blocksize <= (unsigned char) GT_INTWORDSIZE);
  popcount_tab = gt_malloc(sizeof (GtPopcountTab));
  popcount_tab->blocks = gt_GtCompactulongstore_new(num_of_blocks, blocksize);
  popcount_tab->offsets = gt_malloc(sizeof (popcount_tab->offsets) * blocksize);
  init_offset_tab(popcount_tab->offsets, blocksize, num_of_blocks);
  idx_check = init_blocks_tab(popcount_tab->blocks, blocksize);
  gt_assert(idx_check == num_of_blocks);
  return popcount_tab;
}

void gt_popcount_tab_delete(GtPopcountTab *popcount_tab)
{
  if (popcount_tab != NULL) {
    gt_free(popcount_tab->offsets);
    gt_GtCompactulongstore_delete(popcount_tab->blocks);
    gt_free(popcount_tab);
  }
}
