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

#ifdef POPCOUNT_TL
static const uint8_t B1CntBytes[] = {
  0, (uint8_t) 1, (uint8_t) 1, (uint8_t) 2,
  (uint8_t) 1, (uint8_t) 2, (uint8_t) 2, (uint8_t) 3,
  (uint8_t) 1, (uint8_t) 2, (uint8_t) 2, (uint8_t) 3,
  (uint8_t) 2, (uint8_t) 3, (uint8_t) 3, (uint8_t) 4,
  (uint8_t) 1, (uint8_t) 2, (uint8_t) 2, (uint8_t) 3,
  (uint8_t) 2, (uint8_t) 3, (uint8_t) 3, (uint8_t) 4,
  (uint8_t) 2, (uint8_t) 3, (uint8_t) 3, (uint8_t) 4,
  (uint8_t) 3, (uint8_t) 4, (uint8_t) 4, (uint8_t) 5,
  (uint8_t) 1, (uint8_t) 2, (uint8_t) 2, (uint8_t) 3,
  (uint8_t) 2, (uint8_t) 3, (uint8_t) 3, (uint8_t) 4,
  (uint8_t) 2, (uint8_t) 3, (uint8_t) 3, (uint8_t) 4,
  (uint8_t) 3, (uint8_t) 4, (uint8_t) 4, (uint8_t) 5,
  (uint8_t) 2, (uint8_t) 3, (uint8_t) 3, (uint8_t) 4,
  (uint8_t) 3, (uint8_t) 4, (uint8_t) 4, (uint8_t) 5,
  (uint8_t) 3, (uint8_t) 4, (uint8_t) 4, (uint8_t) 5,
  (uint8_t) 4, (uint8_t) 5, (uint8_t) 5, (uint8_t) 6,
  (uint8_t) 1, (uint8_t) 2, (uint8_t) 2, (uint8_t) 3,
  (uint8_t) 2, (uint8_t) 3, (uint8_t) 3, (uint8_t) 4,
  (uint8_t) 2, (uint8_t) 3, (uint8_t) 3, (uint8_t) 4,
  (uint8_t) 3, (uint8_t) 4, (uint8_t) 4, (uint8_t) 5,
  (uint8_t) 2, (uint8_t) 3, (uint8_t) 3, (uint8_t) 4,
  (uint8_t) 3, (uint8_t) 4, (uint8_t) 4, (uint8_t) 5,
  (uint8_t) 3, (uint8_t) 4, (uint8_t) 4, (uint8_t) 5,
  (uint8_t) 4, (uint8_t) 5, (uint8_t) 5, (uint8_t) 6,
  (uint8_t) 2, (uint8_t) 3, (uint8_t) 3, (uint8_t) 4,
  (uint8_t) 3, (uint8_t) 4, (uint8_t) 4, (uint8_t) 5,
  (uint8_t) 3, (uint8_t) 4, (uint8_t) 4, (uint8_t) 5,
  (uint8_t) 4, (uint8_t) 5, (uint8_t) 5, (uint8_t) 6,
  (uint8_t) 3, (uint8_t) 4, (uint8_t) 4, (uint8_t) 5,
  (uint8_t) 4, (uint8_t) 5, (uint8_t) 5, (uint8_t) 6,
  (uint8_t) 4, (uint8_t) 5, (uint8_t) 5, (uint8_t) 6,
  (uint8_t) 5, (uint8_t) 6, (uint8_t) 6, (uint8_t) 7,
  (uint8_t) 1, (uint8_t) 2, (uint8_t) 2, (uint8_t) 3,
  (uint8_t) 2, (uint8_t) 3, (uint8_t) 3, (uint8_t) 4,
  (uint8_t) 2, (uint8_t) 3, (uint8_t) 3, (uint8_t) 4,
  (uint8_t) 3, (uint8_t) 4, (uint8_t) 4, (uint8_t) 5,
  (uint8_t) 2, (uint8_t) 3, (uint8_t) 3, (uint8_t) 4,
  (uint8_t) 3, (uint8_t) 4, (uint8_t) 4, (uint8_t) 5,
  (uint8_t) 3, (uint8_t) 4, (uint8_t) 4, (uint8_t) 5,
  (uint8_t) 4, (uint8_t) 5, (uint8_t) 5, (uint8_t) 6,
  (uint8_t) 2, (uint8_t) 3, (uint8_t) 3, (uint8_t) 4,
  (uint8_t) 3, (uint8_t) 4, (uint8_t) 4, (uint8_t) 5,
  (uint8_t) 3, (uint8_t) 4, (uint8_t) 4, (uint8_t) 5,
  (uint8_t) 4, (uint8_t) 5, (uint8_t) 5, (uint8_t) 6,
  (uint8_t) 3, (uint8_t) 4, (uint8_t) 4, (uint8_t) 5,
  (uint8_t) 4, (uint8_t) 5, (uint8_t) 5, (uint8_t) 6,
  (uint8_t) 4, (uint8_t) 5, (uint8_t) 5, (uint8_t) 6,
  (uint8_t) 5, (uint8_t) 6, (uint8_t) 6, (uint8_t) 7,
  (uint8_t) 2, (uint8_t) 3, (uint8_t) 3, (uint8_t) 4,
  (uint8_t) 3, (uint8_t) 4, (uint8_t) 4, (uint8_t) 5,
  (uint8_t) 3, (uint8_t) 4, (uint8_t) 4, (uint8_t) 5,
  (uint8_t) 4, (uint8_t) 5, (uint8_t) 5, (uint8_t) 6,
  (uint8_t) 3, (uint8_t) 4, (uint8_t) 4, (uint8_t) 5,
  (uint8_t) 4, (uint8_t) 5, (uint8_t) 5, (uint8_t) 6,
  (uint8_t) 4, (uint8_t) 5, (uint8_t) 5, (uint8_t) 6,
  (uint8_t) 5, (uint8_t) 6, (uint8_t) 6, (uint8_t) 7,
  (uint8_t) 3, (uint8_t) 4, (uint8_t) 4, (uint8_t) 5,
  (uint8_t) 4, (uint8_t) 5, (uint8_t) 5, (uint8_t) 6,
  (uint8_t) 4, (uint8_t) 5, (uint8_t) 5, (uint8_t) 6,
  (uint8_t) 5, (uint8_t) 6, (uint8_t) 6, (uint8_t) 7,
  (uint8_t) 4, (uint8_t) 5, (uint8_t) 5, (uint8_t) 6,
  (uint8_t) 5, (uint8_t) 6, (uint8_t) 6, (uint8_t) 7,
  (uint8_t) 5, (uint8_t) 6, (uint8_t) 6, (uint8_t) 7,
  (uint8_t) 6, (uint8_t) 7, (uint8_t) 7, (uint8_t) 8
};
#endif

struct GtPopcountTab
{
  unsigned blocksize;
  unsigned long num_of_blocks;
  GtCompactUlongStore *blocks;
  unsigned long *offsets;
};

static void init_offset_tab(GtPopcountTab *popcount_tab)
{
  unsigned long idx, class_size,
                *offsets = popcount_tab->offsets,
                blocksize = (unsigned long) popcount_tab->blocksize,
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

static unsigned long perm_start(unsigned bits)
{
  return (1UL << bits) - 1;
}

static unsigned long gen_blocks(unsigned popcount_c, unsigned long idx,
                                unsigned long blockmask,
                                GtCompactUlongStore *blocks)
{
  unsigned long v, init;

  v = init = perm_start(popcount_c);
  while (v >= init) {
    gt_compact_ulong_store_update(blocks, idx++, v);
    if (popcount_c == 1U)
      v = (v << 1) & blockmask;
    else
      v = next_perm(v) & blockmask;
  }
  return idx;
}

static unsigned long init_blocks_tab(GtCompactUlongStore *blocks,
                                     unsigned blocksize)
{
  unsigned long idx = 1UL,
                blockmask;
  unsigned popcount_c = 1U;
  blockmask = perm_start(blocksize);
  gt_compact_ulong_store_update(blocks, 0, 0);
  while (popcount_c < blocksize) {
    idx = gen_blocks(popcount_c++, idx, blockmask, blocks);
  }
  gt_compact_ulong_store_update(blocks, idx++, blockmask);
  return idx;
}

GtPopcountTab *gt_popcount_tab_new(unsigned blocksize)
{
  GtPopcountTab *popcount_tab;
  GT_UNUSED unsigned long idx_check;
  gt_assert(blocksize <= (unsigned) GT_INTWORDSIZE);
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
                                  unsigned popcount_c,
                                  unsigned long offset) {
  gt_assert(popcount_c <= popcount_tab->blocksize);
  if (popcount_c == 0) {
    gt_assert(offset == 0);
    return 0;
  }
  if (popcount_c < popcount_tab->blocksize)
    gt_assert(offset < popcount_tab->offsets[popcount_c] -
                       popcount_tab->offsets[popcount_c - 1]);
  else
    gt_assert(offset == 0);
  return gt_compact_ulong_store_get(popcount_tab->blocks,
                                    popcount_tab->offsets[popcount_c - 1] +
                                      offset);
}

size_t gt_popcount_tab_get_size(unsigned blocksize) {
  unsigned long num_of_blocks = 1UL << blocksize;
  size_t size = gt_compact_ulong_store_size(num_of_blocks, blocksize);
  size += sizeof (GtPopcountTab);
  size += sizeof (unsigned long) * blocksize;
  return size;
}

static inline unsigned popcount(unsigned long val)
{
// see page 11, Knuth TAOCP Vol 4 F1A
{
#ifdef __SSE4_2__
  return __builtin_popcountll(val);
#else
  uint64_t x = (uint64_t) val;
#ifdef POPCOUNT_TL
  return (unsigned long)
         B1CntBytes[x         & 0xFFULL] +
         B1CntBytes[(x >>  8) & 0xFFULL] +
         B1CntBytes[(x >> 16) & 0xFFULL] +
         B1CntBytes[(x >> 24) & 0xFFULL] +
         B1CntBytes[(x >> 32) & 0xFFULL] +
         B1CntBytes[(x >> 40) & 0xFFULL] +
         B1CntBytes[(x >> 48) & 0xFFULL] +
         B1CntBytes[(x >> 56) & 0xFFULL];
#else
  x = x - ((x >> 1) & (uint64_t) 0x5555555555555555);
  x = (x & (uint64_t) 0x3333333333333333) +
      ((x >> 2) & (uint64_t) 0x3333333333333333);
  x = (x + (x >> 4)) & (uint64_t) 0x0f0f0f0f0f0f0f0f;
  return (unsigned) ((uint64_t) 0x0101010101010101 * x >> 56);
#endif
#endif
}
}
static inline unsigned rank_1(GtPopcountTab *popcount_tab,
                              unsigned popcount_c,
                              unsigned long offset,
                              unsigned pos)
{
  unsigned long block;
  block = gt_compact_ulong_store_get(popcount_tab->blocks,
                                     popcount_tab->offsets[popcount_c - 1] +
                                       offset);
  block >>= popcount_tab->blocksize - pos - 1;
  return popcount(block);
}

unsigned gt_popcount_tab_rank_1(GtPopcountTab *popcount_tab,
                                     unsigned popcount_c,
                                     unsigned long offset,
                                     unsigned pos) {
  gt_assert(pos < popcount_tab->blocksize);
  gt_assert(popcount_c <= popcount_tab->blocksize);
  gt_assert(popcount_c != 0);
  if (popcount_c < popcount_tab->blocksize)
    gt_assert(offset < popcount_tab->offsets[popcount_c] -
                       popcount_tab->offsets[popcount_c - 1]);
  else {
    gt_assert(offset == 0);
    return popcount_c;
  }
  return rank_1(popcount_tab, popcount_c, offset, pos);
}
unsigned gt_popcount_tab_rank_0(GtPopcountTab *popcount_tab,
                                      unsigned popcount_c,
                                      unsigned long offset,
                                      unsigned pos) {
  gt_assert(pos < popcount_tab->blocksize);
  gt_assert(popcount_c < popcount_tab->blocksize);
  if (popcount_c == 0)
    return pos;
  return pos + 1 - rank_1(popcount_tab, popcount_c, offset, pos);
}

/*@unused@*/
void print_bin(unsigned long x, unsigned bits)
{
  unsigned idx;
  gt_assert(bits <= (unsigned) sizeof (unsigned long) * 8U);

  for (idx = 1U; idx <= bits; idx++) {
    printf("%lu", (x & (1UL << (bits - idx))) >> (bits - idx));
  }
  printf("\n");
}

int gt_popcount_tab_unit_test(GtError *err)
{
  int had_err = 0;
  unsigned long idx, jdx, popc_perm, init,
                blockmask = perm_start(16U);
  unsigned popcount_c;
  static const unsigned blocksize = 4U;
  static const unsigned long blocksize_four[] =
    {0, 1UL, 2UL, 4UL, 8UL, 3UL, 5UL, 6UL,
     9UL, 10UL, 12UL, 7UL, 11UL, 13UL, 14UL, 15UL};
  GtPopcountTab *popcount_t = gt_popcount_tab_new((unsigned) blocksize);

  for (idx = 0; idx < (1UL << blocksize); idx++) {
    gt_ensure(had_err, blocksize_four[idx] ==
              gt_compact_ulong_store_get(popcount_t->blocks, idx));
  }
  for (popcount_c = 0, idx = 0;
       !had_err && popcount_c <= blocksize;
       idx += gt_binomialCoeff_with_ln((unsigned long) blocksize,
                                       (unsigned long) popcount_c),
                                       popcount_c++) {
    for (jdx = 0;
         !had_err && jdx < gt_binomialCoeff_with_ln((unsigned long) blocksize,
                                                    (unsigned long) popcount_c);
         jdx++) {
      gt_ensure(had_err, blocksize_four[idx + jdx] ==
                         gt_popcount_tab_get(popcount_t, popcount_c, jdx));
    }
  }
  gt_ensure(had_err, gt_popcount_tab_rank_1(popcount_t, 2U, 0UL, 1U) == 0);
  gt_ensure(had_err, gt_popcount_tab_rank_0(popcount_t, 2U, 0UL, 1U) == 2U);
  gt_ensure(had_err, gt_popcount_tab_rank_1(popcount_t, 2U, 1UL, 0) == 0);
  gt_ensure(had_err, gt_popcount_tab_rank_0(popcount_t, 2U, 1UL, 0) == 1U);
  gt_ensure(had_err, gt_popcount_tab_rank_1(popcount_t, 2U, 1UL, 1U) == 1U);
  gt_ensure(had_err, gt_popcount_tab_rank_0(popcount_t, 2U, 1UL, 1U) == 1U);
  gt_ensure(had_err, gt_popcount_tab_rank_1(popcount_t, 2U, 1UL, 2U) == 1U);
  gt_ensure(had_err, gt_popcount_tab_rank_0(popcount_t, 2U, 1UL, 2U) == 2U);
  gt_popcount_tab_delete(popcount_t);

  popc_perm = init = perm_start(5U);
  while (!had_err && popc_perm >= init) {
    gt_ensure(had_err, 5U == popcount(popc_perm));
    popc_perm = next_perm(popc_perm) & blockmask;
  }
  return had_err;
}
