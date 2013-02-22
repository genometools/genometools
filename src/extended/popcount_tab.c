/*
  Copyright (c) 2012 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
  Copyright (c) 2012 Center for Bioinformatics, University of Hamburg

  Permission to use, copy, modify, and distribute this software for any
  purpose with or without fee offsets hereby granted, provided that the above
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
#include "core/mathsupport.h"
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
  unsigned int blocksize;
  unsigned long num_of_blocks;
  GtCompactUlongStore *blocks;
  unsigned long *offsets;
  /* this contains a mapping from a block to its offset within its class */
  GtCompactUlongStore *rev_blocks;
};

static void gt_popcount_tab_init_offset_tab(GtPopcountTab *popcount_tab)
{
  unsigned long idx, class_size,
                *offsets = popcount_tab->offsets,
                blocksize = (unsigned long) popcount_tab->blocksize,
                num_of_blocks = popcount_tab->num_of_blocks;

  offsets[0] = 0;
  offsets[1] = 1UL;
  offsets[2] = blocksize + 1UL;
  for (idx = 3UL; idx < blocksize - 1; idx++) {
    class_size = gt_combinatorics_binomial_ln(blocksize, idx - 1);
    offsets[idx] = offsets[idx - 1] + class_size;
  }
  offsets[blocksize - 1] = num_of_blocks - blocksize - 1;
  offsets[blocksize] = num_of_blocks - 1;
}

static unsigned long gt_popcount_tab_next_perm(unsigned long v)
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

static unsigned long gt_popcount_tab_perm_start(unsigned int bits)
{
  return (1UL << bits) - 1;
}

static unsigned long gt_popcount_tab_gen_blocks(unsigned int popcount_c,
                                                unsigned long idx,
                                                unsigned long blockmask,
                                                GtCompactUlongStore *blocks)
{
  unsigned long v, init;

  v = init = gt_popcount_tab_perm_start(popcount_c);
  while (v >= init) {
    gt_compact_ulong_store_update(blocks, idx++, v);
    if (popcount_c == 1U)
      v = (v << 1) & blockmask;
    else
      v = gt_popcount_tab_next_perm(v) & blockmask;
  }
  return idx;
}

static unsigned long gt_popcount_tab_init_blocks_tab(
                                                    GtCompactUlongStore *blocks,
                                                    unsigned int blocksize)
{
  unsigned long idx = 1UL,
                blockmask;
  unsigned int popcount_c = 1U;
  blockmask = gt_popcount_tab_perm_start(blocksize);
  gt_compact_ulong_store_update(blocks, 0, 0);
  while (popcount_c < blocksize) {
    idx = gt_popcount_tab_gen_blocks(popcount_c++, idx, blockmask, blocks);
  }
  gt_compact_ulong_store_update(blocks, idx++, blockmask);
  return idx;
}

GtPopcountTab *gt_popcount_tab_new(unsigned int blocksize)
{
  GtPopcountTab *popcount_tab;
  GT_UNUSED unsigned long idx_check;
  gt_assert(blocksize <= (unsigned) GT_INTWORDSIZE);

  popcount_tab = gt_malloc(sizeof (GtPopcountTab));
  popcount_tab->num_of_blocks = 1UL << blocksize;
  popcount_tab->blocksize = blocksize;
  popcount_tab->blocks = gt_compact_ulong_store_new(popcount_tab->num_of_blocks,
                                                    blocksize);
  popcount_tab->offsets = gt_malloc(sizeof (popcount_tab->offsets) *
                                      (blocksize + 1));
  gt_popcount_tab_init_offset_tab(popcount_tab);
  idx_check = gt_popcount_tab_init_blocks_tab(popcount_tab->blocks, blocksize);
  gt_assert(idx_check == popcount_tab->num_of_blocks);
  popcount_tab->rev_blocks = NULL;
  return popcount_tab;
}

void gt_popcount_tab_delete(GtPopcountTab *popcount_tab)
{
  if (popcount_tab != NULL) {
    gt_free(popcount_tab->offsets);
    gt_compact_ulong_store_delete(popcount_tab->blocks);
    gt_compact_ulong_store_delete(popcount_tab->rev_blocks);
    gt_free(popcount_tab);
  }
}

unsigned long gt_popcount_tab_get(GtPopcountTab *popcount_tab,
                                  unsigned int popcount_c,
                                  unsigned long i)
{
  gt_assert(popcount_c <= popcount_tab->blocksize);
  if (popcount_c == 0) {
    gt_assert(i == 0);
    return 0;
  }
  if (popcount_c < popcount_tab->blocksize)
    gt_assert(i < popcount_tab->offsets[popcount_c + 1] -
                    popcount_tab->offsets[popcount_c]);
  else
    gt_assert(i == 0);
  return gt_compact_ulong_store_get(popcount_tab->blocks,
                                    popcount_tab->offsets[popcount_c] + i);
}

size_t gt_popcount_tab_calculate_size(unsigned int blocksize) {
  unsigned long num_of_blocks = 1UL << blocksize;
  size_t size = gt_compact_ulong_store_size(num_of_blocks, blocksize);
  size += sizeof (GtPopcountTab);
  size += sizeof (unsigned long) * blocksize;
  return size;
}

static inline unsigned int gt_popcount_tab_popcount(unsigned long val)
{
#ifdef __SSE4_2__
  return __builtin_popcountl(val);
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
  /* see page 11, Knuth TAOCP Vol 4 F1A */
  x = x - ((x >> 1) & (uint64_t) 0x5555555555555555ULL);
  x = (x & (uint64_t) 0x3333333333333333ULL) +
      ((x >> 2) & (uint64_t) 0x3333333333333333ULL);
  x = (x + (x >> 4)) & (uint64_t) 0x0f0f0f0f0f0f0f0fULL;
  return (unsigned) ((uint64_t) 0x0101010101010101ULL * x >> 56);
#endif
#endif
}

unsigned int gt_popcount_tab_class(unsigned long block,
                                   GT_UNUSED unsigned int blocksize)
{
  gt_assert(block >> blocksize == 0);
  return gt_popcount_tab_popcount(block);
}

unsigned int gt_popcount_tab_offset_bits(unsigned int blocksize,
                                         unsigned int class)
{
  gt_assert(class <= blocksize);
#ifdef __SSE4__
  return (unsigned int) (sizeof (unsigned long) * CHAR_BIT -
      __builtin_clzl(gt_combinatorics_binomial_ln((unsigned long) blocksize,
                                                  (unsigned long) class)));
#else
  return gt_determinebitspervalue(
      gt_combinatorics_binomial_ln((unsigned long) blocksize,
                                   (unsigned long) class));
#endif
}

static inline unsigned int _gt_popcount_tab_rank_1(GtPopcountTab *popcount_tab,
                                                   unsigned int popcount_c,
                                                   unsigned long i,
                                                   unsigned int pos)
{
  unsigned long block;
  gt_assert(pos < popcount_tab->blocksize);
  gt_assert(popcount_c <= popcount_tab->blocksize);

  if (popcount_c == 0)
    return 0;
  if (popcount_c < popcount_tab->blocksize)
    gt_assert(i < popcount_tab->offsets[popcount_c + 1] -
                       popcount_tab->offsets[popcount_c]);
  else {
    gt_assert(i == 0);
    return pos+1;
  }
  block = gt_compact_ulong_store_get(popcount_tab->blocks,
                                     popcount_tab->offsets[popcount_c] + i);
  block >>= popcount_tab->blocksize - pos - 1;
  return gt_popcount_tab_popcount(block);
}

unsigned int gt_popcount_tab_rank_1(GtPopcountTab *popcount_tab,
                                    unsigned int popcount_c,
                                    unsigned long i,
                                    unsigned int pos)
{
  return _gt_popcount_tab_rank_1(popcount_tab, popcount_c, i, pos);
}

unsigned int gt_popcount_tab_rank_0(GtPopcountTab *popcount_tab,
                                    unsigned int popcount_c,
                                    unsigned long i,
                                    unsigned int pos)
{
  if (popcount_c == 0)
    return pos + 1;
  return pos + 1 - _gt_popcount_tab_rank_1(popcount_tab, popcount_c, i, pos);
}

static inline void gt_popcount_tab_init_rev_offset(GtPopcountTab *popcount_tab)
{
  unsigned long idx, current_offset = 0, current_block;
  unsigned int max_bits = gt_popcount_tab_offset_bits(
                                              popcount_tab->blocksize,
                                              GT_DIV2(popcount_tab->blocksize)),
               current_class = 0;
  max_bits++;
  popcount_tab->rev_blocks = gt_compact_ulong_store_new(
                                                    popcount_tab->num_of_blocks,
                                                    max_bits);
  for (idx = 0; idx < popcount_tab->num_of_blocks; idx++) {
    current_block = gt_compact_ulong_store_get(popcount_tab->blocks, idx);
    if (current_class < popcount_tab->blocksize &&
        idx == popcount_tab->offsets[current_class+1])
      current_class++;
    if (current_class != 0)
      current_offset = idx - popcount_tab->offsets[current_class];
    gt_compact_ulong_store_update(popcount_tab->rev_blocks,
                                  current_block,
                                  current_offset);
  }
}

unsigned long gt_popcount_tab_get_offset_for_block(GtPopcountTab *popcount_tab,
                                                   unsigned long block)
{
  gt_assert(popcount_tab != NULL);
  if (popcount_tab->rev_blocks == NULL)
    gt_popcount_tab_init_rev_offset(popcount_tab);
  return gt_compact_ulong_store_get(popcount_tab->rev_blocks, block);
}

int gt_popcount_tab_unit_test(GtError *err)
{
  int had_err = 0;
  unsigned long idx, jdx, popc_perm, init, offset,
                blockmask = gt_popcount_tab_perm_start(16U);
  unsigned int popcount_c, class_size;
  const unsigned int blocksize = 4U;
  const unsigned long blocksize_four[] =
    {0, 1UL, 2UL, 4UL, 8UL, 3UL, 5UL, 6UL,
     9UL, 10UL, 12UL, 7UL, 11UL, 13UL, 14UL, 15UL};
  GtPopcountTab *popcount_t = gt_popcount_tab_new((unsigned) blocksize);

  gt_error_check(err);

  for (idx = 0; idx < (1UL << blocksize); idx++) {
    gt_ensure(had_err, blocksize_four[idx] ==
              gt_compact_ulong_store_get(popcount_t->blocks, idx));
  }
  for (popcount_c = 0, idx = 0;
       !had_err && popcount_c <= blocksize;
       idx += class_size,
       popcount_c++) {
    class_size = (unsigned int)
      gt_combinatorics_binomial_ln((unsigned long) blocksize,
                                   (unsigned long) popcount_c);
    for (jdx = 0;
         !had_err && jdx < (unsigned long) class_size;
         jdx++) {
      gt_ensure(had_err, blocksize_four[idx + jdx] ==
                         gt_popcount_tab_get(popcount_t, popcount_c, jdx));
    }
  }
  gt_ensure(had_err, gt_popcount_tab_rank_1(popcount_t, 2U, 0UL, 1U) == 0U);
  gt_ensure(had_err, gt_popcount_tab_rank_0(popcount_t, 2U, 0UL, 1U) == 2U);
  gt_ensure(had_err, gt_popcount_tab_rank_1(popcount_t, 2U, 1UL, 0U) == 0U);
  gt_ensure(had_err, gt_popcount_tab_rank_0(popcount_t, 2U, 1UL, 0U) == 1U);
  gt_ensure(had_err, gt_popcount_tab_rank_1(popcount_t, 2U, 1UL, 1U) == 1U);
  gt_ensure(had_err, gt_popcount_tab_rank_0(popcount_t, 2U, 1UL, 1U) == 1U);
  gt_ensure(had_err, gt_popcount_tab_rank_1(popcount_t, 2U, 1UL, 2U) == 1U);
  gt_ensure(had_err, gt_popcount_tab_rank_0(popcount_t, 2U, 1UL, 2U) == 2U);
  gt_ensure(had_err, gt_popcount_tab_rank_1(popcount_t, 4U, 0UL, 1U) == 2U);

  for (idx = 0; !had_err && idx < (1UL << blocksize); idx++) {
    popcount_c = gt_popcount_tab_class(idx, blocksize);
    offset = gt_popcount_tab_get_offset_for_block(popcount_t, idx);

    for (jdx = 0; !had_err && jdx < (unsigned long) blocksize; jdx++) {
      gt_ensure(had_err,
                gt_popcount_tab_rank_1(popcount_t, popcount_c,
                                       offset, (unsigned int) jdx) +
                  gt_popcount_tab_rank_0(popcount_t, popcount_c,
                                         offset, (unsigned int) jdx) ==
                  (unsigned int) jdx + 1U);
    }
  }

  offset = gt_popcount_tab_get_offset_for_block(popcount_t, 4UL);
  gt_ensure(had_err, offset == 2UL);
  offset = gt_popcount_tab_get_offset_for_block(popcount_t, 15UL);
  gt_ensure(had_err, offset == 0);
  gt_popcount_tab_delete(popcount_t);

  popcount_t = gt_popcount_tab_new(10U);
  for (idx = 0; !had_err && idx < 1UL<<10UL; idx++) {
    offset = gt_popcount_tab_get_offset_for_block(popcount_t, idx);
    popcount_c = gt_popcount_tab_popcount(idx);
    jdx = gt_popcount_tab_get(popcount_t, popcount_c, offset);
    gt_ensure(had_err, idx == jdx);
  }
  gt_popcount_tab_delete(popcount_t);
  popc_perm = init = gt_popcount_tab_perm_start(5U);
  while (!had_err && popc_perm >= init) {
    gt_ensure(had_err, gt_popcount_tab_popcount(popc_perm) == 5U);
    popc_perm = gt_popcount_tab_next_perm(popc_perm) & blockmask;
  }
  return had_err;
}
