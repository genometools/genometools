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

#include <limits.h>

#include "core/assert_api.h"
#include "core/byte_popcount_api.h"
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

struct GtPopcountTab
{
  unsigned int        *bit_sizes;
  /* the actual blocks in popcount order */
  GtCompactUlongStore *blocks,
                      *offsets,
  /* this contains a mapping from a block to its offset within its class */
                      *rev_blocks;
  GtUword        num_of_blocks;
  unsigned int         blocksize;
};

static void gt_popcount_tab_init_offset_tab(GtPopcountTab *popcount_tab)
{
  GtUword idx, class_size,
                blocksize = (GtUword) popcount_tab->blocksize,
                num_of_blocks = popcount_tab->num_of_blocks;
  GtCompactUlongStore *offsets = popcount_tab->offsets;

  gt_compact_ulong_store_update(offsets, 0, 0);
  gt_compact_ulong_store_update(offsets, 1UL, 1UL);
  gt_compact_ulong_store_update(offsets, 2UL, blocksize + 1UL);
  for (idx = 3UL; idx < blocksize - 1; idx++) {
    class_size = gt_combinatorics_binomial_ln(blocksize, idx - 1);
    gt_compact_ulong_store_update(offsets, idx,
                                  class_size + gt_compact_ulong_store_get(
                                                                      offsets,
                                                                      idx - 1));

  }
  gt_compact_ulong_store_update(offsets, blocksize - 1,
                                num_of_blocks - blocksize - 1);
  gt_compact_ulong_store_update(offsets, blocksize, num_of_blocks - 1);
}

static GtUword gt_popcount_tab_next_perm(GtUword v)
{
  GtUword head, tail;

  head = (v | (v - 1)) + 1;
  tail = v & (((head & -head) >> 1) -1);
  if (tail != 0) {
#ifdef __SSE__
    int zero_trail = __builtin_ctzl(tail);
    tail >>= zero_trail;
#else
    while (( tail & 1 ) == 0)
      tail = tail >> 1;
#endif
  }
  return head | tail;
}

static GtUword gt_popcount_tab_perm_start(unsigned int bits)
{
  return (1UL << bits) - 1;
}

static GtUword gt_popcount_tab_gen_blocks(unsigned int popcount_c,
                                                GtUword idx,
                                                GtUword blockmask,
                                                GtCompactUlongStore *blocks)
{
  GtUword v, init;

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

static GtUword
gt_popcount_tab_init_blocks_tab(GtCompactUlongStore *blocks,
                                unsigned int blocksize)
{
  GtUword idx = 1UL,
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
  GT_UNUSED GtUword idx_check;
  gt_assert(blocksize <= (unsigned) GT_INTWORDSIZE);

  popcount_tab = gt_malloc(sizeof (GtPopcountTab));
  popcount_tab->rev_blocks = NULL;
  popcount_tab->bit_sizes = NULL;
  popcount_tab->num_of_blocks = 1UL << blocksize;
  popcount_tab->blocksize = blocksize;

  popcount_tab->blocks = gt_compact_ulong_store_new(popcount_tab->num_of_blocks,
                                                    blocksize);
  popcount_tab->offsets =
    gt_compact_ulong_store_new((GtUword) blocksize + 1, blocksize);
  gt_popcount_tab_init_offset_tab(popcount_tab);
  idx_check = gt_popcount_tab_init_blocks_tab(popcount_tab->blocks, blocksize);
  gt_assert(idx_check == popcount_tab->num_of_blocks);
  return popcount_tab;
}

void gt_popcount_tab_delete(GtPopcountTab *popcount_tab)
{
  if (popcount_tab != NULL) {
    gt_compact_ulong_store_delete(popcount_tab->blocks);
    gt_compact_ulong_store_delete(popcount_tab->offsets);
    gt_compact_ulong_store_delete(popcount_tab->rev_blocks);
    gt_free(popcount_tab->bit_sizes);
    gt_free(popcount_tab);
  }
}

GtUword gt_popcount_tab_get(GtPopcountTab *popcount_tab,
                                  unsigned int popcount_c,
                                  GtUword i)
{
  gt_assert(popcount_c <= popcount_tab->blocksize);
  if (popcount_c == 0) {
    gt_assert(i == 0);
    return 0;
  }
  if (popcount_c < popcount_tab->blocksize)
    gt_assert(i <
              gt_compact_ulong_store_get(popcount_tab->offsets,
                                         (GtUword) popcount_c + 1) -
              gt_compact_ulong_store_get(popcount_tab->offsets,
                                         (GtUword) popcount_c));
  else
    gt_assert(i == 0);
  return gt_compact_ulong_store_get(
                    popcount_tab->blocks,
                    gt_compact_ulong_store_get(popcount_tab->offsets,
                                               (GtUword) popcount_c) + i);
}

size_t gt_popcount_tab_calculate_size(unsigned int blocksize)
{
  GtUword num_of_blocks = 1UL << blocksize;
  size_t size = gt_compact_ulong_store_size(num_of_blocks, blocksize);
  size += sizeof (GtPopcountTab);
  size += sizeof (GtUword) * blocksize;
  return size;
}

static inline unsigned int gt_popcount_tab_popcount(GtUword val)
{
#ifdef __SSE4_2__
  return __builtin_popcountl(val);
#else
  uint64_t x = (uint64_t) val;
#ifdef POPCOUNT_TL
  return (GtUword)
         gt_byte_popcount[x         & 0xFFULL] +
         gt_byte_popcount[(x >>  8) & 0xFFULL] +
         gt_byte_popcount[(x >> 16) & 0xFFULL] +
         gt_byte_popcount[(x >> 24) & 0xFFULL] +
         gt_byte_popcount[(x >> 32) & 0xFFULL] +
         gt_byte_popcount[(x >> 40) & 0xFFULL] +
         gt_byte_popcount[(x >> 48) & 0xFFULL] +
         gt_byte_popcount[(x >> 56) & 0xFFULL];
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

unsigned int gt_popcount_tab_class(GT_UNUSED GtPopcountTab *popcount_tab,
                                   GtUword block)
{
  gt_assert(popcount_tab != NULL);
  gt_assert(block >> popcount_tab->blocksize == 0);
  return gt_popcount_tab_popcount(block);
}

static void gt_popcount_tab_init_bit_sizes(unsigned int *bit_sizes,
                                           unsigned int blocksize)
{
  unsigned int class;
  for (class = 0; class <= blocksize; ++class) {
#ifdef __SSE4_2__
    bit_sizes[class] = (unsigned int) (sizeof (GtUword) * CHAR_BIT -
       __builtin_clzl(gt_combinatorics_binomial_dp((GtUword) blocksize,
                                                   (GtUword) class)));
#else
    bit_sizes[class] = gt_determinebitspervalue(
        gt_combinatorics_binomial_dp((GtUword) blocksize,
                                     (GtUword) class));
#endif
  }
}

unsigned int gt_popcount_tab_offset_bits(GtPopcountTab *popcount_tab,
                                         unsigned int class)
{
  gt_assert(popcount_tab != NULL);
  gt_assert(class <= popcount_tab->blocksize);
  if (popcount_tab->bit_sizes == NULL) {
    popcount_tab->bit_sizes = gt_calloc((size_t) popcount_tab->blocksize + 1,
                                        sizeof (*(popcount_tab->bit_sizes)));
    gt_popcount_tab_init_bit_sizes(popcount_tab->bit_sizes,
                                   popcount_tab->blocksize);
  }
  return popcount_tab->bit_sizes[class];
}

unsigned int gt_popcount_tab_rank_1(GtPopcountTab *popcount_tab,
                                    unsigned int popcount_c,
                                    GtUword i,
                                    unsigned int pos)
{
  GtUword block;
  gt_assert(pos < popcount_tab->blocksize);
  gt_assert(popcount_c <= popcount_tab->blocksize);

  if (popcount_c == 0)
    return 0;
  if (popcount_c < popcount_tab->blocksize)
    gt_assert(i <
              gt_compact_ulong_store_get(popcount_tab->offsets,
                                         (GtUword) popcount_c + 1) -
              gt_compact_ulong_store_get(popcount_tab->offsets,
                                         (GtUword) popcount_c));
  else {
    gt_assert(i == 0);
    return pos+1;
  }
  block =
    gt_compact_ulong_store_get(
                    popcount_tab->blocks,
                    gt_compact_ulong_store_get(popcount_tab->offsets,
                                               (GtUword) popcount_c) + i);
  block >>= popcount_tab->blocksize - pos - 1;
  return gt_popcount_tab_popcount(block);
}

unsigned int gt_popcount_tab_rank_0(GtPopcountTab *popcount_tab,
                                    unsigned int popcount_c,
                                    GtUword i,
                                    unsigned int pos)
{
  if (popcount_c == 0)
    return pos + 1;
  return pos + 1 - gt_popcount_tab_rank_1(popcount_tab, popcount_c, i, pos);
}

void gt_popcount_tab_block_to_str(GtPopcountTab *popcount_tab,
                                  GtUword block,
                                  char *buffer)
{
  unsigned int idx;
  buffer[popcount_tab->blocksize] = 0;
  for (idx = 0; idx < popcount_tab->blocksize; ++idx) {
    buffer[idx] =
      ((block >> (GtUword) (popcount_tab->blocksize - idx - 1)) &
       (GtUword) 1) ==
      (GtUword) 1 ?
      '1' : '0';
  }
}

static inline void gt_popcount_tab_init_rev_offset(GtPopcountTab *popcount_tab)
{
  GtUword idx, current_offset = 0, current_block;
  unsigned int max_bits =
    gt_popcount_tab_offset_bits(popcount_tab, GT_DIV2(popcount_tab->blocksize)),
               current_class = 0;
  popcount_tab->rev_blocks =
    gt_compact_ulong_store_new(popcount_tab->num_of_blocks, max_bits);
  for (idx = 0; idx < popcount_tab->num_of_blocks; idx++) {
    current_block = gt_compact_ulong_store_get(popcount_tab->blocks, idx);
    if (current_class < popcount_tab->blocksize &&
        idx ==
          gt_compact_ulong_store_get(popcount_tab->offsets,
                                     (GtUword) current_class + 1))
      current_class++;
    if (current_class != 0)
      current_offset =
        idx - gt_compact_ulong_store_get(popcount_tab->offsets,
                                         (GtUword) current_class);
    gt_compact_ulong_store_update(popcount_tab->rev_blocks,
                                  current_block,
                                  current_offset);
  }
}

GtUword gt_popcount_tab_get_offset_for_block(GtPopcountTab *popcount_tab,
                                                   GtUword block)
{
  gt_assert(popcount_tab != NULL);
  if (popcount_tab->rev_blocks == NULL)
    gt_popcount_tab_init_rev_offset(popcount_tab);
  return gt_compact_ulong_store_get(popcount_tab->rev_blocks, block);
}

int gt_popcount_tab_unit_test(GtError *err)
{
  int had_err = 0;
  GtUword idx, jdx, popc_perm, init, offset,
                blockmask = gt_popcount_tab_perm_start(16U);
  unsigned int popcount_c, class_size;
  const unsigned int blocksize = 4U;
  const GtUword blocksize_four[] =
    {0, 1UL, 2UL, 4UL, 8UL, 3UL, 5UL, 6UL,
     9UL, 10UL, 12UL, 7UL, 11UL, 13UL, 14UL, 15UL};
  const unsigned int offsets_bits[] = {1U, 3U, 3U, 3U, 1U};
  GtPopcountTab *popcount_t = gt_popcount_tab_new((unsigned) blocksize);

  gt_error_check(err);

  for (idx = 0; !had_err && idx <= (GtUword) blocksize; idx++) {
    gt_ensure(
              gt_popcount_tab_offset_bits(popcount_t, (unsigned int) idx) ==
                offsets_bits[idx]);
  }

  for (idx = 0; !had_err && idx < (1UL << blocksize); idx++) {
    gt_ensure(blocksize_four[idx] ==
              gt_compact_ulong_store_get(popcount_t->blocks, idx));
  }
  for (popcount_c = 0, idx = 0;
       !had_err && popcount_c <= blocksize;
       idx += class_size,
       popcount_c++) {
    class_size = (unsigned int)
      gt_combinatorics_binomial_ln((GtUword) blocksize,
                                   (GtUword) popcount_c);
    for (jdx = 0;
         !had_err && jdx < (GtUword) class_size;
         jdx++) {
      gt_ensure(idx ==
                gt_compact_ulong_store_get(popcount_t->offsets,
                                           (GtUword) popcount_c));
      gt_ensure(
                blocksize_four[idx + jdx] ==
                          (GtUword) gt_popcount_tab_get(popcount_t,
                                                              popcount_c, jdx));
    }
  }
  gt_ensure(gt_popcount_tab_rank_1(popcount_t, 2U, 0UL, 1U) == 0U);
  gt_ensure(gt_popcount_tab_rank_0(popcount_t, 2U, 0UL, 1U) == 2U);
  gt_ensure(gt_popcount_tab_rank_1(popcount_t, 2U, 1UL, 0U) == 0U);
  gt_ensure(gt_popcount_tab_rank_0(popcount_t, 2U, 1UL, 0U) == 1U);
  gt_ensure(gt_popcount_tab_rank_1(popcount_t, 2U, 1UL, 1U) == 1U);
  gt_ensure(gt_popcount_tab_rank_0(popcount_t, 2U, 1UL, 1U) == 1U);
  gt_ensure(gt_popcount_tab_rank_1(popcount_t, 2U, 1UL, 2U) == 1U);
  gt_ensure(gt_popcount_tab_rank_0(popcount_t, 2U, 1UL, 2U) == 2U);
  gt_ensure(gt_popcount_tab_rank_1(popcount_t, 4U, 0UL, 1U) == 2U);

  for (idx = 0; !had_err && idx < (1UL << blocksize); idx++) {
    popcount_c = gt_popcount_tab_class(popcount_t, idx);
    offset = gt_popcount_tab_get_offset_for_block(popcount_t, idx);

    for (jdx = 0; !had_err && jdx < (GtUword) blocksize; jdx++) {
      gt_ensure(
                gt_popcount_tab_rank_1(popcount_t, popcount_c,
                                       offset, (unsigned int) jdx) +
                  gt_popcount_tab_rank_0(popcount_t, popcount_c,
                                         offset, (unsigned int) jdx) ==
                  (unsigned int) jdx + 1U);
    }
  }

  if (!had_err) {
    offset = gt_popcount_tab_get_offset_for_block(popcount_t, 4UL);
    gt_ensure(offset == 2UL);
    offset = gt_popcount_tab_get_offset_for_block(popcount_t, 15UL);
    gt_ensure(offset == 0);
  }
  gt_popcount_tab_delete(popcount_t);

  if (!had_err) {
    popcount_t = gt_popcount_tab_new(10U);
    for (idx = 0; !had_err && idx < 1UL<<10UL; idx++) {
      offset = gt_popcount_tab_get_offset_for_block(popcount_t, idx);
      popcount_c = gt_popcount_tab_popcount(idx);
      jdx = gt_popcount_tab_get(popcount_t, popcount_c, offset);
      gt_ensure(idx == jdx);
    }
    gt_popcount_tab_delete(popcount_t);
    popc_perm = init = gt_popcount_tab_perm_start(5U);
    while (!had_err && popc_perm >= init) {
      gt_ensure(gt_popcount_tab_popcount(popc_perm) == 5U);
      popc_perm = gt_popcount_tab_next_perm(popc_perm) & blockmask;
    }
  }
  return had_err;
}
