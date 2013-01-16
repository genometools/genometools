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

#include "extended/compressed_bitsequence.h"
#include "core/ma_api.h"
#include "core/ensure.h"
#include "core/compact_ulong_store.h"
#include "extended/popcount_tab.h"
#include "core/combinatorics.h"
#include "core/safearith.h"
#include "core/unused_api.h"
#include "core/mathsupport.h"

/* this seems to be a good default value. maybe change this in the future */
#define GT_COMP_BITSEQ_BLOCKSIZE 15U

typedef struct {
  unsigned int class;
  unsigned int block_len;
  unsigned long rank_sum;
  unsigned long idx;
  unsigned long block_offset;
} GtCompressedBitsequenceBlockInfo;

/* static void debug_cbs_bi(GtCompressedBitsequenceBlockInfo *cbs_bi)
{
  printf("cls: %u,", cbs_bi->class);
  printf(" blklen: %u,", cbs_bi->block_len);
  printf(" ofst: %lu,", cbs_bi->block_offset);
  printf(" ranks: %lu,", cbs_bi->rank_sum);
  printf(" idx: %lu\n", cbs_bi->idx);
} */

struct GtCompressedBitsequence
{
  GtPopcountTab *popcount_tab;
  GtCompactUlongStore *classes;
  GtBitsequence *c_offsets;
  GtCompressedBitsequenceBlockInfo *cbs_bi;
  unsigned long *superblockranks;
  unsigned long *superblockoffsets;
  unsigned long num_of_blocks;
  unsigned long c_offsets_size;
  unsigned long num_of_bits;
  unsigned long num_of_superblocks;
  unsigned int blocksize;
  unsigned int superblocksize;
  unsigned int last_block_len;
};

/* Assumes most significant bits are filled first in the vector, returns bits in
   <len> least significant positions */
static inline GtBitsequence gt_compressed_bitsequence_get_variable_field(
                                                          GtBitsequence *bitseq,
                                                          unsigned long start,
                                                          unsigned int len)
{
  GtBitsequence result = 0;
  gt_assert(len <= (unsigned int) GT_INTWORDSIZE);
  if (len != 0) {
    const unsigned int bit_offset = (unsigned int) GT_MODWORDSIZE(start);
    const unsigned long word_offset = GT_DIVWORDSIZE(start);

    result = bitseq[word_offset] << (GtBitsequence) bit_offset;
    if ((unsigned int) bit_offset + len > (unsigned int) GT_INTWORDSIZE)
      result |= bitseq[word_offset + 1] >>
                                  (GtBitsequence) (GT_INTWORDSIZE - bit_offset);
    result >>= (GtBitsequence) (GT_INTWORDSIZE-len);
  }
  return result;
}

/* Assumes most significant bits are filled first in the vector and value is
   filled in <len> least significant bits (as usualy numbers are) */
static inline void gt_compressed_bitsequence_set_variable_field(
                                                          GtBitsequence *bitseq,
                                                          unsigned long start,
                                                          unsigned int len,
                                                          GtBitsequence value)
{
  const unsigned int bit_offset = (unsigned int) GT_MODWORDSIZE(start);
  unsigned int bits_left = (unsigned int) (GT_INTWORDSIZE - bit_offset);
  const unsigned long word_offset = GT_DIVWORDSIZE(start);
  GtBitsequence mask = ~((GtBitsequence) 0) << (GtBitsequence) bits_left;

  if (len > bits_left) {
    unsigned int overhang = len - bits_left;

    bitseq[word_offset] &= mask;
    bitseq[word_offset] |= (value >> (GtBitsequence) overhang);

    mask = ~((GtBitsequence) 0) >> (GtBitsequence) overhang;
    bitseq[word_offset + 1] &= mask;
    bitseq[word_offset + 1] |= value << (GtBitsequence)
                                                    (GT_INTWORDSIZE - overhang);
  }
  else {
    mask |= ~((GtBitsequence) 0) >> (GtBitsequence) (bit_offset + len);
    bitseq[word_offset] &= mask;
    bitseq[word_offset] |= value << (GtBitsequence) (bits_left - len);
  }
}

static void gt_compressed_bitsequence_fill_tabs(GtCompressedBitsequence *cbs,
                                                GtBitsequence *bitseq)
{
  unsigned int class,
               samplecount = 0,
               c_offset_bits,
               block_len = cbs->blocksize;
  unsigned long idx,
                rank_sum = 0,
                bitpos = 0,
                samplenum = 0;
  GtBitsequence block, offsets_bitpos;

  for (idx = 0; idx < cbs->num_of_blocks; idx++) {
    if (idx == cbs->num_of_blocks -1) {
      block_len = cbs->last_block_len;
    }
    class = (unsigned int) gt_compact_ulong_store_get(cbs->classes, idx);
    c_offset_bits = gt_popcount_tab_offset_bits(cbs->popcount_tab, class);
    block = gt_compressed_bitsequence_get_variable_field(bitseq,
                                                         idx * cbs->blocksize,
                                                         block_len);
    offsets_bitpos = (GtBitsequence) gt_popcount_tab_get_offset_for_block(
                                                         cbs->popcount_tab,
                                                         (unsigned long) block);
    gt_compressed_bitsequence_set_variable_field(cbs->c_offsets, bitpos,
                                                 c_offset_bits, offsets_bitpos);
    rank_sum += class;
    bitpos += c_offset_bits;
    samplecount++;
    if (cbs->superblocksize == samplecount) {
      samplecount = 0;
      cbs->superblockranks[samplenum] = rank_sum;
      cbs->superblockoffsets[samplenum] = bitpos;
      samplenum++;
    }
  }
  if (samplecount != 0) {
    cbs->superblockranks[samplenum] = rank_sum;
    cbs->superblockoffsets[samplenum] = bitpos;
  }
}

static void gt_compressed_bitsequence_init_s_tabs(GtCompressedBitsequence *cbs)
{
  cbs->num_of_superblocks = cbs->num_of_blocks / cbs->superblocksize;
  if (cbs->num_of_blocks % cbs->superblocksize != 0)
    cbs->num_of_superblocks++;
  cbs->superblockranks = gt_calloc((size_t) cbs->num_of_superblocks,
                                   sizeof (cbs->superblockranks));
  cbs->superblockoffsets = gt_calloc((size_t) cbs->num_of_superblocks,
                                     sizeof (cbs->superblockoffsets));
}

static void gt_compressed_bitsequence_init_c_tab(GtCompressedBitsequence *cbs)
{
  cbs->num_of_blocks = cbs->num_of_bits / cbs->blocksize;
  if (cbs->num_of_bits % cbs->blocksize != 0)
    cbs->num_of_blocks++;
  cbs->classes = gt_compact_ulong_store_new(cbs->num_of_blocks, cbs->blocksize);
}

static void gt_compressed_bitsequence_fill_c_tab_calc_o_size(
                                                   GtCompressedBitsequence *cbs,
                                                   GtBitsequence *bitseq)
{
  unsigned int block_len = cbs->blocksize,
               current_class;
  unsigned long idx, o_size = 0;
  GtBitsequence current_blk;

  for (idx = 0; idx < cbs->num_of_blocks; idx++) {
    if (idx == cbs->num_of_blocks -1) {
      block_len = cbs->last_block_len;
    }
    current_blk = gt_compressed_bitsequence_get_variable_field(
                                                           bitseq,
                                                           idx * cbs->blocksize,
                                                           block_len);
    current_class = gt_popcount_tab_class(cbs->popcount_tab,
                                          (unsigned long) current_blk);
    gt_compact_ulong_store_update(cbs->classes, idx,
                                  (unsigned long) current_class);
    o_size += gt_popcount_tab_offset_bits(cbs->popcount_tab, current_class);
  }
  cbs->c_offsets_size = GT_DIVWORDSIZE(o_size);
  if (GT_MODWORDSIZE(o_size))
    cbs->c_offsets_size++;
}

GtCompressedBitsequence *gt_compressed_bitsequence_new(
                                                      GtBitsequence *bitseq,
                                                      unsigned int samplerate,
                                                      unsigned long num_of_bits)
{
  GtCompressedBitsequence *cbs;

  cbs = gt_malloc(sizeof (*cbs));
  cbs->cbs_bi = NULL;
  cbs->blocksize = GT_COMP_BITSEQ_BLOCKSIZE;
  cbs->superblocksize = samplerate;
  cbs->popcount_tab = gt_popcount_tab_new(cbs->blocksize);
  cbs->c_offsets_size = 0;
  cbs->num_of_bits = num_of_bits;
  cbs->last_block_len = (unsigned int) (cbs->num_of_bits % cbs->blocksize);
  if (cbs->last_block_len == 0)
    cbs->last_block_len = cbs->blocksize;
  gt_compressed_bitsequence_init_c_tab(cbs);
  gt_compressed_bitsequence_fill_c_tab_calc_o_size(cbs, bitseq);
  cbs->c_offsets = gt_calloc((size_t) cbs->c_offsets_size,
                             sizeof (*cbs->c_offsets));
  gt_compressed_bitsequence_init_s_tabs(cbs);
  gt_compressed_bitsequence_fill_tabs(cbs, bitseq);
  return cbs;
}

static inline void gt_compressed_bitsequence_calc_block_info(
                                                   GtCompressedBitsequence *cbs,
                                                   unsigned long position)
{
  unsigned long idx;
  GtCompressedBitsequenceBlockInfo *bi = cbs->cbs_bi;

  idx = position / cbs->blocksize;

  if (cbs->cbs_bi == NULL) {
    cbs->cbs_bi = gt_malloc(sizeof (*cbs->cbs_bi));
    bi = cbs->cbs_bi;
    bi->idx = idx + 1;
  }

  if (idx != bi->idx) {
    unsigned int offset_bits;
    unsigned long jdx, sample, offsets_bitpos;

    bi->idx = idx;
    bi->block_len = cbs->blocksize;
    if (idx == cbs->num_of_blocks -1)
      bi->block_len = cbs->last_block_len;

    sample = idx / cbs->superblocksize;
    if (sample==0) {
      offsets_bitpos = 0;
      bi->rank_sum = 0;
    }
    else {
      offsets_bitpos = cbs->superblockoffsets[sample - 1];
      bi->rank_sum = cbs->superblockranks[sample - 1];
    }
    for (jdx = sample * cbs->superblocksize; jdx < idx; jdx++) {
      bi->class = (unsigned int) gt_compact_ulong_store_get(cbs->classes, jdx);
      bi->rank_sum += bi->class;
      offsets_bitpos += gt_popcount_tab_offset_bits(cbs->popcount_tab,
                                                    bi->class);
    }
    bi->class = (unsigned int) gt_compact_ulong_store_get(cbs->classes, idx);
    offset_bits = gt_popcount_tab_offset_bits(cbs->popcount_tab, bi->class);
    bi->block_offset = (unsigned long)
      gt_compressed_bitsequence_get_variable_field(cbs->c_offsets,
                                                   offsets_bitpos,
                                                   offset_bits);
  }
}

int gt_compressed_bitsequence_access(GtCompressedBitsequence *cbs,
                                     unsigned long position)
{
  int bit;
  unsigned int pos_in_block;
  unsigned long block;
  GtCompressedBitsequenceBlockInfo *cbs_bi;

  gt_assert(cbs != NULL);
  gt_assert(position < cbs->num_of_bits);

  pos_in_block = (unsigned int) position % cbs->blocksize;

  gt_compressed_bitsequence_calc_block_info(cbs, position);
  cbs_bi = cbs->cbs_bi;

  block = gt_popcount_tab_get(cbs->popcount_tab,
                              cbs_bi->class,
                              cbs_bi->block_offset);

  /* account for shorter blocks at the end */
  block <<= (cbs->blocksize - cbs_bi->block_len);

  bit = (int) ((block >> (cbs->blocksize - pos_in_block - 1)) & 1UL);
  return bit;
}

unsigned long gt_compressed_bitsequence_rank_1(GtCompressedBitsequence *cbs,
                                               unsigned long position)
{
  unsigned int pos_in_block;
  GtCompressedBitsequenceBlockInfo *cbs_bi;

  gt_assert(cbs != NULL);
  gt_assert(position < cbs->num_of_bits);

  pos_in_block = (unsigned int) position % cbs->blocksize;

  gt_compressed_bitsequence_calc_block_info(cbs, position);
  cbs_bi = cbs->cbs_bi;

  pos_in_block += cbs->blocksize - cbs_bi->block_len;
  return cbs_bi->rank_sum + gt_popcount_tab_rank_1(cbs->popcount_tab,
                                                   cbs_bi->class,
                                                   cbs_bi->block_offset,
                                                   pos_in_block);
}

unsigned long gt_compressed_bitsequence_rank_0(GtCompressedBitsequence *cbs,
                                               unsigned long position)
{
  unsigned int pos_in_block;
  unsigned long rank_sum0;
  GtCompressedBitsequenceBlockInfo *cbs_bi;

  gt_assert(cbs != NULL);
  gt_assert(position < cbs->num_of_bits);

  pos_in_block = (unsigned int) position % cbs->blocksize;

  gt_compressed_bitsequence_calc_block_info(cbs, position);
  cbs_bi = cbs->cbs_bi;

  pos_in_block += cbs->blocksize - cbs_bi->block_len;
  rank_sum0 = ((cbs_bi->idx) * cbs->blocksize) - cbs_bi->rank_sum;
  rank_sum0 -= cbs->blocksize - cbs_bi->block_len;
  return rank_sum0 + gt_popcount_tab_rank_0(cbs->popcount_tab,
                                            cbs_bi->class,
                                            cbs_bi->block_offset,
                                            pos_in_block);
}

unsigned long gt_compressed_bitsequence_select_1(GtCompressedBitsequence *cbs,
                                                 unsigned long num)
{
  unsigned int class = cbs->blocksize + 1;
  unsigned long start_s_block, middle_s_block, end_s_block,
                containing_s_block,
                block_idx, blocks_offset_pos, block,
                rank_sum = 0,
                select_mask,
                position;

  gt_assert(num != 0);
  gt_assert(cbs != NULL);
  gt_assert(num < cbs->num_of_bits);

  if (num > cbs->superblockranks[cbs->num_of_superblocks - 1])
    return cbs->num_of_bits;
  if (num <= cbs->superblockranks[0]) {
    rank_sum = 0;
    containing_s_block = 0;
    blocks_offset_pos = 0;
  }
  else {
    /* search for superblock */
    start_s_block = num / (cbs->blocksize * cbs->superblocksize);
    end_s_block = cbs->num_of_superblocks;
    middle_s_block = GT_DIV2(start_s_block + end_s_block);
    while (start_s_block < end_s_block) {
      if (cbs->superblockranks[middle_s_block] < num) {
        if (cbs->superblockranks[middle_s_block + 1] < num)
          start_s_block = middle_s_block;
        else
          break;
      }
      else {
        if (cbs->superblockranks[middle_s_block - 1] >= num)
          end_s_block = middle_s_block;
        else {
          middle_s_block--;
          break;
        }
      }
      middle_s_block = GT_DIV2(start_s_block + end_s_block);
    }
    rank_sum = cbs->superblockranks[middle_s_block];
    blocks_offset_pos = cbs->superblockoffsets[middle_s_block];
    containing_s_block = middle_s_block+1;
  }

  /* search within superblock */
  for (block_idx = containing_s_block * cbs->superblocksize;
       block_idx < cbs->num_of_blocks;
       block_idx++) {
    class = (unsigned int) gt_compact_ulong_store_get(cbs->classes, block_idx);
    if (num <= rank_sum + class)
      break;
    blocks_offset_pos += gt_popcount_tab_offset_bits(cbs->popcount_tab, class);
    rank_sum += class;
  }
  position = block_idx * cbs->blocksize;
  gt_assert(class != cbs->blocksize + 1);

  /* search within block */
  if (block_idx != cbs->num_of_blocks - 1)
    select_mask = 1UL << (cbs->blocksize - 1UL);
  else
    select_mask = 1UL << (cbs->last_block_len - 1UL);

  block =
    gt_popcount_tab_get(cbs->popcount_tab, class,
                        (unsigned long)
                          gt_compressed_bitsequence_get_variable_field(
                                  cbs->c_offsets, blocks_offset_pos,
                                  gt_popcount_tab_offset_bits(cbs->popcount_tab,
                                                              class)));
  while (rank_sum < num && select_mask != 0) {
    rank_sum += (block & select_mask) != 0 ? 1 : 0;
    select_mask >>= 1;
    position++;
  }
  gt_assert(rank_sum == num);

  return position - 1;
}

unsigned long gt_compressed_bitsequence_select_0(GtCompressedBitsequence *cbs,
                                                 unsigned long num)
{
  unsigned int class = cbs->blocksize + 1;
  unsigned long start_s_block, middle_s_block, end_s_block,
                containing_s_block,
                block_idx, blocks_offset_pos, block,
                rank_sum = 0,
                select_mask,
                position,
                max_0_rank,
                first_0_superblock_rank,
                s_block_bits;

  gt_assert(num != 0);
  gt_assert(cbs != NULL);
  gt_assert(num < cbs->num_of_bits);

  s_block_bits = (unsigned long) cbs->blocksize * cbs->superblocksize;
  max_0_rank =
    cbs->num_of_bits - cbs->superblockranks[cbs->num_of_superblocks - 1];
  if (num > max_0_rank)
    return cbs->num_of_bits;

  first_0_superblock_rank = s_block_bits - cbs->superblockranks[0];
  if (num <= first_0_superblock_rank) {
    rank_sum = 0;
    containing_s_block = 0;
    blocks_offset_pos = 0;
  }
  else {
    unsigned long middle_s_block_max,
                  middle_s_block_next_max,
                  middle_s_block_prev_max;
    /* search for superblock */
    start_s_block = num / s_block_bits; /* cannot be smaller than this */
    end_s_block = cbs->num_of_superblocks;
    middle_s_block = GT_DIV2(start_s_block + end_s_block);
    while (start_s_block < end_s_block) {
      middle_s_block_max = s_block_bits * (middle_s_block+1);
      middle_s_block_next_max = middle_s_block_max + s_block_bits;
      middle_s_block_prev_max = middle_s_block_max - s_block_bits;
      if ((middle_s_block_max - cbs->superblockranks[middle_s_block]) < num) {
        if ((middle_s_block_next_max -
              cbs->superblockranks[middle_s_block + 1]) < num)
          start_s_block = middle_s_block;
        else
          break;
      }
      else {
        if ((middle_s_block_prev_max -
               cbs->superblockranks[middle_s_block - 1]) >= num)
          end_s_block = middle_s_block;
        else {
          middle_s_block--;
          break;
        }
      }
      middle_s_block = GT_DIV2(start_s_block + end_s_block);
    }
    middle_s_block_max = s_block_bits * (middle_s_block + 1);
    rank_sum = middle_s_block_max - cbs->superblockranks[middle_s_block];
    blocks_offset_pos = cbs->superblockoffsets[middle_s_block];
    containing_s_block = middle_s_block+1;
  }

  /* search within superblock */
  for (block_idx = containing_s_block * cbs->superblocksize;
       block_idx < cbs->num_of_blocks;
       block_idx++) {
    class = (unsigned int) gt_compact_ulong_store_get(cbs->classes, block_idx);
    if (num <= rank_sum + (cbs->blocksize - class))
      break;
    blocks_offset_pos += gt_popcount_tab_offset_bits(cbs->popcount_tab, class);
    rank_sum += cbs->blocksize - class;
  }
  position = block_idx * cbs->blocksize;
  gt_assert(class != cbs->blocksize + 1);

  /* search within block */
  if (block_idx != cbs->num_of_blocks - 1)
    select_mask = 1UL << (cbs->blocksize - 1UL);
  else
    select_mask = 1UL << (cbs->last_block_len - 1UL);

  block =
    gt_popcount_tab_get(cbs->popcount_tab, class,
                        (unsigned long)
                          gt_compressed_bitsequence_get_variable_field(
                                  cbs->c_offsets, blocks_offset_pos,
                                  gt_popcount_tab_offset_bits(cbs->popcount_tab,
                                                              class)));
  block = ~block; /* invert, we are searching for 0 not for 1 */
  while (rank_sum < num && select_mask != 0) {
    rank_sum += (block & select_mask) != 0 ? 1 : 0;
    select_mask >>= 1;
    position++;
  }
  gt_assert(rank_sum == num);

  return position - 1;
}

void gt_compressed_bitsequence_delete(GtCompressedBitsequence *cbs)
{
  if (cbs != NULL) {
    gt_popcount_tab_delete(cbs->popcount_tab);
    gt_compact_ulong_store_delete(cbs->classes);
    gt_free(cbs->c_offsets);
    gt_free(cbs->superblockranks);
    gt_free(cbs->superblockoffsets);
    gt_free(cbs->cbs_bi);
    gt_free(cbs);
  }
}

static int gt_compressed_bitsequence_unit_test_variable_field(
                                            GtError *err,
                                            GtBitsequence *bitseq)
{
  int had_err = 0;
  GtBitsequence result,
                value = ((GtBitsequence) 0xFF) << 28;

  gt_error_check(err);

  result = gt_compressed_bitsequence_get_variable_field(bitseq, 0, 16U);
  gt_ensure(had_err, result == (GtBitsequence) 0xAAAA);
  result = gt_compressed_bitsequence_get_variable_field(bitseq, 16UL, 16U);
  gt_ensure(had_err, result == (GtBitsequence) 0xCCCC);
  result = gt_compressed_bitsequence_get_variable_field(bitseq, 32UL, 64U);
  gt_ensure(had_err, result == (GtBitsequence) 0xAAAACCCCAAAACCCC);

  gt_compressed_bitsequence_set_variable_field(bitseq, 32UL, 64U, value);
  result = gt_compressed_bitsequence_get_variable_field(bitseq, 0, 64U);
  gt_ensure(had_err, result == (GtBitsequence) 0xAAAACCCC0000000F);
  result = gt_compressed_bitsequence_get_variable_field(bitseq, 64UL, 64U);
  gt_ensure(had_err, result == (GtBitsequence) 0xF0000000AAAACCCC);
  gt_compressed_bitsequence_set_variable_field(bitseq, 128UL, 8U, value>>28);
  result = gt_compressed_bitsequence_get_variable_field(bitseq, 128UL, 64U);
  gt_ensure(had_err, result == (GtBitsequence) 0xFFAACCCCAAAACCCC);

  return had_err;
}

static int gt_compressed_bitsequence_unit_test_block_identical(
                                            GtError *err,
                                            GtBitsequence *bitseq,
                                            const unsigned int sample_testratio,
                                            const unsigned long cbs_testsize)
{
  int had_err = 0;
  unsigned int class, len, block_len;
  unsigned long block_offset = 0,
                offsets_bitpos, sample, block,
                idx, jdx;
  GtBitsequence orig_block;
  GtCompressedBitsequence *cbs;

  gt_error_check(err);

  cbs = gt_compressed_bitsequence_new(bitseq, sample_testratio, cbs_testsize);
  gt_ensure(had_err, cbs != NULL);
  if (cbs != NULL) {
    block_len = cbs->blocksize;
  }

  for (idx = 0;
       !had_err && cbs!=NULL && idx < cbs->num_of_blocks;
       idx++) {
    if (idx == cbs->num_of_blocks -1) {
      block_len = cbs->last_block_len;
    }
    sample = idx / cbs->superblocksize;
    if (sample==0)
      offsets_bitpos = 0;
    else {
      offsets_bitpos = cbs->superblockoffsets[sample - 1];
    }
    for (jdx = sample * cbs->superblocksize; jdx < idx; jdx++) {
      class = (unsigned int) gt_compact_ulong_store_get(cbs->classes, jdx);
      offsets_bitpos += gt_popcount_tab_offset_bits(cbs->popcount_tab, class);
    }
    class = (unsigned int) gt_compact_ulong_store_get(cbs->classes, idx);
    len = gt_popcount_tab_offset_bits(cbs->popcount_tab, class);
    block_offset = (unsigned long)
      gt_compressed_bitsequence_get_variable_field(cbs->c_offsets,
                                                   offsets_bitpos, len);
    block = gt_popcount_tab_get(cbs->popcount_tab, class, block_offset);
    orig_block = gt_compressed_bitsequence_get_variable_field(
                                                         bitseq,
                                                         idx * cbs->blocksize,
                                                         block_len);
    gt_ensure(had_err, block == (unsigned long) orig_block);
  }
  gt_compressed_bitsequence_delete(cbs);
  return had_err;
}

int gt_compressed_bitsequence_unit_test(GtError *err)
{
  const unsigned int sample_testratio = 32U;
  int had_err = 0;
  const unsigned long bitseq_testsize = 256UL,
                      cbs_testsize = 16380UL;
  unsigned long idx;
  GtBitsequence *bitseq;
  GtCompressedBitsequence *cbs;

  gt_error_check(err);

  bitseq = gt_malloc((size_t) bitseq_testsize * sizeof (*bitseq));

  for (idx = 0; idx < bitseq_testsize; idx++) {
    bitseq[idx] = (GtBitsequence) 0xAAAACCCCAAAACCCC;
  }

  had_err = gt_compressed_bitsequence_unit_test_variable_field(err, bitseq);

  if (!had_err) {
    for (idx = 0; idx < bitseq_testsize; idx++) {
      bitseq[idx] = (GtBitsequence) 0xAAAACCCCAAAACCCCULL;
    }
    had_err =
      gt_compressed_bitsequence_unit_test_block_identical(err, bitseq,
                                                          sample_testratio,
                                                          cbs_testsize);
  }

  if (!had_err) {
    cbs = gt_compressed_bitsequence_new(bitseq, sample_testratio, cbs_testsize);
    gt_ensure(had_err, cbs != NULL);

    if (cbs != NULL)
    {

      for (idx = 1UL; idx < cbs_testsize; idx += GT_DIV2(GT_INTWORDSIZE))
      {
        gt_ensure(had_err, 1 == gt_compressed_bitsequence_access(cbs, idx - 1));
        gt_ensure(had_err, 0 == gt_compressed_bitsequence_access(cbs, idx));
      }
      gt_compressed_bitsequence_delete(cbs);
    }
  }

  /* set some blocks to 0  or ~0 to see if correct block ist found */
  bitseq[GT_DIV2(bitseq_testsize)] = (GtBitsequence) 0;
  bitseq[GT_DIV2(bitseq_testsize)+1] = (GtBitsequence) ~0ULL;
  if (!had_err) {
    cbs = gt_compressed_bitsequence_new(bitseq, sample_testratio, cbs_testsize);
    gt_ensure(had_err, cbs != NULL);

    if (cbs != NULL)
    {
      unsigned long rank1, rank0;

      rank1 = gt_compressed_bitsequence_rank_1(cbs,
                               (unsigned long) GT_INTWORDSIZE - 1UL);
      gt_ensure(had_err, rank1 == (unsigned long) GT_DIV2(GT_INTWORDSIZE));
      for (idx = (unsigned long) GT_INTWORDSIZE - 1UL;
           !had_err && idx < cbs_testsize;
           idx += idx + 1)
      {
        rank1 = gt_compressed_bitsequence_rank_1(cbs, idx);
        rank0 = gt_compressed_bitsequence_rank_0(cbs, idx);
        gt_ensure(had_err, rank1 + rank0 == idx+1);
      }
      rank1 = gt_compressed_bitsequence_rank_1(cbs, cbs_testsize-1);
      rank0 = gt_compressed_bitsequence_rank_0(cbs, cbs_testsize-1);
      gt_ensure(had_err, rank1 + rank0 == cbs_testsize);
      gt_compressed_bitsequence_delete(cbs);
    }
  }

  if (!had_err) {
    cbs = gt_compressed_bitsequence_new(bitseq, sample_testratio, cbs_testsize);
    if (cbs != NULL) {
      unsigned long ranktotal1, select1, rank1, select0, rank0;
      ranktotal1 = gt_compressed_bitsequence_rank_1(cbs, cbs_testsize - 1);
      for (idx = 1UL; !had_err && idx <= ranktotal1; idx++) {
        select1 = gt_compressed_bitsequence_select_1(cbs, idx);
        gt_ensure(had_err, gt_compressed_bitsequence_access(cbs, select1) == 1);
        rank1 = gt_compressed_bitsequence_rank_1(cbs, select1);
        gt_ensure(had_err, idx == rank1);
      }
      for (idx = 1UL; !had_err && idx <= cbs_testsize - ranktotal1; idx++) {
        select0 = gt_compressed_bitsequence_select_0(cbs, idx);
        gt_ensure(had_err, gt_compressed_bitsequence_access(cbs, select0) == 0);
        rank0 = gt_compressed_bitsequence_rank_0(cbs, select0);
        gt_ensure(had_err, idx == rank0);
      }
    }
    gt_compressed_bitsequence_delete(cbs);
  }

  gt_free(bitseq);

  return had_err;
}
