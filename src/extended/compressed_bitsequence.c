/*
  Copyright (c) 2012-2013 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
  Copyright (c) 2012-2013 Center for Bioinformatics, University of Hamburg

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

#include "core/byte_popcount_api.h"
#include "core/byte_select_api.h"
#include "core/combinatorics.h"
#include "core/ensure.h"
#include "core/fa.h"
#include "core/fileutils_api.h"
#include "core/log_api.h"
#include "core/ma_api.h"
#include "core/mathsupport.h"
#include "core/safearith.h"
#include "core/unused_api.h"
#include "extended/compressed_bitsequence.h"
#include "extended/popcount_tab.h"

/* this seems to be a good default value. maybe change this in the future */
#define GT_COMP_BITSEQ_BLOCKSIZE 15U

/* gt_compressed_bitsequence_ps_overflow contains a bit mask x consisting of 8
   bytes x[7],...,x[0] and each is set to 128-i */
const uint64_t gt_compressed_bitsequence_ps_overflow[] = {
  (uint64_t) 0x8080808080808080ULL, (uint64_t) 0x7f7f7f7f7f7f7f7fULL,
  (uint64_t) 0x7e7e7e7e7e7e7e7eULL, (uint64_t) 0x7d7d7d7d7d7d7d7dULL,
  (uint64_t) 0x7c7c7c7c7c7c7c7cULL, (uint64_t) 0x7b7b7b7b7b7b7b7bULL,
  (uint64_t) 0x7a7a7a7a7a7a7a7aULL, (uint64_t) 0x7979797979797979ULL,
  (uint64_t) 0x7878787878787878ULL, (uint64_t) 0x7777777777777777ULL,
  (uint64_t) 0x7676767676767676ULL, (uint64_t) 0x7575757575757575ULL,
  (uint64_t) 0x7474747474747474ULL, (uint64_t) 0x7373737373737373ULL,
  (uint64_t) 0x7272727272727272ULL, (uint64_t) 0x7171717171717171ULL,
  (uint64_t) 0x7070707070707070ULL, (uint64_t) 0x6f6f6f6f6f6f6f6fULL,
  (uint64_t) 0x6e6e6e6e6e6e6e6eULL, (uint64_t) 0x6d6d6d6d6d6d6d6dULL,
  (uint64_t) 0x6c6c6c6c6c6c6c6cULL, (uint64_t) 0x6b6b6b6b6b6b6b6bULL,
  (uint64_t) 0x6a6a6a6a6a6a6a6aULL, (uint64_t) 0x6969696969696969ULL,
  (uint64_t) 0x6868686868686868ULL, (uint64_t) 0x6767676767676767ULL,
  (uint64_t) 0x6666666666666666ULL, (uint64_t) 0x6565656565656565ULL,
  (uint64_t) 0x6464646464646464ULL, (uint64_t) 0x6363636363636363ULL,
  (uint64_t) 0x6262626262626262ULL, (uint64_t) 0x6161616161616161ULL,
  (uint64_t) 0x6060606060606060ULL, (uint64_t) 0x5f5f5f5f5f5f5f5fULL,
  (uint64_t) 0x5e5e5e5e5e5e5e5eULL, (uint64_t) 0x5d5d5d5d5d5d5d5dULL,
  (uint64_t) 0x5c5c5c5c5c5c5c5cULL, (uint64_t) 0x5b5b5b5b5b5b5b5bULL,
  (uint64_t) 0x5a5a5a5a5a5a5a5aULL, (uint64_t) 0x5959595959595959ULL,
  (uint64_t) 0x5858585858585858ULL, (uint64_t) 0x5757575757575757ULL,
  (uint64_t) 0x5656565656565656ULL, (uint64_t) 0x5555555555555555ULL,
  (uint64_t) 0x5454545454545454ULL, (uint64_t) 0x5353535353535353ULL,
  (uint64_t) 0x5252525252525252ULL, (uint64_t) 0x5151515151515151ULL,
  (uint64_t) 0x5050505050505050ULL, (uint64_t) 0x4f4f4f4f4f4f4f4fULL,
  (uint64_t) 0x4e4e4e4e4e4e4e4eULL, (uint64_t) 0x4d4d4d4d4d4d4d4dULL,
  (uint64_t) 0x4c4c4c4c4c4c4c4cULL, (uint64_t) 0x4b4b4b4b4b4b4b4bULL,
  (uint64_t) 0x4a4a4a4a4a4a4a4aULL, (uint64_t) 0x4949494949494949ULL,
  (uint64_t) 0x4848484848484848ULL, (uint64_t) 0x4747474747474747ULL,
  (uint64_t) 0x4646464646464646ULL, (uint64_t) 0x4545454545454545ULL,
  (uint64_t) 0x4444444444444444ULL, (uint64_t) 0x4343434343434343ULL,
  (uint64_t) 0x4242424242424242ULL, (uint64_t) 0x4141414141414141ULL,
  (uint64_t) 0x4040404040404040ULL
};

typedef struct
{
  GtUword      block_offset,
               idx,
               rank_sum;
  unsigned int class,
               block_len;
} GtCompressedBitsequenceBlockInfo;

typedef struct
{
  GtUword      *c_offsets_size,
               *classes_size,
               *num_of_bits,
               *num_of_blocks,
               *num_of_superblocks,
               *superblockoffsets_size,
               *superblockranks_size;
  unsigned int *blocksize,
               *class_bits,
               *last_block_len,
               *superblockoffsets_bits,
               *superblockranks_bits,
               *superblocksize;
}GtCompressedBitsequenceHeaderPtr;

struct GtCompressedBitsequence
{
  GtCompressedBitsequenceHeaderPtr  header;
  GtPopcountTab                    *popcount_tab;
  GtBitsequence                    *c_offsets,
                                   *classes,
                                   *superblockoffsets,
                                   *superblockranks;
  GtCompressedBitsequenceBlockInfo *cbs_bi;
  void                             *mmapped;
  GtUword                           c_offsets_size,
                                    classes_size,
                                    num_of_bits,
                                    num_of_blocks,
                                    num_of_superblocks,
                                    superblockoffsets_size,
                                    superblockranks_size;
  unsigned int                      blocksize,
                                    class_bits,
                                    last_block_len,
                                    superblockoffsets_bits,
                                    superblockranks_bits,
                                    superblocksize;
  bool                              from_file;
};

static void gt_compressed_bitsequence_header_setup_mapspec(GtMapspec *mapspec,
                                                           void *data,
                                                           bool write)
{
  GtCompressedBitsequence *cbs = data;

  if (write) {
    cbs->header.c_offsets_size         = &(cbs->c_offsets_size);
    cbs->header.classes_size           = &(cbs->classes_size);
    cbs->header.num_of_bits            = &(cbs->num_of_bits);
    cbs->header.num_of_blocks          = &(cbs->num_of_blocks);
    cbs->header.num_of_superblocks     = &(cbs->num_of_superblocks);
    cbs->header.superblockoffsets_size = &(cbs->superblockoffsets_size);
    cbs->header.superblockranks_size   = &(cbs->superblockranks_size);

    cbs->header.blocksize              = &(cbs->blocksize);
    cbs->header.class_bits             = &(cbs->class_bits);
    cbs->header.last_block_len         = &(cbs->last_block_len);
    cbs->header.superblockoffsets_bits = &(cbs->superblockoffsets_bits);
    cbs->header.superblockranks_bits   = &(cbs->superblockranks_bits);
    cbs->header.superblocksize         = &(cbs->superblocksize);
  }
  gt_mapspec_add_ulong(mapspec, cbs->header.c_offsets_size, 1UL);
  gt_mapspec_add_ulong(mapspec, cbs->header.classes_size, 1UL);
  gt_mapspec_add_ulong(mapspec, cbs->header.num_of_bits, 1UL);
  gt_mapspec_add_ulong(mapspec, cbs->header.num_of_blocks, 1UL);
  gt_mapspec_add_ulong(mapspec, cbs->header.num_of_superblocks, 1UL);
  gt_mapspec_add_ulong(mapspec, cbs->header.superblockoffsets_size, 1UL);
  gt_mapspec_add_ulong(mapspec, cbs->header.superblockranks_size, 1UL);

  gt_mapspec_add_uint(mapspec, cbs->header.blocksize, 1UL);
  gt_mapspec_add_uint(mapspec, cbs->header.class_bits, 1UL);
  gt_mapspec_add_uint(mapspec, cbs->header.last_block_len, 1UL);
  gt_mapspec_add_uint(mapspec, cbs->header.superblockoffsets_bits, 1UL);
  gt_mapspec_add_uint(mapspec, cbs->header.superblockranks_bits, 1UL);
  gt_mapspec_add_uint(mapspec, cbs->header.superblocksize, 1UL);
}

static void gt_compressed_bitsequence_data_setup_mapspec(GtMapspec *mapspec,
                                                         void *data,
                                                         GT_UNUSED bool write)
{
  GtCompressedBitsequence *cbs = data;

  /* this makes writing easy, and solves the problem of knowing where the data
     starts within the file/mapped area */
  gt_compressed_bitsequence_header_setup_mapspec(mapspec, data, write);

  gt_mapspec_add_bitsequence(mapspec, cbs->c_offsets, cbs->c_offsets_size);
  gt_mapspec_add_bitsequence(mapspec, cbs->classes, cbs->classes_size);
  gt_mapspec_add_bitsequence(mapspec, cbs->superblockoffsets,
                             cbs->superblockoffsets_size);
  gt_mapspec_add_bitsequence(mapspec, cbs->superblockranks,
                             cbs->superblockranks_size);
}

/* Assumes most significant bits are filled first in the vector, returns bits in
   <len> least significant positions */
static inline GtBitsequence gt_compressed_bitsequence_get_variable_field(
                                                          GtBitsequence *bitseq,
                                                          GtUword start,
                                                          unsigned int len)
{
  GtBitsequence result = 0;
  gt_assert(len <= (unsigned int) GT_INTWORDSIZE);
  if (len != 0) {
    const unsigned int bit_offset = (unsigned int) GT_MODWORDSIZE(start);
    const GtUword word_offset = GT_DIVWORDSIZE(start);

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
static inline void
gt_compressed_bitsequence_set_variable_field(GtBitsequence *bitseq,
                                             GT_UNUSED GtUword size,
                                             GtUword start,
                                             unsigned int len,
                                             GtBitsequence value)
{
  const unsigned int bit_offset = (unsigned int) GT_MODWORDSIZE(start);
  const GtUword word_offset = GT_DIVWORDSIZE(start);
  const unsigned int bits_left = (unsigned int) (GT_INTWORDSIZE - bit_offset);
  GtBitsequence mask = ~((GtBitsequence) 0) << (GtBitsequence) bits_left;
  if (bits_left == (unsigned int) GT_INTWORDSIZE)
    mask = 0;

  if (len > bits_left) {
    unsigned int overhang = len - bits_left;
    gt_assert(word_offset + 1 < size);

    bitseq[word_offset] &= mask;
    bitseq[word_offset] |= (value >> (GtBitsequence) overhang);

    mask = ~((GtBitsequence) 0) >> (GtBitsequence) overhang;
    bitseq[word_offset + 1] &= mask;
    bitseq[word_offset + 1] |=
      value << (GtBitsequence) (GT_INTWORDSIZE - overhang);
  }
  else {
    gt_assert(word_offset < size);
    if (bit_offset + len != (unsigned int) GT_INTWORDSIZE)
      mask |= ~((GtBitsequence) 0) >> (GtBitsequence) (bit_offset + len);

    bitseq[word_offset] &= mask;
    bitseq[word_offset] |= value << (GtBitsequence) (bits_left - len);
  }
}

static inline unsigned int
gt_compressed_bitsequence_get_class(GtCompressedBitsequence *cbs,
                                    GtUword idx)
{
  gt_assert(cbs != NULL);
  gt_assert(cbs->classes != NULL);
  return (unsigned int)
    gt_compressed_bitsequence_get_variable_field(cbs->classes,
                                                 idx * cbs->class_bits,
                                                 cbs->class_bits);
}

static inline void
gt_compressed_bitsequence_set_class(GtCompressedBitsequence *cbs,
                                    unsigned int class,
                                    GtUword idx)
{
  gt_compressed_bitsequence_set_variable_field(cbs->classes,
                                               cbs->classes_size,
                                               idx * cbs->class_bits,
                                               cbs->class_bits,
                                               (GtBitsequence) class);
}

static void gt_compressed_bitsequence_fill_tabs(GtCompressedBitsequence *cbs,
                                                GtBitsequence *bitseq)
{
  unsigned int class,
               samplecount = 0,
               c_offset_bits,
               block_len = cbs->blocksize;
  GtUword idx,
          samplenum = 0;
  GtBitsequence block, offsets_bitpos,
                bitpos = 0,
                rank_sum = 0;

  for (idx = 0; idx < cbs->num_of_blocks; idx++) {
    if (idx == cbs->num_of_blocks -1) {
      block_len = cbs->last_block_len;
    }
    class = gt_compressed_bitsequence_get_class(cbs, idx);
    c_offset_bits = gt_popcount_tab_offset_bits(cbs->popcount_tab, class);
    block = gt_compressed_bitsequence_get_variable_field(bitseq,
                                                         idx * cbs->blocksize,
                                                         block_len);
    offsets_bitpos = (GtBitsequence) gt_popcount_tab_get_offset_for_block(
                                                         cbs->popcount_tab,
                                                         (GtUword) block);
    gt_compressed_bitsequence_set_variable_field(cbs->c_offsets,
                                                 cbs->c_offsets_size,
                                                 (GtUword) bitpos,
                                                 c_offset_bits, offsets_bitpos);
    rank_sum += class;
    bitpos += c_offset_bits;
    samplecount++;
    if (cbs->superblocksize == samplecount) {
      samplecount = 0;
      gt_compressed_bitsequence_set_variable_field(
                                        cbs->superblockoffsets,
                                        cbs->superblockoffsets_size,
                                        samplenum * cbs->superblockoffsets_bits,
                                        cbs->superblockoffsets_bits,
                                        bitpos);
      gt_compressed_bitsequence_set_variable_field(
                                          cbs->superblockranks,
                                          cbs->superblockranks_size,
                                          samplenum * cbs->superblockranks_bits,
                                          cbs->superblockranks_bits,
                                          rank_sum);
      samplenum++;
    }
  }
  if (samplecount != 0) {
    gt_compressed_bitsequence_set_variable_field(
                                        cbs->superblockoffsets,
                                        cbs->superblockoffsets_size,
                                        samplenum * cbs->superblockoffsets_bits,
                                        cbs->superblockoffsets_bits,
                                        bitpos);
    gt_compressed_bitsequence_set_variable_field(
                                          cbs->superblockranks,
                                          cbs->superblockranks_size,
                                          samplenum * cbs->superblockranks_bits,
                                          cbs->superblockranks_bits,
                                          rank_sum);
  }
}

static void gt_compressed_bitsequence_init_s_tabs(GtCompressedBitsequence *cbs)
{
  cbs->num_of_superblocks = cbs->num_of_blocks / cbs->superblocksize;
  if (cbs->num_of_blocks % cbs->superblocksize != 0)
    cbs->num_of_superblocks++;

  cbs->superblockoffsets_size = (GtUword)
    GT_NUMOFINTSFORBITS(cbs->num_of_superblocks * cbs->superblockoffsets_bits);
  cbs->superblockranks_size = (GtUword)
    GT_NUMOFINTSFORBITS(cbs->num_of_superblocks * cbs->superblockranks_bits);
  GT_INITBITTAB(cbs->superblockoffsets,
                cbs->num_of_superblocks * cbs->superblockoffsets_bits);
  GT_INITBITTAB(cbs->superblockranks,
                cbs->num_of_superblocks * cbs->superblockranks_bits);
}

static void gt_compressed_bitsequence_init_c_tab(GtCompressedBitsequence *cbs)
{
  GtUword classes_bits;
  cbs->num_of_blocks = cbs->num_of_bits / cbs->blocksize;
  if (cbs->num_of_bits % cbs->blocksize != 0 || cbs->num_of_blocks == 0)
    cbs->num_of_blocks++;
#ifdef __SSE4_2__
  cbs->class_bits = (unsigned int)
    (sizeof (unsigned int) * CHAR_BIT - __builtin_clz(cbs->blocksize));
#else
  cbs->class_bits = gt_determinebitspervalue((GtUword) cbs->blocksize);
#endif
  classes_bits = cbs->num_of_blocks * cbs->class_bits;
  cbs->classes_size = (GtUword) GT_NUMOFINTSFORBITS(classes_bits);
  GT_INITBITTAB(cbs->classes, classes_bits);
}

static void
gt_compressed_bitsequence_fill_c_tab_init_o_tab(GtCompressedBitsequence *cbs,
                                                GtBitsequence *bitseq)
{
  unsigned int block_len = cbs->blocksize,
               current_class;
  GtUword idx, o_size = 0, ones = 0;
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
                                          (GtUword) current_blk);
    gt_compressed_bitsequence_set_class(cbs, current_class, idx);
    o_size += gt_popcount_tab_offset_bits(cbs->popcount_tab, current_class);
    ones += current_blk;
  }
  cbs->c_offsets_size = (GtUword) GT_NUMOFINTSFORBITS(o_size);
#ifdef __SSE4_2__
  cbs->superblockoffsets_bits = __builtin_clzl(o_size);
  cbs->superblockranks_bits = __builtin_clzl(ones);
#else
  cbs->superblockoffsets_bits = gt_determinebitspervalue(o_size);
  cbs->superblockranks_bits = gt_determinebitspervalue(ones);
#endif
  GT_INITBITTAB(cbs->c_offsets, o_size);
}

static GtCompressedBitsequence* gt_compressed_bitsequence_new_empty(void)
{
  GtCompressedBitsequence *cbs;

  cbs = gt_calloc((size_t) 1, sizeof (*cbs));
  cbs->blocksize = GT_COMP_BITSEQ_BLOCKSIZE;
  cbs->mmapped = NULL;
  return cbs;
}

static inline void gt_compressed_bitsequence_init(GtCompressedBitsequence *cbs,
                                                  unsigned int samplerate,
                                                  GtUword num_of_bits)
{
  cbs->superblocksize = samplerate;
  cbs->popcount_tab = gt_popcount_tab_new(cbs->blocksize);
  cbs->c_offsets_size = 0;
  cbs->num_of_bits = num_of_bits;
  cbs->last_block_len = (unsigned int) (cbs->num_of_bits % cbs->blocksize);
  if (cbs->last_block_len == 0)
    cbs->last_block_len = cbs->blocksize;
  cbs->mmapped = NULL;
}

GtCompressedBitsequence *
gt_compressed_bitsequence_new(GtBitsequence *bitseq,
                              unsigned int samplerate,
                              GtUword num_of_bits)
{
  GtCompressedBitsequence *cbs;

  gt_assert(samplerate != 0);

  cbs = gt_compressed_bitsequence_new_empty();
  gt_compressed_bitsequence_init(cbs, samplerate, num_of_bits);

  gt_compressed_bitsequence_init_c_tab(cbs);
  gt_compressed_bitsequence_fill_c_tab_init_o_tab(cbs, bitseq);
  gt_compressed_bitsequence_init_s_tabs(cbs);
  gt_compressed_bitsequence_fill_tabs(cbs, bitseq);
  cbs->from_file = false;
  gt_log_log("new cbs:\n"
             "blzise: %u\n"
             "offsize: " GT_WU "\n"
             "cbits: %u\n"
             "csize: " GT_WU "\n"
             "lastblen: %u\n",
             cbs->blocksize,
             cbs->c_offsets_size,
             cbs->class_bits,
             cbs->classes_size,
             cbs->last_block_len);
  return cbs;
}

static inline void
gt_compressed_bitsequence_calc_block_info(GtCompressedBitsequence *cbs,
                                          GtUword position)
{
  GtUword idx;
  GtCompressedBitsequenceBlockInfo *bi = cbs->cbs_bi;

  idx = position / cbs->blocksize;

  if (cbs->cbs_bi == NULL) {
    cbs->cbs_bi = gt_malloc(sizeof (*cbs->cbs_bi));
    bi = cbs->cbs_bi;
    bi->idx = idx + 1;
  }

  if (idx != bi->idx) {
    unsigned int offset_bits;
    GtUword jdx, sample, offsets_bitpos;

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
      offsets_bitpos = (GtUword)
        gt_compressed_bitsequence_get_variable_field(
                                     cbs->superblockoffsets,
                                     (sample - 1) * cbs->superblockoffsets_bits,
                                     cbs->superblockoffsets_bits);
      bi->rank_sum = (GtUword)
        gt_compressed_bitsequence_get_variable_field(
                                       cbs->superblockranks,
                                       (sample - 1) * cbs->superblockranks_bits,
                                       cbs->superblockranks_bits);
    }
    for (jdx = sample * cbs->superblocksize; jdx < idx; jdx++) {
      bi->class = gt_compressed_bitsequence_get_class(cbs, jdx);
      bi->rank_sum += bi->class;
      offsets_bitpos += gt_popcount_tab_offset_bits(cbs->popcount_tab,
                                                    bi->class);
    }
    bi->class = gt_compressed_bitsequence_get_class(cbs, idx);
    offset_bits = gt_popcount_tab_offset_bits(cbs->popcount_tab, bi->class);
    bi->block_offset = (GtUword)
      gt_compressed_bitsequence_get_variable_field(cbs->c_offsets,
                                                   offsets_bitpos,
                                                   offset_bits);
  }
}

int gt_compressed_bitsequence_access(GtCompressedBitsequence *cbs,
                                     GtUword position)
{
  int bit;
  unsigned int pos_in_block;
  GtUword block;
  GtCompressedBitsequenceBlockInfo *cbs_bi;

  gt_assert(cbs != NULL);
  gt_assert(position < cbs->num_of_bits);

  pos_in_block = (unsigned int) position % cbs->blocksize;

  gt_compressed_bitsequence_calc_block_info(cbs, position);
  cbs_bi = cbs->cbs_bi;
  if (cbs_bi->class == 0 || cbs_bi->class == cbs->blocksize)
    return cbs_bi->class == 0 ? 0 : 1;

  block = gt_popcount_tab_get(cbs->popcount_tab,
                              cbs_bi->class,
                              cbs_bi->block_offset);

  /* account for shorter blocks at the end */
  block <<= (cbs->blocksize - cbs_bi->block_len);

  bit = (int) ((block >> (cbs->blocksize - pos_in_block - 1)) & 1UL);
  return bit;
}

GtUword gt_compressed_bitsequence_rank_1(GtCompressedBitsequence *cbs,
                                         GtUword position)
{
  unsigned int pos_in_block;
  GtCompressedBitsequenceBlockInfo *cbs_bi;

  gt_assert(cbs != NULL);
  gt_assert(position < cbs->num_of_bits);

  pos_in_block = (unsigned int) position % cbs->blocksize;

  gt_compressed_bitsequence_calc_block_info(cbs, position);
  cbs_bi = cbs->cbs_bi;

  pos_in_block += cbs->blocksize - cbs_bi->block_len;
  if (cbs_bi->class == 0)
    return cbs_bi->rank_sum;
  if (cbs_bi->class == cbs->blocksize)
    return cbs_bi->rank_sum + pos_in_block + 1;

  return cbs_bi->rank_sum + gt_popcount_tab_rank_1(cbs->popcount_tab,
                                                   cbs_bi->class,
                                                   cbs_bi->block_offset,
                                                   pos_in_block);
}

GtUword gt_compressed_bitsequence_rank_0(GtCompressedBitsequence *cbs,
                                         GtUword position)
{
  unsigned int pos_in_block;
  GtUword rank_sum0;
  GtCompressedBitsequenceBlockInfo *cbs_bi;

  gt_assert(cbs != NULL);
  gt_assert(position < cbs->num_of_bits);

  pos_in_block = (unsigned int) position % cbs->blocksize;

  gt_compressed_bitsequence_calc_block_info(cbs, position);
  cbs_bi = cbs->cbs_bi;

  pos_in_block += cbs->blocksize - cbs_bi->block_len;
  rank_sum0 = ((cbs_bi->idx) * cbs->blocksize) - cbs_bi->rank_sum;
  rank_sum0 -= cbs->blocksize - cbs_bi->block_len;
  if (cbs_bi->class == 0)
    return rank_sum0 + pos_in_block + 1;
  if (cbs_bi->class == cbs->blocksize)
    return rank_sum0;
  return rank_sum0 + gt_popcount_tab_rank_0(cbs->popcount_tab,
                                            cbs_bi->class,
                                            cbs_bi->block_offset,
                                            pos_in_block);
}

static inline unsigned int
gt_compressed_bitsequence_select_1_word(uint64_t word, unsigned int i)
{
#ifdef __SSE4_2__
  uint64_t s = word, b;
  unsigned int byte_nr;
  s = s - ((s >> 1) & (uint64_t) 0x5555555555555555ULL);
  s = (s & (uint64_t) 0x3333333333333333ULL) +
      ((s >> 2) & (uint64_t) 0x3333333333333333ULL);
  s = (s + (s >> 4)) & (uint64_t) 0x0F0F0F0F0F0F0F0FULL;
  /* s *= (uint64_t) 0x0101010101010101ULL; [>this will add the running sums to
                                            the most significant byte<] */
  /* analog to multiplication which would use << */
  s = s         + (s >> 8)  + (s >> 16) + (s >> 24) +
      (s >> 32) + (s >> 40) + (s >> 48) + (s >> 52);
  /* now s contains 8 bytes s[0],...,s[7], s[i] contains the cumulative sum
     of (i+1)*8 least significant bits of s */
  b = (s + gt_compressed_bitsequence_ps_overflow[i]) &
      (uint64_t) 0x8080808080808080ULL;
  /* ps_overflow contains a bit mask mask consisting of 8 bytes
     mask[7],...,mask[0] and each set to 128-i
     => a byte b[i] in b is >= 128 if cum sum >= i */
  b >>= 7;

  /* __builtin_clzll returns the number of leading zeros, if b!=0 */
  byte_nr = __builtin_clzll(b) >> 3;   /* byte nr in [0..7] */
  /* subtract the cumulative sum of bits of all bytes before byte_nr */
  s >>= 8; /* remove total sum */
  i -= (s >> ((7 - byte_nr) << 3)) & 0xFFULL;
  return (byte_nr << 3) + (unsigned int)
    gt_byte_select[((i-1) << 8) +
                   ((word >> ((7 - byte_nr) << 3)) & 0xFFULL)];
#else
  unsigned int bytecount,
               idx,
               ranksum = 0,
               shift = (unsigned int) ((CHAR_BIT - 1) * sizeof (word));
  for (idx = 0; idx < (unsigned int) sizeof (word); ++idx, shift -= 8) {
    bytecount =
      (unsigned int) gt_byte_popcount[(word >> shift) & 0xFFULL];
    if (ranksum + bytecount >= i) {
      i -= ranksum;
      return (unsigned int) (idx * CHAR_BIT +
             gt_byte_select[((i - 1) << 8) + ((word >> shift) & 0xFFULL)]);
    }
    else
      ranksum += bytecount;
  }
  /* 64 if bit is not present */
  return (unsigned int) (CHAR_BIT * sizeof (word));
#endif
}

GtUword gt_compressed_bitsequence_select_1(GtCompressedBitsequence *cbs,
                                           GtUword num)
{
  unsigned int block_offset_bits,
               class = cbs->blocksize + 1;
  GtUword block_idx,
          blocks_offset_pos,
          containing_s_block,
          end_s_block, middle_s_block, start_s_block,
          position,
          rank_sum = 0,
          start, start_l, start_r;
  uint64_t block;

  gt_assert(num != 0);
  gt_assert(cbs != NULL);
  gt_assert(num < cbs->num_of_bits);

  /* if larger then max rank1 */
  if (num > (GtUword) gt_compressed_bitsequence_get_variable_field(
                      cbs->superblockranks,
                      (cbs->num_of_superblocks - 1) * cbs->superblockranks_bits,
                      cbs->superblockranks_bits))
    return cbs->num_of_bits;
  if (num <= (GtUword)
      gt_compressed_bitsequence_get_variable_field(cbs->superblockranks, 0,
                                                   cbs->superblockranks_bits)) {
    rank_sum = 0;
    containing_s_block = 0;
    blocks_offset_pos = 0;
  }
  else {
    /* binary search for superblock */
    start_s_block = num / (cbs->blocksize * cbs->superblocksize);
    end_s_block = cbs->num_of_superblocks;
    middle_s_block = GT_DIV2(start_s_block + end_s_block);
    while (start_s_block < end_s_block) {
      start = middle_s_block * cbs->superblockranks_bits;
      if ((GtUword) gt_compressed_bitsequence_get_variable_field(
                                                    cbs->superblockranks, start,
                                                    cbs->superblockranks_bits)
          < num) {
        start_r = start + cbs->superblockranks_bits;
        if ((GtUword) gt_compressed_bitsequence_get_variable_field(
                                                  cbs->superblockranks, start_r,
                                                  cbs->superblockranks_bits)
            < num)
          start_s_block = middle_s_block;
        else
          break;
      }
      else {
        start_l = start - cbs->superblockranks_bits;
        if ((GtUword) gt_compressed_bitsequence_get_variable_field(
                                                  cbs->superblockranks, start_l,
                                                  cbs->superblockranks_bits)
            >= num)
          end_s_block = middle_s_block;
        else {
          middle_s_block--;
          break;
        }
      }
      middle_s_block = GT_DIV2(start_s_block + end_s_block);
    }
    blocks_offset_pos = (GtUword)
      gt_compressed_bitsequence_get_variable_field(
                                   cbs->superblockoffsets,
                                   middle_s_block * cbs->superblockoffsets_bits,
                                   cbs->superblockoffsets_bits);
    rank_sum = (GtUword) gt_compressed_bitsequence_get_variable_field(
                                     cbs->superblockranks,
                                     middle_s_block * cbs->superblockranks_bits,
                                     cbs->superblockranks_bits);
    containing_s_block = middle_s_block+1;
  }

  /* search within superblock */
  for (block_idx = containing_s_block * cbs->superblocksize;
       block_idx < cbs->num_of_blocks;
       block_idx++) {
    class = gt_compressed_bitsequence_get_class(cbs, block_idx);
    if (num <= rank_sum + class)
      break;
    blocks_offset_pos += gt_popcount_tab_offset_bits(cbs->popcount_tab, class);
    rank_sum += class;
  }
  gt_assert(class != cbs->blocksize + 1);
  position = block_idx * cbs->blocksize;
  if (class == cbs->blocksize) {
    position += num - rank_sum - 1;
  }
  else {
    block_offset_bits = gt_popcount_tab_offset_bits(cbs->popcount_tab, class);

    /* search within block */
    block = (uint64_t)
      gt_popcount_tab_get(cbs->popcount_tab, class, (GtUword)
                          gt_compressed_bitsequence_get_variable_field(
                                              cbs->c_offsets, blocks_offset_pos,
                                              block_offset_bits));
    if (block_idx != cbs->num_of_blocks - 1)
      block <<= ((sizeof (block) * CHAR_BIT) - cbs->blocksize);
    else
      block <<= ((sizeof (block) * CHAR_BIT) - cbs->last_block_len);

    position +=
      gt_compressed_bitsequence_select_1_word(block,
                                              (unsigned int) (num - rank_sum));
  }

  return position;
}

GtUword gt_compressed_bitsequence_select_0(GtCompressedBitsequence *cbs,
                                           GtUword num)
{
  unsigned int block_offset_bits,
               class = cbs->blocksize + 1;
  GtUword block_idx,
          blocks_offset_pos,
          containing_s_block,
          end_s_block, middle_s_block, start_s_block,
          first_0_superblock_rank,
          max_0_rank,
          position,
          rank_sum = 0,
          s_block_bits,
          start, start_l, start_r;
  uint64_t block;

  gt_assert(num != 0);
  gt_assert(cbs != NULL);
  gt_assert(num < cbs->num_of_bits);

  s_block_bits = (GtUword) cbs->blocksize * cbs->superblocksize;
  max_0_rank = cbs->num_of_bits - gt_compressed_bitsequence_get_variable_field(
                      cbs->superblockranks,
                      (cbs->num_of_superblocks - 1) * cbs->superblockranks_bits,
                      cbs->superblockranks_bits);
  if (num > max_0_rank)
    return cbs->num_of_bits;

  first_0_superblock_rank =
    s_block_bits - gt_compressed_bitsequence_get_variable_field(
                                                     cbs->superblockranks, 0,
                                                     cbs->superblockranks_bits);
  if (num <= first_0_superblock_rank) {
    rank_sum = 0;
    containing_s_block = 0;
    blocks_offset_pos = 0;
  }
  else {
    GtUword middle_s_block_max,
            middle_s_block_r_max,
            middle_s_block_l_max;
    /* search for superblock */
    start_s_block = num / s_block_bits; /* cannot be smaller than this */
    end_s_block = cbs->num_of_superblocks;
    middle_s_block = GT_DIV2(start_s_block + end_s_block);
    while (start_s_block < end_s_block) {
      middle_s_block_max = s_block_bits * (middle_s_block+1);
      start = middle_s_block * cbs->superblockranks_bits;
      if ((middle_s_block_max -
            gt_compressed_bitsequence_get_variable_field(
                                                  cbs->superblockranks, start,
                                                  cbs->superblockranks_bits))
          < num) {
        middle_s_block_r_max = middle_s_block_max + s_block_bits;
        start_r = start + cbs->superblockranks_bits;
        if ((middle_s_block_r_max -
              gt_compressed_bitsequence_get_variable_field(
                                                  cbs->superblockranks, start_r,
                                                  cbs->superblockranks_bits))
            < num)
          start_s_block = middle_s_block;
        else
          break;
      }
      else {
        middle_s_block_l_max = middle_s_block_max - s_block_bits;
        start_l = start - cbs->superblockranks_bits;
        if ((middle_s_block_l_max -
               gt_compressed_bitsequence_get_variable_field(
                                                  cbs->superblockranks, start_l,
                                                  cbs->superblockranks_bits))
            >= num)
          end_s_block = middle_s_block;
        else {
          middle_s_block--;
          break;
        }
      }
      middle_s_block = GT_DIV2(start_s_block + end_s_block);
    }
    middle_s_block_max = s_block_bits * (middle_s_block + 1);
    blocks_offset_pos = (GtUword)
      gt_compressed_bitsequence_get_variable_field(
                                   cbs->superblockoffsets,
                                   middle_s_block * cbs->superblockoffsets_bits,
                                   cbs->superblockoffsets_bits);
    rank_sum = (GtUword)
      middle_s_block_max - gt_compressed_bitsequence_get_variable_field(
                                     cbs->superblockranks,
                                     middle_s_block * cbs->superblockranks_bits,
                                     cbs->superblockranks_bits);
    containing_s_block = middle_s_block+1;
  }

  /* search within superblock */
  for (block_idx = containing_s_block * cbs->superblocksize;
       block_idx < cbs->num_of_blocks;
       block_idx++) {
    class = gt_compressed_bitsequence_get_class(cbs, block_idx);
    if (num <= rank_sum + (cbs->blocksize - class))
      break;
    blocks_offset_pos += gt_popcount_tab_offset_bits(cbs->popcount_tab, class);
    rank_sum += cbs->blocksize - class;
  }
  position = block_idx * cbs->blocksize;
  gt_assert(class != cbs->blocksize + 1);
  if (class == 0)
    position += num - rank_sum - 1;
  else {
    block_offset_bits = gt_popcount_tab_offset_bits(cbs->popcount_tab, class);

    /* search within block */
    block = (uint64_t)
      gt_popcount_tab_get(cbs->popcount_tab, class, (GtUword)
                          gt_compressed_bitsequence_get_variable_field(
                                              cbs->c_offsets, blocks_offset_pos,
                                              block_offset_bits));
    if (block_idx != cbs->num_of_blocks - 1)
      block <<= ((sizeof (block) * CHAR_BIT) - cbs->blocksize);
    else
      block <<= ((sizeof (block) * CHAR_BIT) - cbs->last_block_len);

    /* invert because we search for 0 */
    position +=
      gt_compressed_bitsequence_select_1_word(~block,
                                              (unsigned int) (num - rank_sum));
  }

  return position;
}

static size_t
gt_compressed_bitsequence_header_size(GtCompressedBitsequence *cbs)
{
  size_t size = 0;
  gt_assert(cbs != NULL);
  size += sizeof (cbs->c_offsets_size);
  size += sizeof (cbs->classes_size);
  size += sizeof (cbs->num_of_bits);
  size += sizeof (cbs->num_of_blocks);
  size += sizeof (cbs->num_of_superblocks);
  size += sizeof (cbs->superblockoffsets_size);
  size += sizeof (cbs->superblockranks_size);

  size += sizeof (cbs->blocksize);
  size += sizeof (cbs->class_bits);
  size += sizeof (cbs->last_block_len);
  size += sizeof (cbs->superblockoffsets_bits);
  size += sizeof (cbs->superblockranks_bits);
  size += sizeof (cbs->superblocksize);

  return size;
}

size_t gt_compressed_bitsequence_file_size(GtCompressedBitsequence *cbs)
{
  size_t size = gt_compressed_bitsequence_header_size(cbs);

  size += sizeof (cbs->c_offsets[0]) * cbs->c_offsets_size;
  size += sizeof (cbs->classes[0]) * cbs->classes_size;
  size += sizeof (cbs->superblockoffsets[0]) * cbs->superblockoffsets_size;
  size += sizeof (cbs->superblockranks[0]) * cbs->superblockranks_size;

  return size;
}

size_t gt_compressed_bitsequence_size(GtCompressedBitsequence *cbs)
{
  size_t size =  sizeof (cbs) +
    gt_popcount_tab_calculate_size(cbs->blocksize) +
    sizeof (cbs->cbs_bi) +
    sizeof (cbs->c_offsets[0]) * cbs->c_offsets_size +
    sizeof (cbs->classes[0]) * cbs->classes_size +
    sizeof (cbs->superblockoffsets[0]) * cbs->superblockoffsets_size +
    sizeof (cbs->superblockranks[0]) * cbs->superblockranks_size;

  return size;
}

int gt_compressed_bitsequence_write(GtCompressedBitsequence *cbs,
                                    char *filename,
                                    GtError *err)
{
  int had_err = 0;
  GtUword expectedsize = 0;
  FILE *fp = NULL;

  fp = gt_fa_fopen(filename, "w", err);
  if (fp == NULL)
    had_err = -1;

  if (!had_err) {
    expectedsize = (GtUword) gt_compressed_bitsequence_file_size(cbs);
    had_err = gt_mapspec_write(gt_compressed_bitsequence_data_setup_mapspec,
                               fp, cbs, expectedsize, err);
  }
  if (!had_err) {
    gt_log_log("blocksize: %u\n"
           "class_offsets_size: " GT_WU "\n"
           "class_bits: %u\n"
           "classes_size: " GT_WU "\n"
           "last_block_len: %u\n"
           "num_of_bits: " GT_WU "\n"
           "num_of_blocks: " GT_WU "\n"
           "num_of_superblocks: " GT_WU "\n"
           "superblocksize: %u\n"
           "size: " GT_WU "\n",
           cbs->blocksize,
           cbs->c_offsets_size,
           cbs->class_bits,
           cbs->classes_size,
           cbs->last_block_len,
           cbs->num_of_bits,
           cbs->num_of_blocks,
           cbs->num_of_superblocks,
           cbs->superblocksize,
           expectedsize);
  }
  gt_fa_fclose(fp);

  return had_err;
}

GtCompressedBitsequence *
gt_compressed_bitsequence_new_from_file(const char *filename,
                                        GtError *err)
{
  int had_err = 0;
  GtUword expectedsize = 0;
  GtCompressedBitsequence *cbs = gt_compressed_bitsequence_new_empty();
  gt_assert(filename != NULL);
  if (gt_file_exists(filename)) {
    expectedsize = (GtUword) gt_compressed_bitsequence_header_size(cbs);
    had_err = gt_mapspec_read_header(
                                 gt_compressed_bitsequence_header_setup_mapspec,
                                 cbs, filename, expectedsize,
                                 &(cbs->mmapped), err);
    cbs->c_offsets_size = cbs->header.c_offsets_size[0];
    cbs->classes_size = cbs->header.classes_size[0];
    cbs->num_of_bits = cbs->header.num_of_bits[0];
    cbs->num_of_blocks = cbs->header.num_of_blocks[0];
    cbs->num_of_superblocks = cbs->header.num_of_superblocks[0];
    cbs->superblockoffsets_size = cbs->header.superblockoffsets_size[0];
    cbs->superblockranks_size = cbs->header.superblockranks_size[0];

    cbs->blocksize = cbs->header.blocksize[0];
    cbs->class_bits = cbs->header.class_bits[0];
    cbs->last_block_len = cbs->header.last_block_len[0];
    cbs->superblockoffsets_bits = cbs->header.superblockoffsets_bits[0];
    cbs->superblockranks_bits = cbs->header.superblockranks_bits[0];
    cbs->superblocksize = cbs->header.superblocksize[0];

    gt_fa_xmunmap(cbs->mmapped);
    cbs->mmapped = NULL;
    if (!had_err) {
      expectedsize = (GtUword) gt_compressed_bitsequence_file_size(cbs);
      gt_log_log("new expected: " GT_WU "\n", expectedsize);
      had_err = gt_mapspec_read(gt_compressed_bitsequence_data_setup_mapspec,
                                cbs, filename, expectedsize,
                                &(cbs->mmapped), err);
    }
    if (!had_err) {
      gt_log_log("blocksize: %u\n"
                 "class_offsets_size: " GT_WU "\n"
                 "class_bits: %u\n"
                 "classes_size: " GT_WU "\n"
                 "last_block_len: %u\n"
                 "num_of_bits: " GT_WU "\n"
                 "num_of_blocks: " GT_WU "\n"
                 "num_of_superblocks: " GT_WU "\n"
                 "superblocksize: %u\n"
                 "size: " GT_WU "\n",
                 cbs->blocksize,
                 cbs->c_offsets_size,
                 cbs->class_bits,
                 cbs->classes_size,
                 cbs->last_block_len,
                 cbs->num_of_bits,
                 cbs->num_of_blocks,
                 cbs->num_of_superblocks,
                 cbs->superblocksize,
                 expectedsize);
    }
  }
  else {
    gt_error_set(err, "file %s does not exist!", filename);
    had_err = -1;
  }
  if (had_err) {
    gt_compressed_bitsequence_delete(cbs);
    return NULL;
  }
  cbs->popcount_tab = gt_popcount_tab_new(cbs->blocksize);
  cbs->from_file = true;
  return cbs;
}

void gt_compressed_bitsequence_delete(GtCompressedBitsequence *cbs)
{
  if (cbs != NULL) {
    gt_popcount_tab_delete(cbs->popcount_tab);
    if (cbs->from_file) {
      gt_fa_xmunmap(cbs->mmapped);
    }
    else {
      gt_free(cbs->classes);
      gt_free(cbs->c_offsets);
      gt_free(cbs->superblockranks);
      gt_free(cbs->superblockoffsets);
    }
    gt_free(cbs->cbs_bi);
    gt_free(cbs);
  }
}

static int gt_compressed_bitsequence_unit_test_variable_field(
                                                          GtError *err,
                                                          GtBitsequence *bitseq,
                                                          GtUword size)
{
  int had_err = 0;
  GtBitsequence result,
                value = ((GtBitsequence) 0xff) << 12;

  gt_error_check(err);

  result = gt_compressed_bitsequence_get_variable_field(bitseq, 0, 16U);
  gt_ensure(result == (GtBitsequence) 0xaaaa);
  result = gt_compressed_bitsequence_get_variable_field(bitseq, 16UL, 16U);
  gt_ensure(result == (GtBitsequence) 0xcccc);
  result = gt_compressed_bitsequence_get_variable_field(bitseq, 32UL, 32U);
  gt_ensure(result == (GtBitsequence) 0xaaaaccccUL);
  if (GT_LOGWORDSIZE == 6) {
    result = gt_compressed_bitsequence_get_variable_field(bitseq, 32UL, 64U);
    gt_ensure(result == (GtBitsequence) 0xaaaaccccaaaaccccULL);
  }

  gt_compressed_bitsequence_set_variable_field(bitseq, size, 16UL, 32U, value);
  result = gt_compressed_bitsequence_get_variable_field(bitseq, 0, 32U);
  gt_ensure(result == (GtBitsequence) 0xaaaa000fUL);
  result = gt_compressed_bitsequence_get_variable_field(bitseq, 32UL, 32U);
  gt_ensure(result == (GtBitsequence) 0xf000ccccUL);
  gt_compressed_bitsequence_set_variable_field(bitseq, size,
                                               16UL, 32U,
                                               (GtBitsequence) 0xccccaaaaUL);
  gt_compressed_bitsequence_set_variable_field(bitseq, size, 128UL, 8U,
                                               (GtBitsequence) 0xFF);
  result = gt_compressed_bitsequence_get_variable_field(bitseq, 128UL, 32U);
  gt_ensure(result == (GtBitsequence) 0xffaaccccUL);

  if (GT_LOGWORDSIZE == 6) {
    value <<= 16; /* total of << 28 = 0x0000000FF0000000 */
    gt_compressed_bitsequence_set_variable_field(bitseq, size,
                                                 32UL, 64U, value);
    result = gt_compressed_bitsequence_get_variable_field(bitseq, 0, 64U);
    gt_ensure(result == (GtBitsequence) 0xaaaacccc0000000fULL);
    result = gt_compressed_bitsequence_get_variable_field(bitseq, 64UL, 64U);
    gt_ensure(result == (GtBitsequence) 0xf0000000aaaaccccULL);
    gt_compressed_bitsequence_set_variable_field(
                bitseq, size, 32UL, 64U, (GtBitsequence) 0xaaaaccccaaaaccccULL);
    result = gt_compressed_bitsequence_get_variable_field(bitseq, 120UL, 64U);
    gt_ensure(result == (GtBitsequence) 0xccffaaccccaaaaccULL);
  }

  return had_err;
}

static int gt_compressed_bitsequence_unit_test_block_identical(
                                            GtError *err,
                                            GtBitsequence *bitseq,
                                            const unsigned int sample_testratio,
                                            const GtUword cbs_testsize)
{
  int had_err = 0;
  unsigned int class, len, block_len;
  GtUword block_offset = 0,
          offsets_bitpos, sample, block,
          idx, jdx;
  GtBitsequence orig_block;
  GtCompressedBitsequence *cbs;

  gt_error_check(err);

  cbs = gt_compressed_bitsequence_new(bitseq, sample_testratio, cbs_testsize);
  gt_ensure(cbs != NULL);
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
      offsets_bitpos = (GtUword)
        gt_compressed_bitsequence_get_variable_field(
                                     cbs->superblockoffsets,
                                     (sample - 1) * cbs->superblockoffsets_bits,
                                     cbs->superblockoffsets_bits);
    }
    for (jdx = sample * cbs->superblocksize; jdx < idx; jdx++) {
      class = gt_compressed_bitsequence_get_class(cbs, jdx);
      offsets_bitpos += gt_popcount_tab_offset_bits(cbs->popcount_tab, class);
    }
    class = gt_compressed_bitsequence_get_class(cbs, idx);
    len = gt_popcount_tab_offset_bits(cbs->popcount_tab, class);
    block_offset = (GtUword)
      gt_compressed_bitsequence_get_variable_field(cbs->c_offsets,
                                                   offsets_bitpos, len);
    block = gt_popcount_tab_get(cbs->popcount_tab, class, block_offset);
    orig_block = gt_compressed_bitsequence_get_variable_field(
                                                         bitseq,
                                                         idx * cbs->blocksize,
                                                         block_len);
    gt_ensure(block == (GtUword) orig_block);
  }
  gt_compressed_bitsequence_delete(cbs);
  return had_err;
}

int gt_compressed_bitsequence_unit_test(GtError *err)
{
  const unsigned int sample_testratio = 32U;
  int had_err = 0;
  const GtUword bitseq_testsize = 256UL,
                cbs_testsize = (GT_LOGWORDSIZE == 6) ?
                               16380UL :
                               8192UL;
  GtUword idx;
  GtBitsequence *bitseq;
  GtCompressedBitsequence *cbs;

  gt_error_check(err);

  bitseq = gt_malloc((size_t) bitseq_testsize * sizeof (*bitseq));

  for (idx = 0; idx < bitseq_testsize; idx++) {
    bitseq[idx] = (GtBitsequence) ((GT_LOGWORDSIZE == 6) ?
                                   (GtBitsequence) 0xaaaaccccaaaaccccULL :
                                   (GtBitsequence) 0xaaaaccccULL);
  }

  had_err =
    gt_compressed_bitsequence_unit_test_variable_field(err, bitseq,
                                                       bitseq_testsize);

  if (!had_err) {
    for (idx = 0; idx < bitseq_testsize; idx++) {
      bitseq[idx] = (GtBitsequence) ((GT_LOGWORDSIZE == 6) ?
                                     (GtBitsequence) 0xaaaaccccaaaaccccULL :
                                     (GtBitsequence) 0xaaaaccccUL);
    }
    had_err =
      gt_compressed_bitsequence_unit_test_block_identical(err, bitseq,
                                                          sample_testratio,
                                                          cbs_testsize);
  }

  if (!had_err) {
    cbs = gt_compressed_bitsequence_new(bitseq, sample_testratio, cbs_testsize);
    gt_ensure(cbs != NULL);

    if (cbs != NULL) {

      for (idx = 1UL; !had_err && idx < cbs_testsize; idx += 32UL) {
        gt_ensure(1 == gt_compressed_bitsequence_access(cbs, idx - 1));
        gt_ensure(0 == gt_compressed_bitsequence_access(cbs, idx));
      }
      gt_compressed_bitsequence_delete(cbs);
    }
  }

  /* set some blocks to 0  or ~0 to see if correct block ist found */
  bitseq[GT_DIV2(bitseq_testsize)] = (GtBitsequence) 0;
  bitseq[GT_DIV2(bitseq_testsize)+1] = (GtBitsequence) ~0ULL;
  if (!had_err) {
    cbs = gt_compressed_bitsequence_new(bitseq, sample_testratio, cbs_testsize);
    gt_ensure(cbs != NULL);

    if (cbs != NULL) {
      GtUword rank1, rank0;

      rank1 = gt_compressed_bitsequence_rank_1(cbs,
                               (GtUword) GT_INTWORDSIZE - 1UL);
      gt_ensure(rank1 == (GtUword) GT_DIV2(GT_INTWORDSIZE));
      for (idx = (GtUword) GT_INTWORDSIZE - 1UL;
           !had_err && idx < cbs_testsize;
           idx += idx + 1) {
        rank1 = gt_compressed_bitsequence_rank_1(cbs, idx);
        rank0 = gt_compressed_bitsequence_rank_0(cbs, idx);
        gt_ensure(rank1 + rank0 == idx+1);
      }
      rank1 = gt_compressed_bitsequence_rank_1(cbs, cbs_testsize-1);
      rank0 = gt_compressed_bitsequence_rank_0(cbs, cbs_testsize-1);
      gt_ensure(rank1 + rank0 == cbs_testsize);
      gt_compressed_bitsequence_delete(cbs);
    }
  }

  if (!had_err) {
    cbs = gt_compressed_bitsequence_new(bitseq, sample_testratio, cbs_testsize);
    if (cbs != NULL) {
      GtUword ranktotal1, select1, rank1, select0, rank0;
      ranktotal1 = gt_compressed_bitsequence_rank_1(cbs, cbs_testsize - 1);
      for (idx = 1UL; !had_err && idx <= ranktotal1; idx++) {
        select1 = gt_compressed_bitsequence_select_1(cbs, idx);
        gt_ensure(gt_compressed_bitsequence_access(cbs, select1) == 1);
        rank1 = gt_compressed_bitsequence_rank_1(cbs, select1);
        gt_ensure(idx == rank1);
      }
      for (idx = 1UL; !had_err && idx <= cbs_testsize - ranktotal1; idx++) {
        select0 = gt_compressed_bitsequence_select_0(cbs, idx);
        gt_ensure(gt_compressed_bitsequence_access(cbs, select0) == 0);
        rank0 = gt_compressed_bitsequence_rank_0(cbs, select0);
        gt_ensure(idx == rank0);
      }
    }
    gt_compressed_bitsequence_delete(cbs);
  }

  gt_free(bitseq);

  return had_err;
}
