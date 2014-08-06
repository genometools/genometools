/*
   Copyright (c) 2013 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
   Copyright (c) 2013 Center for Bioinformatics, University of Hamburg

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

#include <stdbool.h>
#ifndef S_SPLINT_S
#include <ctype.h>
#endif

#include "core/alphabet_api.h"
#include "core/chardef.h"
#include "core/divmodmul.h"
#include "core/encseq_api.h"
#include "core/intbits.h"
#include "core/log_api.h"
#include "core/ma_api.h"
#include "core/mathsupport.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "extended/compressed_bitsequence.h"
#include "extended/wtree_encseq.h"
#include "extended/wtree_rep.h"

typedef struct GtWtreeEncseqFillOffset {
  struct GtWtreeEncseqFillOffset *left,
                                 *right;
  GtUword                         offset,
                                  left_size;
}GtWtreeEncseqFillOffset;

struct GtWtreeEncseq {
  GtWtree                  parent_instance;
  GtEncseq                *encseq;
  GtAlphabet              *alpha;
  GtBitsequence           *bits;
  GtWtreeEncseqFillOffset *root_fo,
                          *current_fo;
  GtCompressedBitsequence *c_bits;
  GtUword                  bits_size,
                           node_start,
                           num_of_bits;
  unsigned int             alpha_size,
                           levels;
};

const GtWtreeClass* gt_wtree_encseq_class(void);

#define gt_wtree_encseq_cast(wtree) \
  gt_wtree_cast(gt_wtree_encseq_class(), wtree)

static GtWtreeSymbol gt_wtree_encseq_access_rec(GtWtreeEncseq *we,
                                                GtUword pos,
                                                GtUword node_start,
                                                GtUword node_size,
                                                unsigned int alpha_start,
                                                unsigned int alpha_end)
{
  unsigned int middle = GT_DIV2(alpha_start + alpha_end);
  int bit;
  GtUword zero_rank_prefix = 0,
          one_rank_prefix = 0,
          left_child_size;
  gt_assert(pos < node_size);

  if (alpha_start < alpha_end) {
    bit = gt_compressed_bitsequence_access(we->c_bits, node_start + pos);
    if (node_start != 0)
      zero_rank_prefix =
        gt_compressed_bitsequence_rank_0(we->c_bits, node_start - 1);
    left_child_size =
      gt_compressed_bitsequence_rank_0(we->c_bits, node_start + node_size - 1) -
      zero_rank_prefix;

    if (bit == 0) {
      pos = gt_compressed_bitsequence_rank_0(we->c_bits, node_start + pos) -
        zero_rank_prefix - 1; /*convert count (rank) to positon */
      alpha_end = middle;
      node_start += we->parent_instance.members->length;
      node_size = left_child_size;
      return gt_wtree_encseq_access_rec(we, pos, node_start,
                                        node_size, alpha_start, alpha_end);
    }
    else {
      if (node_start != 0)
        one_rank_prefix =
          gt_compressed_bitsequence_rank_1(we->c_bits, node_start - 1);
      pos = gt_compressed_bitsequence_rank_1(we->c_bits, node_start + pos) -
        one_rank_prefix - 1; /*convert count (rank) to positon */
      alpha_start = middle + 1;
      node_size =
        gt_compressed_bitsequence_rank_1(we->c_bits,
                                         node_start + node_size - 1) -
        one_rank_prefix;
      node_start +=
        we->parent_instance.members->length + left_child_size;
      return gt_wtree_encseq_access_rec(we, pos, node_start,
                                        node_size, alpha_start, alpha_end);
    }
  }
  return (GtWtreeSymbol) alpha_start;
}

static GtWtreeSymbol gt_wtree_encseq_access(GtWtree *wtree,
                                            GtUword pos)
{
  unsigned int alpha_start = 0,
               alpha_end;
  GtUword node_start = 0,
          node_size;
  GtWtreeEncseq *we;
  gt_assert(wtree != NULL);

  we = gt_wtree_encseq_cast(wtree);
  gt_assert(pos < wtree->members->length);

  alpha_end = we->alpha_size - 1;
  node_size = wtree->members->length;

  return gt_wtree_encseq_access_rec(we, pos, node_start, node_size, alpha_start,
                                    alpha_end);
}

static GtUword gt_wtree_encseq_rank_rec(GtWtreeEncseq *we,
                                        GtUword pos,
                                        GtWtreeSymbol sym,
                                        GtUword node_start,
                                        GtUword node_size,
                                        unsigned int alpha_start,
                                        unsigned int alpha_end)
{
  unsigned int middle = GT_DIV2(alpha_start + alpha_end);
  int bit;
  GtUword zero_rank_prefix = 0,
          one_rank_prefix = 0,
          left_child_size,
          rank;
  gt_log_log("alphabet: %u-%u-%u, sym: " GT_WU,
             alpha_start, middle, alpha_end, (GtUword) sym);
  gt_log_log("pos: "GT_WU"", pos);
  gt_assert(pos < node_size);

  if (alpha_start < alpha_end) {
    bit = middle < (unsigned int) sym ? 1 : 0;
    if (node_start != 0)
      zero_rank_prefix =
        gt_compressed_bitsequence_rank_0(we->c_bits, node_start - 1);
    left_child_size =
      gt_compressed_bitsequence_rank_0(we->c_bits, node_start + node_size - 1) -
      zero_rank_prefix;

    if (bit == 0) {
      rank = gt_compressed_bitsequence_rank_0(we->c_bits, node_start + pos) -
        zero_rank_prefix;
      alpha_end = middle;
      node_start += we->parent_instance.members->length;
      node_size = left_child_size;
    }
    else {
      if (node_start != 0)
        one_rank_prefix =
          gt_compressed_bitsequence_rank_1(we->c_bits, node_start - 1);
      rank = gt_compressed_bitsequence_rank_1(we->c_bits, node_start + pos) -
        one_rank_prefix;
      alpha_start = middle + 1;
      node_size =
        gt_compressed_bitsequence_rank_1(we->c_bits,
                                         node_start + node_size - 1) -
        one_rank_prefix;
      node_start +=
        we->parent_instance.members->length + left_child_size;
    }
    gt_log_log("bit: %d, nodesize: "GT_WU"", bit, node_size);
    if (node_size != 0 && rank != 0) {
      pos = rank - 1;
      return gt_wtree_encseq_rank_rec(we, pos, sym,
                                      node_start, node_size, alpha_start,
                                      alpha_end);
    }
    return 0;
  }
  gt_log_log("found: rank="GT_WU"", pos + 1);
  return pos + 1; /* convert position to count */
}

static GtUword gt_wtree_encseq_rank(GtWtree *wtree,
                                    GtUword pos,
                                    GtWtreeSymbol symbol)
{
  unsigned int alpha_start = 0,
               alpha_end;
  GtUword node_start = 0,
          node_size;
  GtWtreeEncseq *we;
  gt_assert(wtree != NULL);

  we = gt_wtree_encseq_cast(wtree);
  gt_assert(pos < wtree->members->length);

  alpha_end = we->alpha_size - 1;
  node_size = wtree->members->length;

  return gt_wtree_encseq_rank_rec(we, pos, symbol, node_start, node_size,
                                  alpha_start, alpha_end);
}

static GtUword gt_wtree_encseq_select_rec(GtWtreeEncseq *we,
                                          GtUword i,
                                          GtWtreeSymbol sym,
                                          GtUword node_start,
                                          GtUword node_size,
                                          unsigned int alpha_start,
                                          unsigned int alpha_end)
{
  unsigned int middle = GT_DIV2(alpha_start + alpha_end);
  int bit;
  GtUword zero_rank_prefix = 0,
          one_rank_prefix = 0,
          left_child_size, child_start;

  if (alpha_start < alpha_end) {
    bit = middle < (unsigned int) sym ? 1 : 0;
    if (node_start != 0)
      zero_rank_prefix =
        gt_compressed_bitsequence_rank_0(we->c_bits, node_start - 1);
    left_child_size =
      gt_compressed_bitsequence_rank_0(we->c_bits, node_start + node_size - 1) -
      zero_rank_prefix;

    if (bit == 0) {
      alpha_end = middle;
      child_start = node_start + we->parent_instance.members->length;
      node_size = left_child_size;
    }
    else {
      if (node_start != 0)
        one_rank_prefix =
          gt_compressed_bitsequence_rank_1(we->c_bits, node_start - 1);
      alpha_start = middle + 1;
      node_size =
        gt_compressed_bitsequence_rank_1(we->c_bits,
                                         node_start + node_size - 1) -
        one_rank_prefix;
      child_start =
        node_start + we->parent_instance.members->length + left_child_size;
    }
    if (node_size != 0) {
      i = gt_wtree_encseq_select_rec(we, i, sym, child_start, node_size,
                                     alpha_start, alpha_end);
      if (i < node_size) {
        return (bit == 0 ?
                gt_compressed_bitsequence_select_0(we->c_bits,
                                                   zero_rank_prefix + i + 1) :
                gt_compressed_bitsequence_select_1(we->c_bits,
                                                   one_rank_prefix + i + 1)) -
          node_start;
      }
    }
    return ULONG_MAX;
  }
  if (i <= node_size)
    return i - 1;
  return ULONG_MAX;
}

static GtUword gt_wtree_encseq_select(GtWtree *wtree,
                                      GtUword i,
                                      GtWtreeSymbol symbol)
{
  unsigned int alpha_start = 0,
               alpha_end;
  GtUword node_start = 0,
          node_size;
  GtWtreeEncseq *we;
  gt_assert(wtree != NULL);

  we = gt_wtree_encseq_cast(wtree);
  gt_assert(i <= wtree->members->length);
  gt_assert(i != 0);

  alpha_end = we->alpha_size - 1;
  node_size = wtree->members->length;

  return gt_wtree_encseq_select_rec(we, i, symbol, node_start, node_size,
                                    alpha_start, alpha_end);
}

static void gt_wtree_encseq_delete(GT_UNUSED GtWtree *wtree)
{
  if (wtree != NULL) {
    GtWtreeEncseq *wtree_encseq = gt_wtree_encseq_cast(wtree);
    gt_encseq_delete(wtree_encseq->encseq);
    gt_alphabet_delete(wtree_encseq->alpha);
    gt_free(wtree_encseq->bits);
    gt_compressed_bitsequence_delete(wtree_encseq->c_bits);
  }
}

static inline GtWtreeSymbol gt_wtree_encseq_map(GtWtreeEncseq *wtree_encseq,
                                                GtUchar symbol)
{
  if (ISNOTSPECIAL(symbol))
    return (GtWtreeSymbol) symbol;
  else {
    if (symbol == (GtUchar) SEPARATOR) {
      return (GtWtreeSymbol) wtree_encseq->alpha_size - 1;
    }
    if (symbol == (GtUchar) WILDCARD)
      return (GtWtreeSymbol) wtree_encseq->alpha_size - 2;
  }
  gt_assert(symbol == (GtUchar) UNDEFCHAR);
  return (GtWtreeSymbol) wtree_encseq->alpha_size - 3;
}

char gt_wtree_encseq_unmap_decoded(GtWtree *wtree,
                                   GtWtreeSymbol symbol)
{
  GtWtreeEncseq *wtree_encseq;
  GtUchar encseq_sym = (GtUchar) symbol;
  gt_assert((symbol & ((GtUchar) 0)) == 0);
  gt_assert(wtree != NULL);
  wtree_encseq = gt_wtree_encseq_cast(wtree);
  switch (wtree_encseq->alpha_size - encseq_sym) {
    case 1:
      return (char) SEPARATOR;
    case 2:
      return gt_alphabet_decode(wtree_encseq->alpha, (GtUchar) WILDCARD);
    case 3:
      return (char) UNDEFCHAR;
    default:
      return gt_alphabet_decode(wtree_encseq->alpha, encseq_sym);
  }
}

/* map static local methods to interface */
const GtWtreeClass* gt_wtree_encseq_class(void)
{
  static const GtWtreeClass *this_c = NULL;
  if (this_c == NULL) {
    this_c =
      gt_wtree_class_new(sizeof (GtWtreeEncseq), gt_wtree_encseq_access,
                         gt_wtree_encseq_rank, gt_wtree_encseq_select,
                         gt_wtree_encseq_delete);
  }
  return this_c;
}

static inline GtWtreeEncseqFillOffset* gt_wtree_encseq_fill_offset_new(void)
{
  GtWtreeEncseqFillOffset *fo = gt_malloc(sizeof (*fo));
  fo->left = fo->right = NULL;
  fo->offset = 0;
  fo->left_size = 0;
  return fo;
}

static inline void
gt_wtree_encseq_fill_offset_delete(GtWtreeEncseqFillOffset *fo)
{
  if (fo != NULL) {
    gt_wtree_encseq_fill_offset_delete(fo->left);
    gt_wtree_encseq_fill_offset_delete(fo->right);
    gt_free(fo);
  }
}

static bool gt_wtree_encseq_set_nodestart_and_current_fo(GtWtreeEncseq *we,
                                                         unsigned int level,
                                                         GtWtreeSymbol sym)
{
  unsigned int alpha_end = we->alpha_size - 1,
               alpha_start = 0,
               middle = GT_DIV2(alpha_end);
  unsigned int c_level = 0;
  gt_assert(sym <= (GtWtreeSymbol) alpha_end);
  we->current_fo = we->root_fo;
  we->node_start = 0;

  while (c_level < level &&
         alpha_end > alpha_start) {
    if (sym <= (GtWtreeSymbol) middle) {
      alpha_end = middle;
      if (alpha_end > alpha_start &&
          we->current_fo->left == NULL)
        we->current_fo->left = gt_wtree_encseq_fill_offset_new();

      we->node_start =
        we->node_start + we->parent_instance.members->length;
      we->current_fo = we->current_fo->left;
    }
    else {
      alpha_start = middle + 1;
      if (alpha_end > alpha_start &&
          we->current_fo->right == NULL)
        we->current_fo->right = gt_wtree_encseq_fill_offset_new();

      /* start of right child: start of left + size of left */
      we->node_start =
        we->node_start + we->parent_instance.members->length +
        we->current_fo->left_size;
      we->current_fo = we->current_fo->right;
    }
    middle = GT_DIV2(alpha_start + alpha_end);
    c_level++;
  }
  return (sym <= (GtWtreeSymbol) middle);
}

static void gt_wtree_encseq_fill_bits(GtWtreeEncseq *we)
{
  unsigned int level_idx;
  GtUword sym_idx;
  GtEncseqReader *er =
    gt_encseq_create_reader_with_readmode(we->encseq, GT_READMODE_FORWARD, 0);
  gt_assert(we != NULL);

  for (level_idx = 0; level_idx < we->levels; level_idx++) {
    for (sym_idx = 0;
         sym_idx < we->parent_instance.members->length;
         sym_idx++) {
      GtWtreeSymbol c_sym =
        gt_wtree_encseq_map(we, gt_encseq_reader_next_encoded_char(er));
      if (gt_wtree_encseq_set_nodestart_and_current_fo(we, level_idx, c_sym)) {
        /*0*/
        if (we->current_fo != NULL) {
          gt_assert(we->node_start + we->current_fo->offset < we->num_of_bits);
          we->current_fo->offset++;
          we->current_fo->left_size++;
        }
      }
      else {
        if (we->current_fo != NULL) {
          gt_assert(we->node_start + we->current_fo->offset < we->num_of_bits);
          GT_SETIBIT(we->bits, we->node_start + we->current_fo->offset);
          we->current_fo->offset++;
        }
      }
    }
    gt_encseq_reader_reinit_with_readmode(er, we->encseq,
                                          GT_READMODE_FORWARD, 0);
  }
  gt_encseq_reader_delete(er);
  gt_wtree_encseq_fill_offset_delete(we->root_fo);
  we->root_fo = we->current_fo = NULL;
  gt_encseq_delete(we->encseq);
  we->encseq = NULL;
}

GtWtree* gt_wtree_encseq_new(GtEncseq *encseq)
{
  /* sample rate for compressd bitseq */
  const unsigned int samplerate = 32U;
  GtWtree *wtree;
  GtWtreeEncseq *wtree_encseq;
  wtree = gt_wtree_create(gt_wtree_encseq_class());
  wtree_encseq = gt_wtree_encseq_cast(wtree);
  wtree_encseq->encseq = gt_encseq_ref(encseq);
  wtree_encseq->alpha = gt_alphabet_ref(gt_encseq_alphabet(encseq));
  /* encoded chars + WC given by gt_alphabet_size,
     we have to encode UNDEFCHAR and SEPARATOR too */
  wtree_encseq->alpha_size = gt_alphabet_size(wtree_encseq->alpha) + 2;
  wtree->members->num_of_symbols = (GtUword) wtree_encseq->alpha_size;
  /* levels in tree: \lceil log_2(\sigma)\rceil */
  wtree_encseq->levels =
    gt_determinebitspervalue((GtUword) wtree_encseq->alpha_size);
  wtree_encseq->root_fo = gt_wtree_encseq_fill_offset_new();
  wtree_encseq->current_fo = wtree_encseq->root_fo;
  wtree->members->length =
    gt_encseq_total_length(encseq);
  /* each level has number of symbols bits */
  wtree_encseq->num_of_bits =
    wtree_encseq->levels *
    wtree->members->length;
  wtree_encseq->bits_size =
    wtree_encseq->num_of_bits / (sizeof (GtBitsequence) * CHAR_BIT);
  if (wtree_encseq->num_of_bits % (sizeof (GtBitsequence) * CHAR_BIT) != 0)
    wtree_encseq->bits_size++;
  wtree_encseq->bits =
    gt_calloc((size_t) wtree_encseq->bits_size, sizeof (GtBitsequence));
  wtree_encseq->node_start = 0;
  gt_wtree_encseq_fill_bits(wtree_encseq);
  wtree_encseq->c_bits =
    gt_compressed_bitsequence_new(wtree_encseq->bits,
                                  samplerate,
                                  wtree_encseq->num_of_bits);
  gt_free(wtree_encseq->bits);
  wtree_encseq->bits = NULL;
  return wtree;
}
