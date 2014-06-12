/*
  Copyright (c) 2011 Joachim Bonnet <joachim.bonnet@studium.uni-hamburg.de>
  Copyright (c) 2012 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>

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

#include <inttypes.h>
#ifndef S_SPLINT_S
#include <unistd.h>
#include <sys/stat.h>
#include <fcntl.h>
#endif

#include "core/ensure.h"
#include "core/log_api.h"
#include "core/ma_api.h"
#include "core/mathsupport.h"
#include "core/safearith.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "extended/huffcode.h"
#include "extended/rbtree.h"

typedef struct GtHuffmanSymbol {
  GtUint64 freq;
  GtUword      symbol;
} GtHuffmanSymbol;

typedef struct GtHuffmanCode {
  GtBitsequence code;
  unsigned int  numofbits;
} GtHuffmanCode;

typedef struct GtHuffmanTree {
  GtHuffmanCode         code;
  GtHuffmanSymbol       symbol;
  struct GtHuffmanTree *leftchild,
                       *rightchild;
  unsigned int          reference_count;
} GtHuffmanTree;

struct GtHuffman {
  uint64_t       num_of_text_bits,    /* total bits needed to represent the text
                                       */
                 num_of_text_symbols; /* total number of characters in text */
  GtHuffmanTree *root_huffman_tree;   /* stores the final huffmantree */
  GtRBTree      *rbt_root;            /* red black tree */
  GtHuffmanCode *code_tab;            /* table for encoding */
  GtUword  num_of_coded_symbols, /* number of nodes in red black tree, */
                                      /* e.g. symbols with frequency > 0*/
                 num_of_symbols;      /* symbols with frequency >= 0 */
};

struct GtHuffmanDecoder {
  GtBitsequence             *bitsequence;
  GtHuffman                 *huffman;
  GtHuffmanDecoderGetMemFunc mem_func;
  GtHuffmanTree             *cur_node;
  void                      *info;
  GtUword              cur_bitseq,
                             cur_bit,
                             pad_length,
                             length;
  int                        mem_func_stat;
};

struct GtHuffmanBitwiseDecoder {
  GtHuffman     *huffman;
  GtHuffmanTree *cur_node;
};

static int huffman_tree_cmp(const void *tree1, const void *tree2,
                            GT_UNUSED void *unused)
{
  GtHuffmanTree *h1 = (GtHuffmanTree*) tree1,
                *h2 = (GtHuffmanTree*) tree2;

  gt_assert(h1 && h2);

  if (h1->symbol.freq < h2->symbol.freq)
     return (int) -1;
  if (h1->symbol.freq > h2->symbol.freq)
    return (int) 1;
  if (h1->symbol.symbol < h2->symbol.symbol)
    return (int) -1;
  if (h1->symbol.symbol > h2->symbol.symbol)
    return (int) 1;

  return 0;
}

static void huffman_tree_delete(void *tree)
{
  if (tree != NULL) {
    GtHuffmanTree *h_tree = (GtHuffmanTree*) tree;
    if (h_tree->reference_count > 0) {
      h_tree->reference_count--;
    }
    else {
      if (h_tree->leftchild != NULL)
        huffman_tree_delete(h_tree->leftchild);
      if (h_tree->rightchild != NULL)
        huffman_tree_delete(h_tree->rightchild);
      gt_free(h_tree);
    }
  }
}

static GtHuffmanTree *huffman_tree_new(GtUword symbol,
                                       GtUint64 freq)
{
  GtHuffmanTree *huffptr;
  huffptr = gt_malloc(sizeof (GtHuffmanTree));
  huffptr->code.code = 0;
  huffptr->code.numofbits = 0;
  huffptr->symbol.symbol = symbol;
  huffptr->symbol.freq = freq;
  huffptr->leftchild = NULL;
  huffptr->rightchild = NULL;
  huffptr->reference_count = 0;
  return huffptr;
}

static GtHuffmanTree *huffman_tree_ref(GtHuffmanTree *tree)
{
  gt_assert(tree);
  tree->reference_count++;
  return tree;
}

static void initialise_rbt(GtHuffman *huffman,
                           const void *distr,
                           GtDistrFunc distrfunc)
{
  GtHuffmanTree *huffptr;
  GT_UNUSED GtHuffmanTree *huffptr2;
  GtUword i;
  GtUint64 ret;
  bool nodecreated;

  huffman->num_of_coded_symbols = 0;
  huffman->rbt_root = gt_rbtree_new(huffman_tree_cmp,
                                    huffman_tree_delete,
                                    NULL);
  gt_assert(huffman->rbt_root);

  for (i = 0; i < huffman->num_of_symbols; i++) {
    ret = distrfunc(distr, i);
    if (ret > 0) {
      huffptr = huffman_tree_new(i, distrfunc(distr, i));
      nodecreated = false;
      huffptr2 = (GtHuffmanTree*)gt_rbtree_search(huffman->rbt_root,
                                                  huffptr,
                                                  &nodecreated);
      gt_assert(nodecreated && huffptr2 != NULL);
      huffman->num_of_coded_symbols++;
    }
  }
  gt_log_log("added " GT_WU " of " GT_WU " symbols as nodes",
             huffman->num_of_coded_symbols,
             huffman->num_of_symbols);
}

static void make_huffman_tree(GtHuffman *huffman)
{
  GtHuffmanTree *n1 = NULL,
                *n2 = NULL,
                *newnode = NULL;
  GT_UNUSED GtHuffmanTree *n3 = NULL;
  GtUword i,
                symbol;
  GtUint64 freq;
  int GT_UNUSED deleted;
  bool nodecreated = false;
  gt_assert(huffman->num_of_coded_symbols <= huffman->num_of_symbols);

  if (huffman->num_of_coded_symbols == 0)
    huffman->root_huffman_tree = NULL;
  else if (huffman->num_of_coded_symbols == 1UL) {
    huffman->root_huffman_tree =
      gt_rbtree_root_key(huffman->rbt_root);
    huffman->root_huffman_tree->code.code = 0;
    huffman->root_huffman_tree->code.numofbits = 1U;
  }
  else {
    for (i = 0; i < huffman->num_of_coded_symbols - 1; i++) {
      n1 = gt_rbtree_minimum_key(huffman->rbt_root);
      n1 = huffman_tree_ref(n1);
      deleted = gt_rbtree_erase(huffman->rbt_root,
                                n1);
      gt_assert(deleted == 0);
      n2 = gt_rbtree_minimum_key(huffman->rbt_root);
      n2 = huffman_tree_ref(n2);
      deleted = gt_rbtree_erase(huffman->rbt_root,
                                n2);
      gt_assert(deleted == 0);
      symbol = n1->symbol.symbol
               < n2->symbol.symbol ? n2->symbol.symbol : n1->symbol.symbol;
      freq = n1->symbol.freq + n2->symbol.freq;
      newnode = huffman_tree_new(symbol, freq);

      if (n1->symbol.freq < n2->symbol.freq) {
        newnode->leftchild = n2;
        newnode->rightchild = n1;
      }
      else {
        newnode->leftchild = n1;
        newnode->rightchild = n2;
      }
      gt_assert(huffman->rbt_root);
      n3 = gt_rbtree_search(huffman->rbt_root, newnode, &nodecreated);
      gt_assert(nodecreated && n3);
    }
    gt_assert(newnode);
    huffman->root_huffman_tree = (GtHuffmanTree*)newnode;
    huffman->root_huffman_tree->code.code = 0;
    huffman->root_huffman_tree->code.numofbits = 0;
  }
}

static void print_huff_code(unsigned length, GtBitsequence code)
{
  if (length > 0) {
    unsigned leftbit;

    for (leftbit = 1U << (length-1);
        leftbit != 0;
        leftbit >>= 1)
      (void) putchar((code & leftbit) ? '1' : '0');
  }
}

static int store_codes(GtUword symbol,
                       GT_UNUSED GtUint64 freq,
                       const GtBitsequence code,
                       unsigned int code_len,
                       void *huffman)
{
  GtHuffman *huff = (GtHuffman*)huffman;
  huff->code_tab[symbol].code = code;
  gt_safe_assign(huff->code_tab[symbol].numofbits, code_len);
  return 0;
}

static int print_codes(GtUword symbol,
                       GtUint64 freq,
                       const GtBitsequence code,
                       unsigned int code_len,
                       GT_UNUSED void *unused)
{
#ifndef S_SPLINT_S
  printf("control symbol " GT_WU ", freq "GT_LLU", codelength %u: ",
         symbol,
         freq,
         code_len);
#else
  printf("control symbol " GT_WU ", freq " GT_WU ", codelength %u: ",
         symbol,
         (GtUword) freq,
         code_len);
#endif
  print_huff_code(code_len, code);

  printf("\n");
  return 0;
}

static int calc_size(GT_UNUSED GtUword symbol,
                     GtUint64 freq,
                     GT_UNUSED const GtBitsequence code,
                     unsigned int code_len,
                     void *huffman)
{
  GtHuffman *huff = (GtHuffman*)huffman;
  huff->num_of_text_bits += code_len * freq;
  huff->num_of_text_symbols += freq;
  return 0;
}

static void huffman_tree_set_codes_rec(GtHuffmanTree *h_tree)
{
  const char leftbit = 0, rightbit = 1;
  if (h_tree != NULL) {
    if (h_tree->leftchild != NULL) {
        h_tree->leftchild->code.code = h_tree->code.code<< 1 | leftbit;
        h_tree->rightchild->code.code = h_tree->code.code<< 1 | rightbit;
        h_tree->leftchild->code.numofbits =
          h_tree->rightchild->code.numofbits = h_tree->code.numofbits + 1;
        huffman_tree_set_codes_rec(h_tree->leftchild);
        huffman_tree_set_codes_rec(h_tree->rightchild);
    }
  }
}

static inline int huffman_leaf_call_actfunc(GtHuffmanTree *h_tree,
                                            void *actinfo,
                                            GtHuffmanActFunc actfun) {
  return actfun(h_tree->symbol.symbol,
                h_tree->symbol.freq,
                h_tree->code.code,
                h_tree->code.numofbits,
                actinfo);
}

static int visit_huffman_leaves_rec(GtHuffmanTree *h_tree,
                                    void *actinfo,
                                    GtHuffmanActFunc actfun)
{
  int had_err = 0;

  if (h_tree->leftchild == NULL) {
      had_err = huffman_leaf_call_actfunc(h_tree,
                                          actinfo,
                                          actfun);
  } else {
    had_err =  visit_huffman_leaves_rec(h_tree->leftchild,
                                        actinfo,
                                        actfun);
    if (!had_err) {
      had_err =  visit_huffman_leaves_rec(h_tree->rightchild,
                                          actinfo,
                                          actfun);
    }
  }
  return had_err;
}

GtHuffman *gt_huffman_new(const void *distribution,
                          GtDistrFunc distr_func,
                          GtUword num_of_symbols)
{
  GtHuffman *huff;

  gt_assert(distribution && distr_func != NULL);
  gt_assert(num_of_symbols > 0);

  huff = gt_malloc(sizeof (GtHuffman));
  huff->num_of_symbols = num_of_symbols;
  initialise_rbt(huff, distribution, distr_func);
  make_huffman_tree(huff);

  huff->code_tab = gt_calloc((size_t) huff->num_of_symbols,
                             sizeof (GtHuffmanCode));

  huff->num_of_text_symbols = 0;
  huff->num_of_text_bits = 0;
  huffman_tree_set_codes_rec(huff->root_huffman_tree);
  (void) gt_huffman_iterate(huff, calc_size, huff);
  (void) gt_huffman_iterate(huff, store_codes, huff);

  return huff;
}

void gt_huffman_delete(GtHuffman *huffman)
{
  if (huffman != NULL) {
    gt_rbtree_delete(huffman->rbt_root);
    gt_free(huffman->code_tab);
  }
  gt_free(huffman);
}

int gt_huffman_iterate(const GtHuffman *huffman,
                       GtHuffmanActFunc action_func,
                       void *action_info)
{
  gt_assert(huffman != NULL);
  gt_assert(action_func != NULL);

  if (huffman->root_huffman_tree != NULL) {
    return visit_huffman_leaves_rec(huffman->root_huffman_tree,
                                    action_info,
                                    action_func);
  }
  return 0;
}

void gt_huffman_size(const GtHuffman *huffman, uint64_t *bits,
                     uint64_t *chars)
{
  gt_assert(huffman != NULL);
  if (bits != NULL)
    *bits = huffman->num_of_text_bits;
  if (chars != NULL)
    *chars = huffman->num_of_text_symbols;
}

void gt_huffman_print_codes(const GtHuffman *huffman)
{
  gt_assert(huffman != NULL);
  (void) gt_huffman_iterate(huffman, print_codes, NULL);
}

void gt_huffman_encode(const GtHuffman *huffman,
                       GtUword symbol,
                       GtBitsequence *code,
                       unsigned int *codelength)
{
  gt_assert(huffman != NULL);
  gt_assert(symbol < huffman->num_of_symbols);
  *code = huffman->code_tab[symbol].code;
  *codelength = huffman->code_tab[symbol].numofbits;
}

GtUword gt_huffman_numofsymbols(const GtHuffman *huffman)
{
  gt_assert(huffman != NULL);
  return huffman->num_of_coded_symbols;
}

GtUword gt_huffman_totalnumofsymbols(const GtHuffman *huffman)
{
  gt_assert(huffman != NULL);
  return huffman->num_of_symbols;
}

GtHuffmanDecoder *gt_huffman_decoder_new(GtHuffman *huffman,
                                         GtBitsequence *bitsequence,
                                         GtUword length,
                                         GtUword bit_offset,
                                         GtUword pad_length)
{
  GtHuffmanDecoder *huff_decoder = gt_malloc(sizeof (*huff_decoder));

  gt_assert(huffman != NULL);

  huff_decoder->huffman = huffman;
  huff_decoder->cur_node = huff_decoder->huffman->root_huffman_tree;
  huff_decoder->bitsequence = bitsequence;
  huff_decoder->cur_bit = bit_offset;
  huff_decoder->pad_length = pad_length;
  huff_decoder->length = length;
  huff_decoder->mem_func = NULL;
  huff_decoder->info = NULL;
  huff_decoder->cur_bitseq = 0;
  huff_decoder->mem_func_stat = 0;
  return huff_decoder;
}

GtHuffmanDecoder *gt_huffman_decoder_new_from_memory(
                                     GtHuffman *huffman,
                                     GtHuffmanDecoderGetMemFunc mem_func,
                                     void *info,
                                     GtError *err)
{
  GtHuffmanDecoder *huff_decoder = gt_malloc(sizeof (*huff_decoder));

  gt_assert(huffman != NULL);

  huff_decoder->huffman = huffman;
  huff_decoder->cur_node = huff_decoder->huffman->root_huffman_tree;
  huff_decoder->mem_func = mem_func;
  huff_decoder->info = info;
  huff_decoder->cur_bitseq = 0;
  huff_decoder->mem_func_stat =
    huff_decoder->mem_func(&huff_decoder->bitsequence,
                           &huff_decoder->length,
                           &huff_decoder->cur_bit,
                           &huff_decoder->pad_length,
                           huff_decoder->info);
  if (huff_decoder->mem_func_stat == -1) {
    gt_error_set(err, "error calling mem_func");
    return NULL;
  }
  gt_log_log(GT_WU ", " GT_WU ", " GT_WU, huff_decoder->length,
             huff_decoder->cur_bit, huff_decoder->pad_length);
  gt_log_log("got new memchunk, returned %d", huff_decoder->mem_func_stat);
  gt_assert(huff_decoder->mem_func_stat == 1);
  return huff_decoder;
}

int gt_huffman_decoder_get_new_mem_chunk(GtHuffmanDecoder *huff_decoder,
                                         GtError *err)
{
  int had_err = 0;
  gt_assert(huff_decoder != NULL);
  gt_assert(huff_decoder->mem_func != NULL);

  huff_decoder->mem_func_stat =
    huff_decoder->mem_func(&huff_decoder->bitsequence,
                           &huff_decoder->length,
                           &huff_decoder->cur_bit,
                           &huff_decoder->pad_length,
                           huff_decoder->info);
  if (huff_decoder->mem_func_stat == -1) {
    gt_error_set(err, "GtHuffmanDecoder error: error calling callback "
                      "function to get new memory chunk!");
    return huff_decoder->mem_func_stat;
  }
  huff_decoder->cur_bitseq = 0;
  gt_log_log(GT_WU ", " GT_WU ", " GT_WU, huff_decoder->length,
             huff_decoder->cur_bit, huff_decoder->pad_length);
  gt_log_log("got new memchunk, returned %d", huff_decoder->mem_func_stat);
  return had_err;
}

int gt_huffman_decoder_next(GtHuffmanDecoder *huff_decoder,
                            GtArray *symbols,
                            GtUword symbols_to_read,
                            GtError *err)
{
  int had_err = 0,
      bits_to_read = GT_INTWORDSIZE;
  GtUword read_symbols = 0;

  gt_assert((symbols_to_read > 0) && huff_decoder &&
            (gt_array_elem_size(symbols) == sizeof (GtUword)));

  if (huff_decoder->cur_bitseq == huff_decoder->length - 1)
    gt_safe_assign(bits_to_read, (GT_INTWORDSIZE - huff_decoder->pad_length));

  /* reached end of data in a previous call, didn't get reset */
  if (huff_decoder->cur_bitseq == huff_decoder->length) {
    if (!huff_decoder->mem_func)
      return 0;  /*EOF*/
    else if (huff_decoder->mem_func_stat == 0)
      return 0;
    else {
      had_err = -1;
      gt_error_set(err, "huff decoder reached EOF");
    }
  }

  while (!had_err && read_symbols < symbols_to_read) {

    /* huffman was initialized with empty dist */
    gt_assert(huff_decoder->cur_node != NULL);

    if (!had_err && huff_decoder->cur_bit == (GtUword) bits_to_read) {
      huff_decoder->cur_bitseq++;

      if (huff_decoder->cur_bitseq == huff_decoder->length - 1)
        gt_safe_assign(bits_to_read,
                       (GT_INTWORDSIZE - huff_decoder->pad_length));

      /* no next bitseq ask for more if possible */
      if (huff_decoder->cur_bitseq == huff_decoder->length ) {
        if (!huff_decoder->mem_func) {
          return 0;  /*EOF*/
        }

        gt_log_log("reset because end of block");
        huff_decoder->mem_func_stat =
          huff_decoder->mem_func(&huff_decoder->bitsequence,
                                 &huff_decoder->length,
                                 &huff_decoder->cur_bit,
                                 &huff_decoder->pad_length,
                                 huff_decoder->info);
        gt_log_log(GT_WU ", " GT_WU ", " GT_WU, huff_decoder->length,
                   huff_decoder->cur_bit, huff_decoder->pad_length);
        if (huff_decoder->mem_func_stat == -1) {
          gt_error_set(err, "error calling mem_func");
          had_err = huff_decoder->mem_func_stat;
        }
        else if (huff_decoder->mem_func_stat == 0) {
          return huff_decoder->mem_func_stat; /*EOF*/
        }

        huff_decoder->cur_bitseq = 0;
        bits_to_read = GT_INTWORDSIZE;
        if (huff_decoder->cur_bitseq == huff_decoder->length - 1)
          bits_to_read -= huff_decoder->pad_length;
      }
      huff_decoder->cur_bit = 0;
    }

    /* read the next bit */
    else if (!had_err) {
      /* nothing to do if huffman tree has only one node, just advance the bits
       */
      if (huff_decoder->cur_node->leftchild == NULL)
        huff_decoder->cur_bit++;
      else {
        if (GT_ISBITSET(huff_decoder->bitsequence[huff_decoder->cur_bitseq],
                        huff_decoder->cur_bit++))
          huff_decoder->cur_node = huff_decoder->cur_node->rightchild;
        else
          huff_decoder->cur_node = huff_decoder->cur_node->leftchild;
      }
    }

    /* symbol found, reset huffman */
    if (!had_err && huff_decoder->cur_node->leftchild == NULL) {
      gt_array_add(symbols, huff_decoder->cur_node->symbol.symbol);
      read_symbols++;
      huff_decoder->cur_node = huff_decoder->huffman->root_huffman_tree;
    }
  }
  if (!had_err) {
    /* all read, but not eof */
    return 1;
  }
  return had_err;
}

void gt_huffman_decoder_delete(GtHuffmanDecoder *huff_decoder)
{
  gt_free(huff_decoder);
}

GtHuffmanBitwiseDecoder *gt_huffman_bitwise_decoder_new(GtHuffman *huffman,
                                                        GT_UNUSED GtError *err)
{
  GtHuffmanBitwiseDecoder *hbwd;
  gt_assert(huffman != NULL);
  hbwd = gt_calloc((size_t) 1, sizeof (*hbwd));
  hbwd->huffman = huffman;
  hbwd->cur_node = huffman->root_huffman_tree;
  return hbwd;
}

int gt_huffman_bitwise_decoder_next(GtHuffmanBitwiseDecoder *hbwd,
                                    bool bit, GtUword *symbol,
                                    GT_UNUSED GtError *err)
{
  int status = 0;
  gt_assert(hbwd != NULL);

  gt_assert(hbwd->cur_node != NULL);
  if (hbwd->cur_node->leftchild == NULL) {
    *symbol = hbwd->cur_node->symbol.symbol;
    hbwd->cur_node = hbwd->huffman->root_huffman_tree;
  }
  else {
    if (bit)
      hbwd->cur_node = hbwd->cur_node->rightchild;
    else
      hbwd->cur_node = hbwd->cur_node->leftchild;

    if (hbwd->cur_node->leftchild == NULL) {
      *symbol = hbwd->cur_node->symbol.symbol;
      hbwd->cur_node = hbwd->huffman->root_huffman_tree;
    }
    else
      status = 1;
  }
  return status;
}

void gt_huffman_bitwise_decoder_delete(GtHuffmanBitwiseDecoder *hbwd)
{
  gt_free(hbwd);
}

GtUint64 unit_test_distr_func(const void *distr, GtUword symbol)
{
  GtUint64 *distrull = (GtUint64*) distr;
  return distrull[symbol];
}

static int test_bitwise(GtError *err)
{
  int had_err = 0,
      stat;
  GtUword i, j,
                symbol;
  unsigned length1 = 0,
           length2 = 0;
  GtHuffman *huffman;
  GtHuffmanBitwiseDecoder *hbwd;
  GtBitsequence bitseq;
  GtUint64 distr[6] = {45ULL, 16ULL, 13ULL, 12ULL, 9ULL, 5ULL};

  huffman = gt_huffman_new(&distr, unit_test_distr_func, 6UL);

  hbwd = gt_huffman_bitwise_decoder_new(huffman, err);
  gt_ensure(hbwd);
  for (i = 0; !had_err && i < 6UL; i++) {
    gt_huffman_encode(huffman, i, &bitseq, &length2);
    if (i > 0) {
      gt_ensure(length1 <= length2);
    }
    j = 1UL;
    while (1 == (stat = gt_huffman_bitwise_decoder_next(hbwd,
                                           ((bitseq >> (length2 - j)) & 1) != 0,
                                           &symbol, err))) {
      j++;
    }
    gt_ensure(stat == 0);
    gt_ensure(j == (GtUword) length2);
    gt_ensure(symbol == i);
    length1 = length2;
  }

  stat = gt_huffman_bitwise_decoder_next(hbwd, true, &symbol, err);

  gt_ensure(stat == 0);
  gt_ensure(symbol == 0);

  stat = gt_huffman_bitwise_decoder_next(hbwd, false, &symbol, err);
  gt_ensure(stat == 1);
  stat = gt_huffman_bitwise_decoder_next(hbwd, false, &symbol, err);
  gt_ensure(stat == 1);
  stat = gt_huffman_bitwise_decoder_next(hbwd, false, &symbol, err);
  gt_ensure(stat == 0);
  gt_ensure(symbol == 1UL);

  stat = gt_huffman_bitwise_decoder_next(hbwd, false, &symbol, err);
  gt_ensure(stat == 1);
  stat = gt_huffman_bitwise_decoder_next(hbwd, true, &symbol, err);
  gt_ensure(stat == 1);
  stat = gt_huffman_bitwise_decoder_next(hbwd, false, &symbol, err);
  gt_ensure(stat == 0);
  gt_ensure(symbol == 2UL);

  stat = gt_huffman_bitwise_decoder_next(hbwd, false, &symbol, err);
  gt_ensure(stat == 1);
  stat = gt_huffman_bitwise_decoder_next(hbwd, true, &symbol, err);
  gt_ensure(stat == 1);
  stat = gt_huffman_bitwise_decoder_next(hbwd, true, &symbol, err);
  gt_ensure(stat == 0);
  gt_ensure(symbol == 3UL);

  stat = gt_huffman_bitwise_decoder_next(hbwd, false, &symbol, err);
  gt_ensure(stat == 1);
  stat = gt_huffman_bitwise_decoder_next(hbwd, false, &symbol, err);
  gt_ensure(stat == 1);
  stat = gt_huffman_bitwise_decoder_next(hbwd, true, &symbol, err);
  gt_ensure(stat == 1);
  stat = gt_huffman_bitwise_decoder_next(hbwd, false, &symbol, err);
  gt_ensure(stat == 0);
  gt_ensure(symbol == 4UL);

  stat = gt_huffman_bitwise_decoder_next(hbwd, false, &symbol, err);
  gt_ensure(stat == 1);
  stat = gt_huffman_bitwise_decoder_next(hbwd, false, &symbol, err);
  gt_ensure(stat == 1);
  stat = gt_huffman_bitwise_decoder_next(hbwd, true, &symbol, err);
  gt_ensure(stat == 1);
  stat = gt_huffman_bitwise_decoder_next(hbwd, true, &symbol, err);
  gt_ensure(stat == 0);
  gt_ensure(symbol == 5UL);

  gt_huffman_delete(huffman);
  gt_huffman_bitwise_decoder_delete(hbwd);
  return had_err;
}

typedef struct huffman_unit_test_meminfo {
  GtBitsequence *data;
  GtUword  chunk,
                 chunks,
                 size,
                 lastchunk_size,
                 padding;
} HuffmanUnitTestMeminfo;

static int huffman_unit_get_next_block(GtBitsequence **bits,
                                       GtUword *length,
                                       GtUword *offset,
                                       GtUword *pad_length,
                                       void *meminfo)
{
  HuffmanUnitTestMeminfo *info = (HuffmanUnitTestMeminfo*) meminfo;
  if (info->chunk == info->chunks)
    return 0;
  *offset = 0;
  if (info->chunk == info->chunks - 1) {
    *bits = info->data + (info->chunk * info->size);
    *length = info->lastchunk_size;
    *pad_length = info->padding;
  }
  else {
    *bits = info->data + (info->chunk * info->size);
    *length = info->size;
    *pad_length = 0;
  }
  info->chunk++;
  return 1;
}

int test_mem(GtError *err)
{
  unsigned char bits = (unsigned char) (sizeof (GtBitsequence) * 8),
                bits_remain = bits;
  unsigned int code_len;
  int decoder_stat = 1,
      had_err = 0;
  GtUword const max_num = (GtUword) (1000 * gt_rand_0_to_1() + 100),
                      dist_size = (GtUword) (10 + 10 * gt_rand_0_to_1()),
                      step_size = (GtUword) (10 + 5 * gt_rand_0_to_1());
  GtUword range = 0,
                idx, idx_j;
  uint64_t total_bits = 0;
  GtUint64 *distribution = gt_malloc(sizeof (*distribution) *
                                     dist_size);
  GtBitsequence buffer = 0,
                code;
  GtArray *numbers = gt_array_new(sizeof (unsigned char));
  GtArray *encoded = gt_array_new(sizeof (GtUword));
  GtArray *codes = gt_array_new(sizeof (GtBitsequence));
  GtHuffman *huff = NULL;
  GtHuffmanDecoder *huffdec = NULL;
  HuffmanUnitTestMeminfo *meminfo = gt_malloc(sizeof (*meminfo));

  for (idx = 0; idx < dist_size; idx++) {
    distribution[idx] = 0;
  }

  /*produce numbers with more numbers 0-(dist_size/2) */
  while (gt_array_size(numbers) < max_num) {
    unsigned char x;
    if (range) {
      x = (unsigned char) ((dist_size >> 1) * gt_rand_0_to_1());
      range = 0;
    }
    else {
      x = (unsigned char) (dist_size * gt_rand_0_to_1());
      range = 1UL;
    }
    gt_array_add(numbers, x);
  }

  gt_ensure(gt_array_size(numbers) == max_num);

  /*calculate distribution*/
  for (idx = 0; !had_err && idx < max_num; idx++) {
    unsigned char *c = (unsigned char*) gt_array_get(numbers, idx);
    distribution[*c]++;
  }

  if (!had_err) {
    uint64_t chars;
    huff = gt_huffman_new(distribution, unit_test_distr_func, dist_size);
    gt_ensure(huff);
    gt_huffman_size(huff, &total_bits, &chars);
    gt_ensure(chars == (uint64_t) max_num);
  }

  for (idx = 0; !had_err && idx < max_num; idx++) {
    gt_huffman_encode(huff,
                      *(unsigned char*) gt_array_get(numbers, idx),
                      &code,
                      &code_len);
    if (bits_remain < code_len) {
      unsigned char overhang = code_len - bits_remain;
      buffer |= code >> overhang;
      gt_array_add(codes, buffer);
      buffer = 0;
      bits_remain = bits - overhang;
    }
    else
      bits_remain -= code_len;
    buffer |= code << bits_remain;
  }
  if (!had_err)
    gt_array_add(codes, buffer);

  gt_ensure(
            total_bits == (uint64_t) gt_array_size(codes) *
                          gt_array_elem_size(codes) * 8 - bits_remain);

  if (!had_err) {
    meminfo->chunk = 0;
    meminfo->padding = bits_remain;
    meminfo->size = gt_array_size(codes) / 4;
    meminfo->lastchunk_size = gt_array_size(codes) % meminfo->size;
    if (meminfo->lastchunk_size) {
      meminfo->chunks = gt_array_size(codes) / meminfo->size + 1;
    }
    else {
      meminfo->chunks = gt_array_size(codes) / meminfo->size;
      meminfo->lastchunk_size = meminfo->size;
    }
    meminfo->data = gt_array_get_space(codes);

    huffdec = gt_huffman_decoder_new(huff,
                                     gt_array_get_space(codes),
                                     gt_array_size(codes),
                                     0, /*offset*/
                                     bits_remain);
    gt_ensure(huffdec);
  }

  for (idx = 0; !had_err && idx < max_num; idx += step_size) {
    gt_array_reset(encoded);
    gt_ensure(decoder_stat == 1);
    decoder_stat = gt_huffman_decoder_next(huffdec,
                                           encoded,
                                           (GtUword) step_size,
                                           err);
    if (decoder_stat == 1)
      gt_ensure(gt_array_size(encoded) == (GtUword) step_size);
    else
      gt_ensure(decoder_stat == 0);

    gt_ensure(decoder_stat != -1);
    if (!had_err) {
      for (idx_j = 0; !had_err && idx_j < gt_array_size(encoded); idx_j++) {
        GtUword ist = *(GtUword*) gt_array_get(encoded, idx_j);
        unsigned char soll =
          *(unsigned char*)gt_array_get(numbers, idx + idx_j);
        gt_ensure((ist == soll));
      }
    }
  }

  gt_ensure(idx >= max_num);
  if (!had_err)
    decoder_stat = gt_huffman_decoder_next(huffdec,
                                           encoded,
                                           1UL,
                                           err);
  gt_ensure(decoder_stat == 0);

  if (!had_err) {
    gt_huffman_decoder_delete(huffdec);

    huffdec = gt_huffman_decoder_new_from_memory(huff,
                                                 huffman_unit_get_next_block,
                                                 meminfo,
                                                 err);
    gt_ensure(huffdec);
  }

  decoder_stat = 1;
  for (idx = 0; !had_err && idx < max_num; idx += step_size) {
    gt_array_reset(encoded);
    gt_ensure(decoder_stat == 1);
    decoder_stat = gt_huffman_decoder_next(huffdec,
                                           encoded,
                                           step_size,
                                           err);
    if (decoder_stat == 1)
      gt_ensure(gt_array_size(encoded) == step_size);
    else
      gt_ensure(decoder_stat == 0);

    gt_ensure(decoder_stat != -1);
    if (!had_err) {
      for (idx_j = 0; !had_err && idx_j < gt_array_size(encoded); idx_j++) {
        GtUword ist = *(GtUword*)gt_array_get(encoded, idx_j);
        unsigned char soll =
          *(unsigned char*)gt_array_get(numbers, idx + idx_j);
        gt_ensure((ist == soll));
      }
    }
  }

  gt_ensure(idx >= max_num);
  if (!had_err)
    decoder_stat = gt_huffman_decoder_next(huffdec,
                                           encoded,
                                           1UL,
                                           err);
  gt_ensure(decoder_stat == 0);

  gt_huffman_decoder_delete(huffdec);
  gt_free(meminfo);
  gt_huffman_delete(huff);
  gt_array_delete(numbers);
  gt_array_delete(codes);
  gt_array_delete(encoded);
  gt_free(distribution);
  return had_err;
}

int gt_huffman_unit_test(GtError *err)
{
  int had_err = 0;

  had_err = test_bitwise(err);

  if (!had_err)
    had_err = test_mem(err);

  return had_err;
}
