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
#include <unistd.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>

#include "core/ensure.h"
#include "core/ma_api.h"
#include "core/mathsupport.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "core/log_api.h"
#include "extended/huffman.h"
#include "extended/rbtree.h"

typedef struct GtHuffmanSymbol {
  unsigned long symbol;
  unsigned long long freq;
} GtHuffmanSymbol;

typedef struct GtHuffmanCode {
  unsigned long numofbits;
  GtBitsequence code;
} GtHuffmanCode;

typedef struct GtHuffmanTree {
  GtHuffmanSymbol symbol;
  GtHuffmanCode code;
  struct GtHuffmanTree *leftchild,
                       *rightchild;
  unsigned reference_count;
} GtHuffmanTree;

struct GtHuffman {
  GtHuffmanTree *roothuffmantree;  /* stores the final huffmantree */
  unsigned long numofsymbols,      /* number of nodes in red black tree, */
                                   /* e.g. symbols with frequency > 0*/
                totalnumofsymbols; /* symbols with frequency >= 0 */

  uint64_t totalnumofbits,         /* total bits needed to represent the text */
           totalnumofchars;        /* total number of characters in text */

  GtRBTree *rbt_root;              /* red black tree*/
  GtHuffmanCode *code_tab;         /* table for encoding */
};

struct GtHuffmanDecoder {
  unsigned long cur_bitseq,            /*position in bitsequence*/
                cur_bit,               /*position in element of bitsequence*/
                pad_length,            /*number of non data bits in last element
                                        */
                length;                /*number of elements in bitsequence*/
  GtBitsequence *bitsequence;          /*array of bitstrings*/
  GtHuffman *huffman;                  /*encoding/decoding information*/
  GtHuffmanTree *cur_node;             /*position in encoding*/
  GtHuffmanDecoderGetMemFunc mem_func; /*callback to get more data, can be NULL
                                        */
  int mem_func_stat;                   /*status of mem_func, 0, 1, -1*/
  void *info;                          /*info for callback*/
};

struct GtHuffmanBitwiseDecoder {
  GtHuffman *huffman;
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
  if (!tree)
    return;
  GtHuffmanTree *h_tree = (GtHuffmanTree*)tree;
  if (h_tree->reference_count) {
    h_tree->reference_count--;
    return;
  }

  if (h_tree->leftchild != NULL)
    huffman_tree_delete(h_tree->leftchild);
  if (h_tree->rightchild != NULL)
    huffman_tree_delete(h_tree->rightchild);
  gt_free(h_tree);
}

static GtHuffmanTree *huffman_tree_new(unsigned long symbol,
                                       unsigned long long freq)
{
  GtHuffmanTree *huffptr;
  huffptr = gt_malloc(sizeof (GtHuffmanTree));
  huffptr->code.code = GT_UNDEF_UINT;
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
  GtHuffmanTree *huffptr, GT_UNUSED *huffptr2;
  unsigned long i;
  bool nodecreated = false;

  huffman->numofsymbols = 0;
  huffman->rbt_root = gt_rbtree_new(huffman_tree_cmp,
                                    huffman_tree_delete,
                                    NULL);
  gt_assert(huffman->rbt_root);

  for (i = 0; i < huffman->totalnumofsymbols; i++) {
    if (distrfunc(distr, i) > 0) {
      huffptr = huffman_tree_new(i, distrfunc(distr, i));
      huffptr2 = (GtHuffmanTree*)gt_rbtree_search(huffman->rbt_root,
                                                  huffptr,
                                                  &nodecreated);
      gt_assert(nodecreated && huffptr2);
      huffman->numofsymbols++;
    }
  }
}

static int make_huffman_tree(GtHuffman *huffman)
{
  GtHuffmanTree *n1 = NULL,
                *n2 = NULL,
                GT_UNUSED *n3 = NULL,
                *newnode = NULL;
  unsigned long i,
                symbol;
  unsigned long long freq;
  int GT_UNUSED deleted;
  bool nodecreated = false;

  if (huffman->numofsymbols == 0)
    huffman->roothuffmantree = NULL;
  else if (huffman->numofsymbols == 1) {
    huffman->roothuffmantree =
      gt_rbtree_root_key(huffman->rbt_root);
    huffman->roothuffmantree->code.code = 0;
    huffman->roothuffmantree->code.numofbits = 1;
  }
  else {
    for (i = 0; i < huffman->numofsymbols - 1; i++) {
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
    huffman->roothuffmantree = (GtHuffmanTree*)newnode;
    huffman->roothuffmantree->code.code = 0;
    huffman->roothuffmantree->code.numofbits = 0;
  }
  return 0;
}

static void print_huff_code(unsigned long length, GtBitsequence code)
{
  if (length > 0) {
    unsigned long leftbit;

    for (leftbit = 1 << (length-1);
        leftbit != 0;
        leftbit >>= 1)
      (void) putchar((code & leftbit) ? '1' : '0');
  }
}

static int store_codes(unsigned long symbol,
                       GT_UNUSED unsigned long long freq,
                       const GtBitsequence code,
                       unsigned long code_len,
                       void *huffman)
{
  GtHuffman *huff = (GtHuffman*)huffman;
  huff->code_tab[symbol].code = code;
  huff->code_tab[symbol].numofbits = code_len;
  return 0;
}

static int print_codes(unsigned long symbol,
                       unsigned long long freq,
                       const GtBitsequence code,
                       unsigned long code_len,
                       GT_UNUSED void *unused)
{
  printf("control symbol %lu, freq %llu, codelength %lu: ",
         symbol,
         freq,
         code_len);
  print_huff_code(code_len, code);

  printf("\n");
  return 0;
}

static int calc_size(GT_UNUSED unsigned long symbol,
                     unsigned long long freq,
                     GT_UNUSED const GtBitsequence code,
                     unsigned long code_len,
                     void *huffman)
{
  GtHuffman *huff = (GtHuffman*)huffman;
  huff->totalnumofbits += code_len * freq;
  huff->totalnumofchars += freq;
  return 0;
}

static void huffman_tree_set_codes_rec(GtHuffmanTree *h_tree)
{
  const char leftbit = 0, rightbit = 1;
  if (h_tree->leftchild != NULL) {
      h_tree->leftchild->code.code = h_tree->code.code<< 1 | leftbit;
      h_tree->rightchild->code.code = h_tree->code.code<< 1 | rightbit;
      h_tree->leftchild->code.numofbits =
        h_tree->rightchild->code.numofbits = h_tree->code.numofbits + 1;
      huffman_tree_set_codes_rec(h_tree->leftchild);
      huffman_tree_set_codes_rec(h_tree->rightchild);
  }
}

static int huffman_leaf_call_actfunc(GtHuffmanTree *h_tree,
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

GtHuffman *gt_huffman_new(const void *distr,
                          GtDistrFunc distrfunc,
                          unsigned long totalnumofsymbols)
{
  GtHuffman *huff;
  unsigned long i;

  gt_assert(distr && distrfunc);
  gt_assert(totalnumofsymbols > 0);

  huff = gt_malloc(sizeof (GtHuffman));
  huff->totalnumofsymbols = totalnumofsymbols;
  initialise_rbt(huff, distr, distrfunc);
  make_huffman_tree(huff);

  huff->code_tab = gt_calloc(huff->totalnumofsymbols, sizeof (GtHuffmanCode));
  for (i = 0; i < huff->totalnumofsymbols; i++)
    huff->code_tab[i].code = GT_UNDEF_UINT;

  huff->totalnumofchars = 0;
  huff->totalnumofbits = 0;
  huffman_tree_set_codes_rec(huff->roothuffmantree);
  gt_huffman_iterate(huff, calc_size, huff);
  gt_huffman_iterate(huff, store_codes, huff);

  return huff;
}

void gt_huffman_delete(GtHuffman *huffman)
{
  if (huffman) {
    gt_rbtree_delete(huffman->rbt_root);
    gt_free(huffman->code_tab);
  }
  gt_free(huffman);
}

int gt_huffman_iterate(const GtHuffman *huffman,
                       GtHuffmanActFunc actfun,
                       void *actinfo)
{
  gt_assert(huffman);
  gt_assert(actfun);

  if (huffman->roothuffmantree != NULL) {
    return visit_huffman_leaves_rec(huffman->roothuffmantree,
                                    actinfo,
                                    actfun);
  }
  return 0;
}

void gt_huffman_size(const GtHuffman *huffman, uint64_t *bits,
                     unsigned long *chars)
{
  gt_assert(huffman);
  if (bits != NULL)
    *bits = huffman->totalnumofbits;
  if (chars != NULL)
    *chars = huffman->totalnumofchars;
}

void gt_huffman_print_codes(const GtHuffman *huffman)
{
  gt_assert(huffman);
  gt_huffman_iterate(huffman, print_codes, NULL);
}

void gt_huffman_encode(const GtHuffman *huffman,
                       unsigned long symbol,
                       GtBitsequence *code,
                       unsigned long *codelength)
{
  gt_assert(huffman && symbol < huffman->totalnumofsymbols);
  *code = huffman->code_tab[symbol].code;
  *codelength = huffman->code_tab[symbol].numofbits;
}

unsigned long gt_huffman_numofsymbols(const GtHuffman *huffman)
{
  gt_assert(huffman);
  return huffman->numofsymbols;
}

unsigned long gt_huffman_totalnumofsymbols(const GtHuffman *huffman)
{
  gt_assert(huffman);
  return huffman->totalnumofsymbols;
}

GtHuffmanDecoder *gt_huffman_decoder_new(GtHuffman *huffman,
                                         GtBitsequence *bitsequence,
                                         unsigned long length,
                                         unsigned long bit_offset,
                                         unsigned long pad_length)
{
  GtHuffmanDecoder *hd = gt_malloc(sizeof (*hd));

  hd->huffman = huffman;
  hd->cur_node = hd->huffman->roothuffmantree;
  hd->bitsequence = bitsequence;
  hd->cur_bit = bit_offset;
  hd->pad_length = pad_length;
  hd->length = length;
  hd->mem_func = NULL;
  hd->info = NULL;
  hd->cur_bitseq = 0;
  hd->mem_func_stat = 0;
  return hd;
}

GtHuffmanDecoder *gt_huffman_decoder_new_from_memory(
                                     GtHuffman *huffman,
                                     GtHuffmanDecoderGetMemFunc mem_func,
                                     void *info,
                                     GtError *err)
{
  GtHuffmanDecoder *hd = gt_malloc(sizeof (*hd));

  hd->huffman = huffman;
  hd->cur_node = hd->huffman->roothuffmantree;
  hd->mem_func = mem_func;
  hd->info = info;
  hd->cur_bitseq = 0;
  hd->mem_func_stat = hd->mem_func(&hd->bitsequence,
                                   &hd->length,
                                   &hd->cur_bit,
                                   &hd->pad_length,
                                   hd->info);
  if (hd->mem_func_stat == -1) {
    gt_error_set(err, "error calling mem_func");
    return NULL;
  }
  return hd;
}

int gt_huffman_decoder_get_new_mem_chunk(GtHuffmanDecoder *hd, GtError *err)
{
  gt_assert(hd->mem_func);

  hd->mem_func_stat = hd->mem_func(&hd->bitsequence,
                                   &hd->length,
                                   &hd->cur_bit,
                                   &hd->pad_length,
                                   hd->info);
  if (hd->mem_func_stat == -1) {
    gt_error_set(err, "error calling mem_func");
    return hd->mem_func_stat;
  }
  return 0;
}

int gt_huffman_decoder_next(GtHuffmanDecoder *hd,
                            GtArray *symbols,
                            unsigned long symbols2read,
                            GtError *err)
{
  unsigned long read_symbols = 0,
                bits_to_read = GT_INTWORDSIZE;
  int had_err = 0;

  gt_assert((symbols2read > 0) && hd &&
            (gt_array_elem_size(symbols) == sizeof (unsigned long)));

  if (hd->cur_bitseq == hd->length - 1)
    bits_to_read = GT_INTWORDSIZE - hd->pad_length;

  /* reached end of data in a previous call, didn't get reset */
  if (hd->cur_bitseq == hd->length) {
    if (!hd->mem_func)
      return 0;  /*EOF*/
    else if (hd->mem_func_stat == 0)
      return 0;
    else
      return -1;
  }

  while (!had_err && read_symbols < symbols2read) {

    /* symbol found, reset huffman */
    if (!had_err && hd->cur_node->leftchild == NULL) {
      gt_array_add(symbols, hd->cur_node->symbol.symbol);
      read_symbols++;
      hd->cur_node = hd->huffman->roothuffmantree;
    }

    if (hd->cur_bit == bits_to_read) {
      hd->cur_bitseq++;

      if (hd->cur_bitseq == hd->length - 1)
        bits_to_read = GT_INTWORDSIZE - hd->pad_length;

      /* no next bitseq ask for more if possible */
      if (hd->cur_bitseq == hd->length ) {
        if (!hd->mem_func) {
          return 0;  /*EOF*/
        }

        hd->mem_func_stat = hd->mem_func(&hd->bitsequence,
                                         &hd->length,
                                         &hd->cur_bit,
                                         &hd->pad_length,
                                         hd->info);
        if (hd->mem_func_stat == -1) {
          gt_error_set(err, "error calling mem_func");
          had_err = hd->mem_func_stat;
        }
        else if (hd->mem_func_stat == 0) {
          return hd->mem_func_stat; /*EOF*/
        }

        hd->cur_bitseq = 0;
        bits_to_read = GT_INTWORDSIZE;
        if (hd->cur_bitseq == hd->length - 1)
          bits_to_read -= hd->pad_length;
      }
      hd->cur_bit = 0;
    }

    /* read the next bit */
    else if (!had_err) {
      if (GT_ISBITSET(hd->bitsequence[hd->cur_bitseq], hd->cur_bit++))
        hd->cur_node = hd->cur_node->rightchild;
      else
        hd->cur_node = hd->cur_node->leftchild;
    }
  }
  if (!had_err) {
    /* all read, but not eof */
    return 1;
  }
  return had_err;
}

void gt_huffman_decoder_delete(GtHuffmanDecoder *hd)
{
  gt_free(hd);
}

GtHuffmanBitwiseDecoder *gt_huffman_bitwise_decoder_new(GtHuffman *huffman,
                                                        GT_UNUSED GtError *err)
{
  GtHuffmanBitwiseDecoder *hbwd;
  gt_assert(huffman && huffman->roothuffmantree);
  hbwd = gt_calloc(1, sizeof (*hbwd));
  hbwd->huffman = huffman;
  hbwd->cur_node = huffman->roothuffmantree;
  return hbwd;
}

int gt_huffman_bitwise_decoder_next(GtHuffmanBitwiseDecoder *hbwd,
                                    bool bit, unsigned long *symbol)
{
  gt_assert(hbwd);

  if (hbwd->cur_node->leftchild == NULL) {
    *symbol = hbwd->cur_node->symbol.symbol;
    hbwd->cur_node = hbwd->huffman->roothuffmantree;
    return 0;
  }
  else {
    if (bit)
      hbwd->cur_node = hbwd->cur_node->rightchild;
    else
      hbwd->cur_node = hbwd->cur_node->leftchild;

    if (hbwd->cur_node->leftchild == NULL) {
      *symbol = hbwd->cur_node->symbol.symbol;
      hbwd->cur_node = hbwd->huffman->roothuffmantree;
      return 0;
    }
    else
      return 1;
  }
}

void gt_huffman_bitwise_decoder_delete(GtHuffmanBitwiseDecoder *hbwd)
{
  if (!hbwd)
    return;
  gt_free(hbwd);
}

unsigned long long unit_test_distr_func(const void *distr, unsigned long symbol)
{
  unsigned long long *distrull = (unsigned long long*) distr;
  return distrull[symbol];
}

static int test_bitwise(GtError *err)
{
  int had_err = 0,
      stat;
  unsigned long i,j,
                symbol,
                length1 = 0,
                length2 = 0;
  GtHuffman *huffman;
  GtHuffmanBitwiseDecoder *hbwd;
  GtBitsequence bitseq;
  unsigned long long distr[6] = {45, 16, 13, 12, 9, 5};

  huffman = gt_huffman_new(&distr,unit_test_distr_func, 6);

  hbwd = gt_huffman_bitwise_decoder_new(huffman, err);
  gt_ensure(had_err, hbwd);
  for (i = 0; !had_err && i < 6; i++) {
    gt_huffman_encode(huffman, i, &bitseq, &length2);
    if (i > 0) {
      gt_ensure(had_err, length1 <= length2);
    }
    j = 1;
    while (1 == (stat = gt_huffman_bitwise_decoder_next(hbwd,
                                             (bitseq >> (length2 - j)) & 1,
                                             &symbol))) {
      j++;
    }
    gt_ensure(had_err, stat == 0);
    gt_ensure(had_err, j == length2);
    gt_ensure(had_err, symbol == i);
    length1 = length2;
  }

  stat = gt_huffman_bitwise_decoder_next(hbwd, true, &symbol);

  gt_ensure(had_err, stat == 0);
  gt_ensure(had_err, symbol == 0);

  stat = gt_huffman_bitwise_decoder_next(hbwd, false, &symbol);
  gt_ensure(had_err, stat == 1);
  stat = gt_huffman_bitwise_decoder_next(hbwd, false, &symbol);
  gt_ensure(had_err, stat == 1);
  stat = gt_huffman_bitwise_decoder_next(hbwd, false, &symbol);
  gt_ensure(had_err, stat == 0);
  gt_ensure(had_err, symbol == 1);

  stat = gt_huffman_bitwise_decoder_next(hbwd, false, &symbol);
  gt_ensure(had_err, stat == 1);
  stat = gt_huffman_bitwise_decoder_next(hbwd, true, &symbol);
  gt_ensure(had_err, stat == 1);
  stat = gt_huffman_bitwise_decoder_next(hbwd, false, &symbol);
  gt_ensure(had_err, stat == 0);
  gt_ensure(had_err, symbol == 2);

  stat = gt_huffman_bitwise_decoder_next(hbwd, false, &symbol);
  gt_ensure(had_err, stat == 1);
  stat = gt_huffman_bitwise_decoder_next(hbwd, true, &symbol);
  gt_ensure(had_err, stat == 1);
  stat = gt_huffman_bitwise_decoder_next(hbwd, true, &symbol);
  gt_ensure(had_err, stat == 0);
  gt_ensure(had_err, symbol == 3);

  stat = gt_huffman_bitwise_decoder_next(hbwd, false, &symbol);
  gt_ensure(had_err, stat == 1);
  stat = gt_huffman_bitwise_decoder_next(hbwd, false, &symbol);
  gt_ensure(had_err, stat == 1);
  stat = gt_huffman_bitwise_decoder_next(hbwd, true, &symbol);
  gt_ensure(had_err, stat == 1);
  stat = gt_huffman_bitwise_decoder_next(hbwd, false, &symbol);
  gt_ensure(had_err, stat == 0);
  gt_ensure(had_err, symbol == 4);

  stat = gt_huffman_bitwise_decoder_next(hbwd, false, &symbol);
  gt_ensure(had_err, stat == 1);
  stat = gt_huffman_bitwise_decoder_next(hbwd, false, &symbol);
  gt_ensure(had_err, stat == 1);
  stat = gt_huffman_bitwise_decoder_next(hbwd, true, &symbol);
  gt_ensure(had_err, stat == 1);
  stat = gt_huffman_bitwise_decoder_next(hbwd, true, &symbol);
  gt_ensure(had_err, stat == 0);
  gt_ensure(had_err, symbol == 5);

  gt_huffman_delete(huffman);
  gt_huffman_bitwise_decoder_delete(hbwd);
  return had_err;
}

typedef struct huffman_unit_test_meminfo {
  unsigned long chunk, chunks, size, lastchunk_size, padding;
  GtBitsequence *data;
} HuffmanUnitTestMeminfo;

static int huffman_unit_get_next_block(GtBitsequence **bits,
                                       unsigned long *length,
                                       unsigned long *offset,
                                       unsigned long *pad_length,
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
  int had_err = 0,
      range = 0,
      decoder_stat = 1;
  unsigned char bits = sizeof (GtBitsequence) * 8,
                bits_remain = bits;
  int const max_num = 1000 * gt_rand_0_to_1() + 100,
            dist_size = 10 + 10 * gt_rand_0_to_1(),
            step_size = 10 + 5 * gt_rand_0_to_1();
  unsigned long idx, idx_j,
                code_len;
  uint64_t total_bits;
  GtArray *numbers = gt_array_new(sizeof (unsigned char));
  GtArray *encoded = gt_array_new(sizeof (unsigned long));
  GtArray *codes = gt_array_new(sizeof (GtBitsequence));
  GtBitsequence buffer = 0, code;
  unsigned long long *distribution = gt_malloc(sizeof (*distribution) *
                                     dist_size);
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
      x = (dist_size >> 1) * gt_rand_0_to_1();
      range = 0;
    }
    else {
      x = dist_size * gt_rand_0_to_1();
      range = 1;
    }
    gt_array_add(numbers, x);
  }

  gt_ensure(had_err, gt_array_size(numbers) == max_num);

  /*calculate distribution*/
  for (idx = 0; !had_err && idx < max_num; idx++) {
    distribution[*(unsigned char*) gt_array_get(numbers, idx)]++;
  }

  if (!had_err) {
    unsigned long chars;
    huff = gt_huffman_new(distribution, unit_test_distr_func, dist_size);
    gt_ensure(had_err, huff);
    gt_huffman_size(huff, &total_bits, &chars);
    gt_ensure(had_err, chars == max_num);
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

  gt_ensure(had_err,
            total_bits == gt_array_size(codes) *
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
                                     0,/*offset*/
                                     bits_remain);
    gt_ensure(had_err, huffdec);
  }

  for (idx = 0; !had_err && idx < max_num; idx += step_size) {
    gt_array_reset(encoded);
    gt_ensure(had_err, decoder_stat == 1);
    decoder_stat = gt_huffman_decoder_next(huffdec,
                                           encoded,
                                           step_size,
                                           err);
    if (decoder_stat == 1)
      gt_ensure(had_err, gt_array_size(encoded) == step_size);
    else
      gt_ensure(had_err, decoder_stat == 0);

    gt_ensure(had_err, decoder_stat != -1);
    if (!had_err) {
      for (idx_j = 0; !had_err && idx_j < gt_array_size(encoded); idx_j++) {
        unsigned long ist = *(unsigned long*)gt_array_get(encoded, idx_j);
        unsigned char soll =
          *(unsigned char*)gt_array_get(numbers, idx + idx_j);
        gt_ensure(had_err, (ist == soll));
      }
    }
  }

  gt_ensure(had_err, idx >= max_num);
  if (!had_err)
    decoder_stat = gt_huffman_decoder_next(huffdec,
                                           encoded,
                                           1,
                                           err);
  gt_ensure(had_err, decoder_stat == 0);

  if (!had_err) {
    gt_huffman_decoder_delete(huffdec);

    huffdec = gt_huffman_decoder_new_from_memory(huff,
                                                 huffman_unit_get_next_block,
                                                 meminfo,
                                                 err);
    gt_ensure(had_err, huffdec);
  }

  decoder_stat = 1;
  for (idx = 0; !had_err && idx < max_num; idx += step_size) {
    gt_array_reset(encoded);
    gt_ensure(had_err, decoder_stat == 1);
    decoder_stat = gt_huffman_decoder_next(huffdec,
                                           encoded,
                                           step_size,
                                           err);
    if (decoder_stat == 1)
      gt_ensure(had_err, gt_array_size(encoded) == step_size);
    else
      gt_ensure(had_err, decoder_stat == 0);

    gt_ensure(had_err, decoder_stat != -1);
    if (!had_err) {
      for (idx_j = 0; !had_err && idx_j < gt_array_size(encoded); idx_j++) {
        unsigned long ist = *(unsigned long*)gt_array_get(encoded, idx_j);
        unsigned char soll =
          *(unsigned char*)gt_array_get(numbers, idx + idx_j);
        gt_ensure(had_err, (ist == soll));
      }
    }
  }

  gt_ensure(had_err, idx >= max_num);
  if (!had_err)
    decoder_stat = gt_huffman_decoder_next(huffdec,
                                           encoded,
                                           1,
                                           err);
  gt_ensure(had_err, decoder_stat == 0);

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
    test_mem(err);

  return had_err;
}
