/*
  Copyright (c) 2011 Joachim Bonnet <joachim.bonnet@studium.uni-hamburg.de>
  Copyright (c) 2012 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
  Copyright (c) 2010-2012 Center for Bioinformatics, University of Hamburg

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

#ifndef HUFFCODE_H
#define HUFFCODE_H

#include "core/array_api.h"
#include "core/intbits.h"
#include "core/error_api.h"

/* Class <GtHuffman> holds information to encode arbitrary symbols with an
   efficient variable bit length. */
typedef struct GtHuffman GtHuffman;

/* Class <GtHuffmanDecoder> is used to decode consecutive symbols from Huffman
   encoded data */
typedef struct GtHuffmanDecoder GtHuffmanDecoder;

/* Class <GtHuffmanBitwiseDecoder> is used to decode bitstreams bit by bit */
typedef struct GtHuffmanBitwiseDecoder GtHuffmanBitwiseDecoder;

/* Function type, used to get the count of a specific symbol represented by
   <symbol_number> in <distribution>. Needed by <GtHuffman>. */
typedef unsigned long long (*GtDistrFunc) (const void *distribution,
                                           unsigned long symbol_number);

/* Function type, used when iterating over all encoded symbols in a <GtHuffman>
   object. <symbol> is the number of the current symbol, <freq> its frequency in
   the distribution of symbols, <code> points to the bitstring and <code_length>
   is its corresponding length.
   <action_info> points to a structure that might contain information for the
   <GtHuffmanActFunc> for example when calculating the average code length.
   Should return 0 on success and -1 on error. */
typedef int (*GtHuffmanActFunc) (unsigned long symbol,
                                 unsigned long long freq,
                                 GtBitsequence code,
                                 unsigned int code_length,
                                 void *action_info);

/* Returns a new <GtHuffman> object. <distribution> points to a data
   structure containing the distribution of a set of symbols. <distr_func> is a
   function returning the number of occurrences for a given symbol. Symbols are
   represented by numbers in range [0,num_of_symbols - 1]. */
GtHuffman*               gt_huffman_new(const void *distribution,
                                        GtDistrFunc distr_func,
                                        unsigned long num_of_symbols);

/* Traverses the Huffman tree and calculates total number of symbols in text and
   total number of bits needed to represent the text. */
void                     gt_huffman_size(const GtHuffman *huffman,
                                         uint64_t *bits,
                                         uint64_t *chars);

/* Prints the code for each symbol to stdout. One symbol per line. */
void                     gt_huffman_print_codes(const GtHuffman *huffman);

/* Iterate over all encoded symbols and perform <action_func> with its
   statistics. <action_info> can be null or any struct to hold information for
   the action_func*/
int                      gt_huffman_iterate(const GtHuffman *huffman,
                                            GtHuffmanActFunc action_func,
                                            void *action_info);

/* Encodes a single symbol and writes its code to <code> and the number of bits
   for its code to <codelength>. The symbol has to be in the range of the
   symbols with which the Huffman-object was initiated. If the symbol was not
   present in the distribution its code will be 0. */
void                     gt_huffman_encode(const GtHuffman *huffman,
                                           unsigned long symbol,
                                           GtBitsequence *code,
                                           unsigned int *codelength);

/* Returns the number of symbols with frequency > 0. */
unsigned long            gt_huffman_numofsymbols(const GtHuffman *huffman);

/* Returns the number of symbols with frequency >= 0. */
unsigned long            gt_huffman_totalnumofsymbols(const GtHuffman *huffman);

/* Returns a new <GtHuffmanDecoder> object. This decoder is meant to decode a
   Huffman encoded bitstring. <length> is the number of elements in
   <bitsequence>. The <bit_offset> tells the decoder at which bit of the first
   <GtBitsequence> to start, <pad_length> is the number of bits in the last
   element of <bitsequence> that are not part of the encoded data. */
GtHuffmanDecoder*        gt_huffman_decoder_new(GtHuffman *huffman,
                                                GtBitsequence *bitsequence,
                                                unsigned long length,
                                                unsigned long bit_offset,
                                                unsigned long pad_length);

/* Function type, used by <GtHuffmanDecoder> to ask for a pointer to
   <GtBitsequence> to read. Sets <bitsequence> to the new memory, <length> to
   the length of the <GtBitsequence> array, <bit_offset> to the offset in the
   first element of <bitsequence> and <pad_length> to the number of bits in the
   last element of <bitsequence> not part of the encoded data. <mem_info>
   contains info for the callback function.  Returns -1 on error, 1 if
   successfully provided new data, and 0 if all data has been provided */
typedef int (*GtHuffmanDecoderGetMemFunc) (GtBitsequence **bitsequence,
                                           unsigned long *length,
                                           unsigned long *bit_offset,
                                           unsigned long *pad_length,
                                           void *info);

/* <mem_func> and <info> are needed if the decoder reaches the end of the
   <bitsequence> but needs more to read, for example if the <bitsequence> is
   part of mmaped space. <mem_func> will be called initially to get the first
   chunk of data.
   This type of encoder is helpful if the encoded data is sampled or read in
   chunks from mmapped files. */
GtHuffmanDecoder*        gt_huffman_decoder_new_from_memory(
                                            GtHuffman *huffman,
                                            GtHuffmanDecoderGetMemFunc mem_func,
                                            void *info,
                                            GtError *err);

/* Make the decoder <huff_decoder> call the <mem_func> to get new memory.
   Returns -1 on error. This can be used to reset the decoder. */
int                      gt_huffman_decoder_get_new_mem_chunk(
                                                 GtHuffmanDecoder *huff_decoder,
                                                 GtError *err);

/* Decodes at most the next <num_of_symbols> symbols and adds them to
   <symbols> and returns the decoder's status. Returns 0 if EOF was reached, 1
   if there are more symbols to read and -1 if an error occurred. <symbols> must
   be an array with type unsigned long. */
int                      gt_huffman_decoder_next(GtHuffmanDecoder *huff_decoder,
                                                 GtArray *symbols,
                                                 unsigned long symbols_to_read,
                                                 GtError *err);

/* Returns a new GtHuffmanBitwiseDecoder object. This decoder that is
   meant to decode a code bit by bit.*/
GtHuffmanBitwiseDecoder* gt_huffman_bitwise_decoder_new(GtHuffman *huffman,
                                                        GtError *err);

/* Moves inside the Huffman tree exactly by one node according to <bit>
  (true: visit right child, false: visit left child). Returns 0 if a leaf was
  reached and then writes the corresponding symbol to <symbol>. Returns 1 if no
  leaf was reached yet.*/
int                      gt_huffman_bitwise_decoder_next(
                                                  GtHuffmanBitwiseDecoder *hbwd,
                                                  bool bit,
                                                  unsigned long *symbol,
                                                  GtError *err);

/* Deletes <huff_decoder>. */
void                     gt_huffman_decoder_delete(
                                                GtHuffmanDecoder *huff_decoder);

/* Deletes <huffman>. */
void                     gt_huffman_delete(GtHuffman *huffman);

/* Deletes <hbwd>. */
void                     gt_huffman_bitwise_decoder_delete(
                                                GtHuffmanBitwiseDecoder *hbwd);

int                      gt_huffman_unit_test(GtError *err);

#endif
