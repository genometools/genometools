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

#ifndef COMPRESSED_BITSEQUENCE_H
#define COMPRESSED_BITSEQUENCE_H

#include "core/error_api.h"
#include "core/intbits.h"

/* The <GtCompressedBitsequence> class stores a bitvector in a compressed way
   known as an RRR-bitvector like Raman, Raman and Rao described it in 2002. It
   gives constant time access and rank on the bitvector represented. */
typedef struct GtCompressedBitsequence GtCompressedBitsequence;

/* Returns a new <GtCompressedBitsequence> object. <bitseq> points to the bit
   sequence to be compressed, <samplerate> defines the rate of sampling, which
   is the constant factor for access and rank queries, <num_of_bits> is the
   number of bits to be read from <bitseq>. Words in <bitseq> are assumed to be
   filled continuously. That is, the last <GtBitsequence> in <bitseq> has its
   most significant bits filled first. */
GtCompressedBitsequence* gt_compressed_bitsequence_new(
                                                     GtBitsequence *bitseq,
                                                     unsigned int samplerate,
                                                     unsigned long num_of_bits);

/* Returns 0 or 1 according to the bit at <position> in <cbs>. Note that
   <position> has to be smaller than the length  of <cbs>. */
int                      gt_compressed_bitsequence_access(
                                                  GtCompressedBitsequence *cbs,
                                                  unsigned long position);

/* Returns the number of 1 bits in <cbs> upto and including <position>. Note
   that <position> has to be smaller than the length of <cbs>. */
unsigned long            gt_compressed_bitsequence_rank_1(
                                                   GtCompressedBitsequence *cbs,
                                                   unsigned long position);

/* Returns the number of 0 bits in <cbs> upto and including <position>. Note
   that <position> has to be smaller than the length of <cbs>. */
unsigned long            gt_compressed_bitsequence_rank_0(
                                                   GtCompressedBitsequence *cbs,
                                                   unsigned long position);

/* Returns the position of the <num>th bit set to 1 in <cbs>. Returns length of
   <cbs> if there are less than <num> bits set to 1. */
unsigned long            gt_compressed_bitsequence_select_1(
                                                   GtCompressedBitsequence *cbs,
                                                   unsigned long num);

/* Returns the position of the <num>th bit set to 0 in <cbs>. Returns length of
   <cbs> if there are less than <num> bits set to 0. */
unsigned long            gt_compressed_bitsequence_select_0(
                                                   GtCompressedBitsequence *cbs,
                                                   unsigned long num);

/* Frees the memory of <cbs>. */
void                     gt_compressed_bitsequence_delete(
                                                  GtCompressedBitsequence *cbs);

int                      gt_compressed_bitsequence_unit_test(GtError *err);
#endif
