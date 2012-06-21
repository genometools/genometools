/*
  Copyright (c) 2011 Joachim Bonnet <joachim.bonnet@studium.uni-hamburg.de>
  Copyright (c) 2012 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
  Copyright (c) 2003-2012 Center for Bioinformatics, University of Hamburg

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

#ifndef ELIAS_GAMMA_H
#define ELIAS_GAMMA_H

#include "core/bittab_api.h"
#include "core/error_api.h"

/* The <GtEliasGammaBitwiseDecoder> class is used to decode Elias gamma encoded
   integers. For details see Elias, Peter: "Universal codeword sets and
   representations of the integers", IEEE Transactions on Information Theory,
   21(2):194--203 (1975). */
typedef struct GtEliasGammaBitwiseDecoder GtEliasGammaBitwiseDecoder;

/* Encodes the given positive integer <x>, writes it to a new <GtBittab> object
   and returns it. */
GtBittab*                   gt_elias_gamma_encode(unsigned long x);

/* Returns a new <GtEliasGammaBitwiseDecoder> object. This decoder is
   meant to decode a code word bit by bit.*/
GtEliasGammaBitwiseDecoder* gt_elias_gamma_bitwise_decoder_new(void);

/* Appends <bit> to the code word currently being decoded by <egbd>.
   Returns 0 if the code word was completely read and writes its value to the
   location pointed to by <x>. Returns 1 if more bits are needed to decode the
   current code word. */
int                         gt_elias_gamma_bitwise_decoder_next(
                                               GtEliasGammaBitwiseDecoder *egbd,
                                               bool bit, unsigned long *x);

/* Deletes <egbd> and frees all associated memory. */
void                        gt_elias_gamma_bitwise_decoder_delete(
                                              GtEliasGammaBitwiseDecoder *egbd);

int                         gt_elias_gamma_unit_test(GtError *err);

#endif
