/*
  Copyright (c) 2011 Joachim Bonnet <joachim.bonnet@studium.uni-hamburg.de>

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

#ifndef GOLOMB_H
#define GOLOMB_H

#include "core/bittab_api.h"
#include "core/error_api.h"

/* Class <GtGolomb> stores information to encode integers with golomb encoding
   See Golomb, S.W. (1966). , Run-length encodings. */
typedef struct GtGolomb GtGolomb;

/* Class <GtGolombBitwiseDecoder> takes a <GtGolomb> object and can be used to
   decode golomb encoded integers. */
typedef struct GtGolombBitwiseDecoder GtGolombBitwiseDecoder;

/* Returns a new <GtGolomb> object. The parameter <m> influences the length
   of generated code words. It has to be the same for encoding and decoding. */
GtGolomb*               gt_golomb_new(unsigned long m);

/* Encodes the non negative integer <x> for the given <gc>, writes it to a new
   <GtBittab> object and returns it. */
GtBittab*               gt_golomb_encode(GtGolomb *gc, unsigned long x);

/* Returns the parameter <m>, <gc> was initialized with. */
unsigned long           gt_golomb_get_m(GtGolomb *gc);

/* Returns a new <GtGolombBitwiseDecoder> object. This decoder is meant to
   decode a code word bit by bit.*/
GtGolombBitwiseDecoder* gt_golomb_bitwise_decoder_new(GtGolomb *gc);

/* Appends <bit> to the code word currently being decoded. Returns 0 if the
   code word was completely read and writes its decoding to the location
   pointed to by <x>. Returns 1 if more bits are needed to decode the current
   code word. */
int                     gt_golomb_bitwise_decoder_next(
                                                   GtGolombBitwiseDecoder *gcbd,
                                                   bool bit, unsigned long *x);

/* Deletes <gc>. */
void                    gt_golomb_delete(GtGolomb *gc);

/* Deletes <gcbd>. */
void                    gt_golomb_bitwise_decoder_delete(
                                                  GtGolombBitwiseDecoder *gcbd);

int                     gt_golomb_unit_test(GtError *err);

#endif
