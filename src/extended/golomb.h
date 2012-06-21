/*
  Copyright (c) 2011      Joachim Bonnet <joachim.bonnet@studium.uni-hamburg.de>
  Copyright (c)      2012 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
  Copyright (c) 2011-2012 Center for Bioinformatics, University of Hamburg

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

/* Class <GtGolomb> stores information to encode integers with Golomb encoding.
   See Golomb, S.W. (1966), Run-length encodings. */
typedef struct GtGolomb GtGolomb;

/* Class <GtGolombBitwiseDecoder> decode Golomb encoded integers given in a
   <GtGolomb> object. */
typedef struct GtGolombBitwiseDecoder GtGolombBitwiseDecoder;

/* Returns a new <GtGolomb> object. The parameter <median> influences the length
   of generated code words. Identical values are needed for correct encoding and
   decoding. */
GtGolomb*               gt_golomb_new(unsigned long median);

/* Encodes the non negative integer <x> for the given <golomb>, writes it to a
   new <GtBittab> object and returns it. */
GtBittab*               gt_golomb_encode(GtGolomb *golomb, unsigned long x);

/* Returns the parameter <median>, <golomb> was initialized with. */
unsigned long           gt_golomb_get_m(GtGolomb *golomb);

/* Returns a new <GtGolombBitwiseDecoder> object. This decoder is meant to
   decode a code word bit by bit.*/
GtGolombBitwiseDecoder* gt_golomb_bitwise_decoder_new(GtGolomb *golomb);

/* Appends <bit> to the code word currently being decoded. Returns 0 if the
   code word was completely read and writes its decoding to the location
   pointed to by <x>. Returns 1 if more bits are needed to decode the current
   code word. */
int                     gt_golomb_bitwise_decoder_next(
                                                   GtGolombBitwiseDecoder *gbwd,
                                                   bool bit, unsigned long *x);

/* Deletes <golomb>. */
void                    gt_golomb_delete(GtGolomb *golomb);

/* Deletes <gcbd>. */
void                    gt_golomb_bitwise_decoder_delete(
                                                  GtGolombBitwiseDecoder *gbwd);

int                     gt_golomb_unit_test(GtError *err);

#endif
