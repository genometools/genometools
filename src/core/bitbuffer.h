/*
  Copyright (c) 2013 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
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

#ifndef BITBUFFER_H
#define BITBUFFER_H

#include <inttypes.h>
#include <stdio.h>
#include "core/types_api.h"

/* The <GtBitbuffer> class provides means to sequentially write bit-compressed
   integer arrays into a file. */
typedef struct GtBitbuffer GtBitbuffer;

/* Creates a new <GtBitbuffer> for output to <outfp>.
   <bitsperentry> specifies the number of bits per entry if
   greater than 0. In this case a header is written in the input file
   consisting of a <uint64_t>-value specifying the number of elements
   written and a <uint8_t> value specifying the number of bits per
   entry. If <bitsperentry> is 0, then no header is written and
   the number of bits for the value to be written has to be specified
   for each call of <gt_bitbuffer_next_value()>. */
GtBitbuffer* gt_bitbuffer_new(FILE *outfp, unsigned int bitsperentry);

/* Appends GtUword <value> of <bitsforvalue> bits to <bb>. */
void         gt_bitbuffer_next_value(GtBitbuffer *bb, GtUword value,
                                     unsigned int bitsforvalue);

/* Appends GtUword <value> to <bb>. Requires that <bb> has been created
   with a <bitsperentry> value > 0. */
void         gt_bitbuffer_next_fixed_bits_value(GtBitbuffer *bb,
                                                GtUword value);

/* Appends unsigned 32-bit integer array <tab> of length <len> to <bb>. */
void         gt_bitbuffer_next_uint32tab(GtBitbuffer *bb, const uint32_t *tab,
                                         GtUword len);

/* Appends GtUword integer array <tab> of length <len> to <bb>. */
void         gt_bitbuffer_next_ulongtab(GtBitbuffer *bb,
                                        const GtUword *tab,
                                        GtUword len);

/* Deletes <bb> and frees all associated memory. */
void         gt_bitbuffer_delete(GtBitbuffer *bb);

#endif
