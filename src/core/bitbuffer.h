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

#include <stdio.h>

/* The <GtBitbuffer> class provides means to sequentially write bit-compressed
   integer arrays into a file. */
typedef struct GtBitbuffer GtBitbuffer;

/* Creates a new <GtBitbuffer> for output to <outfp>.
   <bitsperentry> specifies the number of bits per entry. */
GtBitbuffer* gt_bitbuffer_new(FILE *outfp, uint8_t bitsperentry);

/* Appends unsigned long <value> to <bb>. */
void         gt_bitbuffer_next_value (GtBitbuffer *bb, unsigned long value);

/* Appends unsigned 32-bit integer array <tab> of length <len> to <bb>. */
void         gt_bitbuffer_next_uint32tab(GtBitbuffer *bb, const uint32_t *tab,
                                         unsigned long len);

/* Appends unsigned long integer array <tab> of length <len> to <bb>. */
void         gt_bitbuffer_next_ulongtab(GtBitbuffer *bb,
                                        const unsigned long *tab,
                                        unsigned long len);

/* Deletes <bb> and frees all associated memory. */
void         gt_bitbuffer_delete(GtBitbuffer *bb);

#endif
