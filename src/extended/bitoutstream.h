/*
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

#ifndef BITOUTSTREAM_H
#define BITOUTSTREAM_H

#include "core/error.h"
#include "core/intbits.h"
#include "core/bittab_api.h"

/* Class <GtBitOutStream> helps writing variable length encoded data to files,
   handling the filling of words of type <GtBitsequence>. */
typedef struct GtBitOutStream GtBitOutStream;

/* Returns a new <GtBitOutStream>, <fp> needs to be valid and opened for
   writing. */
GtBitOutStream* gt_bitoutstream_new(FILE *fp);

/* Append the bitcode <code> to the file associated with <bitstream>.
   <bits_to_write> is the number of bits in <code> that have to be appended.
   Assumes the bits are stored in the least significant bits of <code> like
   standard binary encoding. */
void            gt_bitoutstream_append(GtBitOutStream *bitstream,
                                       GtBitsequence code,
                                       unsigned bits_to_write);

/* Append the bits in <GtBittab> <tab> to the file associated with
   <bitstream>. */
void            gt_bitoutstream_append_bittab(GtBitOutStream *bitstream,
                                              GtBittab *tab);

/* Write all currently appended bitcodes to the file associated with
   <bitstream>. Possibly 'empty' bits in the current word will be set to zero
   and all non empty bits will be shifted to the most significant bits. */
void            gt_bitoutstream_flush(GtBitOutStream *bitstream);

/* Like <gt_bitoutstream_flush()> but moves file pointer to the next start of a
   page. */
void            gt_bitoutstream_flush_advance(GtBitOutStream *bitstream);

/* Returns the position in the file associated with <bitstream>. For reliable
   results gt_bitoutstream_flush needs to be called before. Returns -1 on
   error. */
long            gt_bitoutstream_pos(const GtBitOutStream *bitstream);

void            gt_bitoutstream_delete(GtBitOutStream *bitstream);

#endif
