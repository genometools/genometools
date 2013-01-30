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

#ifndef BITINSTREAM_H
#define BITINSTREAM_H

#ifndef S_SPLINT_S
#include <sys/types.h>
#endif

#include "core/error.h"
#include "core/intbits.h"

/* Class <GtBitInStream> reads variable length encoded data from a mmap()ed
   file. */
typedef struct GtBitInStream GtBitInStream;

/* Returns a new <GtBitInStream>. <path> is a '\0' terminated string with the
   name of the file, <offset> is the page offset where mapping should start (a
   multiple of the page size). <pages_to_map> is the number of pages to map at
   a time. If <path> is not valid exits with <EXIT_FAILURE>. */
GtBitInStream* gt_bitinstream_new(const char *path,
                                  size_t offset,
                                  unsigned long pages_to_map);

/* Tells <bitstream> to remap the file with a new offset. */
void           gt_bitinstream_reinit(GtBitInStream *bitstream,
                                     size_t offset);

/* Reads one more bit and sets <bit> to the read value. Returns 0 if there are
   no more bits to read and 1 if successfully read one bit. */
int            gt_bitinstream_get_next_bit(GtBitInStream *bitstream,
                                           bool *bit);

/* Deletes <bitstream> and frees all associated memory. */
void           gt_bitinstream_delete(GtBitInStream *bitstream);

#endif
