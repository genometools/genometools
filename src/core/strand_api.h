/*
  Copyright (c) 2006-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2008 Center for Bioinformatics, University of Hamburg

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

#ifndef STRAND_API_H
#define STRAND_API_H

typedef enum {
  GT_STRAND_FORWARD, /* '+' */
  GT_STRAND_REVERSE, /* '-' */
  GT_STRAND_BOTH,    /* '.' */
  GT_STRAND_UNKNOWN, /* '?' */
  GT_NUM_OF_STRAND_TYPES
} GtStrand;

#define GT_STRAND_CHARS "+-.?"

/* Returns NUM_OF_GT_STRAND_TYPES if <strand_char> is not a valid one. */
GtStrand gt_strand_get(char strand_char);
GtStrand gt_strand_join(GtStrand, GtStrand);

#endif
