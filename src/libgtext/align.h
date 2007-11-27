/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg

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

#ifndef ALIGN_H
#define ALIGN_H

#include "libgtext/alignment.h"

/* (globally) align <u> and <v> (unit cost) and return one optimal Alignment */
Alignment* align(const char *u, unsigned long ulen,
                 const char *v, unsigned long vlen);

/* align <u> and <v> (unit cost), call proc_alignment for each optimal
   Alignment, and call proc_aligns with the number of optimal alignments*/
void align_all(const char *u, unsigned long ulen,
               const char *v, unsigned long vlen,
               void (*proc_alignment)(const Alignment*, void *data),
               void (*proc_aligns)(unsigned long, void *data), void *data);

#endif
