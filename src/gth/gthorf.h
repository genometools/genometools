/*
  Copyright (c) 2004-2010 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2004-2008 Center for Bioinformatics, University of Hamburg

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

#ifndef GTHORF_H
#define GTHORF_H

#include <stdbool.h>
#include "core/range_api.h"
#include "gth/gthoutput.h"
#include "gth/spliced_seq.h"

void gthshowORFs(char *frame0, char *frame1, char *frame2,
                 unsigned long frame0len, unsigned long frame1len,
                 unsigned long frame2len, bool gen_strand_forward,
                 unsigned long gen_total_length, unsigned long gen_offset,
                 const char *gen_id, unsigned long pglnum,
                 unsigned long agsnum, GthSplicedSeq *splicedseq,
                 unsigned int indentlevel, GthOutput *out);

#endif
