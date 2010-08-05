/*
  Copyright (c) 2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2008 Center for Bioinformatics, University of Hamburg

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

#ifndef GFF3_ESCAPING_H
#define GFF3_ESCAPING_H

#include "core/str.h"

/* Escape <unescaped_seq> of given <length> for GFF3 format and append the
   result to <escaped_seq>. */
void gt_gff3_escape(GtStr *escaped_seq, const char *unescaped_seq,
                    unsigned long length);

/* Unescape GFF3 format <escaped_seq> of given <length> and append the result to
   <unescaped_seq>. */
int  gt_gff3_unescape(GtStr *unescaped_seq, const char *escaped_seq,
                      unsigned long length, GtError*);

/* Perform unit test of GFF3 format escaping module */
int  gt_gff3_escaping_unit_test(GtError*);

#endif
