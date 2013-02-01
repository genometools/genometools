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

#ifndef SPLICEDSEQ_H
#define SPLICEDSEQ_H

#include <stdbool.h>

typedef struct Splicedseq Splicedseq;

Splicedseq*   gt_splicedseq_new(void);
/* adds an ``exon'' to the spliced sequence */
void          gt_splicedseq_add(Splicedseq*, unsigned long start,
                                unsigned long end,
                                const char *original_sequence);
char*         gt_splicedseq_get(const Splicedseq*);
bool          gt_splicedseq_pos_is_border(const Splicedseq*, unsigned long);
/* maps the given position back to the original coordinate system */
unsigned long gt_splicedseq_map(const Splicedseq*, unsigned long);
unsigned long gt_splicedseq_length(const Splicedseq*);
int           gt_splicedseq_reverse(Splicedseq*, GtError*);
void          gt_splicedseq_reset(Splicedseq*);
int           gt_splicedseq_unit_test(GtError*);
void          gt_splicedseq_delete(Splicedseq*);

#endif
