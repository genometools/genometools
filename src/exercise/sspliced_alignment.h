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

#ifndef SSPLICED_ALIGNMENT_H
#define SSPLICED_ALIGNMENT_H

#include <stdbool.h>
#include "core/range.h"

/* the simple spliced alignment class */
typedef struct GtSSplicedAlignment GtSSplicedAlignment;

GtSSplicedAlignment* gt_sspliced_alignment_new(const char *id, bool forward);
void               gt_sspliced_alignment_delete(GtSSplicedAlignment*);
bool               gt_sspliced_alignment_is_forward(const GtSSplicedAlignment*);
void               gt_sspliced_alignment_add_exon(GtSSplicedAlignment*,
                                                  GtRange);
unsigned long      gt_sspliced_alignment_num_of_exons(const
                                                      GtSSplicedAlignment*);
GtRange            gt_sspliced_alignment_get_exon(const GtSSplicedAlignment*,
                                                  unsigned long exon_number);
GtRange            gt_sspliced_alignment_genomic_range(const
                                                       GtSSplicedAlignment*);
int                gt_sspliced_alignment_compare_ptr(const
                                                     GtSSplicedAlignment**,
                                                     const
                                                     GtSSplicedAlignment**);

#endif
