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
#include "libgtcore/range.h"

/* the simple spliced alignment class */
typedef struct SSplicedAlignment SSplicedAlignment;

SSplicedAlignment* sspliced_alignment_new(const char *id, bool forward);
void               sspliced_alignment_delete(SSplicedAlignment*);
bool               sspliced_alignment_is_forward(const SSplicedAlignment*);
void               sspliced_alignment_add_exon(SSplicedAlignment*, Range);
unsigned long      sspliced_alignment_num_of_exons(const SSplicedAlignment*);
Range              sspliced_alignment_get_exon(const SSplicedAlignment*,
                                               unsigned long exon_number);
Range              sspliced_alignment_genomic_range(const SSplicedAlignment*);
int                sspliced_alignment_compare_ptr(const SSplicedAlignment**,
                                                  const SSplicedAlignment**);

#endif
