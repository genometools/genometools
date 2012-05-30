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

#ifndef ALIGNMENT_H
#define ALIGNMENT_H

#include <stdio.h>
#include "core/range.h"
#include "core/types_api.h"
#include "core/error_api.h"

/* the GtAlignment class (an object has to be constructed backwards) */
typedef struct GtAlignment GtAlignment;

GtAlignment*  gt_alignment_new(void);
GtAlignment*  gt_alignment_new_with_seqs(const GtUchar *u, unsigned long ulen,
                                         const GtUchar *v, unsigned long vlen);
void          gt_alignment_set_seqs(GtAlignment*, const GtUchar *u,
                                    unsigned long ulen,
                                    const GtUchar *v, unsigned long vlen);
GtRange       gt_alignment_get_urange(const GtAlignment*);
GtRange       gt_alignment_get_vrange(const GtAlignment*);
void          gt_alignment_set_urange(GtAlignment*, GtRange);
void          gt_alignment_set_vrange(GtAlignment*, GtRange);
void          gt_alignment_add_replacement(GtAlignment*);
void          gt_alignment_add_deletion(GtAlignment*);
void          gt_alignment_add_insertion(GtAlignment*);
unsigned long gt_alignment_get_length(const GtAlignment *a);
/* undo last add operation */
void          gt_alignment_remove_last(GtAlignment*);
/* reset list of edit operations to empty */
void          gt_alignment_reset(GtAlignment *a);
/* returns unit cost */
unsigned long gt_alignment_eval(const GtAlignment*);
long          gt_alignment_eval_with_score(const GtAlignment *a,
                                           long matchscore,
                                           long mismatchscore,
                                           long gapscore);
void          gt_alignment_show(const GtAlignment*, FILE*);
void          gt_alignment_show_with_mapped_chars(const GtAlignment*,
                                                  const GtUchar *characters,
                                                  GtUchar wildcardshow,
                                                  FILE *fp);
void          gt_alignment_show_multieop_list(const GtAlignment*, FILE*);
int           gt_alignment_unit_test(GtError*);
void          gt_alignment_delete(GtAlignment*);

#endif
