/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

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

#ifndef GREEDYEDIST_H
#define GREEDYEDIST_H

#include "core/types_api.h"
#include "core/encseq_api.h"

typedef struct GtGreedyedistSeq GtGreedyedistSeq;

GtGreedyedistSeq *gt_greedyedist_seq_new_empty(void);

void gt_greedyedist_seq_reinit_ptr(GtGreedyedistSeq *greedyedistseq,
                                   const GtUchar *ptr,
                                   unsigned long len,
                                   unsigned long offset);

GtGreedyedistSeq *gt_greedyedist_seq_new_ptr(const GtUchar *ptr,
                                             unsigned long len,
                                             unsigned long offset);

void gt_greedyedist_seq_reinit_encseq(GtGreedyedistSeq *greedyedistseq,
                                      const GtEncseq *encseq,
                                      unsigned long len,
                                      unsigned long offset);

GtGreedyedistSeq *gt_greedyedist_seq_new_encseq(const GtEncseq *encseq,
                                                unsigned long len,
                                                unsigned long offset);

void gt_greedyedist_seq_delete(GtGreedyedistSeq *greedyedistseq);

unsigned long gt_greedyedist_length_get(const GtGreedyedistSeq *greedyedistseq);

unsigned long greedyunitedist(const GtGreedyedistSeq *useq,
                              const GtGreedyedistSeq *vseq);

#endif
