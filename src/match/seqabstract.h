/*
  Copyright (c) 2013 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2013 Center for Bioinformatics, University of Hamburg

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

#ifndef SEQABSTRACT_H
#define SEQABSTRACT_H

#include "core/types_api.h"
#include "core/encseq_api.h"

typedef struct GtSeqabstract GtSeqabstract;

GtSeqabstract* gt_seqabstract_new_empty(void);

GtSeqabstract* gt_seqabstract_new_ptr(const GtUchar *ptr,
                                      unsigned long len,
                                      unsigned long offset);

GtSeqabstract* gt_seqabstract_new_encseq(const GtEncseq *encseq,
                                         unsigned long len,
                                         unsigned long offset);

void           gt_seqabstract_reinit_ptr(GtSeqabstract *sa,
                                         const GtUchar *ptr,
                                         unsigned long len,
                                         unsigned long offset);

void           gt_seqabstract_reinit_encseq(GtSeqabstract *sa,
                                            const GtEncseq *encseq,
                                            unsigned long len,
                                            unsigned long offset);

unsigned long  gt_seqabstract_length_get(const GtSeqabstract *sa);

GtUchar        gt_seqabstract_encoded_char(const GtSeqabstract *sa,
                                          unsigned long idx);

void           gt_seqabstract_delete(GtSeqabstract *sa);

unsigned long gt_seqabstract_lcp(bool forward,
                                 const GtSeqabstract *useq,
                                 const GtSeqabstract *vseq,
                                 unsigned long leftstart,
                                 unsigned long rightstart,
                                 unsigned long minlen);

#endif
