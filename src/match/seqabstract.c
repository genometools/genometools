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

#include <limits.h>
#include "core/chardef.h"
#include "core/types_api.h"
#include "core/ma_api.h"
#include "greedyedist.h"

struct GtSeqabstract
{
  bool isptr;
  unsigned long len, offset;
  union
  {
    const GtUchar *ptr;
    const GtEncseq *encseq;
  } seq;
};

GtSeqabstract *gt_seqabstract_new_empty(void)
{
  GtSeqabstract *greedyedistseq = gt_malloc(sizeof *greedyedistseq);

  greedyedistseq->isptr = true;
  greedyedistseq->len = 0;
  greedyedistseq->offset = 0;
  greedyedistseq->seq.ptr = NULL;
  return greedyedistseq;
}

void gt_seqabstract_reinit_ptr(GtSeqabstract *greedyedistseq,
                                   const GtUchar *ptr,
                                   unsigned long len,
                                   unsigned long offset)
{
  greedyedistseq->isptr = true;
  greedyedistseq->len = len;
  greedyedistseq->offset = offset;
  greedyedistseq->seq.ptr = ptr + offset;
}

GtSeqabstract *gt_seqabstract_new_ptr(const GtUchar *ptr,
                                             unsigned long len,
                                             unsigned long offset)
{
  GtSeqabstract *greedyedistseq = gt_malloc(sizeof *greedyedistseq);

  gt_seqabstract_reinit_ptr(greedyedistseq,ptr,len,offset);
  return greedyedistseq;
}

void gt_seqabstract_reinit_encseq(GtSeqabstract *greedyedistseq,
                                      const GtEncseq *encseq,
                                      unsigned long len,
                                      unsigned long offset)
{
  greedyedistseq->isptr = false;
  greedyedistseq->len = len;
  greedyedistseq->offset = offset;
  greedyedistseq->seq.encseq = encseq;
}

GtSeqabstract *gt_seqabstract_new_encseq(const GtEncseq *encseq,
                                                unsigned long len,
                                                unsigned long offset)
{
  GtSeqabstract *greedyedistseq = gt_malloc(sizeof *greedyedistseq);

  gt_seqabstract_reinit_encseq(greedyedistseq, encseq, len, offset);
  return greedyedistseq;
}

unsigned long gt_seqabstract_length_get(const GtSeqabstract *greedyedistseq)
{
  return greedyedistseq->len;
}

void gt_seqabstract_delete(GtSeqabstract *greedyedistseq)
{
  gt_free(greedyedistseq);
}

GtUchar gt_seqabstract_encoded_char(const GtSeqabstract *greedyedistseq,
                                    unsigned long idx)
{
  return greedyedistseq->isptr
           ? greedyedistseq->seq.ptr[idx]
           : gt_encseq_get_encoded_char(greedyedistseq->seq.encseq,
                                        greedyedistseq->offset + idx,
                                        GT_READMODE_FORWARD);
}
