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
  GtSeqabstract *sa = gt_malloc(sizeof *sa);

  sa->isptr = true;
  sa->len = 0;
  sa->offset = 0;
  sa->seq.ptr = NULL;
  return sa;
}

void gt_seqabstract_reinit_ptr(GtSeqabstract *sa,
                                   const GtUchar *ptr,
                                   unsigned long len,
                                   unsigned long offset)
{
  sa->isptr = true;
  sa->len = len;
  sa->offset = offset;
  sa->seq.ptr = ptr + offset;
}

GtSeqabstract *gt_seqabstract_new_ptr(const GtUchar *ptr,
                                             unsigned long len,
                                             unsigned long offset)
{
  GtSeqabstract *sa = gt_malloc(sizeof *sa);

  gt_seqabstract_reinit_ptr(sa,ptr,len,offset);
  return sa;
}

void gt_seqabstract_reinit_encseq(GtSeqabstract *sa,
                                      const GtEncseq *encseq,
                                      unsigned long len,
                                      unsigned long offset)
{
  sa->isptr = false;
  sa->len = len;
  sa->offset = offset;
  sa->seq.encseq = encseq;
}

GtSeqabstract *gt_seqabstract_new_encseq(const GtEncseq *encseq,
                                                unsigned long len,
                                                unsigned long offset)
{
  GtSeqabstract *sa = gt_malloc(sizeof *sa);

  gt_seqabstract_reinit_encseq(sa, encseq, len, offset);
  return sa;
}

unsigned long gt_seqabstract_length_get(const GtSeqabstract *sa)
{
  return sa->len;
}

void gt_seqabstract_delete(GtSeqabstract *sa)
{
  gt_free(sa);
}

GtUchar gt_seqabstract_encoded_char(const GtSeqabstract *sa,
                                    unsigned long idx)
{
  return sa->isptr
           ? sa->seq.ptr[idx]
           : gt_encseq_get_encoded_char(sa->seq.encseq,
                                        sa->offset + idx,
                                        GT_READMODE_FORWARD);
}

unsigned long gt_seqabstract_lcp(bool *leftsep,
                                 bool *rightsep,
                                 bool forward,
                                 const GtSeqabstract *useq,
                                 const GtSeqabstract *vseq,
                                 unsigned long leftstart,
                                 unsigned long rightstart,
                                 unsigned long minlen)
{
  unsigned long lcp;
  GtUchar a, b;

  *leftsep = false;
  *rightsep = false;
  for (lcp = 0; lcp < minlen; lcp++)
  {
    a = gt_seqabstract_encoded_char(useq,forward ? leftstart + lcp
                                                 : leftstart - lcp);
    if (a == (GtUchar) SEPARATOR)
    {
      *leftsep = true;
      break;
    }
    b = gt_seqabstract_encoded_char(vseq,forward ? rightstart + lcp
                                                 : rightstart - lcp);
    if (b == (GtUchar) SEPARATOR)
    {
      *rightsep = true;
      break;
    }
    if (a != b || a == (GtUchar) WILDCARD)
    {
      break;
    }
  }
  return lcp;
}
