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
#include "core/minmax.h"
#include "core/encseq.h"
#include "stamp.h"
#include "greedyedist.h"
#include "core/xansi_api.h"

typedef enum GtSeqabstractType {
  GT_SEQABSTRACT_UNDEF,
  GT_SEQABSTRACT_STRING,
  GT_SEQABSTRACT_ENCSEQ
} GtSeqabstractType;

struct GtSeqabstract
{
  unsigned long len, offset, totallength;
  GtEncseqReader *esr;
  GtSeqabstractType seqtype;
  bool cmpcharbychar,
       stoppossupport;
  union
  {
    const GtUchar *string;
    const GtEncseq *encseq;
  } seq;
};

GtSeqabstract *gt_seqabstract_new_empty(void)
{
  GtSeqabstract *sa = gt_malloc(sizeof *sa);

  sa->cmpcharbychar = false;
  sa->seqtype = GT_SEQABSTRACT_UNDEF;
  sa->len = 0;
  sa->esr = NULL;
  sa->totallength = 0;
  sa->stoppossupport = false;
  sa->offset = 0;
  sa->seq.string = NULL;
  return sa;
}

void gt_seqabstract_reinit_gtuchar(GtSeqabstract *sa,
                                   const GtUchar *string,
                                   unsigned long len,
                                   unsigned long offset)
{
  if (sa->esr != NULL)
  {
    gt_encseq_reader_delete(sa->esr);
    sa->esr = NULL;
  }
  sa->seqtype = GT_SEQABSTRACT_STRING;
  sa->cmpcharbychar = false;
  sa->stoppossupport = false;
  sa->len = len;
  sa->totallength = 0;
  sa->offset = offset;
  sa->seq.string = string + offset;
}

GtSeqabstract *gt_seqabstract_new_gtuchar(const GtUchar *string,
                                          unsigned long len,
                                          unsigned long offset)
{
  GtSeqabstract *sa = gt_seqabstract_new_empty();

  gt_seqabstract_reinit_gtuchar(sa, string, len, offset);
  return sa;
}

void gt_seqabstract_reinit_encseq(GtSeqabstract *sa,
                                  const GtEncseq *encseq,
                                  unsigned long len,
                                  unsigned long offset)
{
  if (sa->esr != NULL)
  {
    if (encseq != sa->seq.encseq)
    {
      gt_encseq_reader_delete(sa->esr);
      sa->esr = gt_encseq_create_reader_with_readmode(encseq,
                                                      GT_READMODE_FORWARD,0);
    }
  } else
  {
    sa->esr = gt_encseq_create_reader_with_readmode(encseq,
                                                    GT_READMODE_FORWARD,0);
  }
  sa->seqtype = GT_SEQABSTRACT_ENCSEQ;
  sa->cmpcharbychar = gt_encseq_has_twobitencoding(encseq) ? false : true;
  sa->stoppossupport = gt_encseq_has_twobitencoding_stoppos_support(encseq);
  sa->len = len;
  sa->totallength = gt_encseq_total_length(encseq);
  sa->offset = offset;
  sa->seq.encseq = encseq;
}

GtSeqabstract *gt_seqabstract_new_encseq(const GtEncseq *encseq,
                                         unsigned long len,
                                         unsigned long offset)
{
  GtSeqabstract *sa = gt_seqabstract_new_empty();

  gt_seqabstract_reinit_encseq(sa, encseq, len, offset);
  return sa;
}

unsigned long gt_seqabstract_length(const GtSeqabstract *sa)
{
  return sa->len;
}

void gt_seqabstract_delete(GtSeqabstract *sa)
{
  gt_encseq_reader_delete(sa->esr);
  gt_free(sa);
}

GtUchar gt_seqabstract_encoded_char(const GtSeqabstract *sa,
                                    unsigned long idx)
{
  return sa->seqtype == GT_SEQABSTRACT_STRING ?
    sa->seq.string[idx] :
    gt_encseq_get_encoded_char(sa->seq.encseq, sa->offset + idx,
                               GT_READMODE_FORWARD);
}

#define GT_SEQABSTRACT_CHAR_EQUAL_OR_BREAK(VARA,VARB)\
        if ((VARA) != (VARB) || ISSPECIAL(VARA))\
        {\
          break;\
        }

static unsigned long
gt_seqabstract_lcp_gtuchar_gtuchar(bool forward,
                                   const GtSeqabstract *useq,
                                   const GtSeqabstract *vseq,
                                   unsigned long leftstart,
                                   unsigned long rightstart,
                                   unsigned long minlen)
{
  unsigned long lcp;
  GtUchar a, b;
  if (useq->seq.string == vseq->seq.string &&
      useq->offset == vseq->offset &&
      leftstart == rightstart)
  {
    for (lcp = 0; lcp < minlen; lcp++)
    {
      a = useq->seq.string[forward ? leftstart + lcp : leftstart - lcp];
      if (ISSPECIAL(a))
      {
        break;
      }
    }
  } else
  {
    for (lcp = 0; lcp < minlen; lcp++)
    {
      a = useq->seq.string[forward ? leftstart + lcp : leftstart - lcp];
      b = vseq->seq.string[forward ? rightstart + lcp : rightstart - lcp];
      GT_SEQABSTRACT_CHAR_EQUAL_OR_BREAK(a,b);
    }
  }
  return lcp;
}

static unsigned long
gt_seqabstract_lcp_gtuchar_encseq(bool forward,
                                  const GtSeqabstract *useq,
                                  const GtSeqabstract *vseq,
                                  unsigned long leftstart,
                                  unsigned long rightstart,
                                  unsigned long minlen)
{
  unsigned long lcp;
  GtUchar a, b;
  for (lcp = 0; lcp < minlen; lcp++)
  {
    a = useq->seq.string[forward ? leftstart + lcp : leftstart - lcp];
    b = gt_encseq_get_encoded_char(vseq->seq.encseq,
                                   vseq->offset +
                                   (forward ? rightstart + lcp
                                            : rightstart - lcp),
                                   GT_READMODE_FORWARD);
    GT_SEQABSTRACT_CHAR_EQUAL_OR_BREAK(a,b);
  }
  return lcp;
}

static unsigned long gt_seqabstract_lcp_encseq_encseq(bool forward,
                                                      const GtSeqabstract *useq,
                                                      const GtSeqabstract *vseq,
                                                      unsigned long leftstart,
                                                      unsigned long rightstart,
                                                      unsigned long minlen)
{
  unsigned long lcp;
  GtUchar a, b;
  bool is_same_sequence = (useq->seq.encseq == vseq->seq.encseq) &&
    ((useq->offset + leftstart) == (vseq->offset + rightstart));

  if ((useq->cmpcharbychar || vseq->cmpcharbychar) || !is_same_sequence)
  {
    for (lcp = 0; lcp < minlen; lcp++)
    {
      a = gt_encseq_get_encoded_char(useq->seq.encseq,
                                     useq->offset +
                                     (forward ? leftstart + lcp
                                              : leftstart - lcp),
                                     GT_READMODE_FORWARD);
      b = gt_encseq_get_encoded_char(vseq->seq.encseq,
                                     vseq->offset +
                                     (forward ? rightstart + lcp
                                              : rightstart - lcp),
                                     GT_READMODE_FORWARD);
      GT_SEQABSTRACT_CHAR_EQUAL_OR_BREAK(a,b);
    }
  }
  else {
    if (useq->stoppossupport)
    {
      unsigned long stoppos;
      const unsigned long startpos
         = forward ? (useq->offset + leftstart)
                   : GT_REVERSEPOS(useq->totallength,
                                   useq->offset + leftstart);

      gt_encseq_reader_reinit_with_readmode(useq->esr,
                                            useq->seq.encseq,
                                            forward
                                              ? GT_READMODE_FORWARD
                                              : GT_READMODE_REVERSE,
                                            startpos);
      stoppos = gt_getnexttwobitencodingstoppos(forward,useq->esr);
      if (!forward)
      {
        stoppos = GT_REVERSEPOS(useq->totallength+1,stoppos);
      }
      gt_assert(startpos <= stoppos);
      lcp = MIN(minlen,stoppos - startpos);
    } else
    {
      unsigned long startpos = useq->offset + leftstart;
      for (lcp = 0; lcp < minlen; lcp++)
      {
        a = gt_encseq_get_encoded_char(useq->seq.encseq,
                                       forward ? (startpos + lcp)
                                               : (startpos - lcp),
                                       GT_READMODE_FORWARD);
        if (ISSPECIAL(a))
        {
          break;
        }
      }
    }
  }
  return lcp;
}

unsigned long gt_seqabstract_lcp(bool forward,
                                 const GtSeqabstract *useq,
                                 const GtSeqabstract *vseq,
                                 unsigned long leftstart,
                                 unsigned long rightstart,
                                 unsigned long minlen)
{
  unsigned long lcp;

  if (minlen == 0)
  {
    return 0;
  }
  if (useq->seqtype == GT_SEQABSTRACT_STRING)
  {
    if (vseq->seqtype == GT_SEQABSTRACT_STRING)
    {
      lcp = gt_seqabstract_lcp_gtuchar_gtuchar(forward, useq, vseq, leftstart,
                                               rightstart, minlen);
    } else
    {
      lcp = gt_seqabstract_lcp_gtuchar_encseq(forward, useq, vseq, leftstart,
                                              rightstart, minlen);
    }
  } else
  {
    if (vseq->seqtype == GT_SEQABSTRACT_STRING)
    {
      lcp = gt_seqabstract_lcp_gtuchar_encseq(forward, vseq, useq, rightstart,
                                              leftstart, minlen);

    } else
    {
      lcp = gt_seqabstract_lcp_encseq_encseq(forward, vseq, useq, rightstart,
                                             leftstart, minlen);
    }
  }
  return lcp;
}
