/*
  Copyright (c) 2013        Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c)        2014 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
  Copyright (c) 2013 - 2014 Center for Bioinformatics, University of Hamburg

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

#include "core/assert_api.h"
#include "core/chardef.h"
#include "core/encseq.h"
#include "core/ma_api.h"
#include "core/minmax.h"
#include "core/types_api.h"
#include "core/xansi_api.h"
#include "match/greedyedist.h"
#include "match/stamp.h"

typedef enum GtSeqabstractType {
  GT_SEQABSTRACT_UNDEF,
  GT_SEQABSTRACT_STRING,
  GT_SEQABSTRACT_ENCSEQ
} GtSeqabstractType;

struct GtSeqabstract
{
  GtUword len, offset, totallength;
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
                                   GtUword len,
                                   GtUword offset)
{
  gt_assert(sa != NULL);
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
                                          GtUword len,
                                          GtUword offset)
{
  GtSeqabstract *sa = gt_seqabstract_new_empty();

  gt_seqabstract_reinit_gtuchar(sa, string, len, offset);
  return sa;
}

void gt_seqabstract_reinit_encseq(GtSeqabstract *sa,
                                  const GtEncseq *encseq,
                                  GtUword len,
                                  GtUword offset)
{
  gt_assert(sa != NULL);
  if (sa->esr != NULL)
  {
    if (encseq != sa->seq.encseq)
    {
      gt_encseq_reader_delete(sa->esr);
      sa->esr = gt_encseq_create_reader_with_readmode(encseq,
                                                      GT_READMODE_FORWARD,
                                                      offset);
    }
  } else
  {
    sa->esr = gt_encseq_create_reader_with_readmode(encseq,
                                                    GT_READMODE_FORWARD,
                                                    offset);
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
                                         GtUword len,
                                         GtUword offset)
{
  GtSeqabstract *sa = gt_seqabstract_new_empty();

  gt_seqabstract_reinit_encseq(sa, encseq, len, offset);
  return sa;
}

GtUword gt_seqabstract_length(const GtSeqabstract *sa)
{
  gt_assert(sa != NULL);
  return sa->len;
}

void gt_seqabstract_delete(GtSeqabstract *sa)
{
  if (sa != NULL) {
    gt_encseq_reader_delete(sa->esr);
    gt_free(sa);
  }
}

GtUchar gt_seqabstract_encoded_char(const GtSeqabstract *sa,
                                    GtUword idx)
{
  gt_assert(sa != NULL);
  gt_assert(idx < sa->len);
  return sa->seqtype == GT_SEQABSTRACT_STRING ?
         sa->seq.string[sa->offset + idx] :
    gt_encseq_get_encoded_char(sa->seq.encseq, sa->offset + idx,
                               GT_READMODE_FORWARD);
}

#define GT_SEQABSTRACT_CHAR_EQUAL_OR_BREAK(VARA,VARB)\
        if ((VARA) != (VARB) || ISSPECIAL(VARA))\
        {\
          break;\
        }

static GtUword
gt_seqabstract_lcp_gtuchar_gtuchar(bool forward,
                                   const GtSeqabstract *useq,
                                   const GtSeqabstract *vseq,
                                   GtUword ustart,
                                   GtUword vstart,
                                   GtUword maxlen)
{
  GtUword lcp;
  GtUchar a, b;
  if (useq->seq.string == vseq->seq.string &&
      useq->offset == vseq->offset &&
      ustart == vstart)
  {
    for (lcp = 0; lcp < maxlen; lcp++)
    {
      a = useq->seq.string[forward ? ustart + lcp : ustart - lcp];
      if (ISSPECIAL(a))
      {
        break;
      }
    }
  } else
  {
    for (lcp = 0; lcp < maxlen; lcp++)
    {
      a = useq->seq.string[forward ? ustart + lcp : ustart - lcp];
      b = vseq->seq.string[forward ? vstart + lcp : vstart - lcp];
      GT_SEQABSTRACT_CHAR_EQUAL_OR_BREAK(a,b);
    }
  }
  return lcp;
}

static GtUword
gt_seqabstract_lcp_gtuchar_encseq(bool forward,
                                  const GtSeqabstract *useq,
                                  const GtSeqabstract *vseq,
                                  GtUword ustart,
                                  GtUword vstart,
                                  GtUword maxlen)
{
  GtUword lcp;
  GtUchar a, b;
  for (lcp = 0; lcp < maxlen; lcp++)
  {
    a = useq->seq.string[forward ? ustart + lcp : ustart - lcp];
    b = gt_encseq_get_encoded_char(vseq->seq.encseq,
                                   vseq->offset + (forward ?
                                                   vstart + lcp :
                                                   vstart - lcp),
                                   GT_READMODE_FORWARD);
    GT_SEQABSTRACT_CHAR_EQUAL_OR_BREAK(a,b);
  }
  return lcp;
}

static GtUword gt_seqabstract_lcp_encseq_encseq(bool forward,
                                                const GtSeqabstract *useq,
                                                const GtSeqabstract *vseq,
                                                GtUword ustart,
                                                GtUword vstart,
                                                GtUword maxlen)
{
  GtUword lcp;
  GtUchar a, b;
  bool is_same_sequence = (useq->seq.encseq == vseq->seq.encseq) &&
    ((useq->offset + ustart) == (vseq->offset + vstart));

  if ((useq->cmpcharbychar || vseq->cmpcharbychar) || !is_same_sequence)
  {
    for (lcp = 0; lcp < maxlen; lcp++)
    {
      a = gt_encseq_get_encoded_char(useq->seq.encseq,
                                     useq->offset + (forward ?
                                                     ustart + lcp :
                                                     ustart - lcp),
                                     GT_READMODE_FORWARD);
      b = gt_encseq_get_encoded_char(vseq->seq.encseq,
                                     vseq->offset + (forward ?
                                                     vstart + lcp :
                                                     vstart - lcp),
                                     GT_READMODE_FORWARD);
      GT_SEQABSTRACT_CHAR_EQUAL_OR_BREAK(a,b);
    }
  }
  else {
    if (useq->stoppossupport)
    {
      GtUword stoppos;
      const GtUword startpos = forward ?
                               (useq->offset + ustart) :
                               GT_REVERSEPOS(useq->totallength,
                                             useq->offset + ustart);

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
      lcp = MIN(maxlen,stoppos - startpos);
    } else
    {
      GtUword startpos = useq->offset + ustart;
      for (lcp = 0; lcp < maxlen; lcp++)
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

GtUword gt_seqabstract_lcp(bool forward,
                           const GtSeqabstract *useq,
                           const GtSeqabstract *vseq,
                           GtUword ustart,
                           GtUword vstart)
{
  GtUword lcp,
          maxlen;

  gt_assert(useq != NULL && vseq != NULL);
  maxlen = forward ?
           MIN(useq->len - ustart, vseq->len - vstart) :
           MIN(ustart + 1, vstart + 1);

  if (useq->seqtype == GT_SEQABSTRACT_STRING)
  {
    if (vseq->seqtype == GT_SEQABSTRACT_STRING)
    {
      lcp = gt_seqabstract_lcp_gtuchar_gtuchar(forward, useq, vseq, ustart,
                                               vstart, maxlen);
    } else
    {
      lcp = gt_seqabstract_lcp_gtuchar_encseq(forward, useq, vseq, ustart,
                                              vstart, maxlen);
    }
  } else
  {
    if (vseq->seqtype == GT_SEQABSTRACT_STRING)
    {
      lcp = gt_seqabstract_lcp_gtuchar_encseq(forward, vseq, useq, vstart,
                                              ustart, maxlen);

    } else
    {
      lcp = gt_seqabstract_lcp_encseq_encseq(forward, vseq, useq, vstart,
                                             ustart, maxlen);
    }
  }
  return lcp;
}
