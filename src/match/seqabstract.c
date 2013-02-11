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
#include "greedyedist.h"

struct GtSeqabstract
{
  bool isptr, cmpcharbychar;
  unsigned long len, offset, totallength;
  GtEncseqReader *esr;
  bool stoppossupport;
  union
  {
    const GtUchar *ptr;
    const GtEncseq *encseq;
  } seq;
};

GtSeqabstract *gt_seqabstract_new_empty(void)
{
  GtSeqabstract *sa = gt_malloc(sizeof *sa);

  sa->cmpcharbychar = false;
  sa->isptr = true;
  sa->len = 0;
  sa->esr = NULL;
  sa->totallength = 0;
  sa->stoppossupport = false;
  sa->offset = 0;
  sa->seq.ptr = NULL;
  return sa;
}

void gt_seqabstract_reinit_ptr(GtSeqabstract *sa,
                               const GtUchar *ptr,
                               unsigned long len,
                               unsigned long offset)
{
  sa->cmpcharbychar = false;
  sa->isptr = true;
  sa->len = len;
  sa->stoppossupport = false;
  sa->totallength = 0;
  sa->offset = offset;
  if (sa->esr != NULL)
  {
    gt_encseq_reader_delete(sa->esr);
    sa->esr = NULL;
  }
  sa->seq.ptr = ptr + offset;
}

GtSeqabstract *gt_seqabstract_new_ptr(const GtUchar *ptr,
                                      unsigned long len,
                                      unsigned long offset)
{
  GtSeqabstract *sa = gt_seqabstract_new_empty();

  gt_seqabstract_reinit_ptr(sa,ptr,len,offset);
  return sa;
}

void gt_seqabstract_reinit_encseq(GtSeqabstract *sa,
                                  const GtEncseq *encseq,
                                  unsigned long len,
                                  unsigned long offset)
{
  sa->cmpcharbychar = gt_encseq_has_twobitencoding(encseq) ? false : true;
  sa->isptr = false;
  sa->len = len;
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
  sa->stoppossupport = gt_encseq_has_twobitencoding_stoppos_support(encseq);
  sa->offset = offset;
  sa->seq.encseq = encseq;
  sa->totallength = gt_encseq_total_length(encseq);
}

GtSeqabstract *gt_seqabstract_new_encseq(const GtEncseq *encseq,
                                         unsigned long len,
                                         unsigned long offset)
{
  GtSeqabstract *sa = gt_seqabstract_new_empty();

  gt_seqabstract_reinit_encseq(sa, encseq, len, offset);
  return sa;
}

unsigned long gt_seqabstract_length_get(const GtSeqabstract *sa)
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
  return sa->isptr
           ? sa->seq.ptr[idx]
           : gt_encseq_get_encoded_char(sa->seq.encseq,
                                        sa->offset + idx,
                                        GT_READMODE_FORWARD);
}

#define GT_SEQABSTRACT_CHECK_SEP(VAR,SEP)\
        if ((VAR) == (GtUchar) SEPARATOR)\
        {\
          *(SEP) = true;\
          break;\
        }

#define GT_SEQABSTRACT_CMPCHAR(VARA,VARB)\
        if ((VARA) != (VARB) || (VARA) == (GtUchar) WILDCARD)\
        {\
          break;\
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
  if (minlen == 0)
  {
    return 0;
  }
  if (useq->isptr)
  {
    if (vseq->isptr)
    {
      if (useq->seq.ptr == vseq->seq.ptr && useq->offset == vseq->offset &&
          leftstart == rightstart)
      {
        for (lcp = 0; lcp < minlen; lcp++)
        {
          a = useq->seq.ptr[forward ? leftstart + lcp : leftstart - lcp];
          if (ISSPECIAL(a))
          {
            break;
          }
        }
      } else
      {
        for (lcp = 0; lcp < minlen; lcp++)
        {
          a = useq->seq.ptr[forward ? leftstart + lcp : leftstart - lcp];
          GT_SEQABSTRACT_CHECK_SEP(a,leftsep);
          b = vseq->seq.ptr[forward ? rightstart + lcp : rightstart - lcp];
          GT_SEQABSTRACT_CHECK_SEP(b,rightsep);
          GT_SEQABSTRACT_CMPCHAR(a,b);
        }
      }
    } else
    {
      for (lcp = 0; lcp < minlen; lcp++)
      {
        a = useq->seq.ptr[forward ? leftstart + lcp : leftstart - lcp];
        GT_SEQABSTRACT_CHECK_SEP(a,leftsep);
        b = gt_encseq_get_encoded_char(vseq->seq.encseq,
                                       vseq->offset +
                                       (forward ? rightstart + lcp
                                                : rightstart - lcp),
                                       GT_READMODE_FORWARD);
        GT_SEQABSTRACT_CHECK_SEP(b,rightsep);
        GT_SEQABSTRACT_CMPCHAR(a,b);
      }
    }
  } else
  {
    if (vseq->isptr)
    {
      for (lcp = 0; lcp < minlen; lcp++)
      {
        a = gt_encseq_get_encoded_char(useq->seq.encseq,
                                       useq->offset +
                                       (forward ? leftstart + lcp
                                                : leftstart - lcp),
                                       GT_READMODE_FORWARD);
        GT_SEQABSTRACT_CHECK_SEP(a,leftsep);
        b = vseq->seq.ptr[forward ? rightstart + lcp : rightstart - lcp];
        GT_SEQABSTRACT_CHECK_SEP(b,rightsep);
        GT_SEQABSTRACT_CMPCHAR(a,b);
      }
    } else
    {
      if (useq->seq.encseq == vseq->seq.encseq && !useq->cmpcharbychar &&
          forward)
      {
        if (useq->offset + leftstart == vseq->offset + rightstart)
        {
          if (useq->stoppossupport)
          {
            unsigned long stoppos;

            gt_encseq_reader_reinit_with_readmode(useq->esr,
                                                  useq->seq.encseq,
                                                  GT_READMODE_FORWARD,
                                                  useq->offset + leftstart);
            stoppos = gt_getnexttwobitencodingstoppos(true,useq->esr);
            gt_assert(useq->offset + leftstart <= stoppos);
            lcp = stoppos - (useq->offset + leftstart);
            a = gt_encseq_get_encoded_char(useq->seq.encseq,
                                           stoppos,
                                           GT_READMODE_FORWARD);
            if (a == (GtUchar) SEPARATOR)
            {
              *leftsep = true;
            }
          } else
          {
            for (lcp = 0; lcp < minlen; lcp++)
            {
              a = gt_encseq_get_encoded_char(useq->seq.encseq,
                                             useq->offset + leftstart + lcp,
                                             GT_READMODE_FORWARD);
              if (a == (GtUchar) SEPARATOR)
              {
                *leftsep = true;
                break;
              }
              if (a == (GtUchar) WILDCARD)
              {
                break;
              }
            }
          }
        } else
        {
          GtCommonunits commonunits;

          (void) gt_encseq_compare_viatwobitencoding(&commonunits,
                                                     useq->seq.encseq,
                                                     vseq->seq.encseq,
                                                     GT_READMODE_FORWARD,
                                                     useq->esr,
                                                     vseq->esr,
                                                     useq->offset + leftstart,
                                                     vseq->offset + rightstart,
                                                     0,
                                                     minlen);
          lcp = commonunits.finaldepth;
          if (commonunits.leftspecial &&
              (useq->offset + leftstart + lcp == useq->totallength ||
               gt_encseq_position_is_separator(useq->seq.encseq,
                                               useq->offset + leftstart + lcp,
                                               GT_READMODE_FORWARD)))
          {
            *leftsep = true;
          } else
          {
            if (commonunits.rightspecial &&
                (vseq->offset + rightstart + lcp == vseq->totallength ||
                 gt_encseq_position_is_separator(vseq->seq.encseq,
                                                 vseq->offset+ rightstart + lcp,
                                                 GT_READMODE_FORWARD)))
            {
              *rightsep = true;
            }
          }
        }
      } else
      {
        for (lcp = 0; lcp < minlen; lcp++)
        {
          a = gt_encseq_get_encoded_char(useq->seq.encseq,
                                         useq->offset +
                                         (forward ? leftstart + lcp
                                                  : leftstart - lcp),
                                         GT_READMODE_FORWARD);
          GT_SEQABSTRACT_CHECK_SEP(a,leftsep);
          b = gt_encseq_get_encoded_char(vseq->seq.encseq,
                                         vseq->offset +
                                         (forward ? rightstart + lcp
                                                  : rightstart - lcp),
                                         GT_READMODE_FORWARD);
          GT_SEQABSTRACT_CHECK_SEP(b,rightsep);
          GT_SEQABSTRACT_CMPCHAR(a,b);
        }
      }
    }
  }
  return lcp;
}
