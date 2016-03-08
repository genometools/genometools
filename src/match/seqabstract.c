/*
  Copyright (c) 2013, 2015  Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
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
#include "match/seqabstract.h"
#include "match/extend-offset.h"

#define GT_SEQABSTRACT_TOTALLENGTH_UNDEF GT_UWORD_MAX

typedef enum GtSeqabstractType {
  GT_SEQABSTRACT_UNDEF,
  GT_SEQABSTRACT_STRING,
  GT_SEQABSTRACT_ENCSEQ
} GtSeqabstractType;

struct GtSeqabstract
{
  GtUword len, offset, totallength, seqstartpos;
  GtSeqabstractType seqtype;
  bool read_seq_left2right,
       dir_is_complement;
  union
  {
    const GtUchar *string;
    const GtEncseq *encseq;
  } seq;
};

GtSeqabstract *gt_seqabstract_new_empty(void)
{
  GtSeqabstract *sa = gt_malloc(sizeof *sa);

  sa->seqtype = GT_SEQABSTRACT_UNDEF;
  sa->len = sa->offset = 0;
  sa->read_seq_left2right = true;
  sa->dir_is_complement = false;
  sa->totallength = GT_SEQABSTRACT_TOTALLENGTH_UNDEF;
  sa->seqstartpos = 0;
  sa->seq.string = NULL;
  return sa;
}

void gt_seqabstract_reset(GtSeqabstract *sa)
{
  sa->seqtype = GT_SEQABSTRACT_UNDEF;
  sa->len = sa->offset = 0;
  sa->read_seq_left2right = true;
  sa->dir_is_complement = false;
  sa->totallength = GT_SEQABSTRACT_TOTALLENGTH_UNDEF;
  sa->seqstartpos = 0;
  sa->seq.string = NULL;
}

static void gt_seqabstract_init(GtSeqabstract *sa,
                                bool rightextension,
                                GtReadmode readmode,
                                GtUword len,
                                GtUword startpos,
                                GtUword totallength)
{
  sa->len = len;
  gt_assert(startpos >= sa->seqstartpos && (
            !GT_ISDIRREVERSE(readmode) ||
            totallength != GT_SEQABSTRACT_TOTALLENGTH_UNDEF));
  sa->offset = GT_EXTEND_OFFSET(rightextension,
                                readmode,
                                totallength,
                                sa->seqstartpos,
                                startpos,
                                len);
  sa->read_seq_left2right = GT_EXTEND_READ_SEQ_LEFT2RIGHT(rightextension,
                                                          readmode);
  sa->dir_is_complement = GT_ISDIRCOMPLEMENT(readmode) ? true : false;
}

void gt_seqabstract_seqstartpos_set(GtSeqabstract *sa,GtUword seqstartpos)
{
  gt_assert(sa != NULL);
  sa->seqstartpos = seqstartpos;
}

void gt_seqabstract_totallength_set(GtSeqabstract *sa,GtUword totallength)
{
  gt_assert(sa != NULL);
  sa->totallength = totallength;
}

void gt_seqabstract_reinit_gtuchar(bool rightextension,
                                   GtReadmode readmode,
                                   GtSeqabstract *sa,
                                   const GtUchar *string,
                                   GtUword len,
                                   GtUword startpos,
                                   GtUword totallength)
{
  gt_assert(sa != NULL && totallength != GT_SEQABSTRACT_TOTALLENGTH_UNDEF);

  sa->seqtype = GT_SEQABSTRACT_STRING;
  sa->totallength = totallength;
  sa->seq.string = string;
  gt_seqabstract_init(sa,
                      rightextension,
                      readmode,
                      len,
                      startpos,
                      totallength);
}

GtSeqabstract *gt_seqabstract_new_gtuchar(bool rightextension,
                                          GtReadmode readmode,
                                          const GtUchar *string,
                                          GtUword len,
                                          GtUword startpos,
                                          GtUword totallength)
{
  GtSeqabstract *sa = gt_seqabstract_new_empty();

  gt_seqabstract_reinit_gtuchar(rightextension, readmode, sa, string, len,
                                startpos, totallength);
  return sa;
}

void gt_seqabstract_reinit_encseq(bool rightextension,
                                  GtReadmode readmode,
                                  GtSeqabstract *sa,
                                  const GtEncseq *encseq,
                                  GtUword len,
                                  GtUword startpos)
{
  gt_assert(sa != NULL);
  sa->seqtype = GT_SEQABSTRACT_ENCSEQ;
  sa->seq.encseq = encseq;
  gt_seqabstract_init(sa,
                      rightextension,
                      readmode,
                      len,
                      startpos,
                      sa->totallength);
}

GtSeqabstract *gt_seqabstract_new_encseq(bool rightextension,
                                         GtReadmode readmode,
                                         const GtEncseq *encseq,
                                         GtUword len,
                                         GtUword startpos)
{
  GtSeqabstract *sa = gt_seqabstract_new_empty();

  gt_seqabstract_reinit_encseq(rightextension, readmode, sa, encseq, len,
                               startpos);
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
    gt_free(sa);
  }
}

static GtUchar gt_seqabstract_get_encoded_char(GT_UNUSED bool rightextension,
                                               const GtSeqabstract *seq,
                                               GtUword idx)
{
  GtUword accesspos;

  if (seq->read_seq_left2right)
  {
    accesspos = seq->offset + idx;
  } else
  {
    gt_assert(seq->offset >= idx);
    accesspos = seq->offset - idx;
  }
  if (seq->seqtype == GT_SEQABSTRACT_STRING)
  {
    return seq->seq.string[accesspos];
  }
  return gt_encseq_get_encoded_char(seq->seq.encseq,accesspos,
                                    GT_READMODE_FORWARD);
}

GtUword gt_seqabstract_lcp(bool rightextension,
                           const GtSeqabstract *useq,
                           const GtSeqabstract *vseq,
                           GtUword u_start,
                           GtUword v_start)
{
  GtUword lcp, maxlen;

  gt_assert(useq != NULL && vseq != NULL &&
            useq->len >= u_start && vseq->len >= v_start);
  maxlen = MIN(useq->len - u_start, vseq->len - v_start);
  for (lcp = 0; lcp < maxlen; lcp++)
  {
    GtUchar u_cc, v_cc;

    u_cc = gt_seqabstract_get_encoded_char(rightextension,useq,u_start + lcp);
    if (ISSPECIAL(u_cc))
    {
      break;
    }
    if (useq->dir_is_complement)
    {
      u_cc = GT_COMPLEMENTBASE(u_cc);
    }
    v_cc = gt_seqabstract_get_encoded_char(rightextension,vseq,v_start + lcp);
    if (ISSPECIAL(v_cc))
    {
      break;
    }
    if (vseq->dir_is_complement)
    {
      v_cc = GT_COMPLEMENTBASE(v_cc);
    }
    if (u_cc != v_cc)
    {
      break;
    }
  }
  return lcp;
}

char *gt_seqabstract_get(bool rightextension,const GtSeqabstract *seq)
{
  GtUword idx;
  char *buffer = gt_malloc(sizeof *buffer * (seq->len+1));
  char *map = "acgt";

  for (idx = 0; idx < seq->len; idx++)
  {
    GtUchar cc = gt_seqabstract_get_encoded_char(rightextension,seq,idx);

    if (cc == WILDCARD)
    {
      buffer[idx] = '#';
    } else
    {
      if (cc == SEPARATOR)
      {
        buffer[idx] = '$';
      } else
      {
        gt_assert(cc < 4);
        buffer[idx] = map[cc];
      }
    }
  }
  buffer[seq->len] = '\0';
  return buffer;
}
