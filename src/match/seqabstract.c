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
#include "core/xansi_api.h"
#include "match/greedyedist.h"

#define GT_SEQABSTRACT_TOTALLENGTH_UNDEF GT_UWORD_MAX

typedef enum GtSeqabstractType {
  GT_SEQABSTRACT_UNDEF,
  GT_SEQABSTRACT_STRING,
  GT_SEQABSTRACT_ENCSEQ
} GtSeqabstractType;

struct GtSeqabstract
{
  GtUword len, offset, totallength;
  GtSeqabstractType seqtype;
  GtReadmode readmode;
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
  sa->len = 0;
  sa->totallength = GT_SEQABSTRACT_TOTALLENGTH_UNDEF;
  sa->readmode = GT_READMODE_FORWARD; /* default read mode */
  sa->offset = 0;
  sa->seq.string = NULL;
  return sa;
}

void gt_seqabstract_readmode_set(GtSeqabstract *sa,GtReadmode readmode)
{
  gt_assert(sa != NULL);
  sa->readmode = readmode;
}

void gt_seqabstract_reinit_gtuchar(GtSeqabstract *sa,
                                   const GtUchar *string,
                                   GtUword len,
                                   GtUword offset,
                                   GtUword totallength)
{
  gt_assert(sa != NULL && totallength != GT_SEQABSTRACT_TOTALLENGTH_UNDEF);

  sa->seqtype = GT_SEQABSTRACT_STRING;
  sa->len = len;
  sa->totallength = totallength;
  sa->offset = offset;
  sa->seq.string = string;
}

GtSeqabstract *gt_seqabstract_new_gtuchar(const GtUchar *string,
                                          GtUword len,
                                          GtUword offset,
                                          GtUword totallength)
{
  GtSeqabstract *sa = gt_seqabstract_new_empty();

  gt_seqabstract_reinit_gtuchar(sa, string, len, offset, totallength);
  return sa;
}

void gt_seqabstract_reinit_encseq(GtSeqabstract *sa,
                                  const GtEncseq *encseq,
                                  GtUword len,
                                  GtUword offset)
{
  gt_assert(sa != NULL);
  sa->seqtype = GT_SEQABSTRACT_ENCSEQ;
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
    gt_free(sa);
  }
}

static GtUchar gt_seqabstract_get_encoded_char(bool rightextension,
                                               const GtSeqabstract *seq,
                                               GtUword idx)
{
  GtUword accesspos;

  if (rightextension)
  {
    if (GT_ISDIRREVERSE(seq->readmode))
    {
      gt_assert(seq->totallength != GT_SEQABSTRACT_TOTALLENGTH_UNDEF &&
                seq->offset + idx < seq->totallength);
      accesspos = seq->totallength - 1 - seq->offset - idx;
    } else
    {
      accesspos = seq->offset + idx;
    }
  } else
  {
    if (GT_ISDIRREVERSE(seq->readmode))
    {
      gt_assert(seq->totallength >= seq->len);
      accesspos = seq->offset + seq->totallength - seq->len + idx;
    } else
    {
      gt_assert(seq->offset + seq->len > idx);
      accesspos = seq->offset + seq->len - 1 - idx;
    }
  }
  gt_assert(accesspos < seq->totallength);
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
    if (GT_ISDIRCOMPLEMENT(useq->readmode))
    {
      u_cc = GT_COMPLEMENTBASE(u_cc);
    }
    v_cc = gt_seqabstract_get_encoded_char(rightextension,vseq,v_start + lcp);
    if (ISSPECIAL(v_cc))
    {
      break;
    }
    if (GT_ISDIRCOMPLEMENT(vseq->readmode))
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
  char *buffer = malloc(sizeof *buffer * (seq->len+1));
  char *map = "acgt";

  printf("# readmode=%s,rightextension=%s,totallength=" GT_WU ",len=" GT_WU
         ",offset=" GT_WU "\n",
           gt_readmode_show(seq->readmode),rightextension ? "true" : "false",
          seq->totallength,seq->len,seq->offset);
  for (idx = 0; idx < seq->len; idx++)
  {
    GtUchar cc = gt_seqabstract_get_encoded_char(rightextension,seq,idx);

    gt_assert(cc < 4);
    buffer[idx] = map[cc];
  }
  buffer[seq->len] = '\0';
  return buffer;
}
