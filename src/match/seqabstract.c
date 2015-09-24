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
  const GtUchar *origstring;
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
  sa->totallength = GT_UWORD_MAX;
  sa->readmode = GT_READMODE_FORWARD; /* default read mode */
  sa->offset = 0;
  sa->seq.string = NULL;
  sa->origstring = NULL;
  return sa;
}

void gt_seqabstract_readmode_set(GtSeqabstract *sa,GtReadmode readmode)
{
  gt_assert(sa != NULL && (!GT_ISDIRREVERSE(readmode) ||
                           sa->totallength != GT_UWORD_MAX));
  sa->readmode = readmode;
}

void gt_seqabstract_reinit_gtuchar(GtSeqabstract *sa,
                                   const GtUchar *string,
                                   GtUword len,
                                   GtUword offset,
                                   GtUword totallength)
{
  gt_assert(sa != NULL && totallength != GT_UWORD_MAX);

  sa->seqtype = GT_SEQABSTRACT_STRING;
  sa->len = len;
  sa->totallength = totallength;
  sa->offset = offset;
  sa->seq.string = string + sa->offset;
  sa->origstring = string;
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
      gt_assert(seq->offset + idx < seq->totallength);
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
    return seq->origstring[accesspos];
  }
   return gt_encseq_get_encoded_char(seq->seq.encseq,accesspos,
                                     GT_READMODE_FORWARD);
}

#define GT_SEQABSTRACT_ASSIGN_READMODES\
        GtReadmode u_readmode, v_readmode;\
        if (rightextension)\
        {\
          u_readmode = useq->readmode;\
          v_readmode = vseq->readmode;\
        } else\
        {\
          u_readmode = gt_readmode_inverse_dir(useq->readmode);\
          v_readmode = gt_readmode_inverse_dir(vseq->readmode);\
        }

#define GT_SEQABSTRACT_CHAR_EQUAL_OR_BREAK(VARA,RMA,VARB,RMB)\
        if (GT_ISDIRCOMPLEMENT(RMA))\
        {\
          VARA = GT_COMPLEMENTBASE(VARA);\
        }\
        if (GT_ISDIRCOMPLEMENT(v_readmode))\
        {\
          VARB = GT_COMPLEMENTBASE(VARB);\
        }\
        if (VARA != VARB)\
        {\
          break;\
        }

#define GT_SEQABSTRACT_SEQ_ACCESS(VAR,STR,RM,START,LCP)\
        VAR = STR[GT_ISDIRREVERSE(RM) ? ((START) - (LCP)) : ((START) + (LCP))];\
        if (ISSPECIAL(VAR))\
        {\
          break;\
        }

#define GT_SEQABSTRACT_ENC_ACCESS(VAR,ENC,OFF,RM,START,LCP)\
        VAR = gt_encseq_get_encoded_char(ENC,\
                                         (OFF) + (GT_ISDIRREVERSE(RM) \
                                                     ? ((START) - (LCP))\
                                                     : ((START) + (LCP))),\
                                         GT_READMODE_FORWARD);\
        if (ISSPECIAL(VAR))\
        {\
          break;\
        }

static GtUword
gt_seqabstract_lcp_gtuchar_gtuchar(bool rightextension,
                                   const GtSeqabstract *useq,
                                   const GtSeqabstract *vseq,
                                   GtUword ustart,
                                   GtUword vstart,
                                   GtUword maxlen)
{
  GtUword lcp;

  if (useq->seq.string == vseq->seq.string &&
      useq->offset == vseq->offset &&
      ustart == vstart && useq->readmode == vseq->readmode)
  {
    for (lcp = 0; lcp < maxlen; lcp++)
    {
      GtUchar u_cc
        = useq->seq.string[rightextension ? ustart + lcp : ustart - lcp];
      if (ISSPECIAL(u_cc))
      {
        break;
      }
    }
  } else
  {
    GT_SEQABSTRACT_ASSIGN_READMODES
    for (lcp = 0; lcp < maxlen; lcp++)
    {
      GtUchar u_cc, v_cc;

      GT_SEQABSTRACT_SEQ_ACCESS(u_cc,useq->seq.string,u_readmode,ustart,lcp);
      GT_SEQABSTRACT_SEQ_ACCESS(v_cc,vseq->seq.string,v_readmode,vstart,lcp);
      GT_SEQABSTRACT_CHAR_EQUAL_OR_BREAK(u_cc,u_readmode,v_cc,v_readmode);
    }
  }
  return lcp;
}

static GtUword
gt_seqabstract_lcp_gtuchar_encseq(bool rightextension,
                                  const GtSeqabstract *useq,
                                  const GtSeqabstract *vseq,
                                  GtUword ustart,
                                  GtUword vstart,
                                  GtUword maxlen)
{
  GtUword lcp;
  GT_SEQABSTRACT_ASSIGN_READMODES
  for (lcp = 0; lcp < maxlen; lcp++)
  {
    GtUchar u_cc, v_cc;

    GT_SEQABSTRACT_SEQ_ACCESS(u_cc,useq->seq.string,u_readmode,ustart,lcp);
    GT_SEQABSTRACT_ENC_ACCESS(v_cc,vseq->seq.encseq,vseq->offset,v_readmode,
                              vstart,lcp);
    GT_SEQABSTRACT_CHAR_EQUAL_OR_BREAK(u_cc,u_readmode,v_cc,v_readmode);
  }
  return lcp;
}

static GtUword
gt_seqabstract_lcp_encseq_gtuchar(bool rightextension,
                                  const GtSeqabstract *useq,
                                  const GtSeqabstract *vseq,
                                  GtUword ustart,
                                  GtUword vstart,
                                  GtUword maxlen)
{
  GtUword lcp;
  GT_SEQABSTRACT_ASSIGN_READMODES
  for (lcp = 0; lcp < maxlen; lcp++)
  {
    GtUchar u_cc, v_cc;

    GT_SEQABSTRACT_ENC_ACCESS(u_cc,useq->seq.encseq,useq->offset,u_readmode,
                              ustart,lcp);
    GT_SEQABSTRACT_SEQ_ACCESS(v_cc,vseq->seq.string,v_readmode,vstart,lcp);
    GT_SEQABSTRACT_CHAR_EQUAL_OR_BREAK(u_cc,u_readmode,v_cc,v_readmode);
  }
  return lcp;
}

static GtUword gt_seqabstract_lcp_encseq_encseq(bool rightextension,
                                                const GtSeqabstract *useq,
                                                const GtSeqabstract *vseq,
                                                GtUword ustart,
                                                GtUword vstart,
                                                GtUword maxlen)
{
  GtUword lcp;
  bool is_same_sequence = useq->seq.encseq == vseq->seq.encseq &&
                          useq->offset + ustart == vseq->offset + vstart &&
                          useq->readmode == vseq->readmode;
  GT_SEQABSTRACT_ASSIGN_READMODES

  if (!is_same_sequence)
  {
    GtUchar u_cc, v_cc;

    for (lcp = 0; lcp < maxlen; lcp++)
    {
      GT_SEQABSTRACT_ENC_ACCESS(u_cc,useq->seq.encseq,useq->offset,u_readmode,
                                ustart,lcp);
      GT_SEQABSTRACT_ENC_ACCESS(v_cc,vseq->seq.encseq,vseq->offset,v_readmode,
                                vstart,lcp);
      GT_SEQABSTRACT_CHAR_EQUAL_OR_BREAK(u_cc,u_readmode,v_cc,v_readmode);
    }
  } else
  {
    GtUword startpos = useq->offset + ustart;

    for (lcp = 0; lcp < maxlen; lcp++)
    {
      GtUchar u_cc = gt_encseq_get_encoded_char(useq->seq.encseq,
                                                GT_ISDIRREVERSE(u_readmode)
                                                     ? (startpos - lcp)
                                                     : (startpos + lcp),
                                                GT_READMODE_FORWARD);
      if (ISSPECIAL(u_cc))
      {
        break;
      }
    }
  }
  return lcp;
}

GtUword gt_seqabstract_lcp(bool rightextension,
                           const GtSeqabstract *useq,
                           const GtSeqabstract *vseq,
                           GtUword ustart,
                           GtUword vstart)
{
  GtUword lcp,
          maxlen;

  gt_assert(useq != NULL && vseq != NULL &&
            useq->len >= ustart && vseq->len >= vstart);
  maxlen = rightextension ? MIN(useq->len - ustart, vseq->len - vstart)
                          : 1 + MIN(ustart, vstart);

  if (useq->seqtype == GT_SEQABSTRACT_STRING)
  {
    if (vseq->seqtype == GT_SEQABSTRACT_STRING)
    {
      lcp = gt_seqabstract_lcp_gtuchar_gtuchar(rightextension, useq, vseq,
                                               ustart, vstart, maxlen);
    } else
    {
      lcp = gt_seqabstract_lcp_gtuchar_encseq(rightextension, useq, vseq,
                                              ustart, vstart, maxlen);
    }
  } else
  {
    if (vseq->seqtype == GT_SEQABSTRACT_STRING)
    {
      lcp = gt_seqabstract_lcp_encseq_gtuchar(rightextension, useq, vseq,
                                              ustart, vstart, maxlen);

    } else
    {
      lcp = gt_seqabstract_lcp_encseq_encseq(rightextension, useq, vseq, ustart,
                                             vstart, maxlen);
    }
  }
  return lcp;
}

char *gt_seqabstract_get(bool rightextension,const GtSeqabstract *seq)
{
  GtUword idx;
  char *buffer = malloc(sizeof *buffer * (seq->len+1));
  char *map = "acgt";

  printf("readmode=%s,rightextension=%s,totallength=" GT_WU ",len=" GT_WU
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
