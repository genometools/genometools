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
  GtEncseqReader *esr;
  GtSeqabstractType seqtype;
  GtReadmode readmode;
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
  sa->readmode = GT_READMODE_FORWARD;
  sa->offset = 0;
  sa->seq.string = NULL;
  return sa;
}

void gt_seqabstract_readmode_set(GtSeqabstract *sa,GtReadmode readmode)
{
  gt_assert(sa != NULL);
  /*printf("seqabstract: readmode set to %s\n",gt_readmode_show(readmode));*/
  sa->readmode = readmode;
}

void gt_seqabstract_reinit_gtuchar(GtSeqabstract *sa,
                                   const GtUchar *string,
                                   GtUword len,
                                   GtUword offset,
                                   GtUword totallength)
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
  sa->totallength = totallength;
  if (GT_ISDIRREVERSE(sa->readmode))
  {
    gt_assert(offset < totallength);
    sa->offset = GT_REVERSEPOS(totallength,offset);
  } else
  {
    sa->offset = offset;
  }
  sa->seq.string = string + sa->offset;
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

  if (useq->cmpcharbychar || vseq->cmpcharbychar || !is_same_sequence)
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
    if (useq->stoppossupport)
    {
      GtUword stoppos;
      const GtUword startpos
         = GT_ISDIRREVERSE(u_readmode)
              ? GT_REVERSEPOS(useq->totallength,useq->offset + ustart)
              : (useq->offset + ustart);

      gt_encseq_reader_reinit_with_readmode(useq->esr,
                                            useq->seq.encseq,
                                            GT_ISDIRREVERSE(u_readmode)
                                              ? GT_READMODE_REVERSE
                                              : GT_READMODE_FORWARD,
                                            startpos);
      stoppos = gt_getnexttwobitencodingstoppos(GT_ISDIRREVERSE(u_readmode)
                                                  ? false
                                                  : true,
                                                useq->esr);
      if (GT_ISDIRREVERSE(u_readmode))
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
  GtUword idx, startpos;
  GtReadmode readmode;
  char *buffer = malloc(sizeof *buffer * (seq->len+1));
  char *map = "acgt";

  if (rightextension)
  {
    readmode = seq->readmode;
  } else
  {
    readmode = gt_readmode_inverse_dir(seq->readmode);
  }
  startpos = GT_ISDIRREVERSE(readmode) ? seq->len - 1 : 0;
  for (idx = 0; idx < seq->len; idx++)
  {
    GtUchar cc;

    if (seq->seqtype == GT_SEQABSTRACT_STRING)
    {
      GT_SEQABSTRACT_SEQ_ACCESS(cc,seq->seq.string,readmode,startpos,idx);
    } else
    {
      GT_SEQABSTRACT_ENC_ACCESS(cc,seq->seq.encseq,seq->offset,readmode,
                                startpos,idx);
    }
    buffer[idx] = map[cc];
  }
  buffer[seq->len] = '\0';
  return buffer;
}
