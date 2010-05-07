/*
  Copyright (c) 2008 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2008 Center for Bioinformatics, University of Hamburg

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
#include "core/types_api.h"
#include "core/chardef.h"
#include "core/encseq.h"
#include "core/defined-types.h"
#include "initeqsvec.h"
#include "dist-short.h"

#define DECLARELOCALVARS\
        unsigned long Pv = ~0UL,\
                      Mv = 0,\
                      Eq,\
                      Xv,\
                      Xh,\
                      Ph,\
                      Mh,\
                      Ebit = (1UL << (ulen-1)),\
                      distval = ulen

#define COMPUTENEWDIST(CC)\
        gt_assert((CC) != (GtUchar) SEPARATOR);\
        if ((CC) != (GtUchar) WILDCARD)\
        {\
          Eq = eqsvector[(unsigned long) (CC)];\
        } else\
        {\
          Eq = 0;\
        }\
        Xv = Eq | Mv;\
        Xh = (((Eq & Pv) + Pv) ^ Pv) | Eq;\
        Ph = Mv | ~ (Xh | Pv);\
        Mh = Pv & Xh;\
        if (Ph & Ebit)\
        {\
          distval++;\
        } else\
        {\
          if (Mh & Ebit)\
          {\
            distval--;\
          }\
        }\
        /*\
           the original version of Myers included the statement\
           Ph <<= 1. As an effect, the first element of each column\
           was 0. We add a 1 to the shifted Ph vector. As a consequence,\
           we obtain an increment by 1 in the first column.\
        */\
        Ph = (Ph << 1) | 1UL;\
        Pv = (Mh << 1) | ~ (Xv | Ph);\
        Mv = Ph & Xv

unsigned long gt_distanceofshortstringsbytearray(unsigned long *eqsvector,
                                              unsigned int alphasize,
                                              const GtUchar *useq,
                                              unsigned long ulen,
                                              const GtUchar *vseq,
                                              unsigned long vlen)
{
  DECLARELOCALVARS;
  const GtUchar *vptr;

  gt_initeqsvector(eqsvector,(unsigned long) alphasize,useq,ulen);
  for (vptr = vseq; vptr < vseq + vlen; vptr++)
  {
    COMPUTENEWDIST(*vptr);
  }
  return distval;
}

unsigned long gt_distanceofshortstringsencseq(unsigned long *eqsvector,
                                           unsigned int alphasize,
                                           const GtUchar *useq,
                                           unsigned long ulen,
                                           const GtEncseq *encseq,
                                           unsigned long vstartpos,
                                           unsigned long vlen)
{
  DECLARELOCALVARS;
  GtUchar cc;
  unsigned long pos;

  gt_initeqsvector(eqsvector,(unsigned long) alphasize,useq,ulen);
  for (pos = vstartpos; pos < vstartpos + vlen; pos++)
  {
    cc = gt_encseq_get_encoded_char(encseq,pos,GT_READMODE_FORWARD);
    COMPUTENEWDIST(cc);
  }
  return distval;
}

unsigned long gt_reversesuffixmatch(unsigned long *eqsvector,
                                 unsigned int alphasize,
                                 const GtUchar *useq,
                                 unsigned long ulen,
                                 const GtUchar *vseq,
                                 unsigned long vlen,
                                 unsigned long maxdistance)
{
  DECLARELOCALVARS;
  const GtUchar *vptr;

  gt_initeqsvectorrev(eqsvector,(unsigned long) alphasize,useq,ulen);
  gt_assert(maxdistance > 0);
  for (vptr = vseq + vlen - 1; vptr >= vseq; vptr--)
  {
    COMPUTENEWDIST(*vptr);
    if (distval <= maxdistance)
    {
      break;
    }
  }
  /* gt_assert(distval <= maxdistance); */
  return (unsigned long) (vseq + vlen - vptr);
}

Definedunsignedlong gt_forwardprefixmatch(const GtEncseq *encseq,
                                       unsigned int alphasize,
                                       unsigned long startpos,
                                       bool nowildcards,
                                       unsigned long *eqsvector,
                                       const GtUchar *useq,
                                       unsigned long ulen,
                                       unsigned long maxdistance)
{
  DECLARELOCALVARS;
  unsigned long pos, totallength = gt_encseq_total_length(encseq);
  GtUchar cc;
  Definedunsignedlong result;

  gt_initeqsvector(eqsvector,(unsigned long) alphasize,useq,ulen);
  gt_assert(maxdistance > 0);
  for (pos = startpos; /* Nothing */; pos++)
  {
    gt_assert(pos - startpos <= (unsigned long) (ulen + maxdistance));
    cc = gt_encseq_get_encoded_char(encseq,pos,GT_READMODE_FORWARD);
    if (nowildcards && cc == (GtUchar) WILDCARD)
    {
      result.defined = false;
      result.valueunsignedlong = 0;
      return result;
    }
    COMPUTENEWDIST(cc);
    if (distval <= maxdistance || pos == totallength-1)
    {
      break;
    }
  }
  result.defined = true;
  result.valueunsignedlong = (unsigned long) (pos - startpos + 1);
  return result;
}
