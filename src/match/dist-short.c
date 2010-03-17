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
#include "core/symboldef.h"
#include "core/chardef.h"
#include "encodedsequence.h"
#include "defined-types.h"
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

unsigned long distanceofshortstringsbytearray(unsigned long *eqsvector,
                                              unsigned int alphasize,
                                              const GtUchar *useq,
                                              unsigned long ulen,
                                              const GtUchar *vseq,
                                              unsigned long vlen)
{
  DECLARELOCALVARS;
  const GtUchar *vptr;

  initeqsvector(eqsvector,(unsigned long) alphasize,useq,ulen);
  for (vptr = vseq; vptr < vseq + vlen; vptr++)
  {
    COMPUTENEWDIST(*vptr);
  }
  return distval;
}

unsigned long distanceofshortstringsencseq(unsigned long *eqsvector,
                                           unsigned int alphasize,
                                           const GtUchar *useq,
                                           unsigned long ulen,
                                           const GtEncodedsequence *encseq,
                                           Seqpos vstartpos,
                                           Seqpos vlen)
{
  DECLARELOCALVARS;
  GtUchar cc;
  Seqpos pos;

  initeqsvector(eqsvector,(unsigned long) alphasize,useq,ulen);
  for (pos = vstartpos; pos < vstartpos + vlen; pos++)
  {
    cc = gt_encodedsequence_getencodedchar(encseq,pos,GT_READMODE_FORWARD);
    COMPUTENEWDIST(cc);
  }
  return distval;
}

unsigned long reversesuffixmatch(unsigned long *eqsvector,
                                 unsigned int alphasize,
                                 const GtUchar *useq,
                                 unsigned long ulen,
                                 const GtUchar *vseq,
                                 unsigned long vlen,
                                 unsigned long maxdistance)
{
  DECLARELOCALVARS;
  const GtUchar *vptr;

  initeqsvectorrev(eqsvector,(unsigned long) alphasize,useq,ulen);
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

Definedunsignedlong forwardprefixmatch(const GtEncodedsequence *encseq,
                                       unsigned int alphasize,
                                       Seqpos startpos,
                                       bool nowildcards,
                                       unsigned long *eqsvector,
                                       const GtUchar *useq,
                                       unsigned long ulen,
                                       unsigned long maxdistance)
{
  DECLARELOCALVARS;
  Seqpos pos, totallength = gt_encodedsequence_total_length(encseq);
  GtUchar cc;
  Definedunsignedlong result;

  initeqsvector(eqsvector,(unsigned long) alphasize,useq,ulen);
  gt_assert(maxdistance > 0);
  for (pos = startpos; /* Nothing */; pos++)
  {
    gt_assert(pos - startpos <= (Seqpos) (ulen + maxdistance));
    cc = gt_encodedsequence_getencodedchar(encseq,pos,GT_READMODE_FORWARD);
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
