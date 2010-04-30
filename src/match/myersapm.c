/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

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

#include "core/chardef.h"
#include "core/symboldef.h"

#include "spacedef.h"
#include "core/encseq.h"
#include "myersapm.h"
#include "core/defined-types.h"
#include "procmatch.h"
#include "dist-short.h"
#include "initeqsvec.h"

struct Myersonlineresources
{
  GtEncseqReader *esr;
  const GtEncseq *encseq;
  unsigned long totallength;
  unsigned long *eqsvectorrev,
                *eqsvector;
  unsigned int alphasize;
  bool nowildcards;
  Processmatch processmatch;
  void *processmatchinfo;
};

Myersonlineresources *gt_newMyersonlineresources(unsigned int numofchars,
                                                 bool nowildcards,
                                                 const GtEncseq *encseq,
                                                 Processmatch processmatch,
                                                 void *processmatchinfo)
{
  Myersonlineresources *mor;

  ALLOCASSIGNSPACE(mor,NULL,Myersonlineresources,1);
  ALLOCASSIGNSPACE(mor->eqsvectorrev,NULL,unsigned long,numofchars);
  ALLOCASSIGNSPACE(mor->eqsvector,NULL,unsigned long,numofchars);
  mor->encseq = encseq;
  mor->esr = gt_encseq_create_reader_with_readmode(encseq, GT_READMODE_REVERSE,
                                                   0);
  gt_assert(numofchars <= GT_MAXALPHABETCHARACTER);
  mor->alphasize = numofchars;
  mor->totallength = gt_encseq_total_length(encseq);
  mor->nowildcards = nowildcards;
  mor->processmatch = processmatch;
  mor->processmatchinfo = processmatchinfo;
  return mor;
}

void gt_freeMyersonlineresources(Myersonlineresources **ptrmyersonlineresources)
{
  Myersonlineresources *mor = *ptrmyersonlineresources;

  FREESPACE(mor->eqsvectorrev);
  FREESPACE(mor->eqsvector);
  gt_encseq_reader_delete(mor->esr);
  FREESPACE(*ptrmyersonlineresources);
}

void gt_edistmyersbitvectorAPM(Myersonlineresources *mor,
                               const GtUchar *pattern,
                               unsigned long patternlength,
                               unsigned long maxdistance)
{
  unsigned long Pv = ~0UL,
                Mv = 0UL,
                Eq,
                Xv,
                Xh,
                Ph,
                Mh,
                score;
  const unsigned long Ebit = 1UL << (patternlength-1);
  GtUchar cc;
  unsigned long pos;
  const GtReadmode readmode = GT_READMODE_REVERSE;
  GtMatch match;

  gt_initeqsvectorrev(mor->eqsvectorrev,
                   (unsigned long) mor->alphasize,
                   pattern,patternlength);
  score = patternlength;
  gt_encseq_reader_reinit_with_readmode(mor->esr, mor->encseq, readmode, 0);
  match.dbabsolute = NULL;
  match.dbsubstring = NULL;
  match.querystartpos = 0;
  match.querylen = patternlength;
  match.alignment = NULL;
  for (pos = 0; pos < mor->totallength; pos++)
  {
    cc = gt_encseq_reader_next_encoded_char(mor->esr);
    if (cc == (GtUchar) SEPARATOR)
    {
      Pv = ~0UL;
      Mv = 0UL;
      score = patternlength;
    } else
    {
      if (cc == (GtUchar) WILDCARD)
      {
        Eq = 0;
      } else
      {
        Eq = mor->eqsvectorrev[(unsigned long) cc];   /*  6 */
      }
      Xv = Eq | Mv;                                   /*  7 */
      Xh = (((Eq & Pv) + Pv) ^ Pv) | Eq;              /*  8 */

      Ph = Mv | ~ (Xh | Pv);                          /*  9 */
      Mh = Pv & Xh;                                   /* 10 */

      if (Ph & Ebit)
      {
        score++;
      } else
      {
        if (Mh & Ebit)
        {
          gt_assert(score > 0);
          score--;
        }
      }

      Ph <<= 1;                                       /* 15 */
      Pv = (Mh << 1) | ~ (Xv | Ph);                   /* 17 */
      Mv = Ph & Xv;                                   /* 18 */
      if (score <= maxdistance)
      {
        unsigned long dbstartpos = GT_REVERSEPOS(mor->totallength,pos);
        Definedunsignedlong matchlength;

        if (maxdistance > 0)
        {
          matchlength = gt_forwardprefixmatch(mor->encseq,
                                           mor->alphasize,
                                           dbstartpos,
                                           mor->nowildcards,
                                           mor->eqsvector,
                                           pattern,
                                           patternlength,
                                           maxdistance);
        } else
        {
          matchlength.defined = true;
          matchlength.valueunsignedlong = patternlength;
        }
        gt_assert(matchlength.defined || mor->nowildcards);
        if (matchlength.defined)
        {
          match.dbstartpos = dbstartpos;
          match.dblen = (unsigned long) matchlength.valueunsignedlong;
          match.distance = score;
          mor->processmatch(mor->processmatchinfo,&match);
        }
      }
    }
  }
}
