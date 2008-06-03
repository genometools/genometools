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

#include "libgtcore/chardef.h"
#include "libgtcore/symboldef.h"
#include "libgtcore/unused.h"
#include "libgtcore/arraydef.h"
#include "sarr-def.h"
#include "seqpos-def.h"
#include "esa-splititv.h"
#include "spacedef.h"
#include "esa-limdfs.h"

struct Myersonlineresources
{
  Encodedsequencescanstate *esr;
  unsigned long *eqsvectorrev;
  unsigned int alphasize;
};

static void initeqsvectorrev(unsigned long *eqsvector,
                             unsigned long eqslen,
                             const Uchar *u,
                             unsigned long ulen)
{
  unsigned long *vptr, shiftmask;
  const Uchar *uptr;

  for (vptr = eqsvector; vptr < eqsvector + eqslen; vptr++)
  {
    *vptr = 0;
  }
  for (uptr = u, shiftmask = 1UL;
       uptr < u + ulen && shiftmask != 0;
       uptr++, shiftmask <<= 1)
  {
    assert (*uptr != (Uchar) SEPARATOR);
    if (*uptr != (Uchar) WILDCARD)
    {
      eqsvector[(unsigned long) *uptr] |= shiftmask;
    }
  }
}

static void showmatchonline(Seqpos startpos)
{
  printf("match " FormatSeqpos "\n",PRINTSeqposcast(startpos));
}

void edistmyersbitvectorAPM(Myersonlineresources *mor,
                            const Encodedsequence *encseq,
                            const Uchar *pattern,
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
                Ebit,
                score;

  Uchar cc;
  Seqpos pos, totallength;
  const Readmode readmode = Reversemode;

  totallength = getencseqtotallength(encseq);
  initeqsvectorrev(mor->eqsvectorrev,
                   (unsigned long) mor->alphasize,
                   pattern,patternlength);
  Ebit = 1UL << (patternlength-1);
  score = patternlength;
  initEncodedsequencescanstate(mor->esr,
                               encseq,
                               readmode,
                               0);
  for (pos = 0; pos < totallength; pos++)
  {
    cc = sequentialgetencodedchar(encseq,
                                  mor->esr,
                                  pos,
                                  readmode);
    if (cc == (Uchar) SEPARATOR)
    {
      Pv = ~0UL;
      Mv = 0UL;
      score = patternlength;
    } else
    {
      Eq = mor->eqsvectorrev[(unsigned int) cc];      /*  6 */
      Xv = Eq | Mv;                                   /*  7 */
      Xh = (((Eq & Pv) + Pv) ^ Pv) | Eq;              /*  8 */

      Ph = Mv | ~ (Xh | Pv);                          /*  9 */
      Mh = Pv & Xh;                                   /* 10 */

      if (Pv & Ebit)
      {
        score++;
      } else
      {
        if (Mv & Ebit)
        {
          score--;
        }
      }

      Ph <<= 1;                                       /* 15 */
      Pv = (Mh << 1) | ~ (Xv | Ph);                   /* 17 */
      Mv = Ph & Xv;                                   /* 18 */
      if (score <= maxdistance)
      {
        showmatchonline(REVERSEPOS(totallength,pos));
      }
    }
  }
}
