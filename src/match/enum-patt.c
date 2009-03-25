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

#include <stdlib.h>
#include "core/chardef.h"
#include "core/symboldef.h"
#include "encseq-def.h"
#include "spacedef.h"
#include "enum-patt-def.h"

 struct Enumpatterniterator
{
  unsigned long minpatternlen,
                maxpatternlen,
                samplecount,
                *patternstat;
  Uchar *patternspace;
  const Encodedsequence *sampleencseq;
  unsigned int alphasize;
  Seqpos totallength;
  Encodedsequencescanstate *esr;
};

Enumpatterniterator *newenumpatterniterator(unsigned long minpatternlen,
                                            unsigned long maxpatternlen,
                                            const Encodedsequence *encseq,
                                            GtError *err)
{
  Enumpatterniterator *epi = NULL;
  unsigned long i;

  if (maxpatternlen < minpatternlen)
  {
    gt_error_set(err,"maxpatternlen=%lu < %lu\n",
                    maxpatternlen,
                    minpatternlen);
    return NULL;
  }
  ALLOCASSIGNSPACE(epi,NULL,Enumpatterniterator,1);
  epi->totallength = getencseqtotallength(encseq);
  if (epi->totallength <= (Seqpos) maxpatternlen)
  {
    gt_error_set(err,"totallength=" FormatSeqpos " <= maxpatternlen = %lu\n",
                    PRINTSeqposcast(epi->totallength),
                    maxpatternlen);
    FREESPACE(epi);
    return NULL;
  }
  ALLOCASSIGNSPACE(epi->patternspace,NULL,Uchar,maxpatternlen);
  ALLOCASSIGNSPACE(epi->patternstat,NULL,unsigned long,maxpatternlen+1);
  for (i=0; i<=maxpatternlen; i++)
  {
    epi->patternstat[i] = 0;
  }
  epi->minpatternlen = minpatternlen;
  epi->maxpatternlen = maxpatternlen;
  epi->sampleencseq = encseq;
  epi->samplecount = 0;
  epi->alphasize = getencseqAlphabetnumofchars(encseq);
  epi->esr = newEncodedsequencescanstate();
  srand48(42349421);
  return epi;
}

static void reversesequenceinplace(Uchar *s,unsigned long len)
{
  Uchar *front, *back, tmp;

  for (front = s, back = s + len - 1; front < back; front++, back--)
  {
    tmp = *front;
    *front = *back;
    *back = tmp;
  }
}

const Uchar *nextEnumpatterniterator(unsigned long *patternlen,
                                     Enumpatterniterator *epi)
{
  Seqpos start;
  unsigned long j;
  Uchar cc;

  if (epi->minpatternlen == epi->maxpatternlen)
  {
    *patternlen = epi->minpatternlen;
  } else
  {
    *patternlen = (unsigned long) (epi->minpatternlen +
                                   (drand48() *
                                   (double) (epi->maxpatternlen -
                                             epi->minpatternlen+1)));
  }
  start = (Seqpos) (drand48() * (double) (epi->totallength - *patternlen));
  gt_assert(start < (Seqpos) (epi->totallength - *patternlen));
  initEncodedsequencescanstate(epi->esr,epi->sampleencseq,Forwardmode,start);
  for (j=0; j<*patternlen; j++)
  {
    cc = sequentialgetencodedchar(epi->sampleencseq,epi->esr,start+j,
                                  Forwardmode);
    if (ISSPECIAL(cc))
    {
      cc = (Uchar) (drand48() * epi->alphasize);
    }
    epi->patternspace[j] = cc;
  }
  if (epi->samplecount & 1)
  {
    reversesequenceinplace(epi->patternspace,*patternlen);
  }
  epi->samplecount++;
  epi->patternstat[*patternlen]++;
  return epi->patternspace;
}

void showPatterndistribution(const Enumpatterniterator *epi)
{
  unsigned long i;
  double addprob, probsum = 0.0;

  printf("# %lu pattern with the following length distribution:\n",
         epi->samplecount);
  for (i=epi->minpatternlen; i<=epi->maxpatternlen; i++)
  {
    if (epi->patternstat[i] > 0)
    {
      addprob = (double) epi->patternstat[i] / epi->samplecount;
      probsum += addprob;
      printf("# %lu: %lu (prob=%.4f,cumulative=%.4f)\n",
             i,
             epi->patternstat[i],
             addprob,
             probsum);
    }
  }
}

void freeEnumpatterniterator(Enumpatterniterator **epi)
{
  if (!(*epi)) return;
  FREESPACE((*epi)->patternspace);
  FREESPACE((*epi)->patternstat);
  freeEncodedsequencescanstate(&((*epi)->esr));
  FREESPACE(*epi);
}
