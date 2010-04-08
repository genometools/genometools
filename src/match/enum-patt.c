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
#include "core/encodedsequence.h"
#include "spacedef.h"
#include "enum-patt-def.h"

 struct Enumpatterniterator
{
  unsigned long minpatternlen,
                maxpatternlen,
                samplecount,
                *patternstat;
  GtUchar *patternspace;
  const GtEncodedsequence *sampleencseq;
  unsigned int alphasize;
  unsigned long totallength;
  GtEncodedsequenceScanstate *esr;
};

Enumpatterniterator *gt_newenumpatterniterator(unsigned long minpatternlen,
                                            unsigned long maxpatternlen,
                                            const GtEncodedsequence *encseq,
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
  epi->totallength = gt_encodedsequence_totallength(encseq);
  if (epi->totallength <= (unsigned long) maxpatternlen)
  {
    gt_error_set(err,"totallength=%lu <= maxpatternlen = %lu\n",
                    epi->totallength,
                    maxpatternlen);
    FREESPACE(epi);
    return NULL;
  }
  ALLOCASSIGNSPACE(epi->patternspace,NULL,GtUchar,maxpatternlen);
  ALLOCASSIGNSPACE(epi->patternstat,NULL,unsigned long,maxpatternlen+1);
  for (i=0; i<=maxpatternlen; i++)
  {
    epi->patternstat[i] = 0;
  }
  epi->minpatternlen = minpatternlen;
  epi->maxpatternlen = maxpatternlen;
  epi->sampleencseq = encseq;
  epi->samplecount = 0;
  epi->alphasize = gt_alphabet_num_of_chars(
                                           gt_encodedsequence_alphabet(encseq));
  epi->esr = gt_encodedsequence_scanstate_new();
  return epi;
}

static void reversesequenceinplace(GtUchar *s,unsigned long len)
{
  GtUchar *front, *back, tmp;

  for (front = s, back = s + len - 1; front < back; front++, back--)
  {
    tmp = *front;
    *front = *back;
    *back = tmp;
  }
}

const GtUchar *gt_nextEnumpatterniterator(unsigned long *patternlen,
                                       Enumpatterniterator *epi)
{
  unsigned long start;
  unsigned long j;
  GtUchar cc;

  if (epi->minpatternlen == epi->maxpatternlen)
  {
    *patternlen = epi->minpatternlen;
  } else
  {
    *patternlen = (unsigned long) (epi->minpatternlen +
                                   (random() %
                                      (epi->maxpatternlen -
                                       epi->minpatternlen+1)));
  }
  start =
        (unsigned long) (random() % (epi->totallength - *patternlen));
  gt_assert(start < (unsigned long) (epi->totallength - *patternlen));
  gt_encodedsequence_scanstate_init(epi->esr,
                                    epi->sampleencseq,
                                    GT_READMODE_FORWARD,
                                    start);
  for (j=0; j<*patternlen; j++)
  {
    cc = gt_encodedsequence_sequentialgetencodedchar(epi->sampleencseq,
                                                     epi->esr,
                                                     start+j,
                                                     GT_READMODE_FORWARD);
    if (ISSPECIAL(cc))
    {
      cc = (GtUchar) (random() % epi->alphasize);
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

void gt_showPatterndistribution(const Enumpatterniterator *epi)
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

void gt_freeEnumpatterniterator(Enumpatterniterator **epi)
{
  if (!(*epi)) return;
  FREESPACE((*epi)->patternspace);
  FREESPACE((*epi)->patternstat);
  gt_encodedsequence_scanstate_delete(((*epi)->esr));
  FREESPACE(*epi);
}
