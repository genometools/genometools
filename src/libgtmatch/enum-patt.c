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
#include "libgtcore/chardef.h"
#include "libgtcore/symboldef.h"
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
  Seqpos totallength;
};

Enumpatterniterator *newenumpatterniterator(unsigned long minpatternlen,
                                            unsigned long maxpatternlen,
                                            const Encodedsequence *encseq,
                                            Env *env)
{
  Enumpatterniterator *epi = NULL;
  unsigned long i;

  if (maxpatternlen < minpatternlen)
  {
    env_error_set(env,"maxpatternlen=%lu < %lu\n",
                    maxpatternlen,
                    minpatternlen);
    return NULL;
  }
  ALLOCASSIGNSPACE(epi,NULL,Enumpatterniterator,1);
  epi->totallength = getencseqtotallength(encseq);
  if (epi->totallength <= (Seqpos) maxpatternlen)
  {
    env_error_set(env,"totallength=" FormatSeqpos " <= maxpatternlen = %lu\n",
                    PRINTSeqposcast(epi->totallength),
                    maxpatternlen);
    FREESPACE(epi);
    return NULL;
  }
  ALLOCASSIGNSPACE(epi->patternspace,NULL,Uchar,maxpatternlen);
  ALLOCASSIGNSPACE(epi->patternstat,NULL,unsigned long,maxpatternlen+1);
  for (i=0; i<maxpatternlen; i++)
  {
    epi->patternstat[i] = 0;
  }
  epi->minpatternlen = minpatternlen;
  epi->maxpatternlen = maxpatternlen;
  epi->sampleencseq = encseq;
  epi->samplecount = 0;
  srand48(42349421);
  return epi;
}

static void reverseinplace(Uchar *s,unsigned long len)
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
  unsigned long j, requiredpatternlen;

  if (epi->minpatternlen == epi->maxpatternlen)
  {
    requiredpatternlen = epi->minpatternlen;
  } else
  {
    requiredpatternlen = (unsigned long) (epi->minpatternlen +
                                          (drand48() *
                                          (double) (epi->maxpatternlen -
                                                    epi->minpatternlen+1)));
  }
  while (true)
  {
    *patternlen = requiredpatternlen;
    start = (Seqpos) (drand48() * (double) (epi->totallength - *patternlen));
    assert(start < (Seqpos) (epi->totallength - *patternlen));
    for (j=0; j<*patternlen; j++)
    {
      epi->patternspace[j] = getencodedchar(epi->sampleencseq,start+j,
                                            Forwardmode);
      if (ISSPECIAL(epi->patternspace[j]))
      {
        *patternlen = j;
        break;
      }
    }
    if (*patternlen > (unsigned long) 1)
    {
      if (epi->samplecount & 1)
      {
        reverseinplace(epi->patternspace,*patternlen);
      }
      epi->samplecount++;
      epi->patternstat[*patternlen]++;
      break;
    }
  }
  return epi->patternspace;
}

void freeEnumpatterniterator(Enumpatterniterator **epi,Env *env)
{
  if (!(*epi)) return;
  FREESPACE((*epi)->patternspace);
  FREESPACE((*epi)->patternstat);
  FREESPACE(*epi);
}
