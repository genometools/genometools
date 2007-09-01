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
#include "symboldef.h"
#include "encseq-def.h"
#include "spacedef.h"
#include "chardef.h"
#include "enum-patt-def.h"

 struct Enumpatternstate
{
  unsigned long minpatternlen,
                maxpatternlen,
                samplecount,
                *patternstat;
  Uchar *patternspace;
  const Encodedsequence *sampleencseq;
  Seqpos totallength;
};

Enumpatternstate *newenumpattern(unsigned long minpatternlen,
                                 unsigned long maxpatternlen,
                                 const Encodedsequence *encseq,
                                 Env *env)
{
  Enumpatternstate *eps = NULL;
  unsigned long i;

  if (maxpatternlen < minpatternlen)
  {
    env_error_set(env,"maxpatternlen=%lu < %lu\n",
                    maxpatternlen,
                    minpatternlen);
    return NULL;
  }
  ALLOCASSIGNSPACE(eps,NULL,Enumpatternstate,1);
  eps->totallength = getencseqtotallength(encseq);
  if (eps->totallength <= (Seqpos) maxpatternlen)
  {
    env_error_set(env,"totallength=" FormatSeqpos " <= maxpatternlen = %lu\n",
                    PRINTSeqposcast(eps->totallength),
                    maxpatternlen);
    FREESPACE(eps);
    return NULL;
  }
  ALLOCASSIGNSPACE(eps->patternspace,NULL,Uchar,maxpatternlen);
  ALLOCASSIGNSPACE(eps->patternstat,NULL,unsigned long,maxpatternlen+1);
  for (i=0; i<maxpatternlen; i++)
  {
    eps->patternstat[i] = 0;
  }
  eps->minpatternlen = minpatternlen;
  eps->maxpatternlen = maxpatternlen;
  eps->sampleencseq = encseq;
  eps->samplecount = 0;
  srand48(42349421);
  return eps;
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

const Uchar *nextsampledpattern(unsigned long *patternlen,
                                Enumpatternstate *eps)
{
  Seqpos start;
  unsigned long j, requiredpatternlen;

  if (eps->minpatternlen == eps->maxpatternlen)
  {
    requiredpatternlen = eps->minpatternlen;
  } else
  {
    requiredpatternlen = (unsigned long) (eps->minpatternlen +
                                          (drand48() *
                                          (double) (eps->maxpatternlen -
                                                    eps->minpatternlen+1)));
  }
  while (true)
  {
    *patternlen = requiredpatternlen;
    start = (Seqpos) (drand48() * (double) (eps->totallength - *patternlen));
    assert(start < (Seqpos) (eps->totallength - *patternlen));
    for (j=0; j<*patternlen; j++)
    {
      eps->patternspace[j] = getencodedchar(eps->sampleencseq,start+j,
                                            Forwardmode);
      if (ISSPECIAL(eps->patternspace[j]))
      {
        *patternlen = j;
        break;
      }
    }
    if (*patternlen > (unsigned long) 1)
    {
      if (eps->samplecount & 1)
      {
        reverseinplace(eps->patternspace,*patternlen);
      }
      eps->samplecount++;
      eps->patternstat[*patternlen]++;
      break;
    }
  }
  return eps->patternspace;
}

void freeenumpattern(Enumpatternstate *eps,Env *env)
{
  FREESPACE(eps->patternspace);
  FREESPACE(eps->patternstat);
  FREESPACE(eps);
}
