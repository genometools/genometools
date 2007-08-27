#include <stdlib.h>
#include "symboldef.h"
#include "encseq-def.h"
#include "spacedef.h"
#include "chardef.h"

typedef struct
{
  unsigned long minpatternlen,
                maxpatternlen,
                samplecount,
                *patternstat;
  Uchar *patternspace;
  const Encodedsequence *sampleencseq;
  Seqpos totallength; 
} Enumpatternstate;

Enumpatternstate *newenumpattern(unsigned long minpatternlen,
                                 unsigned long maxpatternlen,
                                 const Encodedsequence *encseq,
                                 Env *env)
{
  Enumpatternstate *eps = NULL;
  unsigned long i;

  if(maxpatternlen < minpatternlen)
  {
    env_error_set(env,"maxpatternlen=%lu < %lu\n",
                    maxpatternlen,
                    minpatternlen);
    return NULL;
  }
  ALLOCASSIGNSPACE(eps,NULL,Enumpatternstate,1);
  eps->totallength = getencseqtotallength(encseq);
  if(eps->totallength <= (Seqpos) maxpatternlen)
  {
    env_error_set(env,"totallength=" FormatSeqpos " <= maxpatternlen = %lu\n",
                    PRINTSeqposcast(eps->totallength),
                    maxpatternlen);
    FREESPACE(eps);
    return NULL;
  }
  ALLOCASSIGNSPACE(eps->patternspace,NULL,Uchar,maxpatternlen);
  ALLOCASSIGNSPACE(eps->patternstat,NULL,unsigned long,maxpatternlen+1);
  for(i=0; i<maxpatternlen; i++)
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

  for(front = s, back = s + len - 1; front < back; front++, back--)
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
  unsigned long j;
  bool special;

  if(eps->minpatternlen == eps->maxpatternlen)
  {
    *patternlen = eps->minpatternlen;
  } else
  {
    *patternlen = (unsigned long) (eps->minpatternlen + 
                                   (drand48() * 
                                   (double) (eps->maxpatternlen - 
                                             eps->minpatternlen+1)));
  }
  eps->patternstat[*patternlen]++;
  while(true)
  {
    start = (Seqpos) (drand48() * (double) (eps->totallength - *patternlen));
    assert(start < (Seqpos) (eps->totallength - *patternlen));
    special = false;
    for(j=0; j<*patternlen; j++)
    {
      eps->patternspace[j] = getencodedchar(eps->sampleencseq,start+j,
                                            Forwardmode);
      if(ISSPECIAL(eps->patternspace[j]))
      {
        special = true;
        break;
      }
    }
    if(!special)
    {
      if(eps->samplecount & 1)
      {
        reverseinplace(eps->patternspace,*patternlen);
      }
      eps->samplecount++;
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
