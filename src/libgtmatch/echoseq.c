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

#include <assert.h>
#include "libgtcore/env.h"
#include "spacedef.h"
#include "iterseq.h"
#include "alphadef.h"

unsigned long *calcdescendpositions(const char *destab,
                                    unsigned long destablength,
                                    unsigned long numofsequences,  
                                    Env *env)
{
  unsigned long *descendtab, i, idx = 0;

  ALLOCASSIGNSPACE(descendtab,NULL,unsigned long,numofsequences);
  assert(destab != NULL);
  printf("desclength=%lu\n",destablength);
  for(i=0; i<destablength; i++)
  {
    if(destab[i] == '\n')
    {
      assert(idx < numofsequences - 1);
      descendtab[idx++] = i;
    }
  }
  assert(idx == numofsequences - 1);
  descendtab[idx] = destablength;
  return descendtab;
}

const char *retriesequencedescription(unsigned long *desclen,
                                      const char *destab,
                                      const unsigned long *descendtab,
                                      unsigned long seqnum)
{
  if(seqnum == 0)
  {
    *desclen = descendtab[0];
    return destab;
  }
  assert(descendtab[seqnum-1] < descendtab[seqnum]);
  *desclen = descendtab[seqnum] - descendtab[seqnum-1] - 1;
  return destab + descendtab[seqnum-1] + 1;
}

void checkalldescriptions(const char *destab,unsigned long destablength,
                          unsigned long numofsequences,Env *env)
{
  unsigned long *descendtab, desclen, seqnum, totaldesclength, offset = 0;
  const char *desptr;
  char *copydestab;

  descendtab = calcdescendpositions(destab,
                                    destablength,
                                    numofsequences,
                                    env);
  totaldesclength = numofsequences; /* for each new line */
  for(seqnum = 0; seqnum < numofsequences; seqnum++)
  {
    desptr = retriesequencedescription(&desclen,
                                       destab,
                                       descendtab,
                                       seqnum);
    totaldesclength += desclen;
  }
  ALLOCASSIGNSPACE(copydestab,NULL,char,totaldesclength);
  for(seqnum = 0; seqnum < numofsequences; seqnum++)
  {
    desptr = retriesequencedescription(&desclen,
                                       destab,
                                       descendtab,
                                       seqnum);
    strncpy(copydestab + offset,desptr,(size_t) desclen);
    copydestab[offset+desclen] = '\n';
    offset += (desclen+1);
  }
  if(strncmp(copydestab,destab,(size_t) totaldesclength) != 0)
  {
    fprintf(stderr,"different descriptions\n");
    exit(EXIT_FAILURE);
  }
  FREESPACE(copydestab);
  FREESPACE(descendtab);
}

int echodescriptionandsequence(const StrArray *filenametab,Env *env)
{
  Scansequenceiterator *sseqit;
  char *desc = NULL;
  const Uchar *sequence;
  unsigned long seqlen;
  bool haserr = false;
  int retval;

  sseqit = newScansequenceiterator(filenametab,NULL,true,env);
  while (true)
  {
    retval = nextScansequenceiterator(&sequence,
                                      &seqlen,
                                      &desc,
                                      sseqit,
                                      env);
    if (retval < 0)
    {
      haserr = true;
      break;
    }
    if (retval == 0)
    {
      break;
    }
    fastasymbolstringgeneric(stdout,desc,NULL,sequence,seqlen,
                             (unsigned long) 70);
    FREESPACE(desc);
  }
  freeScansequenceiterator(&sseqit,env);
  return haserr ? -1 : 0;
}
