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
#include "libgtcore/env.h"
#include "encseq-def.h"
#include "arraydef.h"
#include "divmodmul.h"
#include "spacedef.h"

DECLAREARRAYSTRUCT(Seqpos);

static int addmarkpos(ArraySeqpos *asp,
                      const Encodedsequence *encseq,
                      const Sequencerange *seqrange)
{
  Seqpos pos;
  Uchar currentchar;

  for (pos=seqrange->leftpos; pos<seqrange->rightpos; pos++)
  {
    currentchar = getencodedchar(encseq,pos,Forwardmode);
    assert(ISSPECIAL(currentchar));
    if (currentchar == (Uchar) SEPARATOR)
    {
      assert(asp->nextfreeSeqpos < asp->allocatedSeqpos);
      asp->spaceSeqpos[asp->nextfreeSeqpos++] = pos;
    }
  }
  return 0;
}

Seqpos *encseq2markpositions(const Encodedsequence *encseq,
                             unsigned long numofsequences,
                             Env *env)
{
  ArraySeqpos asp;
  Specialrangeiterator *sri;
  Sequencerange range;
  bool haserr = false;

  assert (numofsequences > (unsigned long) 1);
  asp.allocatedSeqpos = numofsequences-1;
  asp.nextfreeSeqpos = 0;
  ALLOCASSIGNSPACE(asp.spaceSeqpos,NULL,Seqpos,asp.allocatedSeqpos);
  sri = newspecialrangeiterator(encseq,true,env);
  while (nextspecialrangeiterator(&range,sri))
  {
    if (addmarkpos(&asp,encseq,&range) != 0)
    {
      haserr = true;
      break;
    }
  }
  freespecialrangeiterator(&sri,env);
  if (haserr)
  {
    FREEARRAY(&asp,Seqpos);
    return NULL;
  }
  return asp.spaceSeqpos;
}

unsigned long *sequence2markpositions(unsigned long *numofsequences,
                                      const Uchar *seq,
                                      unsigned long seqlen,
                                      Env *env)
{
  unsigned long *spacemarkpos, i, allocatedmarkpos, nextfreemarkpos;

  *numofsequences = (unsigned long) 1;
  for (i=0; i<seqlen; i++)
  {
    if (seq[i] == (Uchar) SEPARATOR)
    {
      (*numofsequences)++;
    }
  }
  if (*numofsequences == (unsigned long) 1)
  {
    return NULL;
  }
  allocatedmarkpos = (*numofsequences)-1;
  ALLOCASSIGNSPACE(spacemarkpos,NULL,unsigned long,allocatedmarkpos);
  for (i=0, nextfreemarkpos = 0; i<seqlen; i++)
  {
    if (seq[i] == (Uchar) SEPARATOR)
    {
      spacemarkpos[nextfreemarkpos++] = i;
    }
  }
  return spacemarkpos;
}

unsigned long getrecordnumSeqpos(const Seqpos *recordseps,
                                 unsigned long numofrecords,
                                 Seqpos totalwidth,
                                 Seqpos position,
                                 Env *env)
{
  unsigned long left, mid, right, len;

  assert(numofrecords > 0);
  if (numofrecords == (unsigned long) 1 || position < recordseps[0])
  {
    return 0;
  }
  if (position > recordseps[numofrecords-2])
  {
    if (position < totalwidth)
    {
      return numofrecords - 1;
    }
    env_error_set(env,"getrecordnumSeqpos: cannot find position " FormatSeqpos,
                  PRINTSeqposcast(position));
    return numofrecords; /* failure */
  }
  left = 0;
  right = numofrecords - 2;
  while (left<=right)
  {
    len = (unsigned long) (right-left);
    mid = left + DIV2(len);
#ifdef DEBUG
    printf("left=%lu,right = %lu\n",left,right);
    printf("mid=%lu\n",mid);
#endif
    if (recordseps[mid] < position)
    {
      if (position < recordseps[mid+1])
      {
        return mid + 1;
      }
      left = mid + 1;
    } else
    {
      if (recordseps[mid-1] < position)
      {
        return mid;
      }
      right = mid-1;
    }
  }
  env_error_set(env,"getrecordnumSeqpos: cannot find position " FormatSeqpos,
                PRINTSeqposcast(position));
  return numofrecords; /* failure */
}

unsigned long getrecordnumulong(const unsigned long *recordseps,
                                unsigned long numofrecords,
                                unsigned long totalwidth,
                                unsigned long position,
                                Env *env)
{
  unsigned long left, mid, right, len;

  assert(numofrecords > 0);
  if (numofrecords == (unsigned long) 1 || position < recordseps[0])
  {
    return 0;
  }
  if (position > recordseps[numofrecords-2])
  {
    if (position < totalwidth)
    {
      return numofrecords - 1;
    }
    env_error_set(env,"getrecordnumulong: cannot find position %lu",position);
    return numofrecords; /* failure */
  }
  left = 0;
  right = numofrecords - 2;
  while (left<=right)
  {
    len = (unsigned long) (right-left);
    mid = left + DIV2(len);
#ifdef DEBUG
    printf("left=%lu,right = %lu\n",left,right);
    printf("mid=%lu\n",mid);
#endif
    if (recordseps[mid] < position)
    {
      if (position < recordseps[mid+1])
      {
        return mid + 1;
      }
      left = mid + 1;
    } else
    {
      if (recordseps[mid-1] < position)
      {
        return mid;
      }
      right = mid-1;
    }
  }
  env_error_set(env,"getrecordnumulong: cannot find position %lu",position);
  return numofrecords; /* failure */
}

int checkmarkpos(const Encodedsequence *encseq,
                 unsigned long numofdbsequences,
                 Env *env)
{
  if (numofdbsequences > (unsigned long) 1)
  {
    Seqpos *markpos, totallength, pos;
    unsigned long currentseqnum = 0, seqnum;
    Uchar currentchar;

    markpos = encseq2markpositions(encseq,
                                   numofdbsequences,
                                   env);
    if (markpos == NULL)
    {
      return -1;
    }
    totallength = getencseqtotallength(encseq);
    for (pos=0; pos<totallength; pos++)
    {
      currentchar = getencodedchar(encseq,pos,Forwardmode);
      if (currentchar == (Uchar) SEPARATOR)
      {
        currentseqnum++;
      } else
      {
        seqnum = getrecordnumSeqpos(markpos,
                                    numofdbsequences,
                                    totallength,
                                    pos,
                                    env);
        if (seqnum == numofdbsequences)
        {
          return -1;
        }
        if (seqnum != currentseqnum)
        {
          fprintf(stderr,"pos= " FormatSeqpos
                         " seqnum = %lu != %lu = currentseqnum\n",
                          PRINTSeqposcast(pos),seqnum,currentseqnum);
          exit(EXIT_FAILURE); /* programming error */
        }
      }
    }
    FREESPACE(markpos);
  }
  return 0;
}
