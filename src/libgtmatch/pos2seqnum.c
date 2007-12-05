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

#include "libgtcore/arraydef.h"
#include "libgtcore/chardef.h"
#include "libgtcore/error.h"
#include "encseq-def.h"
#include "divmodmul.h"
#include "spacedef.h"

DECLAREARRAYSTRUCT(Seqpos);

static int addmarkpos(ArraySeqpos *asp,
                      const Encodedsequence *encseq,
                      Encodedsequencescanstate *esr,
                      const Sequencerange *seqrange)
{
  Seqpos pos;
  Uchar currentchar;

  initEncodedsequencescanstate(esr,encseq,Forwardmode,seqrange->leftpos);
  for (pos=seqrange->leftpos; pos<seqrange->rightpos; pos++)
  {
    currentchar = sequentialgetencodedchar(encseq,esr,pos);
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
                             unsigned long numofsequences)
{
  ArraySeqpos asp;
  Specialrangeiterator *sri;
  Sequencerange range;
  bool haserr = false;
  Encodedsequencescanstate *esr;

  assert (numofsequences > 1UL);
  asp.allocatedSeqpos = numofsequences-1;
  asp.nextfreeSeqpos = 0;
  ALLOCASSIGNSPACE(asp.spaceSeqpos,NULL,Seqpos,asp.allocatedSeqpos);
  sri = newspecialrangeiterator(encseq,true);
  esr = newEncodedsequencescanstate();
  while (nextspecialrangeiterator(&range,sri))
  {
    if (addmarkpos(&asp,encseq,esr,&range) != 0)
    {
      haserr = true;
      break;
    }
  }
  freespecialrangeiterator(&sri);
  freeEncodedsequencescanstate(&esr);
  if (haserr)
  {
    FREEARRAY(&asp,Seqpos);
    return NULL;
  }
  return asp.spaceSeqpos;
}

unsigned long *sequence2markpositions(unsigned long *numofsequences,
                                      const Uchar *seq,
                                      unsigned long seqlen)
{
  unsigned long *spacemarkpos, i, allocatedmarkpos, nextfreemarkpos;

  *numofsequences = 1UL;
  for (i=0; i<seqlen; i++)
  {
    if (seq[i] == (Uchar) SEPARATOR)
    {
      (*numofsequences)++;
    }
  }
  if (*numofsequences == 1UL)
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
                                 Error *err)
{
  unsigned long left, mid, right, len;

  assert(numofrecords > 0);
  if (numofrecords == 1UL || position < recordseps[0])
  {
    return 0;
  }
  if (position > recordseps[numofrecords-2])
  {
    if (position < totalwidth)
    {
      return numofrecords - 1;
    }
    error_set(err,"getrecordnumSeqpos: cannot find position " FormatSeqpos,
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
  error_set(err,"getrecordnumSeqpos: cannot find position " FormatSeqpos,
                PRINTSeqposcast(position));
  return numofrecords; /* failure */
}

unsigned long getrecordnumulong(const unsigned long *recordseps,
                                unsigned long numofrecords,
                                unsigned long totalwidth,
                                unsigned long position,
                                Error *err)
{
  unsigned long left, mid, right, len;

  assert(numofrecords > 0);
  if (numofrecords == 1UL || position < recordseps[0])
  {
    return 0;
  }
  if (position > recordseps[numofrecords-2])
  {
    if (position < totalwidth)
    {
      return numofrecords - 1;
    }
    error_set(err,"getrecordnumulong: cannot find position %lu",position);
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
  error_set(err,"getrecordnumulong: cannot find position %lu",position);
  return numofrecords; /* failure */
}

int checkmarkpos(const Encodedsequence *encseq,
                 unsigned long numofdbsequences,
                 Error *err)
{
  bool haserr = false;

  if (numofdbsequences > 1UL)
  {
    Seqpos *markpos, totallength, pos;
    unsigned long currentseqnum = 0, seqnum;
    Uchar currentchar;
    Encodedsequencescanstate *esr;

    markpos = encseq2markpositions(encseq,numofdbsequences);
    if (markpos == NULL)
    {
      return -1;
    }
    totallength = getencseqtotallength(encseq);
    esr = newEncodedsequencescanstate();
    initEncodedsequencescanstate(esr,encseq,Forwardmode,0);
    for (pos=0; pos<totallength; pos++)
    {
      currentchar = sequentialgetencodedchar(encseq,esr,pos);
      if (currentchar == (Uchar) SEPARATOR)
      {
        currentseqnum++;
      } else
      {
        seqnum = getrecordnumSeqpos(markpos,
                                    numofdbsequences,
                                    totallength,
                                    pos,
                                    err);
        if (seqnum == numofdbsequences)
        {
          haserr = true;
          break;
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
    freeEncodedsequencescanstate(&esr);
    FREESPACE(markpos);
  }
  return haserr ? -1 : 0;
}
