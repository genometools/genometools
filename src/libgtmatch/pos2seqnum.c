#include "libgtcore/env.h"
#include "encseq-def.h"
#include "arraydef.h"
#include "chardef.h"
#include "divmodmul.h"

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

Seqpos *calculatemarkpositions(const Encodedsequence *encseq,
                               unsigned long numofdbsequences,
                               Env *env)
{
  ArraySeqpos asp;
  Specialrangeiterator *sri;
  Sequencerange range;
  bool haserr = false;

  assert(numofdbsequences > (unsigned long) 1);
  asp.allocatedSeqpos = numofdbsequences-1;
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

unsigned long getrecordnum(const Seqpos *recordseps,
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
    env_error_set(env,"getrecordnum: cannot find position " FormatSeqpos,
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
  env_error_set(env,"getrecordnum: cannot find position " FormatSeqpos,
                PRINTSeqposcast(position));
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

    markpos = calculatemarkpositions(encseq,
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
        seqnum = getrecordnum(markpos,
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
