#include "libgtcore/env.h"
#include "encseq-def.h"
#include "arraydef.h"
#include "chardef.h"
#include "divmodmul.h"

DECLAREARRAYSTRUCT(Seqpos);

static int addmarkpos(void *info,const Encodedsequence *encseq,
                      const Sequencerange *seqrange,
                      /*@unused@*/ Env *env)
{
  Seqpos pos;
  Uchar currentchar;
  ArraySeqpos *asp = (ArraySeqpos *) info;

  for (pos=seqrange->leftpos; pos<=seqrange->rightpos; pos++)
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

  assert(numofdbsequences > 0);
  if (numofdbsequences == (unsigned long) 1)
  {
    return NULL;
  }
  asp.allocatedSeqpos = numofdbsequences-1;
  asp.nextfreeSeqpos = 0;
  ALLOCASSIGNSPACE(asp.spaceSeqpos,NULL,Seqpos,asp.allocatedSeqpos);
  if (overallspecialranges(encseq,
                          Forwardmode,
                          addmarkpos,
                          &asp,
                          env) != 0)
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
    if (recordseps[mid] < position)
    {
      if (position < recordseps[mid+1])
      {
        return mid - 1;
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
