/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "libgtcore/env.h"
#include "libgtcore/array.h"
#include "spacedef.h"
#include "encseq-def.h"
#include "fbs-def.h"
#include "stamp.h"

#include "sfx-readmode.pr"
#include "fbsadv.pr"
#include "readnextUchar.gen"

int testencodedsequence(const StrArray *filenametab,
                        const Encodedsequence *encseq,
                        Readmode readmode,
                        const Uchar *symbolmap,
                        Env *env)
{
  Seqpos pos, totallength;
  Uchar ccscan = '\0', ccra, ccsr;
  Fastabufferstate fbs;
  int retval;
  bool haserr = false;
  Encodedsequencescanstate *esr;

  printf("# testencodedsequence with readmode = %s\n",showreadmode(readmode));
  env_error_check(env);
  totallength = getencseqtotallength(encseq);
  initformatbufferstate(&fbs,
                        filenametab,
                        symbolmap,
                        false,
                        NULL,
                        NULL,
                        env);
  esr = initEncodedsequencescanstate(encseq,readmode,env);
  for (pos=0; /* Nothing */; pos++)
  {
    if (readmode == Forwardmode)
    {
      retval = readnextUchar(&ccscan,&fbs,env);
      if (retval < 0)
      {
        haserr = true;
        break;
      }
      if (retval == 0)
      {
        break;
      }
    } else
    {
      if (pos >= totallength)
      {
        break;
      }
    }
    ccra = getencodedchar(encseq,pos,readmode);
    if (readmode == Forwardmode)
    {
      if (ccscan != ccra)
      {
        env_error_set(env,"position " FormatSeqpos
                          ": scan (readnextchar) = %u != "
                          "%u = random access (getencodedchar)",
                           pos,
                           (uint32_t) ccscan,
                           (uint32_t) ccra);
        haserr = true;
        break;
      }
    }
    ccsr = sequentialgetencodedchar(encseq,esr,pos);
    if (ccra != ccsr)
    {
      env_error_set(env,"position " FormatSeqpos
                        ": random access (getencodedchar) = %u != "
                        " %u = sequential read (sequentialgetencodedchar)",
                         pos,
                         (uint32_t) ccra,
                         (uint32_t) ccsr);
      haserr = true;
      break;
    }
  }
  if (!haserr)
  {
    if (pos != totallength)
    {
      env_error_set(env,"sequence length must be " FormatSeqpos " but is "
                         FormatSeqpos,totallength,pos);
      haserr = true;
    }
  }
  freeEncodedsequencescanstate(&esr,env);
  return haserr ? -1 : 0;
}

static int addelem(void *processinfo,const Sequencerange *range,Env *env)
{
  env_error_check(env);
  array_add_elem((Array *) processinfo,(void *) range,sizeof (Sequencerange),
                 env);
  return 0;
}

static void reverseSequencerange(Array *a)
{
  unsigned long idx1, idx2;
  Sequencerange tmp, *valptr1, *valptr2;

  for (idx1=0, idx2 = (unsigned long) array_size(a) - 1;
       idx1 < idx2; idx1++, idx2--)
  {
    valptr1 = (Sequencerange *) array_get(a,idx1);
    valptr2 = (Sequencerange *) array_get(a,idx2);
    tmp = *valptr1;
    *valptr1 = *valptr2;
    *valptr2 = tmp;
  }
}

int checkspecialranges(const Encodedsequence *encseq,Env *env)
{
  Array *rangesforward, *rangesbackward;
  Sequencerange *valf, *valb;
  unsigned long idx;
  bool haserr = false;

  env_error_check(env);
  if (!fastspecialranges(encseq))
  {
    return 0;
  }
  rangesforward = array_new(sizeof (Sequencerange),env);
  rangesbackward = array_new(sizeof (Sequencerange),env);

  if (overallspecialrangesfast(encseq,true,Forwardmode,addelem,rangesforward,
                              env) != 0)
  {
    haserr = true;
  }
  if (!haserr)
  {
    if (overallspecialrangesfast(encseq,false,Reversemode,addelem,
                                 rangesbackward,env) != 0)
    {
      haserr = true;
    }
  }
  if (!haserr)
  {
    if (array_size(rangesforward) != array_size(rangesbackward))
    {
      env_error_set(env,"array_size(rangesforward) = %lu != %lu = "
                        "array_size(rangesbackward)\n",
                        (unsigned long) array_size(rangesforward),
                        (unsigned long) array_size(rangesbackward));
      haserr = true;
    }
  }
  if (!haserr)
  {
    reverseSequencerange(rangesbackward);
    printf("# check %lu ranges\n",(unsigned long) array_size(rangesforward));
    for (idx=0; idx<(unsigned long) array_size(rangesforward); idx++)
    {
      valf = (Sequencerange *) array_get(rangesforward,idx);
      valb = (Sequencerange *) array_get(rangesbackward,idx);
      if (valf->leftpos != valb->leftpos || valf->rightpos != valb->rightpos)
      {
        env_error_set(env,
                      "rangesforward[%lu] = (" FormatSeqpos "," FormatSeqpos
                      ") != (" FormatSeqpos "," FormatSeqpos
                      ") = rangesbackward[%lu]\n",
                      idx,
                      PRINTSeqposcast(valf->leftpos),
                      PRINTSeqposcast(valf->rightpos),
                      PRINTSeqposcast(valb->leftpos),
                      PRINTSeqposcast(valb->rightpos),
                      idx);
        haserr = true;
        break;
      }
    }
  }
  array_delete(rangesforward,env);
  array_delete(rangesbackward,env);
  return haserr ? - 1 : 0;
}
