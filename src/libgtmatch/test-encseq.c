/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "libgtcore/env.h"
#include "types.h"
#include "spacedef.h"
#include "encseq-def.h"
#include "fbs-def.h"
#include "stamp.h"

#include "fbsadv.pr"
#include "readnextUchar.gen"

int testencodedsequence(const StrArray *filenametab,
                        const Encodedsequence *encseq,
                        const Uchar *symbolmap,
                        Env *env)
{
  Seqpos pos;
  Uchar cc0, cc1;
  Fastabufferstate fbs;
  int retval;
  PairSeqpos *filelengthtab = NULL;
  bool haserr = false;
  Encodedsequencescanstate *esr;

  env_error_check(env);
  STAMP;
  initfastabufferstate(&fbs,
                       filenametab,
                       symbolmap,
                       &filelengthtab,
                       env);
  STAMP;
  esr = initEncodedsequencescanstate(encseq,env);
  for (pos=0; /* Nothing */; pos++)
  {
    retval = readnextUchar(&cc0,&fbs,env);
    if (retval < 0)
    {
      haserr = true;
      break;
    }
    if (retval == 0)
    {
      break;
    }
    cc1 = getencodedchar(encseq,pos);
    if (cc0 != cc1)
    {
      /*@ignore@*/
      env_error_set(env,"position " FormatSeqpos
                        ": correct = %u != %u = cc1 (getencodedchar64)",
                         pos,
                         (unsigned int) cc0,
                         (unsigned int) cc1);
      /*@end@*/
      haserr = true;
      break;
    }
    cc1 = sequentialgetencodedchar(encseq,esr,pos);
    if (cc0 != cc1)
    {
      /*@ignore@*/
      env_error_set(env,"position " FormatSeqpos
                        ": correct = %u != %u = cc1 "
                        "(sequentialgetencodedchar)",
                         pos,
                         (unsigned int) cc0,
                         (unsigned int) cc1);
      /*@end@*/
      haserr = true;
      break;
    }
  }
  if (!haserr)
  {
    if (pos != getencseqtotallength(encseq))
    {
      env_error_set(env,"sequence length must be " FormatSeqpos " but is "
                         FormatSeqpos,
                         getencseqtotallength(encseq),pos);
      haserr = true;
    }
  }
  freeEncodedsequencescanstate(&esr,env);
  FREESPACE(filelengthtab);
  return haserr ? -1 : 0;
}
