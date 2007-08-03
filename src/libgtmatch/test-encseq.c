/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "libgtcore/env.h"
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
  Uchar ccscan, ccra, ccsr;
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
    if(readmode == Forwardmode)
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
      if(pos >= totallength)
      {
        break;
      }
    }
    ccra = getencodedchar(encseq,pos,readmode);
    if(readmode == Forwardmode)
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
