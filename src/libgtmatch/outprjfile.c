/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <errno.h>
#include "libgtcore/strarray.h"
#include "libgtcore/str.h"
#include "types.h"
#include "spacedef.h"
#include "stamp.h"

#include "endianess.pr"
#include "opensfxfile.pr"

#define PRJSPECIALOUT(VAL)\
        fprintf(outprj,"%s=" FormatSeqpos "\n",#VAL,\
                PRINTSeqposcast(specialcharinfo->VAL))

static void showprjinfo(FILE *outprj,
                        const StrArray *filenametab,
                        const PairSeqpos *filelengthtab,
                        /*@unused@*/ Seqpos totallength,
                        unsigned long numofsequences,
                        const Specialcharinfo *specialcharinfo,
                        uint32_t prefixlength)
{
  unsigned long i;

  assert(filelengthtab != NULL);
  assert(filenametab != NULL);
  STAMP;
  for (i=0; i<strarray_size(filenametab); i++)
  {
    fprintf(outprj,"dbfile=%s " FormatSeqpos " " FormatSeqpos "\n",
                    strarray_get(filenametab,i),
                    PRINTSeqposcast(filelengthtab[i].uint0),
                    PRINTSeqposcast(filelengthtab[i].uint1));
  }
  STAMP;
  fprintf(outprj,"totallength=" FormatSeqpos "\n",
                 PRINTSeqposcast(totallength));
  PRJSPECIALOUT(specialcharacters);
  PRJSPECIALOUT(specialranges);
  PRJSPECIALOUT(lengthofspecialprefix);
  PRJSPECIALOUT(lengthofspecialsuffix);
  fprintf(outprj,"numofsequences=%lu\n",numofsequences);
  fprintf(outprj,"numofdbsequences=%lu\n",numofsequences);
  fprintf(outprj,"numofquerysequences=0\n");
  fprintf(outprj,"prefixlength=%u\n",(unsigned int) prefixlength);
  fprintf(outprj,"integersize=%u\n",
                  (unsigned int) (sizeof (Seqpos) * CHAR_BIT));
  fprintf(outprj,"littleendian=%c\n",islittleendian() ? '1' : '0');
}

int outprjfile(const Str *indexname,
               const StrArray *filenametab,
               const PairSeqpos *filelengthtab,
               Seqpos totallength,
               unsigned long numofsequences,
               const Specialcharinfo *specialcharinfo,
               uint32_t prefixlength,
               Env *env)
{
  FILE *prjfp;
  bool haserr = false;

  env_error_check(env);
  prjfp = opensfxfile(indexname,".prj",env);
  if (prjfp == NULL)
  {
    haserr = true;
  }
  if (!haserr)
  {
    showprjinfo(prjfp,
                filenametab,
                filelengthtab,
                totallength,
                numofsequences,
                specialcharinfo,
                prefixlength);
    env_fa_xfclose(prjfp,env);
  }
  return haserr ? -1 : 0;
}
