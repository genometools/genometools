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
    fprintf(outprj,"dbfile=%s %lu %lu\n",strarray_get(filenametab,i),
                                         (Showuint) filelengthtab[i].uint0,
                                         (Showuint) filelengthtab[i].uint1);
  }
  STAMP;
  fprintf(outprj,"totallength=" FormatSeqpos "\n",PRINTSeqposcast(totallength));
  fprintf(outprj,"specialcharacters=%lu\n",
                  (Showuint) specialcharinfo->specialcharacters);
  fprintf(outprj,"specialranges=%lu\n",
                 (Showuint) specialcharinfo->specialranges);
  fprintf(outprj,"lengthofspecialprefix=%lu\n",
                 (Showuint) specialcharinfo->lengthofspecialprefix);
  fprintf(outprj,"lengthofspecialsuffix=%lu\n",
                 (Showuint) specialcharinfo->lengthofspecialsuffix);
  fprintf(outprj,"numofsequences=%lu\n",(Showuint) numofsequences);
  fprintf(outprj,"numofdbsequences=%lu\n",(Showuint) numofsequences);
  fprintf(outprj,"numofquerysequences=0\n");
  fprintf(outprj,"prefixlength=%lu\n",(Showuint) prefixlength);
  fprintf(outprj,"integersize=%lu\n",(Showuint) (sizeof (Seqpos) * CHAR_BIT));
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
