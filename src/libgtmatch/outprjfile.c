/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <errno.h>
#include "types.h"
#include "spacedef.h"

#include "endianess.pr"
#include "compfilenm.pr"
#include "opensfxfile.pr"

static void showprjinfo(FILE *outprj,
                        const char **filenametab,
                        Uint numoffiles,
                        const PairUint *filelengthtab,
                        Uint64 totallength,
                        Uint numofsequences,
                        const Specialcharinfo *specialcharinfo,
                        unsigned int prefixlength)
{
  Uint i;

  assert(numoffiles > 0);
  assert(filelengthtab != NULL);
  assert(filenametab != NULL);
  for (i=0; i<numoffiles; i++)
  {
    fprintf(outprj,"dbfile=%s %lu %lu\n",filenametab[i],
                                         (Showuint) filelengthtab[i].uint0,
                                         (Showuint) filelengthtab[i].uint1);
  }
  /*@ignore@*/
  fprintf(outprj,"totallength=" FormatUint64 "\n",totallength);
  /*@end@*/
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
  fprintf(outprj,"integersize=%ld\n",(Showsint) (sizeof (Uint) * CHAR_BIT));
  fprintf(outprj,"littleendian=%c\n",islittleendian() ? '1' : '0');
}

int outprjfile(const char *indexname,
               const char **filenamelist,
               unsigned int numoffiles,
               const PairUint *filelengthtab,
               Uint64 totallength,
               Uint numofsequences,
               const Specialcharinfo *specialcharinfo,
               unsigned int prefixlength,
               Env *env)
{
  FILE *prjfp;
  bool haserr = false;

  env_error_check(env);
  prjfp = opensfxfile(indexname,"prj",env);
  if (prjfp == NULL)
  {
    haserr = true;
  }
  if (!haserr)
  {
    showprjinfo(prjfp,
                filenamelist,
                numoffiles,
                filelengthtab,
                totallength,
                numofsequences,
                specialcharinfo,
                prefixlength);
    env_fa_xfclose(prjfp,env);
  }
  return haserr ? -1 : 0;
}
