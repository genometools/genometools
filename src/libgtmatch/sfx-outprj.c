/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <limits.h>
#include "libgtcore/strarray.h"
#include "libgtcore/str.h"
#include "seqpos-def.h"
#include "format64.h"
#include "spacedef.h"
#include "filelength-def.h"
#include "esafileend.h"
#include "readmode-def.h"
#include "stamp.h"

#include "endianess.pr"
#include "opensfxfile.pr"

#define PRJSPECIALOUT(VAL)\
        fprintf(outprj,"%s=" FormatSeqpos "\n",#VAL,\
                PRINTSeqposcast(specialcharinfo->VAL))

static void showprjinfo(FILE *outprj,
                        const StrArray *filenametab,
                        Readmode readmode,
                        const Filelengthvalues *filelengthtab,
                        Seqpos totallength,
                        unsigned long numofsequences,
                        const Specialcharinfo *specialcharinfo,
                        uint32_t prefixlength,
                        Seqpos numoflargelcpvalues,
                        Seqpos maxbranchdepth,
                        const DefinedSeqpos *longest)
{
  unsigned long i;

  assert(filelengthtab != NULL);
  assert(filenametab != NULL);
  for (i=0; i<strarray_size(filenametab); i++)
  {
    fprintf(outprj,"dbfile=%s " Formatuint64_t " " Formatuint64_t "\n",
                    strarray_get(filenametab,i),
                    PRINTuint64_tcast(filelengthtab[i].length),
                    PRINTuint64_tcast(filelengthtab[i].effectivelength));
  }
  fprintf(outprj,"totallength=" FormatSeqpos "\n",PRINTSeqposcast(totallength));
  PRJSPECIALOUT(specialcharacters);
  PRJSPECIALOUT(specialranges);
  PRJSPECIALOUT(lengthofspecialprefix);
  PRJSPECIALOUT(lengthofspecialsuffix);
  fprintf(outprj,"numofsequences=%lu\n",numofsequences);
  fprintf(outprj,"numofdbsequences=%lu\n",numofsequences);
  fprintf(outprj,"numofquerysequences=0\n");
  if (longest->defined)
  {
    fprintf(outprj,"longest=" FormatSeqpos "\n",
            PRINTSeqposcast(longest->valueseqpos));
  }
  fprintf(outprj,"prefixlength=%u\n",(unsigned int) prefixlength);
  fprintf(outprj,"largelcpvalues=" FormatSeqpos "\n",
                   PRINTSeqposcast(numoflargelcpvalues));
  fprintf(outprj,"maxbranchdepth=" FormatSeqpos "\n",
                   PRINTSeqposcast(maxbranchdepth));
  fprintf(outprj,"integersize=%u\n",
                  (unsigned int) (sizeof (Seqpos) * CHAR_BIT));
  fprintf(outprj,"littleendian=%c\n",islittleendian() ? '1' : '0');
  fprintf(outprj,"readmode=%u\n",(unsigned int) readmode);
}

int outprjfile(const Str *indexname,
               const StrArray *filenametab,
               Readmode readmode,
               const Filelengthvalues *filelengthtab,
               Seqpos totallength,
               unsigned long numofsequences,
               const Specialcharinfo *specialcharinfo,
               uint32_t prefixlength,
               Seqpos numoflargelcpvalues,
               Seqpos maxbranchdepth,
               const DefinedSeqpos *longest,
               Env *env)
{
  FILE *prjfp;
  bool haserr = false;

  env_error_check(env);
  prjfp = opensfxfile(indexname,PROJECTFILESUFFIX,"wb",env);
  if (prjfp == NULL)
  {
    haserr = true;
  }
  if (!haserr)
  {
    showprjinfo(prjfp,
                filenametab,
                readmode,
                filelengthtab,
                totallength,
                numofsequences,
                specialcharinfo,
                prefixlength,
                numoflargelcpvalues,
                maxbranchdepth,
                longest);
    env_fa_xfclose(prjfp,env);
  }
  return haserr ? -1 : 0;
}
