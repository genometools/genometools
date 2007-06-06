/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include <string.h>
#include "types.h"
#include "spacedef.h"

/*@notnull@*/ char *composefilenamegeneric(const char *file,
                                           int linenum,
                                           const char *filename,
                                           char sep,
                                           const char *suffix,
                                           Env *env)
{
  Uint i, lenfilename, lensuffix, totalsize;
  char *dest;

  assert(filename != NULL);
  lenfilename = (Uint) strlen(filename);
  assert(suffix != NULL);
  lensuffix = (Uint) strlen(suffix);
  totalsize = lenfilename+lensuffix+1+1;
  ALLOCASSIGNSPACEGENERIC(file,linenum,dest,NULL,char,totalsize);
  assert(dest != NULL);
  for (i=0; i<lenfilename; i++)
  {
    dest[i] = filename[i];
  }
  dest[lenfilename] = sep;
  for (i=0; i<lensuffix; i++)
  {
    dest[lenfilename+1+i] = suffix[i];
  }
  dest[lenfilename+lensuffix+1] = '\0';
  return dest;
}

/*@notnull@*/ char *composefilename(const char *file,
                                    int linenum,
                                    const char *filename,
                                    const char *suffix,
                                    Env *env)
{
  return composefilenamegeneric(file,
                                linenum,
                                filename,
                                '.',
                                suffix,
                                env);
}
