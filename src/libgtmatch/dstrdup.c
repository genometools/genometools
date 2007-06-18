/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include "spacedef.h"

/*EE
  The following function makes a copy of a 0-terminated string pointed to by
  \texttt{source}.
*/

/*@notnull@*/ char *dynamicstrdup(const char *file,int linenum,
                                  const char *source,Env *env)
{
  size_t sourcelength;
  char *dest;

  assert(source != NULL);
  sourcelength = strlen(source);
  ALLOCASSIGNSPACEGENERIC(file,linenum,dest,NULL,char,sourcelength+1);
  strcpy(dest,source);
  return dest;
}
