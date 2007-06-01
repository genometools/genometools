#include <assert.h>
#include "types.h"
#include "spacedef.h"

/*EE
  The following function makes a copy of a 0-terminated string pointed to by
  \texttt{source}.
*/

/*@notnull@*/ char *dynamicstrdup(const char *file,Uint linenum,
                                  const char *source,Env *env)
{
  Uint sourcelength;
  char *dest;

  assert(source != NULL);
  sourcelength = (Uint) strlen(source);
  ALLOCASSIGNSPACEGENERIC(file,linenum,dest,NULL,char,sourcelength+1);
  strcpy(dest,source);
  return dest;
}
