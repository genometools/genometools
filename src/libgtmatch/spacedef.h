/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef SPACEDEF_H
#define SPACEDEF_H
#include <string.h>
#include "libgtcore/env.h"

#ifdef NOSPACEBOOKKEEPING

#define ALLOCASSIGNSPACEGENERIC(FILENAME,LINENUM,V,S,T,N)\
        V = (T *) realloc(S,sizeof (T) * (size_t) (N));\
        if ((V) == NULL)\
        {\
          fprintf(stderr,"file %s, line %d: realloc(%lu) failed\n",\
                  FILENAME,LINENUM,\
                  (unsigned long) (sizeof (T) * (size_t) (N)));\
          exit(EXIT_FAILURE);\
        }

#define ALLOCASSIGNSPACE(V,S,T,N)\
        ALLOCASSIGNSPACEGENERIC(__FILE__,__LINE__,V,S,T,N)\

#define FREESPACE(P)\
        if ((P) != NULL)\
        {\
          free(P);\
          P = NULL;\
        }

#else
#define ALLOCASSIGNSPACEGENERIC(FILENAME,LINENUM,V,S,T,N)\
        V = ma_realloc_mem(env_ma(env), S, sizeof (T) * (N), FILENAME,\
                           LINENUM)

#define ALLOCASSIGNSPACE(V,S,T,N)\
        ALLOCASSIGNSPACEGENERIC(__FILE__,__LINE__,V,S,T,N)

/*
  The macro \texttt{FREESPACE} frees the space pointed to by \texttt{P},
  if this is not \texttt{NULL}. It also sets the
  pointer to \texttt{NULL}.
*/

#define FREESPACE(P)\
        if ((P) != NULL)\
        {\
          env_ma_free(P,env);\
          P = NULL;\
        }

#endif /* NOSPACEBOOKKEEPING */

/*
  The remaining macros call the corresponding function with
  the filename and the line number where the function call
  appears.
*/

/* Some obsolete macros

#define ASSIGNDYNAMICSTRDUP(V,S)\
        V = dynamicstrdup(__FILE__,__LINE__,S,env)

#define COMPOSEFILENAME(FILENAME,SUFFIX)\
        composefilename(__FILE__,__LINE__,FILENAME,SUFFIX,env)

#define COMPOSEFILENAMEGENERIC(FILENAME,SEP,SUFFIX)\
        composefilenamegeneric(__FILE__,__LINE__,FILENAME,SEP,SUFFIX,env)

*/

#endif
