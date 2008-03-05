/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

  Permission to use, copy, modify, and distribute this software for any
  purpose with or without fee is hereby granted, provided that the above
  copyright notice and this permission notice appear in all copies.

  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/

#ifndef SPACEDEF_H
#define SPACEDEF_H
#include <string.h>
#include <assert.h>
#include "libgtcore/error.h"
#include "libgtcore/ma.h"

#ifdef NOSPACEBOOKKEEPING

#define ALLOCASSIGNSPACEGENERIC(FILENAME,LINENUM,V,S,T,N)\
        assert(sizeof(*(V)) == sizeof (T));\
        V = (T *) realloc(S,sizeof (T) * (size_t) (N));\
        if ((V) == NULL)\
        {\
          fprintf(stderr,"file %s, line %d: realloc(%lu) failed\n",\
                  FILENAME,LINENUM,\
                  (unsigned long) (sizeof (T) * (size_t) (N)));\
          exit(EXIT_FAILURE); /* malloc error */ \
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
        assert(sizeof(*(V)) == sizeof (T));\
        V = ma_realloc_mem(S, sizeof (T) * (N), FILENAME,\
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
          ma_free(P);\
          P = NULL;\
        }

#endif /* NOSPACEBOOKKEEPING */

/*
  The remaining macros call the corresponding function with
  the filename and the line number where the function call
  appears.
*/

#define ASSIGNDYNAMICSTRDUP(V,S)\
        V = dynamicstrdup(__FILE__,__LINE__,S,err)

#endif
