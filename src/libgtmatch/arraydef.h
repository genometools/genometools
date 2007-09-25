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

#ifndef ARRAYDEF_H
#define ARRAYDEF_H

#include "libgtcore/symboldef.h"
#include "spacedef.h"

/*
  This file defines macros to conveniently declare and
  manipulate dynamic arrays whose size grow on demand. Each dynamic
  array over some type \texttt{T}
  is implemented by a structure consisting of three components:
  \begin{enumerate}
  \item
  \texttt{space\#\#T} is a pointer to the space block of type \texttt{T}
  allocated for the array.
  \item
  \texttt{allocated\#\#T} is an \texttt{unsigned long}-value storing the number
  of entries in the array currently allocated.
  \item
  \texttt{nextfree\#\#T} holds the smallest index of the array where no
  value is stored.
  \end{enumerate}
  Here \texttt{\#\#} is the concatenation operator of the C-preprocessor.
  The following macro expands to a corresponding type definition over
  some given \texttt{TYPE}.
*/

#define DECLAREARRAYSTRUCT(TYPE)\
        typedef struct\
        {\
          TYPE *space##TYPE;\
          unsigned long allocated##TYPE, nextfree##TYPE;\
        } Array##TYPE

/*
  \texttt{INITARRAY} initializes an empty array.
*/

#define INITARRAY(A,TYPE)\
        (A)->space##TYPE = NULL;\
        (A)->allocated##TYPE = (A)->nextfree##TYPE = 0

/*
  \texttt{COPYARRAY} copies an array.
*/

#define COPYARRAY(A,B)\
        *(A) = *(B)

/*
  \texttt{CHECKARRAYSPACE} checks if the integer \texttt{nextfree\#\#T}
  points to an index for which the space is not allocated yet. If this is
  the case, the number of cells allocated is incremented by \texttt{L}. The
  contents of the previously filled array elements is of course maintained.
*/

#define CHECKARRAYSPACE(A,TYPE,L)\
        if ((A)->nextfree##TYPE >= (A)->allocated##TYPE)\
        {\
          (A)->allocated##TYPE += L;\
          ALLOCASSIGNSPACEGENERIC(__FILE__,__LINE__,(A)->space##TYPE,\
                                  (A)->space##TYPE,\
                                  TYPE,(A)->allocated##TYPE);\
        }

/*
  The next macro is a variation of \texttt{CHECKARRAYSPACE},
  which checks if the next
  \texttt{L} cells have been allocated. If not, then this is done.
*/

#define CHECKARRAYSPACEMULTI(A,TYPE,L)\
        if ((A)->nextfree##TYPE + (L) >= (A)->allocated##TYPE)\
        {\
          (A)->allocated##TYPE += L;\
          ALLOCASSIGNSPACEGENERIC(__FILE__,__LINE__,(A)->space##TYPE,\
                                  (A)->space##TYPE,\
                                  TYPE,(A)->allocated##TYPE);\
        }

/*
  This macro checks the space and delivers a pointer \texttt{P}
  to the next free element in the array.
*/

#define GETNEXTFREEINARRAY(P,A,TYPE,L)\
        CHECKARRAYSPACE(A,TYPE,L);\
        P = (A)->space##TYPE + (A)->nextfree##TYPE++;

/*
  This macro checks the space and stores \texttt{V} in the
  \texttt{nextfree}-component of the array. \texttt{nextfree}
  is incremented.
*/

#define STOREINARRAY(A,TYPE,L,VAL)\
        CHECKARRAYSPACE(A,TYPE,L);\
        (A)->space##TYPE[(A)->nextfree##TYPE++] = VAL

/*
  This macro frees the space for an array if it is not \texttt{NULL}.
*/

#define FREEARRAY(A,TYPE)\
        if ((A)->space##TYPE != NULL)\
        {\
          FREESPACE((A)->space##TYPE);\
        }

/*
  Some declarations for the most common array types.
*/

DECLAREARRAYSTRUCT(Uchar);
DECLAREARRAYSTRUCT(char);

#endif
