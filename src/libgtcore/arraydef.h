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

#include <assert.h>
#include "libgtcore/ma.h"
#include "libgtcore/symboldef.h"

/*
  This file defines macros to conveniently declare and
  manipulate dynamic arrays whose size grow on demand. Each dynamic
  array over some type T
  is implemented by a structure consisting of three components:
  - space##T is a pointer to the space block of type T allocated for the array.
  - allocated##T is an unsigned long value storing the number of entries in the
    array currently allocated.
  - nextfree##T holds the smallest index of the array where no value is stored.
  Here ## is the concatenation operator of the C-preprocessor.
  The following macro expands to a corresponding type definition over
  some given TYPE.
*/

#define DECLAREARRAYSTRUCT(TYPE)\
        typedef struct\
        {\
          TYPE *space##TYPE;\
          unsigned long allocated##TYPE, nextfree##TYPE;\
        } Array##TYPE

/*
  INITARRAY initializes an empty array.
*/

#define INITARRAY(A,TYPE)\
        (A)->space##TYPE = NULL;\
        (A)->allocated##TYPE = (A)->nextfree##TYPE = 0

/*
  COPYARRAY copies an array.
*/

#define COPYARRAY(A,B)\
        *(A) = *(B)

/*
  CHECKARRAYSPACE checks if the integer nextfree##T
  points to an index for which the space is not allocated yet. If this is
  the case, the number of cells allocated is incremented by L. The
  contents of the previously filled array elements is of course maintained.
*/

#define CHECKARRAYSPACE(A,TYPE,L)\
        if ((A)->nextfree##TYPE >= (A)->allocated##TYPE)\
        {\
          (A)->allocated##TYPE += L;\
          (A)->space##TYPE = ma_realloc_mem((A)->space##TYPE,\
                                            sizeof (TYPE) *\
                                            (A)->allocated##TYPE,\
                                            __FILE__, __LINE__);\
        }

/*
  The next macro is a variation of CHECKARRAYSPACE, which checks if the next
  L cells have been allocated. If not, then this is done.
*/

#define CHECKARRAYSPACEMULTI(A,TYPE,L)\
        if ((A)->nextfree##TYPE + (L) >= (A)->allocated##TYPE)\
        {\
          (A)->allocated##TYPE += L;\
          (A)->space##TYPE = ma_realloc_mem((A)->space##TYPE,\
                                            sizeof (TYPE) *\
                                            (A)->allocated##TYPE,\
                                            __FILE__, __LINE__);\
        }

/*
  This macro checks the space and delivers a pointer P
  to the next free element in the array.
*/

#define GETNEXTFREEINARRAY(P,A,TYPE,L)\
        CHECKARRAYSPACE(A,TYPE,L);\
        P = (A)->space##TYPE + (A)->nextfree##TYPE++;

/*
  This macro checks the space and stores V in the
  nextfree-component of the array. nextfree is incremented.
*/

#define STOREINARRAY(A,TYPE,L,VAL)\
        CHECKARRAYSPACE(A,TYPE,L);\
        assert((A)->space##TYPE != NULL);\
        (A)->space##TYPE[(A)->nextfree##TYPE++] = VAL

/*
  This macro frees the space for an array if it is not NULL.
*/

#define FREEARRAY(A,TYPE)\
        if ((A)->space##TYPE != NULL)\
        {\
          ma_free((A)->space##TYPE);\
        }

/*
  Some declarations for the most common array types.
*/

DECLAREARRAYSTRUCT(Uchar);
DECLAREARRAYSTRUCT(char);

#endif
