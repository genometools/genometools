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

#include <inttypes.h>
#include "core/assert_api.h"
#include "core/ma.h"
#include "core/types_api.h"

/*
  This file defines macros to conveniently declare and
  manipulate dynamic arrays whose size grow on demand. Each dynamic
  array over some type T
  is implemented by a structure consisting of three components:
  - space##T is a pointer to the space block of type T allocated for the array.
  - allocated##T is an GtUword value storing the number of entries in the
    array currently allocated.
  - nextfree##T holds the smallest index of the array where no value is stored.
  Here ## is the concatenation operator of the C-preprocessor.
  The following macro expands to a corresponding type definition over
  some given TYPE.
*/
#define GT_DECLAREARRAYSTRUCT(TYPE)     \
        typedef struct                  \
        {                               \
          TYPE *space##TYPE;            \
          GtUword allocated##TYPE,\
          nextfree##TYPE;               \
        } GtArray##TYPE

/*
  GT_INITARRAY initializes an empty array.
*/
#define GT_INITARRAY(A,TYPE)\
        do {\
          (A)->space##TYPE = NULL;\
          (A)->allocated##TYPE = (A)->nextfree##TYPE = 0;\
        } while (false)

/*
  GT_COPYARRAY copies an array.
*/
#define GT_COPYARRAY(A,B)\
        *(A) = *(B)

/*
  GT_CHECKARRAYSPACE checks if the integer nextfree##T
  points to an index for which the space is not allocated yet. If this is
  the case, the number of cells allocated is incremented by L. The
  contents of the previously filled array elements is of course maintained.
*/
#define GT_CHECKARRAYSPACE(A,TYPE,L)\
        do {\
          if ((A)->nextfree##TYPE >= (A)->allocated##TYPE)\
          {\
            (A)->allocated##TYPE += L;\
            (A)->space##TYPE = (TYPE *) gt_realloc_mem((A)->space##TYPE,\
                                              sizeof (TYPE) *\
                                              (A)->allocated##TYPE,\
                                              __FILE__, __LINE__);\
          }\
          gt_assert((A)->space##TYPE != NULL);\
        } while (false)

/*
  The next macro is a variation of GT_CHECKARRAYSPACE, which checks if the next
  L cells have been allocated. If not, then this is done.
*/
#define GT_CHECKARRAYSPACEMULTI(A,TYPE,L)\
        do {\
          if ((A)->nextfree##TYPE + (L) >= (A)->allocated##TYPE)\
          {\
            (A)->allocated##TYPE += L;\
            (A)->space##TYPE = (TYPE *) gt_realloc_mem((A)->space##TYPE,\
                                              sizeof (TYPE) *\
                                              (A)->allocated##TYPE,\
                                              __FILE__, __LINE__);\
          }\
          gt_assert((A)->space##TYPE != NULL);\
        } while (false)

/*
  This macro checks the space and delivers a pointer P
  to the next free element in the array.
*/
#define GT_GETNEXTFREEINARRAY(P,A,TYPE,L)\
        do {\
          GT_CHECKARRAYSPACE(A,TYPE,L);\
          P = (A)->space##TYPE + (A)->nextfree##TYPE++;\
        } while (false)

/*
  This macro checks the space and stores V in the
  nextfree-component of the array. nextfree is incremented.
*/
#define GT_STOREINARRAY(A,TYPE,L,VAL)\
        do {\
          GT_CHECKARRAYSPACE(A,TYPE,L);\
          (A)->space##TYPE[(A)->nextfree##TYPE++] = VAL;\
        } while (false)

/*
  This macro frees the space for an array if it is not NULL.
*/
#define GT_FREEARRAY(A,TYPE)\
        do {\
          if ((A)->space##TYPE != NULL)\
          {\
            gt_free((A)->space##TYPE);\
            (A)->allocated##TYPE = 0;\
            (A)->nextfree##TYPE = 0;\
            (A)->space##TYPE = NULL;\
          }\
        } while (false)

/*
  Some declarations for the most common array types.
*/
GT_DECLAREARRAYSTRUCT(GtUchar);
GT_DECLAREARRAYSTRUCT(GtUlong);
GT_DECLAREARRAYSTRUCT(GtUword);
GT_DECLAREARRAYSTRUCT(char);
GT_DECLAREARRAYSTRUCT(uint32_t);
GT_DECLAREARRAYSTRUCT(uint64_t);

#endif
