/*
  Copyright (c) 2010 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2010 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
  Copyright (c) 2010 Center for Bioinformatics, University of Hamburg

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

#ifndef STACK_INLINED_H
#define STACK_INLINED_H
#include <stdbool.h>
#include <string.h>
#include "core/assert_api.h"
#include "core/ma.h"

/*
  This file defines macros to conveniently declare and
  manipulate dynamic stack whose size grow on demand. Each dynamic
  stack over some type <T> is implemented by a structure consisting of
  six components:
  - <space> is a pointer to the space block of type <T> allocated for the array.
  - <staticspace> is a statically allocated block of memory of type <T>.
    Initially, <space> points to this static array. If the stack grows
    beyond the size of the static array, then dynamically allocated
    memory is used.
  - <allocated> stores the number of entries in the memory referred
    to by <space>.
  - <nextfree> holds the smallest index of the stack where no value is stored.
  - <staticsize> holds the number of elements in <staticspace>
  - <sizeincrement> holds the number by new stack entries in each
    incrementation of the memory
  Here ## is the concatenation operator of the C-preprocessor.
  The following macro expands to a corresponding type definition over
  some given <TYPE>. <STATICSIZE> is the number of elements in the
  static memory area.
*/

#define GT_STACK_DECLARESTRUCT(TYPE,STATICSIZE)\
        typedef struct\
        {\
          unsigned long allocated, maxsize, nextfree,\
                        staticsize, sizeincrement;\
          TYPE staticspace[STATICSIZE], *space;\
          int (*initialiseelement)(void *);\
        } GtStack##TYPE

/*
  Initialize the stack <S> and define the increment.
*/

#define GT_STACK_INIT_WITH_INITFUNC(S,SIZEINCREMENT,INITFUNC)\
        (S)->staticsize = (unsigned long) sizeof ((S)->staticspace)/\
                                          sizeof ((S)->staticspace[0]);\
        (S)->sizeincrement = SIZEINCREMENT;\
        (S)->allocated = (S)->staticsize;\
        (S)->nextfree = 0;\
        (S)->maxsize = 0;\
        (S)->space = &(S)->staticspace[0];\
        (S)->initialiseelement = INITFUNC;\
        if ((S)->initialiseelement != NULL)\
        {\
          unsigned long stackidx;\
          (S)->initialiseelement = INITFUNC;\
          for (stackidx = 0; stackidx < (S)->staticsize; stackidx++)\
          {\
            (void) (S)->initialiseelement((S)->space + stackidx);\
          }\
        }

#define GT_STACK_INIT(S,SIZEINCREMENT)\
        GT_STACK_INIT_WITH_INITFUNC(S,SIZEINCREMENT,NULL)

/*
  Delete the memory allocated for stack contents
*/

#define GT_STACK_DELETE(S)\
        if ((S)->allocated > (S)->staticsize)\
        {\
          gt_free((S)->space);\
          (S)->space = NULL;\
        }

/*
  check if array needs expansion and do so if nessesary
*/

#define GT_STACK_CHECK_SPACE(S)\
        if ((S)->nextfree == (S)->allocated)\
        {\
          size_t sizeoftype = sizeof (*((S)->space));\
          size_t allocsize = sizeoftype * \
                             ((S)->allocated + (S)->sizeincrement);\
          (S)->space = gt_realloc(((S)->allocated == (S)->staticsize)\
                                  ? NULL \
                                  : (S)->space,\
                                  allocsize);\
          if ((S)->allocated == (S)->staticsize)\
          {\
            memcpy((S)->space,&(S)->staticspace[0],sizeoftype*(S)->staticsize);\
          }\
          if ((S)->initialiseelement != NULL)\
          {\
            unsigned long stackidx;\
            for (stackidx = 0; stackidx < (S)->sizeincrement; stackidx++)\
            {\
              (void) (S)->initialiseelement((S)->space + (S)->allocated + \
                                            stackidx);\
            }\
          }\
          (S)->allocated += (S)->sizeincrement;\
        }

/*
  Push a value.
*/

#define GT_STACK_PUSH(S,VALUE)\
        GT_STACK_CHECK_SPACE(S)\
        (S)->space[(S)->nextfree++] = VALUE;\
        if ((S)->maxsize < (S)->nextfree)\
        {\
          (S)->maxsize = (S)->nextfree;\
        }

/*
  check if the stack is empty.
*/

#define GT_STACK_ISEMPTY(S)\
        (((S)->nextfree == 0) ? true : false)

/*
  reduce the stack size by one
*/

#define GT_STACK_DECREMENTTOP(S)\
        (S)->nextfree--

/*
  get the top of the stack without changing stack size
*/

#define GT_STACK_TOP(S)\
        ((S)->space[(S)->nextfree - 1])

/*
  pop element from top of the stack
*/

#define GT_STACK_POP(S)\
        ((S)->space[--(S)->nextfree])

/*
  reduce the stack such that it becomes empty.
*/

#define GT_STACK_MAKEEMPTY(S)\
        (S)->nextfree = 0

/*
  reduce the stack such that it only contains one element.
*/

#define GT_STACK_MAKEALMOSTEMPTY(S)\
        (S)->nextfree = 1UL

/*
  point to next free element, adjust stacksize if nessesary
*/

#define GT_STACK_NEXT_FREE(S,PTR)\
        GT_STACK_CHECK_SPACE(S);\
        PTR = (S)->space + (S)->nextfree++;\
        if ((S)->maxsize < (S)->nextfree)\
        {\
          (S)->maxsize = (S)->nextfree;\
        }

/*
  return the maximum size of a stack
*/

#define GT_STACK_MAXSIZE(S) (S)->maxsize

#endif
