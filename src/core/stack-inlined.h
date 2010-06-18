/*
  Copyright (c) 2010 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
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
          unsigned long allocated, nextfree, staticsize, sizeincrement;\
          TYPE staticspace[STATICSIZE], *space;\
        } GtStack##TYPE

/*
  Initialize the stack <S> and define the increment.
*/

#define GT_STACK_INIT(S,SIZEINCREMENT)\
        (S)->staticsize = (unsigned long) sizeof ((S)->staticspace)/\
                                          sizeof ((S)->staticspace[0]);\
        (S)->sizeincrement = SIZEINCREMENT;\
        (S)->allocated = (S)->staticsize;\
        (S)->nextfree = 0;\
        (S)->space = &(S)->staticspace[0]

/*
<<<<<<< HEAD
  Delete the memory allocated for stack contents
=======
  Delete the stack.
>>>>>>> Added inlined implementation of stacks.
*/

#define GT_STACK_DELETE(S)\
        if ((S)->allocated > (S)->staticsize)\
        {\
          gt_free((S)->space);\
        }

/*
  Push a value.
*/

#define GT_STACK_PUSH(S,VALUE)\
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
          (S)->allocated += (S)->sizeincrement;\
        }\
        (S)->space[(S)->nextfree++] = VALUE

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
  reduce the stack such that it only contains one element.
*/

#define GT_STACK_MAKEALMOSTEMPTY(S)\
        (S)->nextfree = 1UL

#endif
