/*
  Copyright (c) 2007 Gordon Gremme <gordon@gremme.org>
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

#ifndef SAFEARITH_IMP_H
#define SAFEARITH_IMP_H

#include <stdbool.h>
/*
  Low level macros for the safe arithmetics module.

  Many ideas and code fragments have been taken from
  http://www.fefe.de/intof.html
*/

/* generic macros to determine the minimum and maximum value of given <type> */
#define __HALF_MAX_SIGNED(type)  ((type)1 << (type) (sizeof (type) * 8 - 2))
#define __MAX_SIGNED(type)       (__HALF_MAX_SIGNED(type) - 1 +               \
                                  __HALF_MAX_SIGNED(type))
#define __MIN_SIGNED(type)       ((type) -1 - __MAX_SIGNED(type))

#define __MIN(type)              ((type) -1 < (type) 1                        \
                                 ? __MIN_SIGNED(type)                         \
                                 : (type)0)
#define __MAX(type)              ((type) ~__MIN(type))

/* safely assign <src> to <dest>, returns false on success, true otherwise */
#define assign(dest, src)                                                     \
        ({                                                                    \
          __typeof__(src) __x  = (src);                                       \
          __typeof__(dest) __y = (__typeof__(dest)) __x;                      \
          (__x==(__typeof__(src)) __y &&                                      \
           ((__x < (__typeof__(src)) 1) == (__y < (__typeof__(dest)) 1))      \
           ? (void)((dest)=__y), false : true);                               \
         })

/* safely add <a> to <b> and assign the result to <c>,
   returns false on success, true otherwise */
#define add_of(c, a, b)                                                       \
        ({                                                                    \
          __typeof__(a) __a = a;                                              \
          __typeof__(b) __b = b;                                              \
          (__b) < (__typeof__(b)) 1                                           \
          ? ((__MIN(__typeof__(c)) - (__b) <= (__typeof__(c-__b)) (__a))      \
             ? assign(c, __a + __b) : true)                                   \
          : ((__MAX(__typeof__(c)) - (__b) >= (__typeof__(c-__b)) (__a))      \
             ? assign(c, __a + __b) : true);                                  \
        })

/* safely subtract <b> from <a> and assign the result to <c>,
   returns false on success, true otherwise */
#define sub_of(c, a, b)                                                       \
        ({                                                                    \
          __typeof__(a) __a = a;                                              \
          __typeof__(b) __b = b;                                              \
          (__b) < (__typeof__(b)) 1                                           \
          ? ((__MAX(__typeof__(c)) + (__b) >= (__typeof__(c+__b))(__a))       \
             ? assign(c, __a - __b) : true)                                   \
          : ((__MIN(__typeof__(c)) + (__b) <= (__typeof__(c+__b))(__a))       \
             ? assign(c, __a - __b) : true);                                  \
        })

#endif
