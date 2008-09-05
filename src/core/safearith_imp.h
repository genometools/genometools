/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

/*
  Low level macros for the safe arithmetics module.

  Many ideas and code fragments have been taken from
  http://www.fefe.de/intof.html
*/

/* generic macros to determine the minimum and maximum value of given <type> */
#define __HALF_MAX_SIGNED(type)  ((type)1 << (sizeof(type) * 8 - 2))
#define __MAX_SIGNED(type)       (__HALF_MAX_SIGNED(type) - 1 + \
                                  __HALF_MAX_SIGNED(type))
#define __MIN_SIGNED(type)       (-1 - __MAX_SIGNED(type))

#define __MIN(type)              ((type) -1 < 1 ? __MIN_SIGNED(type) :(type)0)
#define __MAX(type)              ((type) ~__MIN(type))

/* safely assign <src> to <dest>, returns 0 on success, 1 otherwise */
#define assign(dest, src)                                                     \
        ({                                                                    \
          typeof(src) __x  = (src);                                           \
          typeof(dest) __y =__x;                                              \
          (__x==__y && ((__x < 1) == (__y < 1)) ? (void)((dest)=__y), 0 : 1); \
         })

/* safely add <a> to <b> and assign the result to <c>,
   returns 0 on success, 1 otherwise */
#define add_of(c, a, b)                                                       \
        ({                                                                    \
          typeof(a) __a = a;                                                  \
          typeof(b) __b = b;                                                  \
          (__b) < 1                                                           \
          ? ((__MIN(typeof(c)) - (__b) <= (__a)) ? assign(c, __a + __b) : 1)  \
          : ((__MAX(typeof(c)) - (__b) >= (__a)) ? assign(c, __a + __b) : 1); \
        })

/* safely subtract <b> from <c> and assign the result to <c>,
   returns 0 on success, 1 otherwise */
#define sub_of(c, a, b)                                                       \
        ({                                                                    \
          typeof(a) __a = a;                                                  \
          typeof(b) __b = b;                                                  \
          (__b) < 1                                                           \
          ? ((__MAX(typeof(c)) - (__b) >= (__a)) ? assign(c, __a - __b) : 1)  \
          : ((__MIN(typeof(c)) + (__b) <= (__a)) ? assign(c, __a - __b) : 1); \
        })

#endif
