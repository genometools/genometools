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

#ifndef SAFEARITH_H
#define SAFEARITH_H

#include <inttypes.h>
#include "libgtcore/env.h"
#include "libgtcore/safearith_imp.h"

/*
  This module allows to do safe arithmetics in C (preventing overflows).
  For an introduction to the theory behind this see

  http://www.fefe.de/intof.html
*/

/* assign <src> to <dest> or exit upon overflow */
#define safe_assign(dest, src)                                               \
        do {                                                                \
          if (assign(dest, src)) {                                          \
            fprintf(stderr, "%s, l.%d: overflow in assignment\n", __FILE__, \
                    __LINE__);                                              \
            exit(EXIT_FAILURE);                                             \
          }                                                                 \
        } while (0)

/* add <a> to <b> and assign the result to <c> or exit upon overflow */
#define safe_add(c, a, b)                                                  \
        do {                                                              \
          if (add_of(c, a, b)) {                                          \
            fprintf(stderr, "%s, l.%d: overflow in addition\n", __FILE__, \
                    __LINE__);                                            \
            exit(EXIT_FAILURE);                                           \
          }                                                               \
        } while (0)

/* subtract <b> from <a> and assign the result to <c> or exit upon overflow */
#define safe_sub(c, a, b)                                                     \
        do {                                                                 \
          if (sub_of(c, a, b)) {                                             \
            fprintf(stderr, "%s, l.%d: overflow in subtraction\n", __FILE__, \
                    __LINE__);                                               \
            exit(EXIT_FAILURE);                                              \
          }                                                                  \
        } while (0)

int           safe_abs(int);
long          safe_labs(long);
long long     safe_llabs(long long);
uint32_t      safe_mult_u32(uint32_t, uint32_t);
uint64_t      safe_mult_u64(uint64_t, uint64_t);
unsigned long safe_mult_ulong(unsigned long, unsigned long);
long          safe_cast2long(unsigned long);
unsigned long safe_cast2ulong(long);
int           safearith_example(Env*);
int           safearith_unit_test(Env*);

#endif
