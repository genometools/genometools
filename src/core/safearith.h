/*
  Copyright (c) 2007-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c)      2010 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2007-2010 Center for Bioinformatics, University of Hamburg

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
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include "core/error.h"
#include "core/safearith_imp.h"

/*
  This module allows one to do safe arithmetics in C (preventing overflows).
  For an introduction to the theory behind this see

  http://www.fefe.de/intof.html
*/

/* assign <src> to <dest> or exit upon overflow */
#define gt_safe_assign(dest, src)                                           \
        do {                                                                \
          if (assign(dest, src)) {                                          \
            fprintf(stderr, "%s, l.%d: overflow in assignment\n", __FILE__, \
                    __LINE__);                                              \
            exit(EXIT_FAILURE);                                             \
          }                                                                 \
        } while (false)

/* add <a> to <b> and assign the result to <c> or exit upon overflow */
#define gt_safe_add(c, a, b)                                              \
        do {                                                              \
          if (add_of(c, a, b)) {                                          \
            fprintf(stderr, "%s, l.%d: overflow in addition\n", __FILE__, \
                    __LINE__);                                            \
            exit(EXIT_FAILURE);                                           \
          }                                                               \
        } while (false)

/* subtract <b> from <a> and assign the result to <c> or exit upon overflow */
#define gt_safe_sub(c, a, b)                                                 \
        do {                                                                 \
          if (sub_of(c, a, b)) {                                             \
            fprintf(stderr, "%s, l.%d: overflow in subtraction\n", __FILE__, \
                    __LINE__);                                               \
            exit(EXIT_FAILURE);                                              \
          }                                                                  \
        } while (false)

typedef void  (GtOverflowHandlerFunc)(const char *src_file, int src_line,
                                      void *data);

void          gt_safe_default_overflow_handler(const char*, int, void*);

int           gt_safe_abs_check_func(int, const char*, int,
                                     GtOverflowHandlerFunc, void*);
#define       gt_safe_abs(j) \
              gt_safe_abs_check_func(j, __FILE__, __LINE__, \
                                     gt_safe_default_overflow_handler, NULL)
#define       gt_safe_abs_check(j, func, data) \
              gt_safe_abs_check_func(j, __FILE__, __LINE__, func, data)

long          gt_safe_labs_check_func(long, const char*, int,
                                      GtOverflowHandlerFunc, void*);
#define       gt_safe_labs(j) \
              gt_safe_labs_check_func(j, __FILE__, __LINE__, \
                                      gt_safe_default_overflow_handler, NULL)
#define       gt_safe_labs_check(j, func, data) \
              gt_safe_labs_check_func(j, __FILE__, __LINE__, func, data)

long long     gt_safe_llabs_check_func(long long, const char*, int,
                                       GtOverflowHandlerFunc, void*);
#define       gt_safe_llabs(j) \
              gt_safe_llabs_check_func(j, __FILE__, __LINE__, \
                                       gt_safe_default_overflow_handler, NULL)
#define       gt_safe_llabs_check(j, func, data) \
              gt_safe_llabs_check_func(j, __FILE__, __LINE__, func, data)

uint32_t      gt_safe_mult_u32_check_func(uint32_t, uint32_t, const char*, int,
                                          GtOverflowHandlerFunc, void*);
#define       gt_safe_mult_u32(i, j) \
              gt_safe_mult_u32_check_func(i, j, __FILE__, __LINE__, \
                                          gt_safe_default_overflow_handler, \
                                          NULL)
#define       gt_safe_mult_u32_check(i, j, func, data) \
              gt_safe_mult_u32_check_func(i, j, __FILE__, __LINE__, func, data)

uint64_t      gt_safe_mult_u64_check_func(uint64_t, uint64_t, const char*, int,
                                          GtOverflowHandlerFunc, void*);
#define       gt_safe_mult_u64(i, j) \
              gt_safe_mult_u64_check_func(i, j, __FILE__, __LINE__, \
                                          gt_safe_default_overflow_handler, \
                                          NULL)
#define       gt_safe_mult_u64_check(i, j, func, data) \
              gt_safe_mult_u64_check_func(i, j, __FILE__, __LINE__, func, data)

unsigned long gt_safe_mult_ulong_check_func(unsigned long, unsigned long,
                                            const char*, int,
                                            GtOverflowHandlerFunc, void*);
#define       gt_safe_mult_ulong(i, j) \
              gt_safe_mult_ulong_check_func(i, j, __FILE__, __LINE__, \
                                          gt_safe_default_overflow_handler, \
                                          NULL)
#define       gt_safe_mult_ulong_check(i, j, func, data) \
              gt_safe_mult_ulong_check_func(i, j, __FILE__, __LINE__, func, \
                                            data)

long          gt_safe_cast2long_check_func(unsigned long, const char*, int,
                                           GtOverflowHandlerFunc, void*);
#define       gt_safe_cast2long(j) \
              gt_safe_cast2long_check_func(j, __FILE__, __LINE__, \
                                      gt_safe_default_overflow_handler, NULL)
#define       gt_safe_cast2long_check(j, func, data) \
              gt_safe_cast2long_check_func(j, __FILE__, __LINE__, func, data)

unsigned long gt_safe_cast2ulong_check_func(long, const char*, int,
                                            GtOverflowHandlerFunc, void*);
#define       gt_safe_cast2ulong(j) \
              gt_safe_cast2ulong_check_func(j, __FILE__, __LINE__, \
                                      gt_safe_default_overflow_handler, NULL)
#define       gt_safe_cast2ulong_check(j, func, data) \
              gt_safe_cast2ulong_check_func(j, __FILE__, __LINE__, func, data)

unsigned long gt_safe_cast2ulong_64_check_func(uint64_t, const char*, int,
                                               GtOverflowHandlerFunc, void*);
#define       gt_safe_cast2ulong_64(j) \
              gt_safe_cast2ulong_64_check_func(j, __FILE__, __LINE__, \
                                      gt_safe_default_overflow_handler, NULL)
#define       gt_safe_cast2ulong_64_check(j, func, data) \
              gt_safe_cast2ulong_64_check_func(j, __FILE__, __LINE__, func, \
                                               data)

int           gt_safearith_example(GtError*);
int           gt_safearith_unit_test(GtError*);

#endif
