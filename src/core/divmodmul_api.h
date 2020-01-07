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

#ifndef DIVMODMUL_API_H
#define DIVMODMUL_API_H

#include <math.h>

/* Divmodmul module */

/* Division by 2. */
#define GT_DIV2(N) \
        ((N) >> (__typeof__(N)) 1)
/* Division by 4. */
#define GT_DIV4(N) \
        ((N) >> (__typeof__(N)) 2)
/* Division by 8. */
#define GT_DIV8(N) \
        ((N) >> (__typeof__(N)) 3)
/* Division by 16. */
#define GT_DIV16(N) \
        ((N) >> (__typeof__(N)) 4)
/* Division by 32. */
#define GT_DIV32(N) \
        ((N) >> (__typeof__(N)) 5)
/* Division by 64. */
#define GT_DIV64(N) \
        ((N) >> (__typeof__(N)) 6)

/* Modulo 2. */
#define GT_MOD2(N) \
        ((N) & (__typeof__(N)) 1)
/* Modulo 4. */
#define GT_MOD4(N) \
        ((N) & (__typeof__(N)) 3)
/* Modulo 8. */
#define GT_MOD8(N) \
        ((N) & (__typeof__(N)) 7)
/* Modulo 16. */
#define GT_MOD16(N) \
        ((N) & (__typeof__(N)) 15)
/* Modulo 32. */
#define GT_MOD32(N) \
        ((N) & (__typeof__(N)) 31)
/* Modulo 64. */
#define GT_MOD64(N) \
        ((N) & (__typeof__(N)) 63)

/* Multiplication by 2. */
#define GT_MULT2(N) \
        ((N) << (__typeof__(N)) 1)
/* Multiplication by 4. */
#define GT_MULT4(N) \
        ((N) << (__typeof__(N)) 2)
/* Multiplication by 8. */
#define GT_MULT8(N) \
        ((N) << (__typeof__(N)) 3)
/* Multiplication by 16. */
#define GT_MULT16(N) \
        ((N) << (__typeof__(N)) 4)
/* Multiplication by 32. */
#define GT_MULT32(N) \
        ((N) << (__typeof__(N)) 5)
/* Multiplication by 64. */
#define GT_MULT64(N) \
        ((N) << (__typeof__(N)) 6)

/* Power of 2 to the N-th. */
#define GT_POW2(N) \
        ((__typeof__(N)) (1) << (__typeof__(N)) (N))
/* Binary logarithm of <v>. */
#define GT_LOG2(v) \
        (log(v) / M_LN2)

#endif
