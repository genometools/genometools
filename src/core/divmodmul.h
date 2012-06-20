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

#ifndef DIVMODMUL_H
#define DIVMODMUL_H

#include <math.h>
/*
  This file defines some macros to define division, multiplication,
  and modulo operations by 2, 4, 8, 16, 32 and 64.
*/

#define GT_DIV2(N)      ((N) >> 1)
#define GT_DIV4(N)      ((N) >> 2)
#define GT_DIV8(N)      ((N) >> 3)
#define GT_DIV16(N)     ((N) >> 4)
#define GT_DIV32(N)     ((N) >> 5)
#define GT_DIV64(N)     ((N) >> 6)

#define GT_MOD2(N)      ((N) & 1)
#define GT_MOD4(N)      ((N) & 3)
#define GT_MOD8(N)      ((N) & 7)
#define GT_MOD16(N)     ((N) & 15)
#define GT_MOD32(N)     ((N) & 31)
#define GT_MOD64(N)     ((N) & 63)

#define GT_MULT2(N)     ((N) << 1)
#define GT_MULT4(N)     ((N) << 2)
#define GT_MULT8(N)     ((N) << 3)
#define GT_MULT16(N)    ((N) << 4)
#define GT_MULT32(N)    ((N) << 5)
#define GT_MULT64(N)    ((N) << 6)

#define GT_POW2(N)      ((1) << (N))
#define GT_LOG2(v)      (log(v) / M_LN2)

#endif
