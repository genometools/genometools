/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef DIVMODMUL_H
#define DIVMODMUL_H

/*
  This file defines some macros to define division, multiplication,
  and modulo operations by \(2, 4, 8, 16, 32\) and \(64\).
*/

#define DIV2(N)      ((N) >> 1)
#define DIV4(N)      ((N) >> 2)
#define DIV8(N)      ((N) >> 3)
#define DIV16(N)     ((N) >> 4)
#define DIV32(N)     ((N) >> 5)
#define DIV64(N)     ((N) >> 6)

#define MOD2(N)      ((N) & 1)
#define MOD4(N)      ((N) & 3)
#define MOD8(N)      ((N) & 7)
#define MOD16(N)     ((N) & 15)
#define MOD32(N)     ((N) & 31)
#define MOD64(N)     ((N) & 63)

#define MULT2(N)     ((N) << 1)
#define MULT4(N)     ((N) << 2)
#define MULT8(N)     ((N) << 3)
#define MULT16(N)    ((N) << 4)
#define MULT32(N)    ((N) << 5)
#define MULT64(N)    ((N) << 6)

#define POW2(N)      ((1) << (N))

#endif
