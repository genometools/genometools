/*
  Copyright (c) 2013 Giorgio Gonnella <gonnella@zbh.uni-hamburg.de>
  Copyright (c) 2013 Center for Bioinformatics, University of Hamburg

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

#ifndef TRIANGULAR_DEF_H
#define TRIANGULAR_DEF_H

/* Macros for a compact representation of triangular and simmetric matrices.

   The macros are in two flavours:
   - with main diagonal (TRIANGULAR)
   - without main diagonal (STRICT_TRIANGULAR).
*/

/* The _SIZE(N) macros give the number of elements which must be actually
   allocated for a matrix of size N*N. */
#define GT_STRICT_TRIANGULAR_SIZE(N) (size_t)(((N)*(N-1))>>1)
#define GT_TRIANGULAR_SIZE(N)        (size_t)(((N+1)*(N))>>1)

/* The _ADDR(N,I,J) macros give the position of cell [i,j] in an (upper)
   triangular matrix. It only works if J >= I. */
#define GT_STRICT_TRIANGULAR_ADDR(N,I,J) (((I)*((N<<1)-3-(I)))>>1)+(J)-1
#define GT_TRIANGULAR_ADDR(N,I,J)        (((I)*((N<<1)-1-(I)))>>1)+(J)

/* The _SADDR(N,I,J) macros give the position of cell [i,j] in a
   simmetrical matrix. I and J may be any positive integer < N. */
#define GT_STRICT_TRIANGULAR_SADDR(N,I,J)\
  (J>=I) ? GT_STRICT_TRIANGULAR_ADDR(N,I,J) : GT_STRICT_TRIANGULAR_ADDR(N,J,I)
#define GT_TRIANGULAR_SADDR(N,I,J)\
  (J>=I) ? GT_TRIANGULAR_ADDR(N,I,J) : GT_TRIANGULAR_ADDR(N,J,I)

#endif
