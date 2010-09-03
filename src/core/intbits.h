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

#ifndef INTBITS_H
#define INTBITS_H

#include <inttypes.h>
#include <string.h>
#include "core/divmodmul.h"
#include "core/ma_api.h"
#include "core/safecast-gen.h"

/*
  This file contains some definitions manipulating bitvectors represented
  by a \texttt{GtBitsequence}. In the comment lines we use $w$ for the word size
  and \texttt{\symbol{94}} for exponentiation of the previous character.
*/

#ifdef _LP64

#define GT_LOGWORDSIZE    6         /* base 2 logarithm of wordsize */
typedef uint64_t GtBitsequence;
#else

#define GT_LOGWORDSIZE   5          /* base 2 logarithm of wordsize */
typedef uint32_t GtBitsequence;

#endif

#define GT_INTWORDSIZE\
        (1 << GT_LOGWORDSIZE) /* # of bits in unsigned long = w */
#define GT_LASTHALVEBITS\
        ((((GtBitsequence) 1) << GT_DIV2(GT_INTWORDSIZE)) - 1)
#define GT_FIRSTHALVEBITS\
        (GT_LASTHALVEBITS << GT_DIV2(GT_INTWORDSIZE))
#define GT_FIRSTBIT\
        (((GtBitsequence) 1) << (GT_INTWORDSIZE-1)) /* \(10^{w-1}\) */
#define GT_ISBITSET(S,I)\
        (((S) << (I)) & GT_FIRSTBIT)         /* is \(i\)th bit set? */
#define GT_ITHBIT(I)\
        (GT_FIRSTBIT >> (I))                 /* \(0^{i}10^{w-i-1}\)  */
#define GT_SECONDBIT\
        (GT_FIRSTBIT >> 1)                   /* \(010^{w-2}\) */
#define GT_THIRDBIT\
        (GT_FIRSTBIT >> 2)                   /* \(0010^{w-3}\) */
#define GT_FOURTHBIT\
        (GT_FIRSTBIT >> 3)                   /* \(00010^{w-4}\) */
#define GT_FIFTHBIT\
        (GT_FIRSTBIT >> 4)                   /* \(000010^{w-3}\) */
#define GT_FIRSTTWOBITS\
        (((GtBitsequence) 3) << (GT_INTWORDSIZE-2)) /* \(11^{w-2}\) */
#define GT_EXCEPTFIRSTBIT\
        (~GT_FIRSTBIT)                       /* \(01^{w-1}\) */
#define GT_EXCEPTFIRSTTWOBITS\
        (GT_EXCEPTFIRSTBIT >> 1)             /* \(001^{w-2}\) */
#define GT_EXCEPTFIRSTTHREEBITS\
        (GT_EXCEPTFIRSTBIT >> 2)             /* \(0001^{w-3}\) */
#define GT_EXCEPTFIRSTFOURBITS\
        (GT_EXCEPTFIRSTBIT >> 3)             /* \(00001^{w-4}\) */

typedef GtBitsequence GtTwobitencoding;
#define GT_UNITSIN2BITENC              GT_DIV2(GT_INTWORDSIZE)
#define GT_DIVBYUNITSIN2BITENC(V)      ((V) >> (GT_LOGWORDSIZE-1))
#define GT_MODBYUNITSIN2BITENC(V)      ((V) & ((1 << (GT_LOGWORDSIZE-1))-1))

#define GT_DIVWORDSIZE(I)\
        ((I) >> GT_LOGWORDSIZE)              /* \((I) div w\) */

#define GT_MODWORDSIZE(I)\
        ((I) & (GT_INTWORDSIZE-1))           /* \((I) mod w\) */

#define GT_MULWORDSIZE(I)\
        ((I) << GT_LOGWORDSIZE)              /* \((I) * w\) */

/*
  The following defines the number of integers for a bitvector with N bits.
*/

#define GT_NUMOFINTSFORBITS(N)\
        ((GT_DIVWORDSIZE(N) == 0)\
           ? (size_t) 1 \
           : ((size_t) 1 + (size_t) GT_DIVWORDSIZE((N) - 1)))

/*
  The following macro allocates a bitarray of \texttt{N} bits. All bits
  are off.
*/

#define GT_INITBITTABGENERIC(TAB,OLDTAB,NUMOFBITS)\
        {\
          size_t tabsize = GT_NUMOFINTSFORBITS(NUMOFBITS);\
          TAB = gt_realloc(OLDTAB,sizeof (GtBitsequence) * tabsize);\
          (void) memset(TAB,0,sizeof (GtBitsequence) * tabsize);\
        }

#define GT_INITBITTAB(TAB,N) GT_INITBITTABGENERIC(TAB,NULL,N)

/*
  The following macro inititalizes a bitarray such that all bits
  are off.
*/

#define GT_CLEARBITTAB(TAB,N)\
        (void) memset(TAB,0,sizeof (GtBitsequence) * GT_NUMOFINTSFORBITS(N))

/*
  \texttt{SETIBIT(TAB,I)} sets the \texttt{I}-th bit in bitarray
  \texttt{TAB} to 1.
*/

#define GT_SETIBIT(TAB,I)    (TAB)[GT_DIVWORDSIZE(I)] |= \
                                    GT_ITHBIT(GT_MODWORDSIZE(I))

/*
  \texttt{UNSSETIBIT(TAB,I)} sets the \texttt{I}-th bit in bitarray
  \texttt{TAB} to 0.
*/

#define GT_UNSETIBIT(TAB,I)  (TAB)[GT_DIVWORDSIZE(I)] &= \
                                     ~(GT_ITHBIT(GT_MODWORDSIZE(I)))

/*
  \texttt{ISIBITSET(TAB,I)} checks if the \texttt{I}-th bit in bitarray
  \texttt{TAB} is 1.
*/

#define GT_ISIBITSET(TAB,I)  ((TAB)[GT_DIVWORDSIZE(I)] & \
                                    GT_ITHBIT(GT_MODWORDSIZE(I)))

/*
  \texttt{BITNUM2WORD(TAB,I)} delivers the integer containing
  the \texttt{I}-th bit.
*/

#define GT_BITNUM2WORD(TAB,I)  (TAB)[GT_DIVWORDSIZE(I)]

/*@unused@*/ static inline void gt_byte2string(char *buffer, unsigned char bs)
{
  unsigned int i;
  unsigned char mask;

  for (i=0, mask = (unsigned char) 128; i < 8U; i++, mask >>= 1)
  {
    buffer[i] = (bs & mask) ? '1' : '0';
  }
  buffer[8] = '\0';
}

/*@unused@*/ static inline void gt_bitsequence_tostring(char *buffer,
                                                        GtBitsequence bs)
{
  unsigned int i;
  GtBitsequence mask;

  for (i=0, mask = GT_FIRSTBIT;
       i < (unsigned int) GT_INTWORDSIZE;
       i++, mask >>= 1)
  {
    buffer[i] = (bs & mask) ? '1' : '0';
  }
  buffer[GT_INTWORDSIZE] = '\0';
}

static inline unsigned long gt_unitsoftwobitencoding(unsigned long totallength)
{
  uint64_t unitsoftwobitencoding;

  if (totallength < (unsigned long) GT_UNITSIN2BITENC)
  {
    return 2UL;
  }
  unitsoftwobitencoding = (uint64_t) (2 +
                          GT_DIVBYUNITSIN2BITENC(totallength - 1));
  return CALLCASTFUNC(uint64_t,unsigned_long,unitsoftwobitencoding);
}

#endif
