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

/*
  This file contains some definitions manipulating bitvectors represented
  by a \texttt{Bitstring}. In the comment lines we use $w$ for the word size
  and \texttt{\symbol{94}} for exponentiation of the previous character.
*/

#ifdef _LP64

#define LOGWORDSIZE    6         /* base 2 logarithm of wordsize */
typedef uint64_t Bitstring;
#else

#define LOGWORDSIZE   5               /* base 2 logarithm of wordsize */
typedef uint32_t Bitstring;

#endif

#define INTWORDSIZE\
        (((Bitstring) 1) << LOGWORDSIZE) /* # of bits in unsigned long = w */
#define FIRSTBIT\
        (((Bitstring) 1) << (INTWORDSIZE-1)) /* \(10^{w-1}\) */
#define ISBITSET(S,I)\
        (((S) << (I)) & FIRSTBIT)         /* is \(i\)th bit set? */
#define ITHBIT(I)\
        (FIRSTBIT >> (I))                 /* \(0^{i}10^{w-i-1}\)  */
#define SECONDBIT\
        (FIRSTBIT >> 1)                   /* \(010^{w-2}\) */
#define THIRDBIT\
        (FIRSTBIT >> 2)                   /* \(0010^{w-3}\) */
#define FOURTHBIT\
        (FIRSTBIT >> 3)                   /* \(00010^{w-4}\) */
#define FIFTHBIT\
        (FIRSTBIT >> 4)                   /* \(000010^{w-3}\) */
#define FIRSTTWOBITS\
        (((Bitstring) 3) << (INTWORDSIZE-2)) /* \(11^{w-2}\) */
#define EXCEPTFIRSTBIT\
        (~FIRSTBIT)                       /* \(01^{w-1}\) */
#define EXCEPTFIRSTTWOBITS\
        (EXCEPTFIRSTBIT >> 1)             /* \(001^{w-2}\) */
#define EXCEPTFIRSTTHREEBITS\
        (EXCEPTFIRSTBIT >> 2)             /* \(0001^{w-3}\) */
#define EXCEPTFIRSTFOURBITS\
        (EXCEPTFIRSTBIT >> 3)             /* \(00001^{w-4}\) */

#define NEWTWOBITENCODING
#ifdef NEWTWOBITENCODING
typedef uint32_t Twobitencoding;
#define UNITSIN2BITENC              16
#define DIVBYUNITSIN2BITENC(V)      DIV16(V)
#define MODBYUNITSIN2BITENC(V)      MOD16(V)

#else
typedef unsigned char Twobitencoding;
#define UNITSIN2BITENC              4
#define DIVBYUNITSIN2BITENC(V)      DIV4(V)
#define MODBYUNITSIN2BITENC(V)      MOD4(V)
#endif

/*@unused@*/ static inline void byte2string(char *buffer,unsigned char bs)
{
  unsigned int i;
  unsigned char mask;

  for (i=0, mask = (unsigned char) 128; i < 8U; i++, mask >>= 1)
  {
    buffer[i] = (bs & mask) ? '1' : '0';
  }
  buffer[8] = '\0';
}

/*@unused@*/ static inline void uint32_t2string(char *buffer,uint32_t bs)
{
  unsigned int i;
  uint32_t mask;

  for (i=0, mask = (uint32_t) 2147483648U; i < 32U; i++, mask >>= 1)
  {
    buffer[i] = (bs & mask) ? '1' : '0';
  }
  buffer[32] = '\0';
}

#endif
