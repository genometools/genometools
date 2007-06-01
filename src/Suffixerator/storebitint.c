#include "stdint.h"

/*
  store bitvectorsi with q bits each in integer array.
  The first vector is stored at 0, the second at bit q,
  i.e. the i-th vector, i>= is stored at bit q*i. This is
  bit (q*i) mod 32 in integer (q*i) div 32
*/

#if _WORDSIZE == 64
#define LOGWORDSIZE 6
#else
#define LOGWORDSIZE 5
#endif

typedef uint64_t Seqlength;

#define ONEVAL           ((uint64_t) 1)
#define VECTORBITS       (ONEVAL << LOGWORDSIZE)
#define FRONTMASK(NBITS) (((ONEVAL << (NBITS)) - 1) << (_WORDSIZE - (NBITS)))
#define BACKMASK(NBITS)  ((ONEVAL << (NBITS)) - 1)

void storeinteger(uintptr_t *vector,uintptr_t bits,
                  Seqlength idx,Seqlength value)
{
  uint64_t bitprod = idx * bits;
  uint64_t firstidx = bitprod >> LOGWORDSIZE;
           firstbit = bitprod & ((ONEVAL << LOGWORDSIZE) - 1),
           lastidx = (bitprod + bits - 1) >> LOGWORDSIZE;

  if(firstindex == lastindex)
  {
    if(firstbit == 0)
    {
      vector[firstindex] = value << (VECTORBITS - bits) | 
                           (vector[firstindex] & BACKMASK(VECTORBITS - bits));
    } else
    {
      vector[firstindex] = value << ((VECTORBITS - firstbit) - bits) | 
                           (vector[firstindex] & FRONTMASK(VECTORBITS - bits))
                                               & BACKMA  ;
    }
  }
}
