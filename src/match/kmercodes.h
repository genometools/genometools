#ifndef KMERCODES_H
#define KMERCODES_H
#include "core/intbits.h"
#include "core/codetype.h"
#include "core/divmodmul.h"

static inline GtCodetype gt_kmercode_at_position(
                                   const GtTwobitencoding *twobitencoding,
                                   unsigned long pos,
                                   unsigned int kmersize)
{
  const unsigned int unitoffset = (unsigned int) GT_MODBYUNITSIN2BITENC(pos);
  const unsigned long unitindex = GT_DIVBYUNITSIN2BITENC(pos);
  const GtCodetype maskright = GT_MASKRIGHT(kmersize);

  if (unitoffset <= (unsigned int) GT_UNITSIN2BITENC - kmersize)
  {
    return (GtCodetype) (twobitencoding[unitindex]
            >> GT_MULT2(GT_UNITSIN2BITENC - kmersize - unitoffset))
           & maskright;
  } else
  {
    unsigned int shiftleft = GT_MULT2(unitoffset+kmersize-GT_UNITSIN2BITENC);
    return (GtCodetype)
           ((twobitencoding[unitindex] << shiftleft) |
            (twobitencoding[unitindex+1] >> (GT_MULT2(GT_UNITSIN2BITENC) -
                                             shiftleft)))
           & maskright;
  }
}

static inline GtCodetype gt_kmercode_complement(GtCodetype kmer,
                                                GtCodetype maskright)
{
  return kmer ^ maskright;
}

#define GT_SWAPBITPAIRS(KMER,L1,L2,D) (((KMER) & (3UL << L1)) >> D) |\
                                      (((KMER) & (3UL << L2)) << D)

/*@unused@*/ static inline GtCodetype gt_kmercode_reverse32(GtCodetype kmer)
{
  return GT_SWAPBITPAIRS(kmer,62,0,62) |
         GT_SWAPBITPAIRS(kmer,60,2,58) |
         GT_SWAPBITPAIRS(kmer,58,4,54) |
         GT_SWAPBITPAIRS(kmer,56,6,50) |
         GT_SWAPBITPAIRS(kmer,54,8,46) |
         GT_SWAPBITPAIRS(kmer,52,10,42) |
         GT_SWAPBITPAIRS(kmer,50,12,38) |
         GT_SWAPBITPAIRS(kmer,48,14,34) |
         GT_SWAPBITPAIRS(kmer,46,16,30) |
         GT_SWAPBITPAIRS(kmer,44,18,26) |
         GT_SWAPBITPAIRS(kmer,42,20,22) |
         GT_SWAPBITPAIRS(kmer,40,22,18) |
         GT_SWAPBITPAIRS(kmer,38,24,14) |
         GT_SWAPBITPAIRS(kmer,36,26,10) |
         GT_SWAPBITPAIRS(kmer,34,28,6) |
         GT_SWAPBITPAIRS(kmer,32,30,2);
}
#endif
