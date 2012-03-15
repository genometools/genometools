/*
  Copyright (c) 2011 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2011 Center for Bioinformatics, University of Hamburg

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

static inline GtCodetype gt_kmercode_reverse(GtCodetype kmer,
                                             unsigned int kmersize)
{
  switch (kmersize)
  {
    case 2:
      return GT_SWAPBITPAIRS(kmer,2,0,2);
    case 3:
      return GT_SWAPBITPAIRS(kmer,4,0,4) |
             (kmer & (3U << 2));
    case 4:
      /* 4 shifts, 4 &, 2 |, 1 assignment = 11 ops */
      kmer = (kmer & 0xF0) >> 4 | (kmer & 0x0F) << 4;
      return (kmer & 0xCC) >> 2 | (kmer & 0x33) << 2;
    case 5:
      return GT_SWAPBITPAIRS(kmer,8,0,8) |
             GT_SWAPBITPAIRS(kmer,6,2,4) |
             (kmer & (3U << 4));
    case 6:
      return GT_SWAPBITPAIRS(kmer,10,0,10) |
             GT_SWAPBITPAIRS(kmer,8,2,6) |
             GT_SWAPBITPAIRS(kmer,6,4,2);
    case 7:
      return GT_SWAPBITPAIRS(kmer,12,0,12) |
             GT_SWAPBITPAIRS(kmer,10,2,8) |
             GT_SWAPBITPAIRS(kmer,8,4,4) |
             (kmer & (3U << 6));
    case 8:
      /* 6 shifts, 6 &, 3 |, 2 assignments = 17 */
      kmer = (kmer & 0xFF00) >> 8 | (kmer & 0xFF) << 8;
      kmer = (kmer & 0XF0F0) >> 4 | (kmer & 0x0F0F) << 4;
      return (kmer & 0xCCCC) >> 2 | (kmer & 0x3333) << 2;
    case 9:
      return GT_SWAPBITPAIRS(kmer,16,0,16) |
             GT_SWAPBITPAIRS(kmer,14,2,12) |
             GT_SWAPBITPAIRS(kmer,12,4,8) |
             GT_SWAPBITPAIRS(kmer,10,6,4) |
             (kmer & (3U << 8));
    case 10:
      return GT_SWAPBITPAIRS(kmer,18,0,18) |
             GT_SWAPBITPAIRS(kmer,16,2,14) |
             GT_SWAPBITPAIRS(kmer,14,4,10) |
             GT_SWAPBITPAIRS(kmer,12,6,6) |
             GT_SWAPBITPAIRS(kmer,10,8,2);
    case 11:
      return GT_SWAPBITPAIRS(kmer,20,0,20) |
             GT_SWAPBITPAIRS(kmer,18,2,16) |
             GT_SWAPBITPAIRS(kmer,16,4,12) |
             GT_SWAPBITPAIRS(kmer,14,6,8) |
             GT_SWAPBITPAIRS(kmer,12,8,4) |
             (kmer & (3U << 10));
    case 12:
      return GT_SWAPBITPAIRS(kmer,22,0,22) |
             GT_SWAPBITPAIRS(kmer,20,2,18) |
             GT_SWAPBITPAIRS(kmer,18,4,14) |
             GT_SWAPBITPAIRS(kmer,16,6,10) |
             GT_SWAPBITPAIRS(kmer,14,8,6) |
             GT_SWAPBITPAIRS(kmer,12,10,2);
    case 13:
      return GT_SWAPBITPAIRS(kmer,24,0,24) |
             GT_SWAPBITPAIRS(kmer,22,2,20) |
             GT_SWAPBITPAIRS(kmer,20,4,16) |
             GT_SWAPBITPAIRS(kmer,18,6,12) |
             GT_SWAPBITPAIRS(kmer,16,8,8) |
             GT_SWAPBITPAIRS(kmer,14,10,4) |
             (kmer & (3U << 12));
    case 14:
      return GT_SWAPBITPAIRS(kmer,26,0,26) |
             GT_SWAPBITPAIRS(kmer,24,2,22) |
             GT_SWAPBITPAIRS(kmer,22,4,18) |
             GT_SWAPBITPAIRS(kmer,20,6,14) |
             GT_SWAPBITPAIRS(kmer,18,8,10) |
             GT_SWAPBITPAIRS(kmer,16,10,6) |
             GT_SWAPBITPAIRS(kmer,14,12,2);
    case 15:
      return GT_SWAPBITPAIRS(kmer,28,0,28) |
             GT_SWAPBITPAIRS(kmer,26,2,24) |
             GT_SWAPBITPAIRS(kmer,24,4,20) |
             GT_SWAPBITPAIRS(kmer,22,6,16) |
             GT_SWAPBITPAIRS(kmer,20,8,12) |
             GT_SWAPBITPAIRS(kmer,18,10,8) |
             GT_SWAPBITPAIRS(kmer,16,12,4) |
             (kmer & (3U << 14));
    case 16:
             /* 8 shifts, 8 &, 4 |, 3 assignments = 23 ops */
             kmer = (kmer & 0xFFFF0000) >> 16 | (kmer & 0x0000FFFF) << 16;
             kmer = (kmer & 0xFF00FF00) >> 8 |  (kmer & 0x00FF00FF) << 8;
             kmer = (kmer & 0XF0F0F0F0) >> 4 |  (kmer & 0x0F0F0F0F) << 4;
             return (kmer & 0xCCCCCCCC) >> 2 |  (kmer & 0x33333333) << 2;
#ifdef _LP64
    case 17:
      return GT_SWAPBITPAIRS(kmer,32,0,32) |
             GT_SWAPBITPAIRS(kmer,30,2,28) |
             GT_SWAPBITPAIRS(kmer,28,4,24) |
             GT_SWAPBITPAIRS(kmer,26,6,20) |
             GT_SWAPBITPAIRS(kmer,24,8,16) |
             GT_SWAPBITPAIRS(kmer,22,10,12) |
             GT_SWAPBITPAIRS(kmer,20,12,8) |
             GT_SWAPBITPAIRS(kmer,18,14,4) |
             (kmer & (3U << 16));
    case 18:
      return GT_SWAPBITPAIRS(kmer,34,0,34) |
             GT_SWAPBITPAIRS(kmer,32,2,30) |
             GT_SWAPBITPAIRS(kmer,30,4,26) |
             GT_SWAPBITPAIRS(kmer,28,6,22) |
             GT_SWAPBITPAIRS(kmer,26,8,18) |
             GT_SWAPBITPAIRS(kmer,24,10,14) |
             GT_SWAPBITPAIRS(kmer,22,12,10) |
             GT_SWAPBITPAIRS(kmer,20,14,6) |
             GT_SWAPBITPAIRS(kmer,18,16,2);
    case 19:
      return GT_SWAPBITPAIRS(kmer,36,0,36) |
             GT_SWAPBITPAIRS(kmer,34,2,32) |
             GT_SWAPBITPAIRS(kmer,32,4,28) |
             GT_SWAPBITPAIRS(kmer,30,6,24) |
             GT_SWAPBITPAIRS(kmer,28,8,20) |
             GT_SWAPBITPAIRS(kmer,26,10,16) |
             GT_SWAPBITPAIRS(kmer,24,12,12) |
             GT_SWAPBITPAIRS(kmer,22,14,8) |
             GT_SWAPBITPAIRS(kmer,20,16,4) |
             (kmer & (3U << 18));
    case 20:
      return GT_SWAPBITPAIRS(kmer,38,0,38) |
             GT_SWAPBITPAIRS(kmer,36,2,34) |
             GT_SWAPBITPAIRS(kmer,34,4,30) |
             GT_SWAPBITPAIRS(kmer,32,6,26) |
             GT_SWAPBITPAIRS(kmer,30,8,22) |
             GT_SWAPBITPAIRS(kmer,28,10,18) |
             GT_SWAPBITPAIRS(kmer,26,12,14) |
             GT_SWAPBITPAIRS(kmer,24,14,10) |
             GT_SWAPBITPAIRS(kmer,22,16,6) |
             GT_SWAPBITPAIRS(kmer,20,18,2);
    case 21:
      return GT_SWAPBITPAIRS(kmer,40,0,40) |
             GT_SWAPBITPAIRS(kmer,38,2,36) |
             GT_SWAPBITPAIRS(kmer,36,4,32) |
             GT_SWAPBITPAIRS(kmer,34,6,28) |
             GT_SWAPBITPAIRS(kmer,32,8,24) |
             GT_SWAPBITPAIRS(kmer,30,10,20) |
             GT_SWAPBITPAIRS(kmer,28,12,16) |
             GT_SWAPBITPAIRS(kmer,26,14,12) |
             GT_SWAPBITPAIRS(kmer,24,16,8) |
             GT_SWAPBITPAIRS(kmer,22,18,4) |
             (kmer & (3U << 20));
    case 22:
      return GT_SWAPBITPAIRS(kmer,42,0,42) |
             GT_SWAPBITPAIRS(kmer,40,2,38) |
             GT_SWAPBITPAIRS(kmer,38,4,34) |
             GT_SWAPBITPAIRS(kmer,36,6,30) |
             GT_SWAPBITPAIRS(kmer,34,8,26) |
             GT_SWAPBITPAIRS(kmer,32,10,22) |
             GT_SWAPBITPAIRS(kmer,30,12,18) |
             GT_SWAPBITPAIRS(kmer,28,14,14) |
             GT_SWAPBITPAIRS(kmer,26,16,10) |
             GT_SWAPBITPAIRS(kmer,24,18,6) |
             GT_SWAPBITPAIRS(kmer,22,20,2);
    case 23:
      return GT_SWAPBITPAIRS(kmer,44,0,44) |
             GT_SWAPBITPAIRS(kmer,42,2,40) |
             GT_SWAPBITPAIRS(kmer,40,4,36) |
             GT_SWAPBITPAIRS(kmer,38,6,32) |
             GT_SWAPBITPAIRS(kmer,36,8,28) |
             GT_SWAPBITPAIRS(kmer,34,10,24) |
             GT_SWAPBITPAIRS(kmer,32,12,20) |
             GT_SWAPBITPAIRS(kmer,30,14,16) |
             GT_SWAPBITPAIRS(kmer,28,16,12) |
             GT_SWAPBITPAIRS(kmer,26,18,8) |
             GT_SWAPBITPAIRS(kmer,24,20,4) |
             (kmer & (3U << 22));
    case 24:
      return GT_SWAPBITPAIRS(kmer,46,0,46) |
             GT_SWAPBITPAIRS(kmer,44,2,42) |
             GT_SWAPBITPAIRS(kmer,42,4,38) |
             GT_SWAPBITPAIRS(kmer,40,6,34) |
             GT_SWAPBITPAIRS(kmer,38,8,30) |
             GT_SWAPBITPAIRS(kmer,36,10,26) |
             GT_SWAPBITPAIRS(kmer,34,12,22) |
             GT_SWAPBITPAIRS(kmer,32,14,18) |
             GT_SWAPBITPAIRS(kmer,30,16,14) |
             GT_SWAPBITPAIRS(kmer,28,18,10) |
             GT_SWAPBITPAIRS(kmer,26,20,6) |
             GT_SWAPBITPAIRS(kmer,24,22,2);
    case 25:
      return GT_SWAPBITPAIRS(kmer,48,0,48) |
             GT_SWAPBITPAIRS(kmer,46,2,44) |
             GT_SWAPBITPAIRS(kmer,44,4,40) |
             GT_SWAPBITPAIRS(kmer,42,6,36) |
             GT_SWAPBITPAIRS(kmer,40,8,32) |
             GT_SWAPBITPAIRS(kmer,38,10,28) |
             GT_SWAPBITPAIRS(kmer,36,12,24) |
             GT_SWAPBITPAIRS(kmer,34,14,20) |
             GT_SWAPBITPAIRS(kmer,32,16,16) |
             GT_SWAPBITPAIRS(kmer,30,18,12) |
             GT_SWAPBITPAIRS(kmer,28,20,8) |
             GT_SWAPBITPAIRS(kmer,26,22,4) |
             (kmer & (3U << 24));
    case 26:
      return GT_SWAPBITPAIRS(kmer,50,0,50) |
             GT_SWAPBITPAIRS(kmer,48,2,46) |
             GT_SWAPBITPAIRS(kmer,46,4,42) |
             GT_SWAPBITPAIRS(kmer,44,6,38) |
             GT_SWAPBITPAIRS(kmer,42,8,34) |
             GT_SWAPBITPAIRS(kmer,40,10,30) |
             GT_SWAPBITPAIRS(kmer,38,12,26) |
             GT_SWAPBITPAIRS(kmer,36,14,22) |
             GT_SWAPBITPAIRS(kmer,34,16,18) |
             GT_SWAPBITPAIRS(kmer,32,18,14) |
             GT_SWAPBITPAIRS(kmer,30,20,10) |
             GT_SWAPBITPAIRS(kmer,28,22,6) |
             GT_SWAPBITPAIRS(kmer,26,24,2);
    case 27:
      return GT_SWAPBITPAIRS(kmer,52,0,52) |
             GT_SWAPBITPAIRS(kmer,50,2,48) |
             GT_SWAPBITPAIRS(kmer,48,4,44) |
             GT_SWAPBITPAIRS(kmer,46,6,40) |
             GT_SWAPBITPAIRS(kmer,44,8,36) |
             GT_SWAPBITPAIRS(kmer,42,10,32) |
             GT_SWAPBITPAIRS(kmer,40,12,28) |
             GT_SWAPBITPAIRS(kmer,38,14,24) |
             GT_SWAPBITPAIRS(kmer,36,16,20) |
             GT_SWAPBITPAIRS(kmer,34,18,16) |
             GT_SWAPBITPAIRS(kmer,32,20,12) |
             GT_SWAPBITPAIRS(kmer,30,22,8) |
             GT_SWAPBITPAIRS(kmer,28,24,4) |
             (kmer & (3U << 26));
    case 28:
      return GT_SWAPBITPAIRS(kmer,54,0,54) |
             GT_SWAPBITPAIRS(kmer,52,2,50) |
             GT_SWAPBITPAIRS(kmer,50,4,46) |
             GT_SWAPBITPAIRS(kmer,48,6,42) |
             GT_SWAPBITPAIRS(kmer,46,8,38) |
             GT_SWAPBITPAIRS(kmer,44,10,34) |
             GT_SWAPBITPAIRS(kmer,42,12,30) |
             GT_SWAPBITPAIRS(kmer,40,14,26) |
             GT_SWAPBITPAIRS(kmer,38,16,22) |
             GT_SWAPBITPAIRS(kmer,36,18,18) |
             GT_SWAPBITPAIRS(kmer,34,20,14) |
             GT_SWAPBITPAIRS(kmer,32,22,10) |
             GT_SWAPBITPAIRS(kmer,30,24,6) |
             GT_SWAPBITPAIRS(kmer,28,26,2);
    case 29:
      return GT_SWAPBITPAIRS(kmer,56,0,56) |
             GT_SWAPBITPAIRS(kmer,54,2,52) |
             GT_SWAPBITPAIRS(kmer,52,4,48) |
             GT_SWAPBITPAIRS(kmer,50,6,44) |
             GT_SWAPBITPAIRS(kmer,48,8,40) |
             GT_SWAPBITPAIRS(kmer,46,10,36) |
             GT_SWAPBITPAIRS(kmer,44,12,32) |
             GT_SWAPBITPAIRS(kmer,42,14,28) |
             GT_SWAPBITPAIRS(kmer,40,16,24) |
             GT_SWAPBITPAIRS(kmer,38,18,20) |
             GT_SWAPBITPAIRS(kmer,36,20,16) |
             GT_SWAPBITPAIRS(kmer,34,22,12) |
             GT_SWAPBITPAIRS(kmer,32,24,8) |
             GT_SWAPBITPAIRS(kmer,30,26,4) |
             (kmer & (3U << 28));
    case 30:
      return GT_SWAPBITPAIRS(kmer,58,0,58) |
             GT_SWAPBITPAIRS(kmer,56,2,54) |
             GT_SWAPBITPAIRS(kmer,54,4,50) |
             GT_SWAPBITPAIRS(kmer,52,6,46) |
             GT_SWAPBITPAIRS(kmer,50,8,42) |
             GT_SWAPBITPAIRS(kmer,48,10,38) |
             GT_SWAPBITPAIRS(kmer,46,12,34) |
             GT_SWAPBITPAIRS(kmer,44,14,30) |
             GT_SWAPBITPAIRS(kmer,42,16,26) |
             GT_SWAPBITPAIRS(kmer,40,18,22) |
             GT_SWAPBITPAIRS(kmer,38,20,18) |
             GT_SWAPBITPAIRS(kmer,36,22,14) |
             GT_SWAPBITPAIRS(kmer,34,24,10) |
             GT_SWAPBITPAIRS(kmer,32,26,6) |
             GT_SWAPBITPAIRS(kmer,30,28,2);
    case 31:
      return GT_SWAPBITPAIRS(kmer,60,0,60) |
             GT_SWAPBITPAIRS(kmer,58,2,56) |
             GT_SWAPBITPAIRS(kmer,56,4,52) |
             GT_SWAPBITPAIRS(kmer,54,6,48) |
             GT_SWAPBITPAIRS(kmer,52,8,44) |
             GT_SWAPBITPAIRS(kmer,50,10,40) |
             GT_SWAPBITPAIRS(kmer,48,12,36) |
             GT_SWAPBITPAIRS(kmer,46,14,32) |
             GT_SWAPBITPAIRS(kmer,44,16,28) |
             GT_SWAPBITPAIRS(kmer,42,18,24) |
             GT_SWAPBITPAIRS(kmer,40,20,20) |
             GT_SWAPBITPAIRS(kmer,38,22,16) |
             GT_SWAPBITPAIRS(kmer,36,24,12) |
             GT_SWAPBITPAIRS(kmer,34,26,8) |
             GT_SWAPBITPAIRS(kmer,32,28,4) |
             (kmer & (3U << 30));
    case 32:
             kmer = (kmer & 0xFFFFFFFF00000000) >> 32 |
                    (kmer & 0x00000000FFFFFFFF) << 32;
             kmer = (kmer & 0xFFFF0000FFFF0000) >> 16 |
                    (kmer & 0x0000FFFF0000FFFF) << 16;
             kmer = (kmer & 0xFF00FF00FF00FF00) >> 8 |
                    (kmer & 0x00FF00FF00FF00FF) << 8;
             kmer = (kmer & 0XF0F0F0F0F0F0F0F0) >> 4 |
                    (kmer & 0x0F0F0F0F0F0F0F0F) << 4;
             return (kmer & 0xCCCCCCCCCCCCCCCC) >> 2 |
                    (kmer & 0x3333333333333333) << 2;
#endif
    default: fprintf(stderr,"illegal kmersize=%u\n",kmersize);
             exit(GT_EXIT_PROGRAMMING_ERROR);
  }
}

#endif
