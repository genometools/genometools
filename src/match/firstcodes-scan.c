#include "core/intbits.h"
#include "core/types_api.h"
#include "core/codetype.h"
#include "core/divmodmul.h"
#include "kmercodes.h"
#include "firstcodes-scan.h"

#define FIRSTCODES_VERIFY
#ifdef FIRSTCODES_VERIFY
static void gt_firstcode_verifycodes(const GtBitsequence *twobitencoding,
                                     unsigned long position,
                                     unsigned int kmersize,
                                     GtCodetype fcode,
                                     GtCodetype rccode)
{
  GtCodetype bfcode = gt_kmercode_at_position(twobitencoding,position,kmersize),
             bfrccode = gt_kmercode_complement(gt_kmercode_reverse32(bfcode),
                                               GT_MASKRIGHT(kmersize));

  gt_assert(fcode == bfcode);
  if (rccode != bfrccode)
  {
    printf("position %lu: fcode=%lu,rccode=%lu != %lu\n",position,fcode,
            rccode,bfrccode);
    exit(EXIT_FAILURE);
  }
}
#endif

static void gt_firstcodes_kmerscan_range(const GtBitsequence *twobitencoding,
                                         unsigned int kmersize,
                                         unsigned long startpos,
                                         unsigned long endpos,
                                         unsigned long maxunitindex,
                                         void (*processcode)(GtCodetype,
                                                             GtCodetype,
                                                             unsigned long,
                                                             void *),
                                         void *data)
{
  unsigned long position, unitindex;
  GtTwobitencoding currentencoding;
  unsigned int shiftright;
  const unsigned int shiftleft = GT_MULT2(kmersize-1);
  GtCodetype fcode, rccode;
  const unsigned long maskright = GT_MASKRIGHT(kmersize);
  GtCodetype cc;

  gt_assert(kmersize == 32U);
  position = startpos;
  fcode = gt_kmercode_at_position(twobitencoding, position, kmersize);
  rccode = gt_kmercode_complement(gt_kmercode_reverse32(fcode),maskright);
  if (processcode != NULL)
  {
    processcode(fcode,rccode,position,data);
  }
  unitindex = GT_DIVBYUNITSIN2BITENC(startpos + kmersize);
  currentencoding = twobitencoding[unitindex];
  shiftright = (unsigned int)
               GT_MULT2(GT_UNITSIN2BITENC - 1 -
                        GT_MODBYUNITSIN2BITENC(startpos + kmersize));
  gt_assert(endpos >= (unsigned long) kmersize);
  endpos -= kmersize;
  while (position < endpos)
  {
    position++;
    cc = (GtCodetype) (currentencoding >> shiftright) & 3;
    fcode = ((fcode << 2) | cc) & maskright;
    rccode = (rccode >> 2) | ((cc ^ 3UL) << shiftleft);
    if (processcode != NULL)
    {
      processcode(fcode,rccode,position,data);
    }
#ifdef FIRSTCODES_VERIFY
    gt_firstcode_verifycodes(twobitencoding,
                             position,
                             kmersize,
                             fcode,
                             rccode);
#endif
    if (shiftright > 0)
    {
      shiftright -= 2;
    } else
    {
      gt_assert(unitindex < maxunitindex-1 || position == endpos);
      if (unitindex < maxunitindex-1)
      {
        currentencoding = twobitencoding[++unitindex];
        shiftright = (unsigned int) (GT_INTWORDSIZE-2);
      }
    }
  }
}

void gt_firstcodes_kmerscan(const GtBitsequence *twobitencoding,
                            unsigned long equallength,
                            unsigned long totallength,
                            unsigned long maxunitindex,
                            unsigned int kmersize,
                            void (*processcode)(GtCodetype,GtCodetype,
                                                unsigned long,void *),
                            void *data)
{
  unsigned long startpos;

  for (startpos = 0; startpos < totallength; startpos += equallength+1)
  {
    /*printf("startpos=%lu,endpos=%lu\n",startpos,startpos+equallength);*/
    gt_firstcodes_kmerscan_range(twobitencoding,
                                 kmersize,
                                 startpos,
                                 startpos + equallength,
                                 maxunitindex,
                                 processcode,
                                 data);
  }
}
