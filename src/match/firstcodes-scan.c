#include "core/intbits.h"
#include "core/types_api.h"
#include "core/codetype.h"
#include "core/divmodmul.h"
#include "kmercodes.h"
#include "firstcodes-scan.h"

static void gt_firstcode_verifycodes(const GtBitsequence *twobitencoding,
                                     unsigned long position,
                                     unsigned int kmersize,
                                     GtCodetype fcode,
                                     GtCodetype rccode)
{
  GtCodetype code = gt_kmercode_at_position(twobitencoding,position,kmersize);

  gt_assert(fcode == code);
  gt_assert(rccode == gt_kmercode_complement(gt_kmercode_reverse32(code),
                                             GT_MASKRIGHT(kmersize)));
}

static void gt_firstcodes_kmerscan_range(const GtBitsequence *twobitencoding,
                                         unsigned int kmersize,
                                         unsigned long startpos,
                                         unsigned long endpos,
                                         unsigned long maxunitindex,
                                         void (*processcode)(bool,GtCodetype,
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
  GtUchar cc;

  gt_assert(kmersize == 32U);
  position = startpos;
  fcode = gt_kmercode_at_position(twobitencoding, position, kmersize);
  processcode(true,fcode,position,data);
  rccode = gt_kmercode_complement(gt_kmercode_reverse32(fcode),maskright);
  processcode(false,rccode,position,data);
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
    cc = (GtUchar) (currentencoding >> shiftright) & 3;
    fcode = ((fcode << 2) | cc) & maskright;
    processcode(true,fcode,position,data);
    rccode = (rccode >> 2) | (unsigned long) ((cc ^ 3) << shiftleft);
    processcode(false,fcode,position,data);
    gt_firstcode_verifycodes(twobitencoding,
                             position,
                             kmersize,
                             fcode,
                             rccode);
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
                            void (*processcode)(bool,GtCodetype,
                                                unsigned long,void *),
                            void *data)
{
  unsigned long startpos;

  for (startpos = 0; startpos < totallength; startpos += equallength+1)
  {
    gt_firstcodes_kmerscan_range(twobitencoding,
                                 kmersize,
                                 startpos,
                                 startpos + equallength,
                                 maxunitindex,
                                 processcode,
                                 data);
  }
}
