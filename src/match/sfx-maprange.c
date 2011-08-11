#include "core/intbits.h"
#include "sfx-maprange.h"

static unsigned long gt_multipleofpagesize(unsigned long code,
                                           bool smaller,
                                           size_t sizeofbasetype,
                                           unsigned long pagesize)
{
  if ((code * sizeofbasetype) % pagesize == 0)
  {
    return code * sizeofbasetype;
  }
  if (smaller)
  {
    return ((code * sizeofbasetype)/pagesize) * pagesize;
  }
  return ((code * sizeofbasetype)/pagesize) * pagesize + pagesize;
}

void gt_mapped_lbrange_get(GtMappedrange *range,
                           size_t sizeofbasetype,
                           unsigned long pagesize,
                           unsigned long mincode,
                           unsigned long maxcode)
{
  range->mapoffset = gt_multipleofpagesize(mincode,true,sizeofbasetype,
                                           pagesize);
  range->mapend = gt_multipleofpagesize(maxcode,false,sizeofbasetype,pagesize);
}

unsigned int gt_mapped_csrange_get(GtMappedrange *range,
                                   size_t sizeofbasetype,
                                   unsigned long numofallcodes,
                                   unsigned int numofchars,
                                   unsigned long pagesize,
                                   GtCodetype mincode,
                                   GtCodetype maxcode)
{
  GtCodetype firstcode, lastcode;
  unsigned int padoffset = 0;

  firstcode = numofallcodes + 1;
  if (sizeofbasetype < GT_WORDSIZE_INBYTES &&
      (firstcode * sizeofbasetype) % GT_WORDSIZE_INBYTES > 0)
  {
    padoffset = 1U;
  }
  if (mincode >= (GtCodetype) (numofchars - 1))
  {
    firstcode += FROMCODE2SPECIALCODE(mincode,numofchars);
  }
  if (maxcode >= (GtCodetype) (numofchars - 1))
  {
    lastcode = numofallcodes + 1 +
               FROMCODE2SPECIALCODE(maxcode,numofchars);
  } else
  {
    lastcode = numofallcodes + 1;
  }
  range->mapoffset
    = gt_multipleofpagesize(firstcode+padoffset,true,sizeofbasetype,pagesize);
  range->mapend = gt_multipleofpagesize(lastcode+padoffset,false,sizeofbasetype,
                                        pagesize);
  return padoffset;
}
