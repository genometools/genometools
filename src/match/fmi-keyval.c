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

#include <math.h>
#include "core/divmodmul.h"
#include "core/mathsupport.h"
#include "core/safecast-gen.h"
#include "fmindex.h"

unsigned long gt_determinenumberofspecialstostore(const GtSpecialcharinfo
                                                        *specialcharinfo)
{
  unsigned long addprefixsuffix = 0;

  if (specialcharinfo->lengthofspecialprefix > 0)
  {
    addprefixsuffix++;
  }
  if (specialcharinfo->lengthofspecialsuffix > 0)
  {
    addprefixsuffix++;
  }
  return specialcharinfo->realspecialranges + 1 - addprefixsuffix;
}

static unsigned long determinefmindexsize (const Fmindex *fm,
                                           const GtSpecialcharinfo
                                              *specialcharinfo,
                                           unsigned int suffixlength,
                                           bool storeindexpos)
{
  uint64_t sumsize = 0;

  sumsize +=
          (uint64_t) sizeof (unsigned long) * (uint64_t) TFREQSIZE(fm->mapsize);
  sumsize += (uint64_t) sizeof (unsigned long) *
             (uint64_t) SUPERBFREQSIZE(fm->mapsize,fm->nofsuperblocks);
  if (storeindexpos)
  {
    sumsize += (uint64_t) sizeof (unsigned long) *
               (uint64_t) MARKPOSTABLELENGTH(fm->bwtlength,fm->markdist);
  }
  if (suffixlength > 0)
  {
    sumsize += (uint64_t) sizeof (GtUlongBound) * (uint64_t) fm->numofcodes;
  }
  if (storeindexpos)
  {
    sumsize += (uint64_t) sizeof (GtPairBwtidx) *
               (uint64_t) gt_determinenumberofspecialstostore(specialcharinfo);
  }
  sumsize += (uint64_t) sizeof (GtUchar) *
             (uint64_t) BFREQSIZE(fm->mapsize,fm->nofblocks);
  return CALLCASTFUNC(uint64_t,unsigned_long,sumsize);
}

void gt_computefmkeyvalues (Fmindex *fm,
                            const GtSpecialcharinfo *specialcharinfo,
                            unsigned long bwtlength,
                            unsigned int log2bsize,
                            unsigned int log2markdist,
                            unsigned int numofchars,
                            unsigned int suffixlength,
                            bool storeindexpos)
{
  fm->mappedptr = NULL;
  fm->log2bsize = log2bsize;
  fm->log2markdist = log2markdist;
  fm->bwtlength = bwtlength;
  fm->log2superbsize = GT_MULT2 (fm->log2bsize);
  fm->bsize = (unsigned int) GT_POW2 (fm->log2bsize);
  fm->bsizehalve = GT_DIV2(fm->bsize);
  fm->superbsize = (unsigned int) GT_POW2 (fm->log2superbsize);
  fm->nofblocks = (unsigned long) (fm->bwtlength / fm->bsize) + 1;
  fm->nofsuperblocks = (unsigned long) (fm->bwtlength / fm->superbsize) + 2;
  fm->markdist = (unsigned long) GT_POW2 (fm->log2markdist);
  fm->markdistminus1 = (unsigned long) (fm->markdist - 1);
  fm->negatebsizeones = ~ (unsigned long) (fm->bsize - 1);
  fm->negatesuperbsizeones = ~ (unsigned long) (fm->superbsize - 1);
  fm->log2superbsizeminuslog2bsize = fm->log2superbsize - fm->log2bsize;
  fm->mapsize = numofchars+1;
  fm->suffixlength = suffixlength;
  if (fm->suffixlength > 0)
  {
    fm->numofcodes = gt_power_for_small_exponents(fm->mapsize-1,
                                                  fm->suffixlength);
  } else
  {
    fm->numofcodes = 0;
  }
  fm->sizeofindex = determinefmindexsize (fm,
                                          specialcharinfo,
                                          suffixlength,
                                          storeindexpos);
}
