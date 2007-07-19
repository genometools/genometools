#include <math.h>
#include "types.h"
#include "divmodmul.h"
#include "fmindex.h"
#include "alphadef.h"
#include "safecast-gen.h"

static Seqpos determinenumofcodes(uint32_t numofchars,uint32_t prefixlength)
{
  return (Seqpos) pow((double) numofchars,(double) prefixlength);
}

Seqpos determinenumberofspecialstostore(
                    const Specialcharinfo *specialcharinfo)
{
  Seqpos addprefixsuffix = 0;

  if(specialcharinfo->lengthofspecialprefix > 0)
  {
    addprefixsuffix++;
  }
  if(specialcharinfo->lengthofspecialsuffix > 0)
  {
    addprefixsuffix++;
  }
  return specialcharinfo->specialranges + 1 - addprefixsuffix;
}

 DECLARESAFECASTFUNCTION(uint64_t,uint64_t,unsigned long,unsigned_long)

static unsigned long determinefmindexsize (const Fmindex *fm,
                                           Seqpos suffixlength,
                                           bool storeindexpos)
{
  uint64_t sumsize = 0;

  sumsize += (uint64_t) sizeof(Seqpos) * (uint64_t) TFREQSIZE(fm->mapsize);
  sumsize += (uint64_t) sizeof(Seqpos) * 
             (uint64_t) SUPERBFREQSIZE(fm->mapsize,fm->nofsuperblocks);
  if(storeindexpos)
  {
    sumsize += (uint64_t) sizeof(Seqpos) * 
               (uint64_t) MARKPOSTABLELENGTH(fm->bwtlength,fm->markdist);
  }
  if(suffixlength > 0)
  {
    sumsize += (uint64_t) sizeof(Bwtbound) * (uint64_t) fm->numofcodes;
  }
  if(storeindexpos)
  {
    sumsize += (uint64_t) sizeof(PairBwtidx) * 
               (uint64_t) determinenumberofspecialstostore(
                                 &fm->specialcharinfo);
  }
  sumsize += (uint64_t) sizeof(Uchar) * 
             (uint64_t) BFREQSIZE(fm->mapsize,fm->nofblocks);
  return CALLCASTFUNC(uint64_t,unsigned_long,sumsize);
}

void computefmkeyvalues (Fmindex *fm,
                         Seqpos bwtlength,
                         uint32_t log2bsize,
                         uint32_t log2markdist,
                         uint32_t mapsize,
                         uint32_t suffixlength,
                         bool storeindexpos,
                         const Specialcharinfo *specialcharinfo)
{
  fm->log2bsize = log2bsize;
  fm->log2markdist = log2markdist;
  fm->bwtlength = bwtlength;
  fm->log2superbsize = MULT2 (fm->log2bsize);
  fm->bsize = (uint32_t) POW2 (fm->log2bsize);
  fm->bsizehalve = DIV2(fm->bsize);
  fm->superbsize = (uint32_t) POW2 (fm->log2superbsize);
  fm->nofblocks = (Seqpos) (fm->bwtlength / fm->bsize) + 1;
  fm->nofsuperblocks = (Seqpos) (fm->bwtlength / fm->superbsize) + 2;
  fm->markdist = (Seqpos) POW2 (fm->log2markdist);
  fm->markdistminus1 = (Seqpos) (fm->markdist - 1);
  fm->negatebsizeones = ~ (Seqpos) (fm->bsize - 1);
  fm->negatesuperbsizeones = ~ (Seqpos) (fm->superbsize - 1);
  fm->log2superbsizeminuslog2bsize = fm->log2superbsize - fm->log2bsize;
  fm->mapsize = mapsize;
  fm->suffixlength = suffixlength;
  if(fm->suffixlength > 0)
  {
    fm->numofcodes = determinenumofcodes(fm->mapsize-1,fm->suffixlength);
  } else
  {
    fm->numofcodes = 0;
  }
  fm->specialcharinfo = *specialcharinfo;
  fm->sizeofindex = determinefmindexsize (fm,
                                          suffixlength,
                                          storeindexpos);
}
