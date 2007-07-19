/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "types.h"
#include "fmindex.h"
#include "mapspec-def.h"

#include "mapspec-gen.pr"

static void assignfmmapspecification(ArrayMapspecification *mapspectable,
                                     void *voidinfo,
                                     Env *env)
{
  Fmindexwithoptions *fmwithoptions = (Fmindexwithoptions *) voidinfo;
  Fmindex *fm;
  Mapspecification *mapspecptr;

  fm = fmwithoptions->fmptr;
  NEWMAPSPEC(fm->tfreq,Seqpos,(unsigned long) TFREQSIZE(fm->mapsize));
  NEWMAPSPEC(fm->superbfreq,Seqpos,
             (unsigned long) SUPERBFREQSIZE(fm->mapsize,fm->nofsuperblocks));
  NEWMAPSPEC(fm->markpostable,Seqpos,
             fmwithoptions->storeindexpos
             ? (unsigned long) MARKPOSTABLELENGTH(fm->bwtlength,fm->markdist)
             : 0);
  NEWMAPSPEC(fm->boundarray,Bwtbound,(unsigned long) fm->numofcodes);
  NEWMAPSPEC(fm->specpos.spacePairBwtidx,PairBwtidx,
             fmwithoptions->storeindexpos 
             ? fm->specpos.nextfreePairBwtidx
             : 0);
  NEWMAPSPEC(fm->bfreq,Uchar,
             (unsigned long) BFREQSIZE(fm->mapsize,fm->nofblocks));
}

int flushfmindex2file(FILE *fp,
                       Fmindex *fm,
                       bool storeindexpos,
                       Env *env)
{
  Fmindexwithoptions fmwithoptions;

  fmwithoptions.fmptr = fm;
  fmwithoptions.storeindexpos = storeindexpos;
  return flushtheindex2file(fp,assignfmmapspecification,
                            (void *) &fmwithoptions,fm->sizeofindex,env);
}

int fillfmmapspecstartptr(Fmindex *fm,bool storeindexpos,
                          const Str *tmpfilename,
                          Env *env)
{
  Fmindexwithoptions fmwithoptions;
  void *mappedusrptr;

  fmwithoptions.fmptr = fm;
  fmwithoptions.storeindexpos = storeindexpos;
  return fillmapspecstartptr(assignfmmapspecification,
                             &mappedusrptr,
                             (void *) &fmwithoptions,
                             tmpfilename,
                             fm->sizeofindex,
                             env);
}
