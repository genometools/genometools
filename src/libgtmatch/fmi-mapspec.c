/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "fmindex.h"
#include "mapspec-def.h"

#include "mapspec-gen.pr"

static void assignfmmapspecification(ArrayMapspecification *mapspectable,
                                     void *voidinfo,
                                     Env *env)
{
  Fmindexwithoptions *fmwithoptions = (Fmindexwithoptions *) voidinfo;
  Fmindex *fmindex;
  Mapspecification *mapspecptr;

  fmindex = fmwithoptions->fmptr;
  NEWMAPSPEC(fmindex->tfreq,Seqpos,(unsigned long) TFREQSIZE(fmindex->mapsize));
  NEWMAPSPEC(fmindex->superbfreq,Seqpos,
             (unsigned long) SUPERBFREQSIZE(fmindex->mapsize,
                                            fmindex->nofsuperblocks));
  NEWMAPSPEC(fmindex->markpostable,Seqpos,
             fmwithoptions->storeindexpos
             ? (unsigned long) MARKPOSTABLELENGTH(fmindex->bwtlength,
                                                  fmindex->markdist)
             : 0);
  NEWMAPSPEC(fmindex->boundarray,Bwtbound,(unsigned long) fmindex->numofcodes);
  NEWMAPSPEC(fmindex->specpos.spacePairBwtidx,PairBwtidx,
             fmwithoptions->storeindexpos 
             ? fmindex->specpos.nextfreePairBwtidx
             : 0);
  NEWMAPSPEC(fmindex->bfreq,Uchar,
             (unsigned long) BFREQSIZE(fmindex->mapsize,fmindex->nofblocks));
}

int flushfmindex2file(FILE *fp,
                      Fmindex *fmindex,
                      bool storeindexpos,
                      Env *env)
{
  Fmindexwithoptions fmwithoptions;

  fmwithoptions.fmptr = fmindex;
  fmwithoptions.storeindexpos = storeindexpos;
  return flushtheindex2file(fp,assignfmmapspecification,
                            (void *) &fmwithoptions,fmindex->sizeofindex,env);
}

int fillfmmapspecstartptr(Fmindex *fmindex,
                          bool storeindexpos,
                          const Str *tmpfilename,
                          Env *env)
{
  Fmindexwithoptions fmwithoptions;

  fmwithoptions.fmptr = fmindex;
  fmwithoptions.storeindexpos = storeindexpos;
  return fillmapspecstartptr(assignfmmapspecification,
                             &fmindex->mappedptr,
                             (void *) &fmwithoptions,
                             tmpfilename,
                             fmindex->sizeofindex,
                             env);
}
