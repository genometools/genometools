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

#include "core/unused_api.h"
#include "core/mapspec.h"
#include "fmindex.h"

static void assignfmmapspecification(GtMapspec *mapspec,
                                     void *data,
                                     GT_UNUSED bool writemode)
{
  Fmindexwithoptions *fmwithoptions = (Fmindexwithoptions *) data;
  Fmindex *fmindex;

  fmindex = fmwithoptions->fmptr;
  gt_mapspec_add_ulong(mapspec, fmindex->tfreq,
                       (unsigned long) TFREQSIZE(fmindex->mapsize));
  gt_mapspec_add_ulong(mapspec, fmindex->superbfreq,
                       (unsigned long) SUPERBFREQSIZE(fmindex->mapsize,
                                                      fmindex->nofsuperblocks));
  gt_mapspec_add_ulong(mapspec, fmindex->markpostable,
                       fmwithoptions->storeindexpos
                       ? (unsigned long) MARKPOSTABLELENGTH(fmindex->bwtlength,
                                                            fmindex->markdist)
                       : 0);
  gt_mapspec_add_ulongbound(mapspec, fmindex->boundarray,
                            (unsigned long) fmindex->numofcodes);
  gt_mapspec_add_pairbwtindex(mapspec, fmindex->specpos.spaceGtPairBwtidx,
                              fmwithoptions->storeindexpos
                              ? fmindex->specpos.nextfreeGtPairBwtidx
                              : 0);
  gt_mapspec_add_uchar(mapspec, fmindex->bfreq,
                (unsigned long) BFREQSIZE(fmindex->mapsize,fmindex->nofblocks));
}

int gt_flushfmindex2file(FILE *fp,
                         Fmindex *fmindex,
                         bool storeindexpos,
                         GtError *err)
{
  Fmindexwithoptions fmwithoptions;

  gt_error_check(err);
  fmwithoptions.fmptr = fmindex;
  fmwithoptions.storeindexpos = storeindexpos;
  return gt_mapspec_write(assignfmmapspecification,fp,
                          (void *) &fmwithoptions,fmindex->sizeofindex,err);
}

int gt_fillfmmapspecstartptr(Fmindex *fmindex,
                             bool storeindexpos,
                             const GtStr *tmpfilename,
                             GtError *err)
{
  Fmindexwithoptions fmwithoptions;

  gt_error_check(err);
  fmwithoptions.fmptr = fmindex;
  fmwithoptions.storeindexpos = storeindexpos;
  return gt_mapspec_read(assignfmmapspecification,
                         (void *) &fmwithoptions,
                         tmpfilename,
                         fmindex->sizeofindex,
                         &fmindex->mappedptr,
                         err);
}
