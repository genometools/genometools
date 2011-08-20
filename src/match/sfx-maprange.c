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

#include <errno.h>
#include <unistd.h>
#include "core/fa.h"
#include "core/xposix.h"
#include "core/intbits.h"
#include "sfx-maprange.h"
#include "stamp.h"

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

typedef struct
{
  unsigned long mapoffset, mapend;
} GtMappedrange;

static void gt_mapped_lbrange_get(GtMappedrange *range,
                                  size_t sizeofbasetype,
                                  unsigned long pagesize,
                                  unsigned long mincode,
                                  unsigned long maxcode)
{
  range->mapoffset = gt_multipleofpagesize(mincode,true,sizeofbasetype,
                                           pagesize);
  range->mapend = gt_multipleofpagesize(maxcode,false,sizeofbasetype,pagesize);
}

struct GtSfxmappedrange
{
  void *ptr, *entire, **usedptrptr;
  GtStr *filename;
  const char *tablename;
  unsigned long pagesize;
  size_t numofunits, sizeofunit;
  GtSfxmappedrangetype type;
  unsigned long(*transformfunc)(unsigned long,unsigned int);
  unsigned transformfunc_data;
  bool writable;
};

void *gt_sfxmappedrange_map_entire(GtSfxmappedrange *sfxmappedrange,
                                   GtError *err)
{
  size_t mappedsize;

  sfxmappedrange->entire = gt_fa_mmap_read(gt_str_get(sfxmappedrange->filename),
                                           &mappedsize,err);
  if (sfxmappedrange->entire == NULL)
  {
    return NULL;
  }
  if (mappedsize != sfxmappedrange->sizeofunit * sfxmappedrange->numofunits)
  {
    gt_error_set(err,"map file %s: mapped size = %lu != %lu = "
                     "expected size",
                      gt_str_get(sfxmappedrange->filename),
                      mappedsize,
                      sfxmappedrange->sizeofunit *
                      sfxmappedrange->numofunits);
    gt_fa_xmunmap(sfxmappedrange->entire);
    sfxmappedrange->entire = NULL;
    return NULL;
  }
  return sfxmappedrange->entire;
}

GtSfxmappedrange *gt_Sfxmappedrange_new_basic(unsigned long numofentries,
                                              GtSfxmappedrangetype type,
                                              unsigned long(*transformfunc)(
                                                  unsigned long,unsigned int),
                                              unsigned int transformfunc_data)
{
  GtSfxmappedrange *sfxmappedrange;

  sfxmappedrange = gt_malloc(sizeof (*sfxmappedrange));
  sfxmappedrange->ptr = NULL;
  sfxmappedrange->pagesize = (unsigned long) sysconf((int) _SC_PAGESIZE);
  sfxmappedrange->usedptrptr = NULL;
  sfxmappedrange->filename = NULL;
  sfxmappedrange->writable = false;
  sfxmappedrange->entire = NULL;
  sfxmappedrange->transformfunc = transformfunc;
  sfxmappedrange->transformfunc_data = transformfunc_data;
  sfxmappedrange->type = type;
  sfxmappedrange->tablename = NULL;
  switch (type)
  {
    case GtSfxGtBitsequence:
      sfxmappedrange->sizeofunit = sizeof (GtBitsequence);
      sfxmappedrange->numofunits = GT_NUMOFINTSFORBITS(numofentries);
      break;
    case GtSfxuint32_t:
      sfxmappedrange->sizeofunit = sizeof (uint32_t);
      sfxmappedrange->numofunits = (size_t) numofentries;
      break;
    case GtSfxunsignedlong:
      sfxmappedrange->sizeofunit = sizeof (unsigned long);
      sfxmappedrange->numofunits = (size_t) numofentries;
      break;
    default:
      gt_assert(false);
      break;
  }
  return sfxmappedrange;
}

int gt_Sfxmappedrange_enhance(GtSfxmappedrange *sfxmappedrange,
                              void **usedptrptr,
                              bool writable,
                              const char *tablename,
                              GtLogger *logger,
                              GtError *err)
{
  bool haserr = false;
  FILE *outfp;

  sfxmappedrange->ptr = NULL;
  sfxmappedrange->usedptrptr = usedptrptr;
  sfxmappedrange->filename = gt_str_new();
  sfxmappedrange->writable = writable;
  outfp = gt_xtmpfp(sfxmappedrange->filename);
  gt_assert(outfp != NULL);
  sfxmappedrange->tablename = tablename;
  gt_logger_log(logger,"write %s to file %s (%lu units of %lu bytes)",
                tablename,gt_str_get(sfxmappedrange->filename),
                sfxmappedrange->numofunits,sfxmappedrange->sizeofunit);
  if (fwrite(*sfxmappedrange->usedptrptr,sfxmappedrange->sizeofunit,
             sfxmappedrange->numofunits,outfp) != sfxmappedrange->numofunits)
  {
    gt_error_set(err,"table %s: cannot write %lu items of size %u: "
                     "errormsg=\"%s\"",
                 tablename,
                 (unsigned long) sfxmappedrange->numofunits,
                 (unsigned int) sfxmappedrange->sizeofunit,
                 strerror(errno));
    haserr = true;
  }
  gt_fa_fclose(outfp);
  gt_free(*sfxmappedrange->usedptrptr);
  *sfxmappedrange->usedptrptr = NULL;
  if (haserr)
  {
    gt_str_delete(sfxmappedrange->filename);
    gt_free(sfxmappedrange);
    return -1;
  }
  return 0;
}

GtSfxmappedrange *gt_Sfxmappedrange_new(void **usedptrptr,
                                        bool writable,
                                        unsigned long numofentries,
                                        /* for GtBitsequence numofentries
                                           means the number of bits, otherwise
                                           the number of array indices */
                                        GtSfxmappedrangetype type,
                                        const char *tablename,
                                        unsigned long(*transformfunc)(
                                             unsigned long,unsigned int),
                                        unsigned int transformfunc_data,
                                        GtLogger *logger,
                                        GtError *err)
{
  GtSfxmappedrange *sfxmappedrange;

  sfxmappedrange = gt_Sfxmappedrange_new_basic(numofentries,
                                               type,
                                               transformfunc,
                                               transformfunc_data);
  if (gt_Sfxmappedrange_enhance(sfxmappedrange,
                                usedptrptr,
                                writable,
                                tablename,
                                logger,
                                err) != 0)
  {
    return NULL;
  }
  return sfxmappedrange;
}

unsigned long gt_Sfxmappedrange_mappedsize(GtSfxmappedrange *sfxmappedrange,
                                           unsigned long minindex,
                                           unsigned long maxindex)
{
  GtMappedrange lbrange;

  if (sfxmappedrange->transformfunc != NULL)
  {
    minindex = sfxmappedrange->transformfunc(minindex,
                                            sfxmappedrange->transformfunc_data);
    maxindex = sfxmappedrange->transformfunc(maxindex,
                                            sfxmappedrange->transformfunc_data);
  }
  gt_mapped_lbrange_get(&lbrange,
                        sfxmappedrange->sizeofunit,
                        sfxmappedrange->pagesize,
                        minindex,
                        maxindex);
  return lbrange.mapend - lbrange.mapoffset + 1;
}

void *gt_Sfxmappedrange_map(GtSfxmappedrange *sfxmappedrange,
                            unsigned int part,
                            unsigned long minindex,
                            unsigned long maxindex,
                            GtLogger *logger)
{
  GtMappedrange lbrange;
  unsigned long unitoffset;

  if (sfxmappedrange->ptr != NULL)
  {
    gt_fa_xmunmap(sfxmappedrange->ptr);
  }
  if (sfxmappedrange->transformfunc != NULL)
  {
    minindex = sfxmappedrange->transformfunc(minindex,
                                            sfxmappedrange->transformfunc_data);
    maxindex = sfxmappedrange->transformfunc(maxindex,
                                            sfxmappedrange->transformfunc_data);
  }
  gt_mapped_lbrange_get(&lbrange,
                        sfxmappedrange->sizeofunit,
                        sfxmappedrange->pagesize,
                        minindex,
                        maxindex);
  if (logger != NULL)
  {
    size_t sizeoftable = sfxmappedrange->sizeofunit *
                         sfxmappedrange->numofunits;
    gt_logger_log(logger,
                  "part %u: mapped %s from %lu to %lu for %s (%.1f%% of all)",
                  part,sfxmappedrange->tablename,lbrange.mapoffset,
                  lbrange.mapend,
                  sfxmappedrange->writable ? "writing" : "reading",
                  (lbrange.mapend - lbrange.mapoffset + 1
                    >= (unsigned long) sizeoftable)
                      ? 100.0
                      : 100.0 * (lbrange.mapend - lbrange.mapoffset + 1)/
                        sizeoftable);
  }
  gt_assert(lbrange.mapoffset <= lbrange.mapend);
  gt_assert(lbrange.mapoffset <= minindex * sfxmappedrange->sizeofunit);
  gt_assert(maxindex * sfxmappedrange->sizeofunit <= lbrange.mapend);
  gt_assert(lbrange.mapoffset % sfxmappedrange->pagesize == 0);
  if (sfxmappedrange->writable)
  {
    sfxmappedrange->ptr
      = gt_fa_xmmap_write_range (gt_str_get(sfxmappedrange->filename),
                                 (size_t) (lbrange.mapend-lbrange.mapoffset+1),
                                 (size_t) lbrange.mapoffset);
  } else
  {
    sfxmappedrange->ptr
      = gt_fa_xmmap_read_range (gt_str_get(sfxmappedrange->filename),
                                (size_t) (lbrange.mapend-lbrange.mapoffset+1),
                                (size_t) lbrange.mapoffset);
  }
  unitoffset = lbrange.mapoffset / sfxmappedrange->sizeofunit;
  switch (sfxmappedrange->type)
  {
    case GtSfxGtBitsequence:
      return ((GtBitsequence *) sfxmappedrange->ptr) - unitoffset;
    case GtSfxuint32_t:
      return ((uint32_t *) sfxmappedrange->ptr) - unitoffset;
    case GtSfxunsignedlong:
      return ((unsigned long *) sfxmappedrange->ptr) - unitoffset;
    default:
      gt_assert(false);
      break;
  }
  return NULL;
}

void gt_Sfxmappedrange_delete(GtSfxmappedrange *sfxmappedrange,GtLogger *logger)
{
  gt_assert(sfxmappedrange != NULL);
  gt_fa_xmunmap(sfxmappedrange->ptr);
  sfxmappedrange->ptr = NULL;
  gt_fa_xmunmap(sfxmappedrange->entire);
  sfxmappedrange->entire = NULL;
  if (sfxmappedrange->usedptrptr != NULL)
  {
    *sfxmappedrange->usedptrptr = NULL;
  }
  if (sfxmappedrange->filename != NULL)
  {
    gt_logger_log(logger,"remove \"%s\"",gt_str_get(sfxmappedrange->filename));
    gt_xunlink(gt_str_get(sfxmappedrange->filename));
  }
  gt_str_delete(sfxmappedrange->filename);
  gt_free(sfxmappedrange);
}
