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

static GtCodetype gt_mapped_transformcode(unsigned long offset,
                                          unsigned int numofchars,
                                          GtCodetype code)
{
  if (code >= (GtCodetype) (numofchars - 1))
  {
    return offset + FROMCODE2SPECIALCODE(code,numofchars);
  } else
  {
    return offset;
  }
}

void gt_mapped_csrange_get(GtMappedrange *range,
                           unsigned long offset,
                           unsigned int numofchars,
                           size_t sizeofbasetype,
                           unsigned long pagesize,
                           GtCodetype mincode,
                           GtCodetype maxcode)
{
  gt_mapped_lbrange_get(range,
                        sizeofbasetype,
                        pagesize,
                        gt_mapped_transformcode(offset,numofchars,mincode),
                        gt_mapped_transformcode(offset,numofchars,maxcode));
}

unsigned int gt_Sfxmappedrange_padoffset(size_t sizeofbasetype,
                                         unsigned long offset)
{
  if (sizeofbasetype < GT_WORDSIZE_INBYTES &&
      ((offset+1) * sizeofbasetype) % GT_WORDSIZE_INBYTES > 0)
  {
    return 1U;
  } else
  {
    return 0;
  }
}

struct GtSfxmappedrange
{
  void *ptr, **usedptrptr;
  GtStr *filename;
  const char *tablename;
  unsigned long pagesize;
  size_t numofunits, sizeofunit;
  GtSfxmappedrangetype type;
  bool writable;
};

GtSfxmappedrange *gt_Sfxmappedrange_new(void **usedptrptr,
                                        unsigned long numofentries,
                                        GtSfxmappedrangetype type,
                                        const char *tablename,
                                        GtLogger *logger,
                                        GtError *err)
{
  GtSfxmappedrange *sfxmappedrange;
  FILE *outfp;
  bool haserr = false;

  sfxmappedrange = gt_malloc(sizeof (*sfxmappedrange));
  sfxmappedrange->ptr = NULL;
  sfxmappedrange->pagesize = (unsigned long) sysconf((int) _SC_PAGESIZE);
  sfxmappedrange->usedptrptr = usedptrptr;
  sfxmappedrange->filename = gt_str_new();
  outfp = gt_xtmpfp(sfxmappedrange->filename);
  sfxmappedrange->writable = false;
  sfxmappedrange->type = type;
  sfxmappedrange->tablename = tablename;
  switch (type)
  {
    case GtSfxGtBitsequence:
      sfxmappedrange->sizeofunit = sizeof (GtBitsequence);
      break;
    case GtSfxuint32_t:
      sfxmappedrange->sizeofunit = sizeof (uint32_t);
      break;
    case GtSfxunsignedlong:
      sfxmappedrange->sizeofunit = sizeof (unsigned long);
      break;
    default:
      gt_assert(false);
      break;
  }
  if (type == GtSfxGtBitsequence)
  {
    sfxmappedrange->numofunits = GT_NUMOFINTSFORBITS(numofentries);
  } else
  {
    sfxmappedrange->numofunits = (size_t) numofentries;
  }
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
    return NULL;
  }
  return sfxmappedrange;
}

void gt_Sfxmappedrange_make_writable(GtSfxmappedrange *sfxmappedrange)
{
  sfxmappedrange->writable = true;
}

void *gt_Sfxmappedrange_map(GtSfxmappedrange *sfxmappedrange,
                            unsigned long offset,
                            unsigned int numofchars,
                            unsigned int part,
                            unsigned long minindex,
                            unsigned long maxindex,
                            GtLogger *logger)
{
  GtMappedrange lbrange;
  unsigned long unitoffset;
  unsigned int padoffset;

  if (sfxmappedrange->ptr != NULL)
  {
    gt_fa_xmunmap(sfxmappedrange->ptr);
  }
  if (offset == 0)
  {
    padoffset = 0;
    gt_mapped_lbrange_get(&lbrange,
                          sfxmappedrange->sizeofunit,
                          sfxmappedrange->pagesize,
                          minindex,
                          maxindex);
    unitoffset = lbrange.mapoffset / sfxmappedrange->sizeofunit;
  } else
  {
    padoffset = gt_Sfxmappedrange_padoffset(sfxmappedrange->sizeofunit,offset);
    gt_mapped_csrange_get(&lbrange,
                          offset + padoffset + 1,
                          numofchars,
                          sfxmappedrange->sizeofunit,
                          sfxmappedrange->pagesize,
                          minindex,
                          maxindex);
    unitoffset = lbrange.mapoffset / sfxmappedrange->sizeofunit +
                 offset + 1 + padoffset;
  }
  if (logger != NULL)
  {
    size_t sizeoftable = sfxmappedrange->sizeofunit *
                         sfxmappedrange->numofunits;
    gt_logger_log(logger,
                  "part %u: mapped %s from %lu to %lu (%.1f%% of all)",
                  part,sfxmappedrange->tablename,lbrange.mapoffset,
                  lbrange.mapend,(lbrange.mapend - lbrange.mapoffset + 1
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

int gt_unlink_possibly_with_error(const char *filename,GtLogger *logger,
                                  GtError *err)
{
  bool haserr = false;

  gt_logger_log(logger,"remove \"%s\"",filename);
  if (unlink(filename) != 0)
  {
    if (err != NULL)
    {
      gt_error_set(err,"Cannot unlink file \"%s\": %s",
                      filename,strerror(errno));
      haserr = true;
    } else
    {
      fprintf(stderr,"Cannot unlink file \"%s\": %s",
                      filename,strerror(errno));
      exit(EXIT_FAILURE);
    }
  }
  return haserr ? -1 : 0;
}

int gt_Sfxmappedrange_delete(GtSfxmappedrange *sfxmappedrange,
                             GtLogger *logger,GtError *err)
{
  bool haserr = false;

  gt_assert(sfxmappedrange != NULL);
  if (sfxmappedrange->ptr != NULL)
  {
    gt_fa_xmunmap(sfxmappedrange->ptr);
  } else
  {
    gt_free(*sfxmappedrange->usedptrptr);
  }
  *sfxmappedrange->usedptrptr = NULL;
  if (sfxmappedrange->filename != NULL)
  {
    if (gt_unlink_possibly_with_error(gt_str_get(sfxmappedrange->filename),
                                      logger,err) != 0)
    {
      haserr = true;
    }
    gt_str_delete(sfxmappedrange->filename);
  }
  gt_free(sfxmappedrange);
  return haserr ? -1 : 0;
}
