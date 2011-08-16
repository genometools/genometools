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

GtSfxmappedrange *gt_Sfxmappedrange_new(void **usedptrptr,
                                        unsigned long numofindexes,
                                        size_t sizeofunit,
                                        const char *tablename,
                                        GtLogger *logger,
                                        GtError *err)
{
  GtSfxmappedrange *sfxmappedrange;
  size_t intsforbits;
  FILE *outfp;
  bool haserr = false;

  sfxmappedrange = gt_malloc(sizeof (*sfxmappedrange));
  sfxmappedrange->ptr = NULL;
  sfxmappedrange->pagesize = (unsigned long) sysconf((int) _SC_PAGESIZE);
  sfxmappedrange->usedptrptr = usedptrptr;
  sfxmappedrange->filename = gt_str_new();
  outfp = gt_xtmpfp(sfxmappedrange->filename);
  sfxmappedrange->tablename = tablename;
  sfxmappedrange->numofindexes = numofindexes;
  sfxmappedrange->sizeofunit = sizeofunit;
  intsforbits = GT_NUMOFINTSFORBITS(sfxmappedrange->numofindexes);
  gt_logger_log(logger,"write %s to file %s",
                tablename,gt_str_get(sfxmappedrange->filename));
  if (fwrite(*sfxmappedrange->usedptrptr,sizeofunit,intsforbits,outfp)
       != intsforbits)
  {
    gt_error_set(err,"table %s: cannot write %lu items of size %u: "
                     "errormsg=\"%s\"",
                 tablename,
                 (unsigned long) intsforbits,
                 (unsigned int) sizeofunit,
                 strerror(errno));
    haserr = true;
  }
  gt_fa_fclose(outfp);
  gt_free(*sfxmappedrange->usedptrptr);
  if (haserr)
  {
    gt_str_delete(sfxmappedrange->filename);
    gt_free(sfxmappedrange);
    return NULL;
  }
  return sfxmappedrange;
}

void *gt_Sfxmappedrange_map(GtSfxmappedrange *sfxmappedrange,
                            unsigned int part,
                            unsigned long minindex,
                            unsigned long maxindex,
                            GtLogger *logger)
{
  GtMappedrange lbrange;

  if (sfxmappedrange->ptr != NULL)
  {
    gt_fa_xmunmap(sfxmappedrange->ptr);
  }
  gt_mapped_lbrange_get(&lbrange, sfxmappedrange->sizeofunit,
                        sfxmappedrange->pagesize, minindex, maxindex);
  if (logger != NULL)
  {
    size_t sizeoftable = sfxmappedrange->sizeofunit *
                         GT_NUMOFINTSFORBITS(sfxmappedrange->numofindexes);
    gt_logger_log(logger,
               "part %u: mapped prefixbuckets from %lu to %lu (%.1f%% of all)",
                 part,lbrange.mapoffset,lbrange.mapend,
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
  sfxmappedrange->ptr
    = gt_fa_xmmap_read_range(gt_str_get(sfxmappedrange->filename),
                              (size_t) (lbrange.mapend - lbrange.mapoffset + 1),
                              (size_t) lbrange.mapoffset);
  if (sfxmappedrange->sizeofunit == sizeof (GtBitsequence))
  {
    return ((GtBitsequence *) sfxmappedrange->ptr) -
            (lbrange.mapoffset / sfxmappedrange->sizeofunit);
  } else
  {
    gt_assert(false);
    return NULL;
  }
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
