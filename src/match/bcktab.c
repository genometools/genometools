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

#include <string.h>
#include <errno.h>
#include <math.h>
#include <unistd.h>
#include "core/error.h"
#include "core/str.h"
#include "core/fa.h"
#include "core/types_api.h"
#include "core/chardef.h"
#include "core/mathsupport.h"
#include "core/unused_api.h"
#include "core/format64.h"
#include "core/mapspec-gen.h"
#include "core/unused_api.h"
#include "core/minmax.h"
#include "core/safecast-gen.h"
#include "intcode-def.h"
#include "esa-fileend.h"
#include "bcktab.h"
#include "initbasepower.h"
#include "stamp.h"

typedef struct
{
  unsigned long nonspecialsmaxbucketsize,
                specialsmaxbucketsize,
                maxbucketsize;
} GtMaxbucketinfo;

struct GtLeftborder
{
  uint32_t *uintbounds;
  uint32_t *uintboundsforpart;
  unsigned long *ulongbounds;
  unsigned long *ulongboundsforpart;
};

struct GtBcktab
{
  GtLeftborder leftborder;
  uint32_t *uintcountspecialcodes,
           *uintcountspecialcodesforpart,
           **uintdistpfxidx;
  unsigned long *ulongcountspecialcodes,
                *ulongcountspecialcodesforpart,
                **ulongdistpfxidx,
                sizeofrep,
                pagesize;
  unsigned int prefixlength,
               optimalnumofbits;
  GtCodetype numofallcodes,
             numofspecialcodes,
             **multimappower,
             *basepower,
             *filltable;
  GtUchar *qgrambuffer;
  GtMaxbucketinfo maxbucketinfo;
  bool allocated,
       withspecialsuffixes,
       useulong;
  void *mappedleftborder,
       *mappedcountspecialcodes,
       *mappedptr;
};

void gt_bcktab_leftborder_addcode(GtLeftborder *lb,GtCodetype code)
{
  gt_assert(lb != NULL);
  if (lb->ulongbounds != NULL)
  {
    lb->ulongbounds[code]++;
  } else
  {
    lb->uintbounds[code]++;
  }
}

static unsigned long checkboundsinsert = 0, checkboundsget = 0,
                     checkspecialcodesaccess = 0;

unsigned long gt_bcktab_leftborder_insertionindex(GtLeftborder *lb,
                                                  GtCodetype code)
{
  gt_assert(lb != NULL);
  if (lb->ulongbounds != NULL)
  {
    return --lb->ulongbounds[code];
  }
  if (lb->uintboundsforpart != NULL)
  {
    gt_assert(lb->uintboundsforpart[code] == lb->uintbounds[code]);
    checkboundsinsert++;
    --lb->uintboundsforpart[code];
  }
  return (unsigned long) --lb->uintbounds[code];
}

void gt_bcktab_leftborder_assign(GtLeftborder *lb,GtCodetype code,
                                 unsigned long value)
{
  gt_assert(lb != NULL);
  if (lb->ulongbounds != NULL)
  {
    lb->ulongbounds[code] = value;
  } else
  {
    gt_assert(value <= (unsigned long) UINT_MAX);
    if (lb->uintboundsforpart != NULL)
    {
      lb->uintboundsforpart[code] = (uint32_t) value;
    }
    lb->uintbounds[code] = (uint32_t) value;
  }
}

unsigned long gt_bcktab_get(const GtBcktab *bcktab,GtCodetype code)
{
  gt_assert(bcktab != NULL);
  if (bcktab->leftborder.ulongbounds != NULL)
  {
    return bcktab->leftborder.ulongbounds[code];
  }
  if (bcktab->leftborder.uintboundsforpart != NULL)
  {
    gt_assert(bcktab->leftborder.uintbounds[code] ==
              bcktab->leftborder.uintboundsforpart[code]);
    checkboundsget++;
  }
  return (unsigned long) bcktab->leftborder.uintbounds[code];
}

static unsigned long gt_bcktab_distpfxidx_get(const GtBcktab *bcktab,
                                              unsigned int prefixindex,
                                              GtCodetype ordercode)
{
  if (bcktab->ulongdistpfxidx != NULL)
  {
    return bcktab->ulongdistpfxidx[prefixindex][ordercode];
  }
  gt_assert(bcktab->uintdistpfxidx != NULL);
  return (unsigned long) bcktab->uintdistpfxidx[prefixindex][ordercode];
}

static void gt_bcktab_distpfxidx_increment(const GtBcktab *bcktab,
                                           unsigned int prefixindex,
                                           GtCodetype ordercode)
{
  if (bcktab->ulongdistpfxidx != NULL)
  {
    bcktab->ulongdistpfxidx[prefixindex][ordercode]++;
  } else
  {
    gt_assert(bcktab->uintdistpfxidx != NULL);
    bcktab->uintdistpfxidx[prefixindex][ordercode]++;
  }
}

size_t gt_bcktab_sizeofbasetype(const GtBcktab *bcktab)
{
  return bcktab->useulong ? sizeof (unsigned long) : sizeof (uint32_t);
}

static unsigned long multipleofpagesize(GtCodetype code,
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

#define FROMCODE2SPECIALCODE(CODE,NUMOFCHARS)\
                            (((NUMOFCHARS) == 4U)\
                            ? ((CODE) >> 2)\
                            : (((CODE) - ((NUMOFCHARS)-1)) / (NUMOFCHARS)))

void gt_bcktab_assignboundsforpart(GtBcktab *bcktab,
                                   const char *bcktmpfilename,
                                   unsigned int part,
                                   unsigned int numofchars,
                                   GtCodetype mincode,
                                   GtCodetype maxcode,
                                   GtLogger *logger)
{
  unsigned long mapoffset, mapend, totalsizeofcodes;

  if (bcktab->mappedleftborder != NULL)
  {
    gt_fa_xmunmap(bcktab->mappedleftborder);
  }
  if (bcktab->mappedcountspecialcodes != NULL)
  {
    gt_fa_xmunmap(bcktab->mappedcountspecialcodes);
  }
  mapoffset = multipleofpagesize(mincode,
                                 true,
                                 gt_bcktab_sizeofbasetype(bcktab),
                                 bcktab->pagesize);
  mapend = multipleofpagesize(maxcode,
                              false,
                              gt_bcktab_sizeofbasetype(bcktab),
                              bcktab->pagesize);
  totalsizeofcodes
    = (bcktab->numofallcodes+1) * gt_bcktab_sizeofbasetype(bcktab);
  gt_logger_log(logger,
             "part %u: mapped leftborder from %lu to %lu (%.2f of all)",
               part,mapoffset,mapend,
                 (mapend - mapoffset + 1 >= totalsizeofcodes)
                   ? 100.0
                   : 100.0 * (mapend - mapoffset + 1)/totalsizeofcodes);
  gt_assert(mapoffset <= mapend);
  gt_assert(mapoffset <= mincode * gt_bcktab_sizeofbasetype(bcktab));
  gt_assert(maxcode * gt_bcktab_sizeofbasetype(bcktab) <= mapend);
  gt_assert(mapoffset % bcktab->pagesize == 0);
  bcktab->mappedleftborder
    = gt_fa_xmmap_write_range(bcktmpfilename,(size_t) (mapend - mapoffset + 1),
                             (size_t) mapoffset);
  if (bcktab->useulong)
  {
    bcktab->leftborder.ulongboundsforpart
      = ((unsigned long *) bcktab->mappedleftborder) -
        (mapoffset / sizeof (unsigned long));
  } else
  {
    bcktab->leftborder.uintboundsforpart
      = ((uint32_t *) bcktab->mappedleftborder) -
         (mapoffset / sizeof (uint32_t));
  }
  if (bcktab->withspecialsuffixes)
  {
    GtCodetype firstcode, lastcode;

    if (mincode >= (GtCodetype) (numofchars - 1))
    {
      firstcode = bcktab->numofallcodes + 1 +
                  FROMCODE2SPECIALCODE(mincode,numofchars);
    } else
    {
      firstcode = bcktab->numofallcodes + 1;
    }
    if (maxcode >= (GtCodetype) (numofchars - 1))
    {
      lastcode = bcktab->numofallcodes + 1 +
                 FROMCODE2SPECIALCODE(maxcode,numofchars);
    } else
    {
      lastcode = bcktab->numofallcodes + 1;
    }
    mapoffset = multipleofpagesize(firstcode,
                                   true,
                                   gt_bcktab_sizeofbasetype(bcktab),
                                   bcktab->pagesize);
    mapend = multipleofpagesize(lastcode,
                                false,
                                gt_bcktab_sizeofbasetype(bcktab),
                                bcktab->pagesize);
    gt_assert(mapoffset % bcktab->pagesize == 0);
    totalsizeofcodes
      = bcktab->numofspecialcodes * gt_bcktab_sizeofbasetype(bcktab);
    gt_logger_log(logger,
                  "part %u: mapped countspecialcodes from %lu to %lu "
                  "(%.2f of all)",part,
                  mapoffset,mapend,
                  (mapend - mapoffset + 1 >= totalsizeofcodes)
                    ? 100.0
                    : 100.0 * (mapend - mapoffset + 1)/totalsizeofcodes);
    bcktab->mappedcountspecialcodes
    = gt_fa_xmmap_write_range(bcktmpfilename,(size_t) (mapend - mapoffset + 1),
                             (size_t) mapoffset);
    if (bcktab->useulong)
    {
      bcktab->ulongcountspecialcodesforpart
        = ((unsigned long *) bcktab->mappedcountspecialcodes) -
          (mapoffset / sizeof (unsigned long)) +
          bcktab->numofallcodes + 1;
    } else
    {
      bcktab->uintcountspecialcodesforpart
        = ((uint32_t *) bcktab->mappedcountspecialcodes) -
           (mapoffset / sizeof (uint32_t)) +
           bcktab->numofallcodes + 1;
    }
  }
}

static unsigned long numofdistpfxidxcounters(const GtCodetype *basepower,
                                             unsigned int prefixlength)
{
  if (prefixlength > 2U)
  {
    unsigned long numofcounters = 0;
    unsigned int idx;

    for (idx=1U; idx < prefixlength-1; idx++)
    {
      numofcounters += basepower[idx];
    }
    return numofcounters;
  }
  return 0;
}

static uint64_t gt_bcktab_sizeoftable_generic(unsigned int prefixlength,
                                              uint64_t numofallcodes,
                                              uint64_t numofspecialcodes,
                                              unsigned long maxvalue,
                                              const GtCodetype *basepower,
                                              bool withspecialsuffixes)
{
  uint64_t sizeoftable;
  size_t sizeofbasetype;

  sizeofbasetype = maxvalue <= (unsigned long) UINT_MAX
                                 ? sizeof (uint32_t)
                                 : sizeof (unsigned long);
  sizeoftable = (uint64_t) sizeofbasetype * (numofallcodes + 1);
  if (withspecialsuffixes)
  {
    sizeoftable += sizeofbasetype * numofspecialcodes;
    sizeoftable += (uint64_t) sizeofbasetype *
                              numofdistpfxidxcounters(basepower,prefixlength);
  }
  return sizeoftable;
}

uint64_t gt_bcktab_sizeoftable(unsigned int numofchars,
                               unsigned int prefixlength,
                               unsigned long maxvalue,
                               bool withspecialsuffixes)
{
  uint64_t numofallcodes, numofspecialcodes, sizeofbuckettable;
  GtCodetype *basepower;

  numofallcodes = (uint64_t) pow((double) numofchars,(double) prefixlength);
  if (withspecialsuffixes)
  {
    if (prefixlength >= 2U)
    {
      basepower = gt_initbasepower(numofchars,prefixlength-2);
    } else
    {
      basepower = NULL;
    }
    numofspecialcodes = (uint64_t) pow((double) numofchars,
                                       (double) (prefixlength-1));
  } else
  {
    basepower = NULL;
    numofspecialcodes = 0;
  }
  sizeofbuckettable
    = gt_bcktab_sizeoftable_generic(prefixlength,
                                    numofallcodes,
                                    numofspecialcodes,
                                    maxvalue,
                                    basepower,
                                    withspecialsuffixes);
  gt_free(basepower);
  return sizeofbuckettable;
}

unsigned long gt_bcktab_sizeofworkspace(unsigned int prefixlength)
{
  size_t size = sizeof (GtCodetype) * (prefixlength+1);
  size += sizeof (GtCodetype) * prefixlength;
  size += sizeof (GtBcktab);
  size += prefixlength;
  return (unsigned long) size;
}

static void ulong_setdistpfxidxptrs(unsigned long **ulongdistpfxidx,
                                    unsigned long *ptr,
                                    const GtCodetype *basepower,
                                    unsigned int prefixlength)
{
  unsigned int idx;

  ulongdistpfxidx[0] = ptr;
  for (idx=1U; idx<prefixlength-1; idx++)
  {
    ulongdistpfxidx[idx] = ulongdistpfxidx[idx-1] + basepower[idx];
  }
}

static void uint_setdistpfxidxptrs(uint32_t **uintdistpfxidx,
                                   uint32_t *ptr,
                                   const GtCodetype *basepower,
                                   unsigned int prefixlength)
{
  unsigned int idx;

  uintdistpfxidx[0] = ptr;
  for (idx=1U; idx<prefixlength-1; idx++)
  {
    uintdistpfxidx[idx] = uintdistpfxidx[idx-1] + basepower[idx];
  }
}

static void allocdistpfxidxcounts(GtBcktab *bcktab,unsigned int prefixlength,
                                  GtLogger *logger)
{
  gt_assert(bcktab->withspecialsuffixes);
  if (prefixlength > 2U)
  {
    unsigned long numofcounters;

    numofcounters = numofdistpfxidxcounters(bcktab->basepower,prefixlength);
    if (numofcounters > 0)
    {
      unsigned long *ulongcounters;
      uint32_t *uintcounters;
      size_t allocsize;

      if (bcktab->useulong)
      {
        bcktab->ulongdistpfxidx = gt_malloc(sizeof (void *) * (prefixlength-1));
        allocsize = sizeof (*ulongcounters) * numofcounters;
        ulongcounters = gt_malloc(allocsize);
        memset(ulongcounters,0,allocsize);
        ulong_setdistpfxidxptrs(bcktab->ulongdistpfxidx,ulongcounters,
                                bcktab->basepower,prefixlength);
      } else
      {
        bcktab->uintdistpfxidx = gt_malloc(sizeof (void *) * (prefixlength-1));
        allocsize = sizeof (*uintcounters) * numofcounters;
        uintcounters = gt_malloc(allocsize);
        memset(uintcounters,0,allocsize);
        uint_setdistpfxidxptrs(bcktab->uintdistpfxidx,uintcounters,
                               bcktab->basepower,prefixlength);
      }
      gt_logger_log(logger,"sizeof (distpfxidx)=%lu bytes",allocsize);
    }
  }
}

static GtBcktab *gt_bcktab_new_withinit(unsigned int numofchars,
                                        unsigned int prefixlength,
                                        unsigned long maxvalue,
                                        bool withspecialsuffixes)
{
  GtBcktab *bcktab;
  uint64_t sizeofrep_uint64_t;

  bcktab = gt_malloc(sizeof *bcktab);
  bcktab->mappedptr = NULL;
  bcktab->mappedleftborder = NULL;
  bcktab->leftborder.ulongbounds = NULL;
  bcktab->leftborder.ulongboundsforpart = NULL;
  bcktab->leftborder.uintbounds = NULL;
  bcktab->leftborder.uintboundsforpart = NULL;
  bcktab->mappedcountspecialcodes = NULL;
  bcktab->ulongcountspecialcodes = NULL;
  bcktab->ulongcountspecialcodesforpart = NULL;
  bcktab->uintcountspecialcodes = NULL;
  bcktab->uintcountspecialcodesforpart = NULL;
  bcktab->ulongdistpfxidx = NULL;
  bcktab->uintdistpfxidx = NULL;
  bcktab->prefixlength = prefixlength;
  bcktab->withspecialsuffixes = withspecialsuffixes;
  bcktab->basepower = gt_initbasepower(numofchars,prefixlength);
  bcktab->filltable = gt_initfilltable(numofchars,prefixlength);
  bcktab->numofallcodes = bcktab->basepower[prefixlength];
  bcktab->numofspecialcodes = bcktab->basepower[prefixlength-1];
  bcktab->multimappower = gt_initmultimappower(numofchars,prefixlength);
  bcktab->pagesize = (unsigned long) sysconf((int) _SC_PAGESIZE);
  gt_assert(bcktab->pagesize % sizeof (unsigned long) == 0);
  bcktab->allocated = false;
  bcktab->useulong = false;
  bcktab->qgrambuffer = gt_malloc(sizeof (*bcktab->qgrambuffer) * prefixlength);
  sizeofrep_uint64_t = gt_bcktab_sizeoftable_generic(
                                          prefixlength,
                                          (uint64_t) bcktab->numofallcodes,
                                          (uint64_t) bcktab->numofspecialcodes,
                                          maxvalue,
                                          bcktab->basepower,
                                          withspecialsuffixes);
  bcktab->sizeofrep = CALLCASTFUNC(uint64_t,unsigned_long,sizeofrep_uint64_t);
  return bcktab;
}

GtBcktab *gt_bcktab_new(unsigned int numofchars,
                        unsigned int prefixlength,
                        unsigned long maxvalue,
                        bool storespecialcodes,
                        bool withspecialsuffixes,
                        GtLogger *logger,
                        GtError *err)
{
  GtBcktab *bcktab;
  bool haserr = false;

  bcktab = gt_bcktab_new_withinit(numofchars,prefixlength,maxvalue,
                                  withspecialsuffixes);
  bcktab->allocated = true;
  if (storespecialcodes && bcktab->numofallcodes > 0 &&
      bcktab->numofallcodes-1 > (unsigned long) MAXCODEVALUE)
  {
    gt_error_set(err,"alphasize^prefixlength-1 = " FormatGtCodetype
                  " does not fit into %u"
                  " bits: choose smaller value for prefixlength",
                  bcktab->numofallcodes-1,
                  CODEBITS);
    haserr = true;
  }
  if (!haserr)
  {
    size_t allocsize_bounds, allocsize_countspecialcodes;
    if (maxvalue <= (unsigned long) UINT_MAX)
    {
      allocsize_bounds = sizeof (*bcktab->leftborder.uintbounds) *
                         (bcktab->numofallcodes+1);
      bcktab->leftborder.uintbounds = gt_malloc(allocsize_bounds);
      memset(bcktab->leftborder.uintbounds,0,allocsize_bounds);
      if (withspecialsuffixes)
      {
        allocsize_countspecialcodes = sizeof (*bcktab->uintcountspecialcodes) *
                                      bcktab->numofspecialcodes;
        bcktab->uintcountspecialcodes = gt_malloc(allocsize_countspecialcodes);
        memset(bcktab->uintcountspecialcodes,0,allocsize_countspecialcodes);
      }
      bcktab->useulong = false;
    } else
    {
      allocsize_bounds = sizeof (*bcktab->leftborder.ulongbounds) *
                         (bcktab->numofallcodes+1);
      bcktab->leftborder.ulongbounds = gt_malloc(allocsize_bounds);
      memset(bcktab->leftborder.ulongbounds,0,allocsize_bounds);
      if (withspecialsuffixes)
      {
        allocsize_countspecialcodes = sizeof (*bcktab->ulongcountspecialcodes) *
                                      bcktab->numofspecialcodes;
        bcktab->ulongcountspecialcodes = gt_malloc(allocsize_countspecialcodes);
        memset(bcktab->ulongcountspecialcodes,0,allocsize_countspecialcodes);
      }
      bcktab->useulong = true;
    }
    gt_logger_log(logger,"sizeof (leftborder)=%lu bytes",
                (unsigned long) allocsize_bounds);
    if (withspecialsuffixes)
    {
      gt_logger_log(logger,"sizeof (countspecialcodes)=%lu bytes",
                (unsigned long) allocsize_countspecialcodes);
      allocdistpfxidxcounts(bcktab,prefixlength,logger);
    }
    gt_logger_log(logger,"sizeof (bcktab)=" Formatuint64_t " bytes",
              PRINTuint64_tcast(gt_bcktab_sizeoftable(numofchars,
                                                      prefixlength,
                                                      maxvalue,
                                                      withspecialsuffixes)));
  }
  if (haserr)
  {
    gt_bcktab_delete(bcktab);
    return NULL;
  }
  return bcktab;
}

static void assignbcktabmapspecification(
                                        GtArrayGtMapspecification *mapspectable,
                                        void *voidinfo,
                                        bool writemode)
{
  GtBcktab *bcktab = (GtBcktab *) voidinfo;
  GtMapspecification *mapspecptr;
  unsigned long numofcounters;

  if (bcktab->useulong)
  {
    NEWMAPSPEC(bcktab->leftborder.ulongbounds,GtUlong,
               (unsigned long) (bcktab->numofallcodes+1));
  } else
  {
    NEWMAPSPEC(bcktab->leftborder.uintbounds,Uint32,
               (unsigned long) (bcktab->numofallcodes+1));
  }
  if (bcktab->withspecialsuffixes)
  {
    if (bcktab->useulong)
    {
      NEWMAPSPEC(bcktab->ulongcountspecialcodes,GtUlong,
                 (unsigned long) bcktab->numofspecialcodes);
    } else
    {
      NEWMAPSPEC(bcktab->uintcountspecialcodes,Uint32,
                 (unsigned long) bcktab->numofspecialcodes);
    }
    numofcounters
      = numofdistpfxidxcounters((const GtCodetype *) bcktab->basepower,
                                bcktab->prefixlength);
    if (numofcounters > 0)
    {
      if (!writemode)
      {
        if (bcktab->useulong)
        {
          bcktab->ulongdistpfxidx = gt_malloc(sizeof (void *) *
                                              (bcktab->prefixlength-1));
        } else
        {
          bcktab->uintdistpfxidx = gt_malloc(sizeof (void *) *
                                             (bcktab->prefixlength-1));
        }
      }
      if (bcktab->useulong)
      {
        NEWMAPSPEC(bcktab->ulongdistpfxidx[0],GtUlong,numofcounters);
      } else
      {
        NEWMAPSPEC(bcktab->uintdistpfxidx[0],Uint32,numofcounters);
      }
    }
  }
}

int gt_bcktab_flush_to_file(FILE *fp,const GtBcktab *bcktab,GtError *err)
{
  gt_error_check(err);
  return gt_mapspec_flushtheindex2file(fp,
                            assignbcktabmapspecification,
                            (GtBcktab *) bcktab,
                            bcktab->sizeofrep,
                            err);
}

static int fillbcktabmapspecstartptr(GtBcktab *bcktab,
                                     const char *indexname,
                                     GtError *err)
{
  bool haserr = false;
  GtStr *tmpfilename;

  gt_error_check(err);
  tmpfilename = gt_str_new_cstr(indexname);
  gt_str_append_cstr(tmpfilename,BCKTABSUFFIX);
  if (gt_mapspec_fillmapspecstartptr(assignbcktabmapspecification,
                          &bcktab->mappedptr,
                          bcktab,
                          tmpfilename,
                          bcktab->sizeofrep,
                          err) != 0)
  {
    haserr = true;
  }
  gt_str_delete(tmpfilename);
  return haserr ? -1 : 0;
}

#define CHECKCOUNTSPECIALCODES
#ifdef CHECKCOUNTSPECIALCODES
static unsigned long fromcode2countspecialcodes(GtCodetype code,
                                                const GtBcktab *bcktab)
{
  if (code >= bcktab->filltable[bcktab->prefixlength-1])
  {
    GtCodetype ordercode = code - bcktab->filltable[bcktab->prefixlength-1];
    GtCodetype divisor = bcktab->filltable[bcktab->prefixlength-1] + 1;
    if (ordercode % divisor == 0)
    {
      ordercode /= divisor;
      if (bcktab->ulongcountspecialcodes != NULL)
      {
        return bcktab->ulongcountspecialcodes[ordercode];
      } else
      {
        gt_assert(bcktab->uintcountspecialcodes != NULL);
        return (unsigned long) bcktab->uintcountspecialcodes[ordercode];
      }
    }
  }
  return 0;
}

static void pfxidxpartialsums(unsigned long *count,
                              GtCodetype code,
                              const GtBcktab *bcktab)
{
  unsigned int prefixindex;
  unsigned long sum = 0, specialsinbucket;
  GtCodetype ordercode, divisor;

  memset(count,0,sizeof (*count) * (size_t) bcktab->prefixlength);
  for (prefixindex=bcktab->prefixlength-2; prefixindex>=1U; prefixindex--)
  {
    if (code >= bcktab->filltable[prefixindex])
    {
      ordercode = code - bcktab->filltable[prefixindex];
      divisor = bcktab->filltable[prefixindex] + 1;
      if (ordercode % divisor == 0)
      {
        ordercode /= divisor;
        count[prefixindex]
          = gt_bcktab_distpfxidx_get(bcktab,prefixindex-1,ordercode);
        sum += count[prefixindex];
      }
    } else
    {
      break;
    }
  }
  specialsinbucket = fromcode2countspecialcodes(code,bcktab);
  gt_assert(sum <= specialsinbucket);
  count[bcktab->prefixlength-1] = specialsinbucket - sum;
  if (bcktab->prefixlength > 2U)
  {
    for (prefixindex = bcktab->prefixlength-2; prefixindex>=1U; prefixindex--)
    {
      count[prefixindex] += count[prefixindex+1];
    }
#ifndef NDEBUG
    if (specialsinbucket != count[1])
    {
      fprintf(stderr,"code " FormatGtCodetype ": sum = %lu != %lu = count[1]\n",
              code,sum,count[1]);
      exit(GT_EXIT_PROGRAMMING_ERROR);
    }
#endif
  }
}

void gt_bcktab_checkcountspecialcodes(const GtBcktab *bcktab)
{
  GtCodetype code;
  unsigned long *count;

  if (bcktab->prefixlength >= 2U)
  {
    count = gt_malloc(sizeof (*count) * bcktab->prefixlength);
    for (code=0; code<bcktab->numofallcodes; code++)
    {
      pfxidxpartialsums(count, code, bcktab);
    }
    gt_free(count);
  }
}
#endif

GtBcktab *gt_bcktab_map(const char *indexname,
                        unsigned int numofchars,
                        unsigned int prefixlength,
                        unsigned long maxvalue,
                        bool withspecialsuffixes,
                        GtError *err)
{
  GtBcktab *bcktab;

  bcktab = gt_bcktab_new_withinit(numofchars,prefixlength,maxvalue,
                                  withspecialsuffixes);
  bcktab->allocated = false;
  if (fillbcktabmapspecstartptr(bcktab, indexname, err) != 0)
  {
    gt_bcktab_delete(bcktab);
    return NULL;
  }
  if (withspecialsuffixes)
  {
    if (bcktab->ulongdistpfxidx != NULL)
    {
      ulong_setdistpfxidxptrs(bcktab->ulongdistpfxidx,
                              bcktab->ulongdistpfxidx[0],
                              bcktab->basepower,bcktab->prefixlength);
    } else
    {
      if (bcktab->uintdistpfxidx != NULL)
      {
        uint_setdistpfxidxptrs(bcktab->uintdistpfxidx,
                               bcktab->uintdistpfxidx[0],
                               bcktab->basepower,bcktab->prefixlength);
      }
    }
#ifdef CHECKCOUNTSPECIALCODES
    gt_bcktab_checkcountspecialcodes(bcktab);
#endif
  }
  return bcktab;
}

void gt_showbcktab(const GtBcktab *bcktab)
{
  unsigned int prefixindex;

  for (prefixindex=1U; prefixindex < bcktab->prefixlength-1; prefixindex++)
  {
    GtCodetype code;
    unsigned long sum = 0, value;
    for (code = 0; code < bcktab->basepower[prefixindex]; code++)
    {
      value = gt_bcktab_distpfxidx_get(bcktab,prefixindex-1,code);
      sum += value;
      printf("distpfxidx[%u][%lu]=%lu\n",prefixindex,code,value);
    }
    printf("sum %lu\n",sum);
  }
}

void gt_bcktab_showleftborder(const GtBcktab *bcktab)
{
  GtCodetype idx;

  for (idx=0; idx<bcktab->numofallcodes; idx++)
  {
    printf("leftborder[" FormatGtCodetype "]=%lu\n",idx,
                                                    gt_bcktab_get(bcktab,idx));
  }
}

void gt_bcktab_delete(GtBcktab *bcktab)
{
  if (bcktab == NULL)
  {
    return;
  }
  /*showbcktab(bcktabptr);*/
  if (bcktab->allocated)
  {
    gt_free(bcktab->leftborder.ulongbounds);
    gt_free(bcktab->leftborder.uintbounds);
    gt_free(bcktab->ulongcountspecialcodes);
    gt_free(bcktab->uintcountspecialcodes);
    if (bcktab->ulongdistpfxidx != NULL)
    {
      gt_free(bcktab->ulongdistpfxidx[0]);
    }
    if (bcktab->uintdistpfxidx != NULL)
    {
      gt_free(bcktab->uintdistpfxidx[0]);
    }
  } else
  {
    if (bcktab->mappedptr != NULL)
    {
      gt_fa_xmunmap(bcktab->mappedptr);
    }
    bcktab->mappedptr = NULL;
    if (bcktab->ulongdistpfxidx != NULL)
    {
      bcktab->ulongdistpfxidx[0] = NULL;
    }
    if (bcktab->uintdistpfxidx != NULL)
    {
      bcktab->uintdistpfxidx[0] = NULL;
    }
  }
  bcktab->leftborder.ulongbounds = NULL;
  bcktab->leftborder.uintbounds = NULL;
  bcktab->ulongcountspecialcodes = NULL;
  bcktab->uintcountspecialcodes = NULL;
  if (bcktab->mappedleftborder != NULL)
  {
    gt_fa_xmunmap(bcktab->mappedleftborder);
  }
  if (bcktab->mappedcountspecialcodes != NULL)
  {
    gt_fa_xmunmap(bcktab->mappedcountspecialcodes);
  }
  gt_free(bcktab->ulongdistpfxidx);
  bcktab->ulongdistpfxidx = NULL;
  gt_free(bcktab->uintdistpfxidx);
  bcktab->uintdistpfxidx = NULL;
  if (bcktab->multimappower != NULL)
  {
    gt_multimappowerfree(&bcktab->multimappower);
  }
  gt_free(bcktab->filltable);
  bcktab->filltable = NULL;
  gt_free(bcktab->basepower);
  bcktab->basepower = NULL;
  gt_free(bcktab->qgrambuffer);
  bcktab->qgrambuffer = NULL;
  gt_free(bcktab);
  /*
  printf("# checkboundsinsert=%lu\n",checkboundsinsert);
  printf("# checkboundsget=%lu\n",checkboundsget);
  printf("# checkspecialcodesaccess=%lu\n",checkspecialcodesaccess);
  */
}

void gt_bcktab_updatespecials(GtBcktab *bcktab,
                              GtCodetype code,
                              unsigned int numofchars,
                              unsigned int prefixindex)
{
  GtCodetype specialcode;

  gt_assert(prefixindex > 0);
  gt_assert(prefixindex <= bcktab->prefixlength);
  gt_assert(code < bcktab->numofallcodes);
  if (prefixindex < bcktab->prefixlength-1)
  {
    GtCodetype ordercode = (code - bcktab->filltable[prefixindex])/
                           (bcktab->filltable[prefixindex]+1);
    gt_bcktab_distpfxidx_increment(bcktab,prefixindex-1,ordercode);
  }
  gt_assert(code >= (GtCodetype) (numofchars-1));
  specialcode = FROMCODE2SPECIALCODE(code,numofchars);
  if (bcktab->ulongcountspecialcodes != NULL)
  {
    bcktab->ulongcountspecialcodes[specialcode]++;
  } else
  {
    gt_assert(bcktab->uintcountspecialcodes != NULL);
    bcktab->uintcountspecialcodes[specialcode]++;
    if (bcktab->uintcountspecialcodesforpart != NULL)
    {
      bcktab->uintcountspecialcodesforpart[specialcode]++;
    }
  }
}

void gt_bcktab_addfinalspecials(GtBcktab *bcktab,unsigned int numofchars,
                                unsigned long specialcharacters)
{
  if (bcktab->withspecialsuffixes)
  {
    GtCodetype specialcode;

    gt_assert(bcktab->filltable[0] >= (GtCodetype) (numofchars-1));
    specialcode = FROMCODE2SPECIALCODE(bcktab->filltable[0],numofchars);
    if (bcktab->ulongcountspecialcodes != NULL)
    {
      bcktab->ulongcountspecialcodes[specialcode]
        += (unsigned long) (specialcharacters + 1);
    } else
    {
      gt_assert(bcktab->uintcountspecialcodes != NULL);
      bcktab->uintcountspecialcodes[specialcode]
        += (uint32_t) (specialcharacters + 1);
      if (bcktab->uintcountspecialcodesforpart != NULL)
      {
        bcktab->uintcountspecialcodesforpart[specialcode]
          += (uint32_t) (specialcharacters + 1);
      }
    }
  }
}

unsigned int gt_bcktab_calcboundsparts(GtBucketspecification *bucketspec,
                                       const GtBcktab *bcktab,
                                       GtCodetype code,
                                       GtCodetype maxcode,
                                       unsigned long totalwidth,
                                       unsigned int rightchar,
                                       unsigned int numofchars)
{
  bucketspec->left = gt_bcktab_get(bcktab,code);
  if (code < maxcode)
  {
    unsigned long nextleftborder = gt_bcktab_get(bcktab,code+1);

    if (nextleftborder > 0)
    {
      bucketspec->nonspecialsinbucket = nextleftborder - bucketspec->left;
    } else
    {
      bucketspec->nonspecialsinbucket = 0;
    }
  } else
  {
    gt_assert(totalwidth >= bucketspec->left);
    bucketspec->nonspecialsinbucket
      = (unsigned long) (totalwidth - bucketspec->left);
  }
  if (bcktab->withspecialsuffixes && rightchar == numofchars - 1)
  {
    GtCodetype specialcode = FROMCODE2SPECIALCODE(code,numofchars);

    gt_assert(code >= (GtCodetype) (numofchars-1));
    if (bcktab->ulongcountspecialcodes != NULL)
    {
      bucketspec->specialsinbucket
        = bcktab->ulongcountspecialcodes[specialcode];
    } else
    {
      gt_assert(bcktab->uintcountspecialcodes != NULL);
      bucketspec->specialsinbucket
        = (unsigned long) bcktab->uintcountspecialcodes[specialcode];
      if (bcktab->uintcountspecialcodesforpart != NULL)
      {
        gt_assert(bucketspec->specialsinbucket ==
                  (unsigned long)
                  bcktab->uintcountspecialcodesforpart[specialcode]);
        checkspecialcodesaccess++;
      }
    }
    if (bucketspec->nonspecialsinbucket >= bucketspec->specialsinbucket)
    {
      bucketspec->nonspecialsinbucket -= bucketspec->specialsinbucket;
    } else
    {
      bucketspec->nonspecialsinbucket = 0;
    }
  } else
  {
    bucketspec->specialsinbucket = 0;
  }
  return (rightchar == numofchars - 1) ? 0 : (rightchar + 1);
}

void gt_bcktab_calcboundaries(GtBucketspecification *bucketspec,
                              const GtBcktab *bcktab,
                              GtCodetype code)
{
  unsigned int numofchars = (unsigned int) bcktab->basepower[1];
  gt_assert(code != bcktab->numofallcodes);
  (void) gt_bcktab_calcboundsparts(bucketspec,
                                   bcktab,
                                   code,
                                   bcktab->numofallcodes, /*code<numofallcodes*/
                                   0,  /* not necessary */
                                   bcktab->withspecialsuffixes
                                     ? (unsigned int) (code % numofchars)
                                     : 0,
                                   numofchars);
}

unsigned long gt_bcktab_calcrightbounds(const GtBcktab *bcktab,
                                        GtCodetype code,
                                        GtCodetype maxcode,
                                        unsigned long totalwidth)
{
  return code < maxcode ? gt_bcktab_get(bcktab,code+1) : totalwidth;
}

GtCodetype gt_bcktab_codedownscale(const GtBcktab *bcktab,
                                   GtCodetype code,
                                   unsigned int prefixindex,
                                   unsigned int maxprefixlen)
{
  unsigned int remain;

  code -= bcktab->filltable[maxprefixlen];
  remain = maxprefixlen-prefixindex;
  code %= (bcktab->filltable[remain]+1);
  code *= bcktab->basepower[remain];
  code += bcktab->filltable[prefixindex];
  return code;
}

void gt_bcktab_determinemaxsize(GtBcktab *bcktab,
                                const GtCodetype mincode,
                                const GtCodetype maxcode,
                                unsigned long partwidth,
                                unsigned int numofchars)
{
  unsigned int rightchar = (unsigned int) (mincode % numofchars);
  GtBucketspecification bucketspec;
  GtCodetype code;

#ifdef SKDEBUG
  printf("mincode=%lu,maxcode=%lu,partwidth=%lu,totallength=%lu\n",
          (unsigned long) mincode,(unsigned long) maxcode,
          partwidth,totallength);
#endif
  bcktab->maxbucketinfo.specialsmaxbucketsize = 1UL;
  bcktab->maxbucketinfo.nonspecialsmaxbucketsize = 1UL;
  bcktab->maxbucketinfo.maxbucketsize = 1UL;
  for (code = mincode; code <= maxcode; code++)
  {
    rightchar = gt_bcktab_calcboundsparts(&bucketspec,
                                          bcktab,
                                          code,
                                          maxcode,
                                          partwidth,
                                          rightchar,
                                          numofchars);
    if (bucketspec.nonspecialsinbucket >
        bcktab->maxbucketinfo.nonspecialsmaxbucketsize)
    {
      bcktab->maxbucketinfo.nonspecialsmaxbucketsize
        = bucketspec.nonspecialsinbucket;
    }
    if (bucketspec.specialsinbucket >
        bcktab->maxbucketinfo.specialsmaxbucketsize)
    {
      bcktab->maxbucketinfo.specialsmaxbucketsize = bucketspec.specialsinbucket;
    }
    if (bucketspec.nonspecialsinbucket + bucketspec.specialsinbucket
        > bcktab->maxbucketinfo.maxbucketsize)
    {
      bcktab->maxbucketinfo.maxbucketsize = bucketspec.nonspecialsinbucket +
                                            bucketspec.specialsinbucket;
    }
  }
  /*
  gt_logger_log(logger,"maxbucket (specials)=%lu",
              bcktab->maxbucketinfo.specialsmaxbucketsize);
  gt_logger_log(logger,"maxbucket (nonspecials)=%lu",
              bcktab->maxbucketinfo.nonspecialsmaxbucketsize);
  gt_logger_log(logger,"maxbucket (all)=%lu",
              bcktab->maxbucketinfo.maxbucketsize);
  */
}

static unsigned long gt_bcktab_nonspecialsmaxsize(const GtBcktab *bcktab)
{
  return bcktab->maxbucketinfo.nonspecialsmaxbucketsize;
}

unsigned int gt_bcktab_prefixlength(const GtBcktab *bcktab)
{
  return bcktab->prefixlength;
}

unsigned int gt_bcktab_singletonmaxprefixindex(const GtBcktab *bcktab,
                                               GtCodetype code)
{
  if (bcktab->prefixlength > 2U)
  {
    GtCodetype ordercode, divisor;
    unsigned int prefixindex;

    for (prefixindex=bcktab->prefixlength-2; prefixindex>=1U; prefixindex--)
    {
      if (code >= bcktab->filltable[prefixindex])
      {
        ordercode = code - bcktab->filltable[prefixindex];
        divisor = bcktab->filltable[prefixindex] + 1;
        if (ordercode % divisor == 0)
        {
          ordercode /= divisor;
          if (gt_bcktab_distpfxidx_get(bcktab,prefixindex-1,ordercode) > 0)
          {
            return prefixindex;
          }
        }
      } else
      {
        break;
      }
    }
  }
  return bcktab->prefixlength-1;
}

/* The following function is not used */

unsigned long gt_bcktab_distpfxidxpartialsums(const GtBcktab *bcktab,
                                              GtCodetype code,
                                              unsigned int lowerbound)
{
  GtCodetype ordercode, divisor;
  unsigned long sum = 0;
  unsigned int prefixindex;

  for (prefixindex=bcktab->prefixlength-2; prefixindex>lowerbound;
       prefixindex--)
  {
    if (code >= bcktab->filltable[prefixindex])
    {
      ordercode = code - bcktab->filltable[prefixindex];
      divisor = bcktab->filltable[prefixindex] + 1;
      if (ordercode % divisor == 0)
      {
        ordercode /= divisor;
        sum += gt_bcktab_distpfxidx_get(bcktab,prefixindex-1,ordercode);
      }
    } else
    {
      break;
    }
  }
  return sum;
}

unsigned int gt_bcktab_pfxidx2lcpvalues_uint8(unsigned int *minprefixindex,
                                              uint8_t *smalllcpvalues,
                                              unsigned long specialsinbucket,
                                              const GtBcktab *bcktab,
                                              GtCodetype code)
{
  unsigned int prefixindex, maxprefixindex = 0;
  unsigned long idx, value, insertpos;
  GtCodetype ordercode, divisor;

  gt_assert(smalllcpvalues != NULL);
  *minprefixindex = bcktab->prefixlength;
  insertpos = specialsinbucket;
  for (prefixindex=1U; prefixindex<bcktab->prefixlength-1; prefixindex++)
  {
    if (code >= bcktab->filltable[prefixindex])
    {
      ordercode = code - bcktab->filltable[prefixindex];
      divisor = bcktab->filltable[prefixindex] + 1;
      if (ordercode % divisor == 0)
      {
        ordercode /= divisor;
        value = gt_bcktab_distpfxidx_get(bcktab,prefixindex-1,ordercode);
        if (value > 0)
        {
          maxprefixindex = prefixindex;
          if (*minprefixindex > prefixindex)
          {
            *minprefixindex = prefixindex;
          }
          for (idx=0; idx < value; idx++)
          {
            gt_assert(insertpos > 0);
            smalllcpvalues[--insertpos] = (uint8_t) prefixindex;
          }
        }
      }
    }
  }
  if (insertpos > 0)
  {
    maxprefixindex = bcktab->prefixlength-1;
    if (*minprefixindex == bcktab->prefixlength)
    {
      *minprefixindex = bcktab->prefixlength-1;
    }
    while (insertpos > 0)
    {
      smalllcpvalues[--insertpos] = (uint8_t) (bcktab->prefixlength-1);
    }
  }
  return maxprefixindex;
}

unsigned int gt_bcktab_pfxidx2lcpvalues_ulong(unsigned int *minprefixindex,
                                              unsigned long *bucketoflcpvalues,
                                              unsigned long specialsinbucket,
                                              const GtBcktab *bcktab,
                                              GtCodetype code)
{
  unsigned int prefixindex, maxprefixindex = 0;
  unsigned long idx, value, insertpos;
  GtCodetype ordercode, divisor;

  gt_assert(bucketoflcpvalues != NULL);
  *minprefixindex = bcktab->prefixlength;
  insertpos = specialsinbucket;
  for (prefixindex=1U; prefixindex<bcktab->prefixlength-1; prefixindex++)
  {
    if (code >= bcktab->filltable[prefixindex])
    {
      ordercode = code - bcktab->filltable[prefixindex];
      divisor = bcktab->filltable[prefixindex] + 1;
      if (ordercode % divisor == 0)
      {
        ordercode /= divisor;
        value = gt_bcktab_distpfxidx_get(bcktab,prefixindex-1,ordercode);
        if (value > 0)
        {
          maxprefixindex = prefixindex;
          if (*minprefixindex > prefixindex)
          {
            *minprefixindex = prefixindex;
          }
          for (idx=0; idx < value; idx++)
          {
            gt_assert(insertpos > 0);
            bucketoflcpvalues[--insertpos] = (unsigned long) prefixindex;
          }
        }
      }
    }
  }
  if (insertpos > 0)
  {
    maxprefixindex = bcktab->prefixlength-1;
    if (*minprefixindex == bcktab->prefixlength)
    {
      *minprefixindex = bcktab->prefixlength-1;
    }
    while (insertpos > 0)
    {
      bucketoflcpvalues[--insertpos] = (unsigned long) (bcktab->prefixlength-1);
    }
  }
  return maxprefixindex;
}

const GtCodetype **gt_bcktab_multimappower(const GtBcktab *bcktab)
{
  return (const GtCodetype **) bcktab->multimappower;
}

GtCodetype gt_bcktab_filltable(const GtBcktab *bcktab,unsigned int idx)
{
  return bcktab->filltable[idx];
}

GtLeftborder *gt_bcktab_leftborder(GtBcktab *bcktab)
{
  return &bcktab->leftborder;
}

GtCodetype gt_bcktab_numofallcodes(const GtBcktab *bcktab)
{
  return bcktab->numofallcodes;
}

unsigned long gt_bcktab_leftborderpartialsums(
                             unsigned long *saved_bucketswithoutwholeleaf,
                             unsigned long *numofsuffixestosort,
                             GtBcktab *bcktab,
                             const GtBitsequence *markwholeleafbuckets)
{
  unsigned long code, largestbucketsize, sumbuckets, saved = 0, currentsize;

  gt_assert(bcktab->numofallcodes > 0);
  if (markwholeleafbuckets != NULL && !GT_ISIBITSET(markwholeleafbuckets,0))
  {
    saved += gt_bcktab_get(bcktab,0);
    gt_bcktab_leftborder_assign(&bcktab->leftborder,0,0);
    sumbuckets = 0;
  } else
  {
    sumbuckets = gt_bcktab_get(bcktab,0);
  }
  largestbucketsize = sumbuckets;
  for (code = 1UL; code < bcktab->numofallcodes; code++)
  {
    currentsize = gt_bcktab_get(bcktab,code);
    if (markwholeleafbuckets != NULL &&
        !GT_ISIBITSET(markwholeleafbuckets,code))
    {
      saved += currentsize;
      gt_bcktab_leftborder_assign(&bcktab->leftborder,code,sumbuckets);
    } else
    {
      sumbuckets += currentsize;
      if (largestbucketsize < currentsize)
      {
        largestbucketsize = currentsize;
      }
      gt_bcktab_leftborder_assign(&bcktab->leftborder,code,sumbuckets);
    }
  }
  gt_bcktab_leftborder_assign(&bcktab->leftborder,bcktab->numofallcodes,
                              sumbuckets);
  if (saved_bucketswithoutwholeleaf != NULL)
  {
    *saved_bucketswithoutwholeleaf = saved;
  }
  if (numofsuffixestosort != NULL)
  {
    *numofsuffixestosort = sumbuckets;
  }
  return largestbucketsize;
}

size_t gt_bcktab_sizeforlcpvalues(const GtBcktab *bcktab)
{
  size_t sizespeciallcps, sizelcps;

  sizespeciallcps = sizeof (uint8_t) *
                    bcktab->maxbucketinfo.specialsmaxbucketsize;
  sizelcps
    = sizeof (unsigned long) * gt_bcktab_nonspecialsmaxsize(bcktab);
  return MAX(sizelcps,sizespeciallcps);
}

GtCodetype gt_bcktab_findfirstlarger(const GtBcktab *bcktab,
                                     unsigned long suftaboffset)
{
  GtCodetype left = 0, right = bcktab->numofallcodes, mid,
             found = bcktab->numofallcodes;
  unsigned long midval;

  while (left+1 < right)
  {
    mid = GT_DIV2(left+right);
    midval = gt_bcktab_get(bcktab,mid);
    if (suftaboffset == midval)
    {
      return mid;
    }
    if (suftaboffset < midval)
    {
      found = mid;
      right = mid - 1;
    } else
    {
      left = mid + 1;
    }
  }
  return found;
}

#ifdef SKDEBUG
#include "qgram2code.h"
void gt_bcktab_consistencyofsuffix(int line,
                                   const GtEncseq *encseq,
                                   GtReadmode readmode,
                                   const GtBcktab *bcktab,
                                   unsigned int numofchars,
                                   const Suffixwithcode *suffix)
{
  unsigned int idx, firstspecial = bcktab->prefixlength, gramfirstspecial;
  GtCodetype qgramcode = 0;
  unsigned long totallength;
  GtUchar cc = 0;

  totallength = gt_encseq_total_length(encseq);
  for (idx=0; idx<bcktab->prefixlength; idx++)
  {
    if (suffix->startpos + idx >= totallength)
    {
      firstspecial = idx;
      break;
    }
    cc = gt_encseq_get_encoded_char(encseq,suffix->startpos + idx,
                                           readmode);
    if (ISSPECIAL(cc))
    {
      firstspecial = idx;
      break;
    }
    bcktab->qgrambuffer[idx] = cc;
  }
  for (idx=firstspecial; idx<bcktab->prefixlength; idx++)
  {
    bcktab->qgrambuffer[idx] = (GtUchar) (numofchars-1);
  }
  gramfirstspecial = qgram2code(&qgramcode,
                                (const GtCodetype **) bcktab->multimappower,
                                bcktab->prefixlength,
                                bcktab->qgrambuffer);
  gt_assert(gramfirstspecial == bcktab->prefixlength);
  gt_assert(qgramcode == suffix->code);
  if (firstspecial != suffix->prefixindex)
  {
    fprintf(stderr,"line %d: code=%u: ",line,suffix->code);
    fprintf(stderr,"firstspecial = %u != %u = suffix->prefixindex\n",
                    firstspecial,suffix->prefixindex);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }
}
#endif
