/*
  Copyright (c) 2007      Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c)      2011 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2007-2011 Center for Bioinformatics, University of Hamburg

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

#include "core/error.h"
#include "core/fa.h"
#include "core/filelengthvalues.h"
#include "core/format64.h"
#include "core/intbits.h"
#include "core/types_api.h"
#include "core/mapspec.h"
#include "core/pairbwtidx.h"
#include "core/safecast-gen.h"
#include "core/str.h"
#include "core/ulongbound.h"
#include "core/xansi_api.h"

typedef enum
{
  GtCharType, /* \0 terminated string */
  GtFilelengthvaluesType,
  GtUcharType,
  Uint16Type,
  Uint32Type,
  Uint64Type,
  GtUlongType,
  GtUlongBoundType,
  GtBitsequenceType,
  GtPairBwtidxType,
  GtTwobitencodingType,
  GtSpecialcharinfoType,
  GtBitElemType
} GtTypespec;

typedef struct
{
  GtTypespec typespec;
  const char *function;
  void *startptr;
  size_t sizeofunit;
  unsigned long numofunits;
} GtMapspecification;

GT_DECLAREARRAYSTRUCT(GtMapspecification);

#define ASSIGNPTR2STARTPTR(TYPE)\
        if (mapspec->numofunits == 0)\
        {\
          *((TYPE **) mapspec->startptr) = NULL;\
        } else\
        {\
          voidptr = (((char *) ptr) + byteoffset);\
          *((TYPE **) mapspec->startptr) = voidptr;\
        }

#define WRITEACTIONWITHTYPE(TYPE)\
        gt_xfwrite(*((TYPE **) mapspecptr->startptr),\
                   mapspecptr->sizeofunit,\
                   (size_t) mapspecptr->numofunits, fp);

static uint64_t detexpectedaccordingtomapspec(const GtArrayGtMapspecification
                                              *mapspectable)
{
  uint64_t sumup = 0;
  GtMapspecification *mapspecptr;

  for (mapspecptr = mapspectable->spaceGtMapspecification;
       mapspecptr < mapspectable->spaceGtMapspecification +
                    mapspectable->nextfreeGtMapspecification; mapspecptr++)
  {
    sumup += (uint64_t) mapspecptr->sizeofunit *
             (uint64_t) mapspecptr->numofunits;
    if (sumup % GT_WORDSIZE_INBYTES > 0)
    {
      sumup += (GT_WORDSIZE_INBYTES - (sumup % GT_WORDSIZE_INBYTES));
    }
  }
  return sumup;
}

#undef SKDEBUG
#ifdef SKDEBUG
static void showmapspec(const GtMapspecification *mapspec)
{
  printf("(%s,size=%lu,elems=%lu)",
           mapspec->function,
           (unsigned long) mapspec->sizeofunit,
           mapspec->numofunits);
}
#endif

static int assigncorrecttype(GtMapspecification *mapspec,
                             void *ptr,
                             unsigned long byteoffset,
                             GtError *err)
{
  void *voidptr;
  int had_err = 0;

  gt_error_check(err);
  switch (mapspec->typespec)
  {
    case GtCharType:
      ASSIGNPTR2STARTPTR(char);
      break;
    case GtFilelengthvaluesType:
      ASSIGNPTR2STARTPTR(GtFilelengthvalues);
      break;
    case GtUcharType:
      ASSIGNPTR2STARTPTR(GtUchar);
      break;
    case Uint16Type:
      ASSIGNPTR2STARTPTR(uint16_t);
      break;
    case Uint32Type:
      ASSIGNPTR2STARTPTR(uint32_t);
      break;
    case GtUlongType:
      ASSIGNPTR2STARTPTR(unsigned long);
      break;
    case Uint64Type:
      ASSIGNPTR2STARTPTR(uint64_t);
      break;
    case GtBitsequenceType:
      ASSIGNPTR2STARTPTR(GtBitsequence);
      break;
    case GtUlongBoundType:
      ASSIGNPTR2STARTPTR(GtUlongBound);
      break;
    case GtPairBwtidxType:
      ASSIGNPTR2STARTPTR(GtPairBwtidx);
      break;
    case GtTwobitencodingType:
      ASSIGNPTR2STARTPTR(GtTwobitencoding);
      break;
    case GtSpecialcharinfoType:
      ASSIGNPTR2STARTPTR(GtSpecialcharinfo);
      break;
    case GtBitElemType:
      ASSIGNPTR2STARTPTR(BitElem);
      break;
    default:
      gt_error_set(err,"no assignment specification for size %lu",
                    (unsigned long) mapspec->sizeofunit);
      had_err = -1;
  }
  return had_err;
}

struct GtMapspec {
  GtArrayGtMapspecification mapspectable;
};

int  gt_mapspec_read(GtMapspecSetupFunc setup, void *data,
                     const GtStr *filename, unsigned long expectedsize,
                     void **mapped, GtError *err)
{
  void *mapptr;
  uint64_t expectedaccordingtomapspec;
  unsigned long byteoffset = 0;
  size_t numofbytes;
  GtMapspec *ms = gt_malloc(sizeof (GtMapspec));
  GtMapspecification *mapspecptr;
  int had_err = 0;
  unsigned long totalpadunits = 0;

  gt_error_check(err);
  GT_INITARRAY(&ms->mapspectable, GtMapspecification);
  setup(ms, data, false);

  mapptr = gt_fa_mmap_read(gt_str_get(filename), &numofbytes, err);
  if (mapptr == NULL)
  {
    had_err = -1;
  }
  *mapped = mapptr;
  if (!had_err)
  {
    if (assigncorrecttype(ms->mapspectable.spaceGtMapspecification,
                          mapptr,0,err) != 0)
    {
      had_err = -1;
    }
  }
  if (!had_err)
  {
    expectedaccordingtomapspec =
                               detexpectedaccordingtomapspec(&ms->mapspectable);
    if (expectedaccordingtomapspec != (uint64_t) numofbytes)
    {
      gt_error_set(err,"%lu bytes read from %s, but " Formatuint64_t
                         " expected",
                         (unsigned long) numofbytes,
                         gt_str_get(filename),
                         PRINTuint64_tcast(expectedaccordingtomapspec));
      had_err = -1;
    }
  }
  if (!had_err)
  {
    mapspecptr = ms->mapspectable.spaceGtMapspecification;
    gt_assert(mapspecptr != NULL);
    byteoffset = CALLCASTFUNC(uint64_t,unsigned_long,
                              (uint64_t) (mapspecptr->sizeofunit *
                                          mapspecptr->numofunits));
    if (byteoffset % (unsigned long) GT_WORDSIZE_INBYTES > 0)
    {
      size_t padunits
        = GT_WORDSIZE_INBYTES - (byteoffset % GT_WORDSIZE_INBYTES);
      byteoffset += (unsigned long) padunits;
      totalpadunits += (unsigned long) padunits;
    }
    for (mapspecptr++;
         mapspecptr < ms->mapspectable.spaceGtMapspecification +
                      ms->mapspectable.nextfreeGtMapspecification; mapspecptr++)
    {
      if (assigncorrecttype(mapspecptr,mapptr,byteoffset,err) != 0)
      {
        had_err = -1;
        break;
      }
      byteoffset = CALLCASTFUNC(uint64_t,unsigned_long,
                                (uint64_t) (byteoffset +
                                            mapspecptr->sizeofunit *
                                            mapspecptr->numofunits));
      if (byteoffset % (unsigned long) GT_WORDSIZE_INBYTES > 0)
      {
        size_t padunits
          = GT_WORDSIZE_INBYTES - (byteoffset % GT_WORDSIZE_INBYTES);
        byteoffset += (unsigned long) padunits;
        totalpadunits += (unsigned long) padunits;
      }
    }
  }
  if (!had_err)
  {
    if (expectedsize + totalpadunits != byteoffset)
    {
      gt_error_set(err,"mapping: expected file size is %lu bytes, "
                       "but file has %lu bytes",
                       expectedsize,byteoffset);
      had_err = -1;
    }
  }
  GT_FREEARRAY(&ms->mapspectable,GtMapspecification);
  gt_free(ms);
  return had_err;
}

int gt_mapspec_pad(FILE *fp, unsigned long *bytes_written,
                   unsigned long byteoffset, GT_UNUSED GtError *err)
{
  if (byteoffset % (unsigned long) GT_WORDSIZE_INBYTES > 0)
  {
    GtUchar padbuffer[GT_WORDSIZE_INBYTES-1] = {0};

    size_t padunits = GT_WORDSIZE_INBYTES - (byteoffset % GT_WORDSIZE_INBYTES);
    gt_xfwrite(padbuffer,sizeof (GtUchar),padunits,fp);
    *bytes_written = (unsigned long) padunits;
  } else
  {
    *bytes_written = 0;
  }
  return 0;
}

int gt_mapspec_write(GtMapspecSetupFunc setup, FILE *fp,
                     void *data, unsigned long expectedsize, GtError *err)
{
  GtMapspecification *mapspecptr;
  unsigned long byteoffset = 0;
  int had_err = 0;
  unsigned long totalpadunits = 0;
  unsigned long byteswritten;
  GtMapspec *ms = gt_malloc(sizeof (GtMapspec));

  gt_error_check(err);
  GT_INITARRAY(&ms->mapspectable,GtMapspecification);
  setup(ms, data, true);
  gt_assert(ms->mapspectable.spaceGtMapspecification != NULL);
  for (mapspecptr = ms->mapspectable.spaceGtMapspecification;
       mapspecptr < ms->mapspectable.spaceGtMapspecification +
                    ms->mapspectable.nextfreeGtMapspecification;
       mapspecptr++)
  {
#ifdef SKDEBUG
    printf("# %s",__func__);
    showmapspec(mapspecptr);
    printf(" at byteoffset %lu\n",byteoffset);
#endif
    if (mapspecptr->numofunits > 0)
    {
      switch (mapspecptr->typespec)
      {
        case GtCharType:
          WRITEACTIONWITHTYPE(char);
          break;
        case GtFilelengthvaluesType:
          WRITEACTIONWITHTYPE(GtFilelengthvalues);
          break;
        case GtUcharType:
          WRITEACTIONWITHTYPE(GtUchar);
          break;
        case Uint16Type:
          WRITEACTIONWITHTYPE(uint16_t);
          break;
        case Uint32Type:
          WRITEACTIONWITHTYPE(uint32_t);
          break;
        case GtUlongType:
          WRITEACTIONWITHTYPE(GtUlong);
          break;
        case Uint64Type:
          WRITEACTIONWITHTYPE(uint64_t);
          break;
        case GtBitsequenceType:
          WRITEACTIONWITHTYPE(GtBitsequence);
          break;
        case GtUlongBoundType:
          WRITEACTIONWITHTYPE(GtUlongBound);
          break;
        case GtPairBwtidxType:
          WRITEACTIONWITHTYPE(GtPairBwtidx);
          break;
        case GtTwobitencodingType:
          WRITEACTIONWITHTYPE(GtTwobitencoding);
          break;
        case GtSpecialcharinfoType:
          WRITEACTIONWITHTYPE(GtSpecialcharinfo);
          break;
        case GtBitElemType:
          WRITEACTIONWITHTYPE(BitElem);
          break;
        default:
           gt_error_set(err,"no map specification for size %lu",
                         (unsigned long) mapspecptr->sizeofunit);
           had_err = -1;
      }
    }
    if (had_err)
    {
      break;
    }
    byteoffset = CALLCASTFUNC(uint64_t,unsigned_long,
                              (uint64_t) (byteoffset +
                                          mapspecptr->sizeofunit *
                                          mapspecptr->numofunits));
    if (gt_mapspec_pad(fp,&byteswritten,byteoffset,err) != 0)
    {
      had_err = -1;
    }
    byteoffset += byteswritten;
    totalpadunits += byteswritten;
  }
  if (!had_err)
  {
    if (expectedsize + totalpadunits != byteoffset)
    {
      gt_error_set(err,"expected file size is %lu bytes, "
                       "but file has %lu bytes",
                       expectedsize,
                       byteoffset);
      had_err = -1;
    }
  }
  GT_FREEARRAY(&ms->mapspectable,GtMapspecification);
  gt_free(ms);
  return had_err;
}

#define NEWMAPSPEC(MS,PTR,TYPE,SIZE,ELEMS)\
        GT_GETNEXTFREEINARRAY(mapspecptr,&MS->mapspectable, \
                              GtMapspecification,10);\
        mapspecptr->typespec = TYPE;\
        mapspecptr->startptr = PTR;\
        mapspecptr->sizeofunit = SIZE;\
        mapspecptr->numofunits = ELEMS;\
        mapspecptr->function = __func__

void gt_mapspec_add_char_ptr(GtMapspec *mapspec, char **ptr, unsigned long n)
{
  GtMapspecification *mapspecptr;
  gt_assert(mapspec && ptr);
  NEWMAPSPEC(mapspec, ptr, GtCharType, sizeof (char), n);
}

void gt_mapspec_add_uchar_ptr(GtMapspec *mapspec, GtUchar **ptr,
                              unsigned long n)
{
  GtMapspecification *mapspecptr;
  gt_assert(mapspec && ptr);
  NEWMAPSPEC(mapspec, ptr, GtUcharType, sizeof (GtUchar), n);
}

void gt_mapspec_add_uint16_ptr(GtMapspec *mapspec, uint16_t **ptr,
                               unsigned long n)
{
  GtMapspecification *mapspecptr;
  gt_assert(mapspec && ptr);
  NEWMAPSPEC(mapspec, ptr, Uint16Type, sizeof (uint16_t), n);
}

void gt_mapspec_add_ulong_ptr(GtMapspec *mapspec, unsigned long **ptr,
                              unsigned long n)
{
  GtMapspecification *mapspecptr;
  gt_assert(mapspec && ptr);
  NEWMAPSPEC(mapspec, ptr, GtUlongType, sizeof (unsigned long), n);
}

void gt_mapspec_add_ulongbound_ptr(GtMapspec *mapspec, GtUlongBound **ptr,
                                   unsigned long n)
{
  GtMapspecification *mapspecptr;
  gt_assert(mapspec && ptr);
  NEWMAPSPEC(mapspec, ptr, GtUlongBoundType, sizeof (GtUlongBound), n);
}

void gt_mapspec_add_uint32_ptr(GtMapspec *mapspec, uint32_t **ptr,
                               unsigned long n)
{
  GtMapspecification *mapspecptr;
  gt_assert(mapspec && ptr);
  NEWMAPSPEC(mapspec, ptr, Uint32Type, sizeof (uint32_t), n);
}

void gt_mapspec_add_uint64_ptr(GtMapspec *mapspec, uint64_t **ptr,
                               unsigned long n)
{
  GtMapspecification *mapspecptr;
  gt_assert(mapspec && ptr);
  NEWMAPSPEC(mapspec, ptr, Uint64Type, sizeof (uint64_t), n);
}

void gt_mapspec_add_bitsequence_ptr(GtMapspec *mapspec, GtBitsequence **ptr,
                                    unsigned long n)
{
  GtMapspecification *mapspecptr;
  gt_assert(mapspec && ptr);
  NEWMAPSPEC(mapspec, ptr, GtBitsequenceType, sizeof (GtBitsequence), n);
}
void gt_mapspec_add_twobitencoding_ptr(GtMapspec *mapspec,
                                       GtTwobitencoding **ptr, unsigned long n)
{
  GtMapspecification *mapspecptr;
  gt_assert(mapspec && ptr);
  NEWMAPSPEC(mapspec, ptr, GtTwobitencodingType, sizeof (GtTwobitencoding), n);
}

void gt_mapspec_add_specialcharinfo_ptr(GtMapspec *mapspec,
                                        GtSpecialcharinfo **ptr,
                                        unsigned long n)
{
  GtMapspecification *mapspecptr;
  gt_assert(mapspec && ptr);
  NEWMAPSPEC(mapspec, ptr, GtSpecialcharinfoType, sizeof (GtSpecialcharinfo),
             n);
}

void gt_mapspec_add_bitelem_ptr(GtMapspec *mapspec, BitElem **ptr,
                                unsigned long n)
{
  GtMapspecification *mapspecptr;
  gt_assert(mapspec && ptr);
  NEWMAPSPEC(mapspec, ptr, GtBitElemType, sizeof (BitElem), n);
}

void gt_mapspec_add_filelengthvalues_ptr(GtMapspec *mapspec,
                                         GtFilelengthvalues **ptr,
                                         unsigned long n)
{
  GtMapspecification *mapspecptr;
  gt_assert(mapspec && ptr);
  NEWMAPSPEC(mapspec, ptr, GtFilelengthvaluesType,
             sizeof (GtFilelengthvalues), n);
}

void gt_mapspec_add_pairbwtindex_ptr(GtMapspec *mapspec, GtPairBwtidx **ptr,
                                     unsigned long n)
{
  GtMapspecification *mapspecptr;
  gt_assert(mapspec && ptr);
  NEWMAPSPEC(mapspec, ptr, GtPairBwtidxType, sizeof (GtPairBwtidx), n);
}
