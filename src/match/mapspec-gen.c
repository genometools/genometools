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

#include <errno.h>
#include <string.h>
#include "core/error.h"
#include "core/fa.h"
#include "core/str.h"
#include "bitpack-itf.h"
#include "ushort-def.h"
#include "fmi-bwtbound.h"
#include "intbits.h"
#include "safecast-gen.h"
#include "mapspec-def.h"
#include "format64.h"
#include "stamp.h"

#define ASSIGNPTR2STARTPTR(TYPE)\
        if (mapspec->numofunits == 0)\
        {\
          *((TYPE **) mapspec->startptr) = NULL;\
        } else\
        {\
          voidptr = (((char *) ptr) + byteoffset);\
          *((TYPE **) mapspec->startptr) = voidptr;\
        }

#define ALIGNSIZE sizeof (void *)

#define WRITEACTIONWITHTYPE(TYPE)\
        if (fwrite(*((TYPE **) mapspecptr->startptr),\
                   mapspecptr->sizeofunit,\
                   (size_t) mapspecptr->numofunits, fp) !=\
                    (size_t) mapspecptr->numofunits)\
        {\
          gt_error_set(err,"cannot write %lu items of size %u: "\
                            "errormsg=\"%s\"",\
                        (unsigned long) mapspecptr->numofunits,\
                        (unsigned int) mapspecptr->sizeofunit,\
                        strerror(errno));\
          haserr = true;\
        }

static uint64_t detexpectedaccordingtomapspec(const GtArrayMapspecification
                                              *mapspectable)
{
  uint64_t sumup = 0;
  Mapspecification *mapspecptr;

  for (mapspecptr = mapspectable->spaceMapspecification;
       mapspecptr < mapspectable->spaceMapspecification +
                    mapspectable->nextfreeMapspecification; mapspecptr++)
  {
    sumup += (uint64_t) mapspecptr->sizeofunit *
             (uint64_t) mapspecptr->numofunits;
    if (sumup % ALIGNSIZE > 0)
    {
      sumup += (ALIGNSIZE - (sumup % ALIGNSIZE));
    }
  }
  return sumup;
}

#undef SKDEBUG
#ifdef SKDEBUG
static void showmapspec(const Mapspecification *mapspec)
{
  printf("(%s,size=%lu,elems=%lu)",
           mapspec->name,
           (unsigned long) mapspec->sizeofunit,
           mapspec->numofunits);
}
#endif

static int assigncorrecttype(Mapspecification *mapspec,
                             void *ptr,
                             unsigned long byteoffset,
                             GtError *err)
{
  void *voidptr;
  bool haserr = false;

  gt_error_check(err);
  switch (mapspec->typespec)
  {
    case CharType:
      ASSIGNPTR2STARTPTR(char);
      break;
    case GtUcharType:
      ASSIGNPTR2STARTPTR(GtUchar);
      break;
    case UshortType:
      ASSIGNPTR2STARTPTR(Ushort);
      break;
    case Uint32Type:
      ASSIGNPTR2STARTPTR(uint32_t);
      break;
    case UnsignedlongType:
      ASSIGNPTR2STARTPTR(Unsignedlong);
      break;
    case Uint64Type:
      ASSIGNPTR2STARTPTR(uint64_t);
      break;
    case BitsequenceType:
      ASSIGNPTR2STARTPTR(Bitsequence);
      break;
    case SeqposType:
      ASSIGNPTR2STARTPTR(Seqpos);
      break;
    case BwtboundType:
      ASSIGNPTR2STARTPTR(Bwtbound);
      break;
    case PairBwtidxType:
      ASSIGNPTR2STARTPTR(PairBwtidx);
      break;
    case TwobitencodingType:
      ASSIGNPTR2STARTPTR(Twobitencoding);
      break;
    case SpecialcharinfoType:
      ASSIGNPTR2STARTPTR(Specialcharinfo);
      break;
    case BitElemType:
      ASSIGNPTR2STARTPTR(BitElem);
      break;
    default:
      gt_error_set(err,"no assignment specification for size %lu",
                    (unsigned long) mapspec->sizeofunit);
      haserr = true;
  }
  return haserr ? -1 : 0;
}

 DECLARESAFECASTFUNCTION(uint64_t,uint64_t,unsigned long,unsigned_long)

int fillmapspecstartptr(Assignmapspec assignmapspec,
                        void **mappeduserptr,
                        void *assignmapinfo,
                        const GtStr *tmpfilename,
                        unsigned long expectedsize,
                        GtError *err)
{
  void *mapptr;
  uint64_t expectedaccordingtomapspec;
  unsigned long byteoffset = 0;
  size_t numofbytes;
  GtArrayMapspecification mapspectable;
  Mapspecification *mapspecptr;
  bool haserr = false;
  unsigned long totalpadunits = 0;

  gt_error_check(err);
  GT_INITARRAY(&mapspectable,Mapspecification);
  assignmapspec(&mapspectable,assignmapinfo,false);
  mapptr = gt_fa_mmap_read(gt_str_get(tmpfilename), &numofbytes, err);
  if (mapptr == NULL)
  {
    haserr = true;
  }
  *mappeduserptr = mapptr;
  if (!haserr)
  {
    if (assigncorrecttype(mapspectable.spaceMapspecification,
                          mapptr,0,err) != 0)
    {
      haserr = true;
    }
  }
  if (!haserr)
  {
    expectedaccordingtomapspec = detexpectedaccordingtomapspec(&mapspectable);
    if (expectedaccordingtomapspec != (uint64_t) numofbytes)
    {
      gt_error_set(err,"%lu bytes read from %s, but " Formatuint64_t
                         " expected",
                         (unsigned long) numofbytes,
                         gt_str_get(tmpfilename),
                         PRINTuint64_tcast(expectedaccordingtomapspec));
      haserr = true;
    }
  }
  if (!haserr)
  {
    mapspecptr = mapspectable.spaceMapspecification;
    gt_assert(mapspecptr != NULL);
    byteoffset = CALLCASTFUNC(uint64_t,unsigned_long,
                              (uint64_t) (mapspecptr->sizeofunit *
                                          mapspecptr->numofunits));
    if (byteoffset % (unsigned long) ALIGNSIZE > 0)
    {
      size_t padunits = ALIGNSIZE - (byteoffset % ALIGNSIZE);
      byteoffset += (unsigned long) padunits;
      totalpadunits += (unsigned long) padunits;
    }
    for (mapspecptr++;
         mapspecptr < mapspectable.spaceMapspecification +
                      mapspectable.nextfreeMapspecification; mapspecptr++)
    {
      if (assigncorrecttype(mapspecptr,mapptr,byteoffset,err) != 0)
      {
        haserr = true;
        break;
      }
      byteoffset = CALLCASTFUNC(uint64_t,unsigned_long,
                                (uint64_t) (byteoffset +
                                            mapspecptr->sizeofunit *
                                            mapspecptr->numofunits));
      if (byteoffset % (unsigned long) ALIGNSIZE > 0)
      {
        size_t padunits = ALIGNSIZE - (byteoffset % ALIGNSIZE);
        byteoffset += (unsigned long) padunits;
        totalpadunits += (unsigned long) padunits;
      }
    }
  }
  if (!haserr)
  {
    if (expectedsize + totalpadunits != byteoffset)
    {
      gt_error_set(err,"mapping: expected size of index is %lu bytes, "
                        "but index has %lu bytes",
                        expectedsize,byteoffset);
      haserr = true;
    }
  }
  GT_FREEARRAY(&mapspectable,Mapspecification);
  return haserr ? -1 : 0;
}

int flushtheindex2file(FILE *fp,
                       Assignmapspec assignmapspec,
                       void *assignmapinfo,
                       unsigned long expectedsize,
                       GtError *err)
{
  GtArrayMapspecification mapspectable;
  Mapspecification *mapspecptr;
  unsigned long byteoffset = 0;
  bool haserr = false;
  GtUchar padbuffer[ALIGNSIZE-1] = {0};
  unsigned long totalpadunits = 0;

  gt_error_check(err);
  GT_INITARRAY(&mapspectable,Mapspecification);
  assignmapspec(&mapspectable,assignmapinfo,true);
  gt_assert(mapspectable.spaceMapspecification != NULL);
  for (mapspecptr = mapspectable.spaceMapspecification;
       mapspecptr < mapspectable.spaceMapspecification +
                    mapspectable.nextfreeMapspecification;
       mapspecptr++)
  {
#ifdef SKDEBUG
    printf("# flushtheindex2file");
    showmapspec(mapspecptr);
    printf(" at byteoffset %lu\n",byteoffset);
#endif
    if (mapspecptr->numofunits > 0)
    {
      switch (mapspecptr->typespec)
      {
        case CharType:
          WRITEACTIONWITHTYPE(char);
          break;
        case GtUcharType:
          WRITEACTIONWITHTYPE(GtUchar);
          break;
        case UshortType:
          WRITEACTIONWITHTYPE(Ushort);
          break;
        case Uint32Type:
          WRITEACTIONWITHTYPE(uint32_t);
          break;
        case UnsignedlongType:
          WRITEACTIONWITHTYPE(Unsignedlong);
          break;
        case Uint64Type:
          WRITEACTIONWITHTYPE(uint64_t);
          break;
        case BitsequenceType:
          WRITEACTIONWITHTYPE(Bitsequence);
          break;
        case SeqposType:
          WRITEACTIONWITHTYPE(Seqpos);
          break;
        case BwtboundType:
          WRITEACTIONWITHTYPE(Bwtbound);
          break;
        case PairBwtidxType:
          WRITEACTIONWITHTYPE(PairBwtidx);
          break;
        case TwobitencodingType:
          WRITEACTIONWITHTYPE(Twobitencoding);
          break;
        case SpecialcharinfoType:
          WRITEACTIONWITHTYPE(Specialcharinfo);
          break;
        case BitElemType:
          WRITEACTIONWITHTYPE(BitElem);
          break;
        default:
           gt_error_set(err,"no map specification for size %lu",
                         (unsigned long) mapspecptr->sizeofunit);
           haserr = true;
      }
    }
    if (haserr)
    {
      break;
    }
    byteoffset = CALLCASTFUNC(uint64_t,unsigned_long,
                              (uint64_t) (byteoffset +
                                          mapspecptr->sizeofunit *
                                          mapspecptr->numofunits));
    if (byteoffset % (unsigned long) ALIGNSIZE > 0)
    {
      size_t padunits = ALIGNSIZE - (byteoffset % ALIGNSIZE);
      if (fwrite(padbuffer,
                sizeof (GtUchar),padunits,fp) != padunits)
      {
        gt_error_set(err,"cannot write %lu items of size %u: "
                          "errormsg=\"%s\"",
                           (unsigned long) padunits,
                           (unsigned int) sizeof (GtUchar),
                           strerror(errno));
        haserr = true;
      }
      byteoffset += (unsigned long) padunits;
      totalpadunits += (unsigned long) padunits;
    }
  }
  if (!haserr)
  {
    if (expectedsize + totalpadunits != byteoffset)
    {
      gt_error_set(err,"flushindex: expected size of index is %lu bytes, "
                        "but index has %lu bytes",
                        expectedsize,
                        byteoffset);
      haserr = true;
    }
  }
  GT_FREEARRAY(&mapspectable,Mapspecification);
  return haserr ? -1 : 0;
}
