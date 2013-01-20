/*
  Copyright (c) 2008 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2008 Center for Bioinformatics, University of Hamburg

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

#include <math.h>
#include <errno.h>
#include "core/defined-types.h"
#include "core/divmodmul.h"
#include "core/fa.h"
#include "core/intbits.h"
#include "core/mathsupport.h"
#include "core/minmax.h"
#include "core/str_api.h"
#include "core/types_api.h"
#include "core/unused_api.h"
#include "core/xansi_api.h"
#include "core/ma_api.h"
#include "tyr-map.h"
#include "tyr-mersplit.h"

#define BUCKETSUFFIX                     ".mbd"
#define MAXUCHARVALUEWITHBITS(BITNUM)    ((1 << (BITNUM)) - 1)
#define ISBOUNDDEFINED(UDB,IDX)          GT_ISIBITSET(UDB,IDX)
#define SETDEFINEDBOUND(UDB,IDX)         GT_SETIBIT(UDB,IDX)

struct Tyrbckinfo
{
  void *mappedmbdfileptr;
  unsigned int prefixlength;
  unsigned long numofcodes,
                *boundisdefined,
                *bounds;
  GtUchar remainmask;
};

typedef struct
{
  const GtUchar *leftmer,
              *rightmer;
} Merbounds;

static unsigned long extractprefixbytecode(unsigned long merbytes,
                                           unsigned int prefixlength,
                                           const GtUchar *bytecode)
{
  unsigned long idx, code = 0;

  for (idx=0; idx < MIN((unsigned long) sizeof (unsigned long),merbytes); idx++)
  {
    code = (code << 8) | bytecode[idx];
    if (GT_MULT4(idx+1) == (unsigned long) prefixlength)
    {
      break;
    }
    if (GT_MULT4(idx+1) > (unsigned long) prefixlength)
    {
      code >>= GT_MULT2(GT_MULT4(idx+1) - prefixlength);
      break;
    }
  }
  return code;
}

static GtUchar extractremainingbytes(const GtUchar remainmask,
                                   unsigned long byteoffset,
                                   const GtUchar *bytecode)
{
  return bytecode[byteoffset] & remainmask;
}

static const GtUchar *remainingleftmost(unsigned long merbytes,
                                      unsigned long byteoffset,
                                      GtUchar remainmask,
                                      GtUchar code,
                                      const GtUchar *leftptr,
                                      const GtUchar *rightptr)
{
  unsigned long len;
  GtUchar midcode;
  const GtUchar *midptr;

  while (leftptr + merbytes < rightptr)
  {
    len = (unsigned long) (rightptr-leftptr)/GT_MULT2(merbytes);
    midptr = leftptr + merbytes * len;
    midcode = extractremainingbytes(remainmask,byteoffset,midptr);
    if (code <= midcode)
    {
      rightptr = midptr;
    } else
    {
      leftptr = midptr;
    }
  }
  return rightptr;
}

static const GtUchar *remainingrightmost(unsigned long merbytes,
                                       unsigned long byteoffset,
                                       GtUchar remainmask,
                                       GtUchar code,
                                       const GtUchar *leftptr,
                                       const GtUchar *rightptr)
{
  unsigned long len;
  GtUchar midcode;
  const GtUchar *midptr;

  while (leftptr + merbytes < rightptr)
  {
    len = (unsigned long) (rightptr-leftptr)/GT_MULT2(merbytes);
    midptr = leftptr + merbytes * len;
    midcode = extractremainingbytes(remainmask,byteoffset,midptr);
    if (code >= midcode)
    {
      leftptr = midptr;
    } else
    {
      rightptr = midptr;
    }
  }
  return leftptr;
}

static bool remainadvance(Merbounds *merbounds,
                          unsigned long merbytes,
                          unsigned long byteoffset,
                          GtUchar remainmask,
                          const GtUchar *searchbytecode,
                          const GtUchar *leftptr,
                          const GtUchar *rightptr)
{
  GtUchar scode, scodeleft, scoderight;

  scode = extractremainingbytes(remainmask,byteoffset,searchbytecode);
  scodeleft = extractremainingbytes(remainmask,byteoffset,leftptr);
  if (scode > scodeleft)
  {
    scoderight = extractremainingbytes(remainmask,byteoffset,rightptr);
    if (scode > scoderight)
    {
      return false;
    }
    merbounds->leftmer = remainingleftmost(merbytes,
                                           byteoffset,
                                           remainmask,
                                           scode,
                                           leftptr,
                                           rightptr);
  }
  if (scode < scodeleft)
  {
    return false;
  }
  if (scode == scodeleft)
  {
    merbounds->leftmer = leftptr;
  }
  scoderight = extractremainingbytes(remainmask,byteoffset,rightptr);
  if (scode >= scoderight)
  {
    merbounds->rightmer = rightptr;
  } else
  {
    merbounds->rightmer = remainingrightmost(merbytes,
                                             byteoffset,
                                             remainmask,
                                             scode,
                                             leftptr,
                                             rightptr);
  }
  return true;
}

const GtUchar *gt_searchinbuckets(const Tyrindex *tyrindex,
                             const Tyrbckinfo *tyrbckinfo,
                             const GtUchar *bytecode)
{
  const GtUchar *result;
  unsigned long prefixcode, leftbound, merbytes;

  gt_assert(tyrbckinfo != NULL);
  merbytes = gt_tyrindex_merbytes(tyrindex);
  prefixcode = extractprefixbytecode(merbytes,
                                     tyrbckinfo->prefixlength,
                                     bytecode);
  leftbound = tyrbckinfo->bounds[prefixcode];
  if (ISBOUNDDEFINED(tyrbckinfo->boundisdefined,prefixcode))
  {
    const GtUchar *mertable;
    unsigned long rightbound;

    mertable = gt_tyrindex_mertable(tyrindex);
    rightbound = tyrbckinfo->bounds[prefixcode+1] - merbytes;
    if (GT_MOD4(tyrbckinfo->prefixlength) == 0)
    {
      result = gt_tyrindex_binmersearch(tyrindex,
                                     (unsigned long)
                                     GT_DIV4(tyrbckinfo->prefixlength),
                                     bytecode,
                                     mertable + leftbound,
                                     mertable + rightbound);
    } else
    {
      Merbounds merbounds;

      merbounds.leftmer = merbounds.rightmer = NULL;
      if (!remainadvance(&merbounds,
                         merbytes,
                         (unsigned long) GT_DIV4(tyrbckinfo->prefixlength),
                         tyrbckinfo->remainmask,
                         bytecode,
                         mertable + leftbound,
                         mertable + rightbound) ||
         merbounds.leftmer == NULL ||
         merbounds.leftmer > merbounds.rightmer)
      {
        result = NULL;
      } else
      {
        result = gt_tyrindex_binmersearch(tyrindex,
                                       1UL + (unsigned long)
                                             GT_DIV4(tyrbckinfo->prefixlength),
                                       bytecode,
                                       merbounds.leftmer,
                                       merbounds.rightmer);
      }
    }
  } else
  {
    result = NULL;
  }
  return result;
}

static const GtUchar *findrightmostmer(unsigned long merbytes,
                                     unsigned int prefixlength,
                                     unsigned long code,
                                     const GtUchar *leftptr,
                                     const GtUchar *rightptr)
{
  unsigned long len, midcode;
  const GtUchar *midptr;

  while (leftptr + merbytes < rightptr)
  {
    len = (unsigned long) (rightptr-leftptr)/GT_MULT2(merbytes);
    midptr = leftptr + merbytes * len;
    midcode = extractprefixbytecode(merbytes,prefixlength,midptr);
    if (midcode > code)
    {
      rightptr = midptr;
    } else
    {
      leftptr = midptr;
    }
  }
  return leftptr;
}

static void splitmerinterval(Tyrbckinfo *tyrbckinfo,
                             const Tyrindex *tyrindex)

{
  const GtUchar *rightbound, *leftptr, *rightptr, *mertable, *lastmer;
  unsigned long code, leftcode, rightcode, merbytes;

  mertable = gt_tyrindex_mertable(tyrindex);
  lastmer = gt_tyrindex_lastmer(tyrindex);
  merbytes = gt_tyrindex_merbytes(tyrindex);
  leftptr = mertable;
  rightptr = lastmer;
  rightcode = extractprefixbytecode(merbytes,
                                    tyrbckinfo->prefixlength,
                                    rightptr);
  while (true)
  {
    leftcode = extractprefixbytecode(merbytes,
                                     tyrbckinfo->prefixlength,
                                     leftptr);
    tyrbckinfo->bounds[leftcode] = (unsigned long) (leftptr - mertable);
    SETDEFINEDBOUND(tyrbckinfo->boundisdefined,leftcode);
    if (leftcode == rightcode)
    {
      break;
    }
    rightbound = findrightmostmer(merbytes,
                                  tyrbckinfo->prefixlength,
                                  leftcode,leftptr,rightptr);
    leftptr = rightbound + merbytes;
  }
  tyrbckinfo->bounds[tyrbckinfo->numofcodes]
    = (unsigned long) (lastmer + merbytes - mertable);
  SETDEFINEDBOUND(tyrbckinfo->boundisdefined,tyrbckinfo->numofcodes);
  for (code = tyrbckinfo->numofcodes - 1; /* Nothing */ ; code--)
  {
    if (!ISBOUNDDEFINED(tyrbckinfo->boundisdefined,code))
    {
      tyrbckinfo->bounds[code] = tyrbckinfo->bounds[code+1];
    }
    if (code == 0)
    {
      break;
    }
  }
}

int gt_constructmerbuckets(const char *inputindex,
                           const Definedunsignedint *callprefixlength,
                           GtError *err)
{
  Tyrindex *tyrindex;
  Tyrbckinfo tyrbckinfo;
  FILE *bucketfp = NULL;
  bool haserr = false;

  gt_error_check(err);
  tyrbckinfo.bounds = NULL;
  tyrbckinfo.boundisdefined = NULL;
  tyrindex = gt_tyrindex_new(inputindex,err);
  if (tyrindex == NULL)
  {
    haserr = true;
  }
  if (!haserr && tyrindex != NULL && !gt_tyrindex_isempty(tyrindex))
  {
    if (gt_determinetyrbckpfxlen(&tyrbckinfo.prefixlength,
                              tyrindex,
                              callprefixlength,
                              err) != 0)
    {
      haserr = true;
    }
  }
  if (!haserr && tyrindex != NULL && !gt_tyrindex_isempty(tyrindex))
  {
    printf("# construct mer buckets for prefixlength %u\n",
            tyrbckinfo.prefixlength);
    tyrbckinfo.numofcodes
      = gt_power_for_small_exponents(gt_tyrindex_alphasize(tyrindex),
                                     tyrbckinfo.prefixlength);
    tyrbckinfo.mappedmbdfileptr = NULL;
    printf("# numofcodes = %lu\n",tyrbckinfo.numofcodes);
    gt_tyrindex_show(tyrindex);
    tyrbckinfo.bounds = gt_malloc(sizeof *tyrbckinfo.bounds
                                  * (tyrbckinfo.numofcodes+1));
    GT_INITBITTAB(tyrbckinfo.boundisdefined,tyrbckinfo.numofcodes+1);
    splitmerinterval(&tyrbckinfo,tyrindex);
    bucketfp = gt_fa_fopen_with_suffix(inputindex,BUCKETSUFFIX,"wb",err);
    if (bucketfp == NULL)
    {
      haserr = true;
    }
  }
  if (!haserr && tyrindex != NULL && !gt_tyrindex_isempty(tyrindex))
  {
    unsigned long pl_long = (unsigned long) tyrbckinfo.prefixlength;
    gt_assert(bucketfp != NULL);
    gt_xfwrite(&pl_long, sizeof (pl_long), (size_t) 1, bucketfp);
  }
  if (!haserr && tyrindex != NULL && !gt_tyrindex_isempty(tyrindex))
  {
    gt_assert(bucketfp != NULL);
    gt_xfwrite(tyrbckinfo.bounds, sizeof (*tyrbckinfo.bounds),
               (size_t) (tyrbckinfo.numofcodes+1), bucketfp);
  }
  if (!haserr && tyrindex != NULL && !gt_tyrindex_isempty(tyrindex))
  {
    gt_assert(bucketfp != NULL);
    gt_xfwrite(tyrbckinfo.boundisdefined, sizeof (*tyrbckinfo.boundisdefined),
               GT_NUMOFINTSFORBITS(tyrbckinfo.numofcodes+1), bucketfp);
  }
  gt_fa_xfclose(bucketfp);
  if (tyrindex != NULL)
  {
    gt_tyrindex_delete(&tyrindex);
  }
  gt_free(tyrbckinfo.bounds);
  gt_free(tyrbckinfo.boundisdefined);
  return haserr ? -1 : 0;
}

Tyrbckinfo *gt_tyrbckinfo_new(const char *tyrindexname,unsigned int alphasize,
                              GtError *err)
{
  size_t numofbytes;
  Tyrbckinfo *tyrbckinfo;
  bool haserr = false;

  tyrbckinfo = gt_malloc(sizeof *tyrbckinfo);
  tyrbckinfo->mappedmbdfileptr = gt_fa_mmap_read_with_suffix(tyrindexname,
                                                          BUCKETSUFFIX,
                                                          &numofbytes,err);
  if (tyrbckinfo->mappedmbdfileptr == NULL)
  {
    haserr = true;
  }
  if (!haserr)
  {
    unsigned long pl_long;

    gt_assert(tyrbckinfo->mappedmbdfileptr != NULL);
    pl_long = *((unsigned long *) tyrbckinfo->mappedmbdfileptr);
    tyrbckinfo->prefixlength = (unsigned int) pl_long;
    tyrbckinfo->numofcodes
      = gt_power_for_small_exponents(alphasize,tyrbckinfo->prefixlength);
    /*check if numofbytes == expected size*/
    gt_assert(numofbytes == sizeof (unsigned long) *
                            (1UL +
                             (tyrbckinfo->numofcodes+1) +
                             GT_NUMOFINTSFORBITS(tyrbckinfo->numofcodes + 1)));
    tyrbckinfo->bounds = ((unsigned long *) tyrbckinfo->mappedmbdfileptr) + 1;
    tyrbckinfo->boundisdefined
      = tyrbckinfo->bounds + tyrbckinfo->numofcodes + 1;
    if (tyrbckinfo->prefixlength > 0 && GT_MOD4(tyrbckinfo->prefixlength) > 0)
    {
      tyrbckinfo->remainmask
        = (GtUchar) MAXUCHARVALUEWITHBITS(GT_MULT2(
                                       4U - GT_MOD4(tyrbckinfo->prefixlength)));
    }
  }
  if (haserr)
  {
    gt_free(tyrbckinfo);
    return NULL;
  }
  return tyrbckinfo;
}

void gt_tyrbckinfo_delete(Tyrbckinfo **tyrbckinfoptr)
{
  Tyrbckinfo *tyrbckinfo = *tyrbckinfoptr;

  gt_fa_xmunmap(tyrbckinfo->mappedmbdfileptr);
  tyrbckinfo->mappedmbdfileptr = NULL;
  gt_free(tyrbckinfo);
  *tyrbckinfoptr = NULL;
}
