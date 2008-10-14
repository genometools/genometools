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
#include "core/str_api.h"
#include "core/minmax.h"
#include "core/symboldef.h"
#include "core/fa.h"
#include "divmodmul.h"
#include "defined-types.h"
#include "intbits.h"
#include "intbits-tab.h"
#include "tyr-map.h"
#include "tyr-mersplit.h"
#include "opensfxfile.h"

#define BUCKETSUFFIX                     ".mbd"
#define ISBOUNDDEFINED(UDB,IDX)          ISIBITSET(UDB,IDX)
#define SETDEFINEDBOUND(UDB,IDX)         SETIBIT(UDB,IDX)

struct Tyrbckinfo
{
  void *mappedmbdfileptr;
  unsigned int prefixlength;
  unsigned long numofcodes,
                *boundisdefined,
                *bounds;
};

static unsigned long extractprefixbytecode(unsigned long merbytes,
                                           unsigned int prefixlength,
                                           const Uchar *bytecode)
{
  unsigned long idx, code = 0;

  for (idx=0; idx < MIN((unsigned long) sizeof (unsigned long),merbytes); idx++)
  {
    code = (code << 8) | bytecode[idx];
    if (MULT4(idx+1) == (unsigned long) prefixlength)
    {
      break;
    }
    if (MULT4(idx+1) > (unsigned long) prefixlength)
    {
      code >>= MULT2(MULT4(idx+1) - prefixlength);
      break;
    }
  }
  return code;
}

static const Uchar *findrightmostmer(unsigned long merbytes,
                                     unsigned int prefixlength,
                                     unsigned long code,
                                     const Uchar *leftptr,
                                     const Uchar *rightptr)
{
  unsigned long len, midcode;
  const Uchar *midptr;

  while (leftptr + merbytes < rightptr)
  {
    len = (unsigned long) (rightptr-leftptr)/MULT2(merbytes);
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
  const Uchar *rightbound, *leftptr, *rightptr;
  unsigned long code, leftcode, rightcode;
  const Uchar *mertable, *lastmer;
  unsigned long merbytes;

  mertable = tyrindex_mertable(tyrindex);
  lastmer = tyrindex_lastmer(tyrindex);
  merbytes = tyrindex_merbytes(tyrindex);
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

int constructmerbuckets(const GtStr *inputindex,
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
  tyrindex = tyrindex_new(inputindex,err);
  if (tyrindex == NULL)
  {
    haserr = true;
  }
  if (!haserr && tyrindex != NULL && !tyrindex_isempty(tyrindex))
  {
    if (determinetyrbckpfxlen(&tyrbckinfo.prefixlength,
                              tyrindex,
                              callprefixlength,
                              err) != 0)
    {
      haserr = true;
    }
  }
  if (!haserr && tyrindex != NULL && !tyrindex_isempty(tyrindex))
  {
    printf("# construct mer buckets for prefixlength %u\n",
            tyrbckinfo.prefixlength);
    tyrbckinfo.numofcodes
      = (unsigned long) pow((double) tyrindex_alphasize(tyrindex),
                            (double) tyrbckinfo.prefixlength);
    tyrbckinfo.mappedmbdfileptr = NULL;
    printf("# numofcodes = %lu\n",tyrbckinfo.numofcodes);
    tyrindex_show(tyrindex);
    ALLOCASSIGNSPACE(tyrbckinfo.bounds,NULL,unsigned long,
                     tyrbckinfo.numofcodes+1);
    INITBITTAB(tyrbckinfo.boundisdefined,tyrbckinfo.numofcodes+1);
    splitmerinterval(&tyrbckinfo,tyrindex);
    bucketfp = opensfxfile(inputindex,BUCKETSUFFIX,"wb",err);
    if (bucketfp == NULL)
    {
      haserr = true;
    }
  }
  if (!haserr && tyrindex != NULL && !tyrindex_isempty(tyrindex))
  {
    unsigned long pl_long = (unsigned long) tyrbckinfo.prefixlength;
    gt_assert(bucketfp != NULL);
    if (fwrite(&pl_long,
               sizeof (pl_long),
               (size_t) 1,
               bucketfp) != (size_t) 1)
    {
      gt_error_set(err,"cannot write 1 item of size %u: errormsg=\"%s\"",
                   (unsigned int) sizeof (pl_long),
                   strerror(errno));
      haserr = true;
    }
  }
  if (!haserr && tyrindex != NULL && !tyrindex_isempty(tyrindex))
  {
    gt_assert(bucketfp != NULL);
    if (fwrite(tyrbckinfo.bounds,
               sizeof (*tyrbckinfo.bounds),
               (size_t) (tyrbckinfo.numofcodes+1),
               bucketfp) != (size_t) (tyrbckinfo.numofcodes+1))
    {
      haserr = true;
    }
  }
  if (!haserr && tyrindex != NULL && !tyrindex_isempty(tyrindex))
  {
    gt_assert(bucketfp != NULL);
    if (fwrite(tyrbckinfo.boundisdefined,
               sizeof (*tyrbckinfo.boundisdefined),
               NUMOFINTSFORBITS(tyrbckinfo.numofcodes+1),
               bucketfp) != NUMOFINTSFORBITS(tyrbckinfo.numofcodes+1))
    {
      haserr = true;
    }
  }
  gt_fa_xfclose(bucketfp);
  if (tyrindex != NULL)
  {
    tyrindex_delete(&tyrindex);
  }
  FREESPACE(tyrbckinfo.bounds);
  FREESPACE(tyrbckinfo.boundisdefined);
  return haserr ? -1 : 0;
}

Tyrbckinfo *tyrbckinfo_new(const GtStr *tyrindexname,unsigned int alphasize,
                           GtError *err)
{
  size_t numofbytes, expectedsize;
  Tyrbckinfo *tyrbckinfo;
  bool haserr = false;

  ALLOCASSIGNSPACE(tyrbckinfo,NULL,Tyrbckinfo,1);
  tyrbckinfo->mappedmbdfileptr = genericmaponlytable(tyrindexname,
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
      = (unsigned long) pow((double) alphasize,
                            (double) tyrbckinfo->prefixlength);
    expectedsize
      = sizeof (unsigned long) *
        (1UL +
         (tyrbckinfo->numofcodes+1) +
         NUMOFINTSFORBITS(tyrbckinfo->numofcodes + 1));
    gt_assert(expectedsize == numofbytes);
    tyrbckinfo->bounds = ((unsigned long *) tyrbckinfo->mappedmbdfileptr) + 1;
    tyrbckinfo->boundisdefined
      = tyrbckinfo->bounds + tyrbckinfo->numofcodes + 1;
  }
  if (haserr)
  {
    FREESPACE(tyrbckinfo);
    return NULL;
  }
  return tyrbckinfo;
}

void tyrbckinfo_delete(Tyrbckinfo **tyrbckinfoptr)
{
  Tyrbckinfo *tyrbckinfo = *tyrbckinfoptr;

  gt_fa_xmunmap(tyrbckinfo->mappedmbdfileptr);
  tyrbckinfo->mappedmbdfileptr = NULL;
  FREESPACE(tyrbckinfo);
  *tyrbckinfoptr = NULL;
}
