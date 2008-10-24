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

#include "core/str_api.h"
#include "core/unused_api.h"
#include "core/minmax.h"
#include "core/symboldef.h"
#include "divmodmul.h"
#include "defined-types.h"
#include "intbits.h"
#include "tyr-search.h"

typedef struct
{
  void *mappedmbdfileptr;
  const GtStr *indexfilename;
  unsigned int prefixlength;
  unsigned long numofcodes,
                *boundisdefined,
                *bounds;
} Tyrbckinfo;

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

void splitmerinterval(Tyrbckinfo *tyrbckinfo,
                             const Uchar *mertable,
                             const Uchar *lastmer,
                             unsigned long merbytes)
{
  const Uchar *rightbound, *leftptr, *rightptr;
  unsigned long code, leftcode, rightcode;

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
    tyrbckinfo->bounds[leftcode]
      = (unsigned long) (leftptr - mertable);
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

int constructmerbuckets(GT_UNUSED const GtStr *inputindex,
                        GT_UNUSED const Definedunsignedint *prefixlength)
{
  return 0;
}
