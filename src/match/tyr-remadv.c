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

#include "core/symboldef.h"
#include "divmodmul.h"
#include "tyr-basic.h"

static Uchar extractremainingbytes(const Uchar remainmask,
                                   unsigned long byteoffset,
                                   const Uchar *bytecode)
{
  return bytecode[byteoffset] & remainmask;
}

static const Uchar *remainingleftmost(unsigned long merbytes,
                                      unsigned long byteoffset,
                                      Uchar remainmask,
                                      Uchar code,
                                      const Uchar *leftptr,
                                      const Uchar *rightptr)
{
  unsigned long len;
  Uchar midcode;
  const Uchar *midptr;

  while (leftptr + merbytes < rightptr)
  {
    len = (unsigned long) (rightptr-leftptr)/MULT2(merbytes);
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

static const Uchar *remainingrightmost(unsigned long merbytes,
                                       unsigned long byteoffset,
                                       Uchar remainmask,
                                       Uchar code,
                                       const Uchar *leftptr,
                                       const Uchar *rightptr)
{
  unsigned long len;
  Uchar midcode;
  const Uchar *midptr;

  while (leftptr + merbytes < rightptr)
  {
    len = (unsigned long) (rightptr-leftptr)/MULT2(merbytes);
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

bool remainadvance(Merbounds *merbounds,
                   unsigned long merbytes,
                   unsigned long byteoffset,
                   Uchar remainmask,
                   const Uchar *searchbytecode,
                   const Uchar *leftptr,
                   const Uchar *rightptr)
{
  Uchar scode, scodeleft, scoderight;

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
