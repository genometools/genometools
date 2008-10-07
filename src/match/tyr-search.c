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

#include "core/fa.h"
#include "tyr-search.h"
#include "divmodmul.h"
#include "opensfxfile.h"

unsigned long decodesingleinteger(const Uchar *start)
{
  unsigned long idx, value;

  value = (unsigned long) start[0];
  for (idx=1UL; idx < (unsigned long) sizeof (unsigned long); idx++)
  {
    value |= (((unsigned long) start[idx]) << MULT8(idx));
  }
  return value;
}

int mapvmerindex(Vmerindex *vmerindex,const GtStr *vmerindexname,GtError *err)
{
  bool haserr = false;
  size_t numofbytes, rest;

  vmerindex->indexfilename = vmerindexname;
  vmerindex->mappedfileptr = genericmaponlytable(vmerindexname,
                                                 MERSUFFIX,&numofbytes,err);
  if (vmerindex->mappedfileptr == NULL)
  {
    haserr = true;
  }
  if (!haserr)
  {
    vmerindex->mertable = (Uchar *) vmerindex->mappedfileptr;
    rest = sizeof (unsigned long) * EXTRAINTEGERS;
    if (rest > numofbytes)
    {
      gt_error_set(err,"index must contain at least %lu bytes",
                        (unsigned long) rest);
      haserr = true;
    }
  }
  if (!haserr)
  {
    assert(vmerindex->mertable != NULL);
    vmerindex->mersize
      = decodesingleinteger(vmerindex->mertable + numofbytes - rest);
    vmerindex->alphasize
      = (unsigned int) decodesingleinteger(vmerindex->mertable +
                                           numofbytes - rest +
                                           sizeof (unsigned long));
    vmerindex->merbytes = MERBYTES(vmerindex->mersize);
    if ((numofbytes - rest) % vmerindex->merbytes != 0)
    {
      gt_error_set(err,"size of index is %lu which is not a multiple of %lu",
                   numofbytes - rest,
                   vmerindex->merbytes);
      haserr = true;
    }
  }
  if (!haserr)
  {
    vmerindex->numofmers
      = (unsigned long) (numofbytes - rest) / vmerindex->merbytes;
    assert(vmerindex->mertable != NULL);
    if (vmerindex->numofmers == 0)
    {
      vmerindex->lastmer = vmerindex->mertable - 1;
    } else
    {
      vmerindex->lastmer
        = vmerindex->mertable + numofbytes - rest - vmerindex->merbytes;
    }
  }
  if (haserr && vmerindex->mappedfileptr != NULL)
  {
    gt_fa_xmunmap(vmerindex->mappedfileptr);
    vmerindex->mappedfileptr = NULL;
  }
  return haserr ? -1 : 0;
}
