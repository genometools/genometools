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
#include "divmodmul.h"
#include "opensfxfile.h"
#include "tyr-search.h"
#include "spacedef.h"

struct MCTinfo
{
  void *mappedmctfileptr;
  const GtStr *indexfilename;
  Uchar *smallcounts;
  Largecount *largecounts;
  unsigned long numoflargecounts;
};

struct Tallymerindex
{
  void *mappedfileptr;
  const GtStr *indexfilename;
  unsigned int alphasize;
  unsigned long numofmers,
                mersize,
                merbytes;
  Uchar *mertable,
        *lastmer;
};

typedef struct
{
  Uchar *leftmer,
        *rightmer;
} Merbounds;

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

Tallymerindex *tallymerindex_new(const GtStr *tallymerindexname,GtError *err)
{
  bool haserr = false;
  size_t numofbytes, rest;
  Tallymerindex *tallymerindex;

  ALLOCASSIGNSPACE(tallymerindex,NULL,Tallymerindex,1);
  tallymerindex->indexfilename = tallymerindexname;
  tallymerindex->mappedfileptr = genericmaponlytable(tallymerindexname,
                                                 MERSUFFIX,&numofbytes,err);
  if (tallymerindex->mappedfileptr == NULL)
  {
    haserr = true;
  }
  if (!haserr)
  {
    tallymerindex->mertable = (Uchar *) tallymerindex->mappedfileptr;
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
    assert(tallymerindex->mertable != NULL);
    tallymerindex->mersize
      = decodesingleinteger(tallymerindex->mertable + numofbytes - rest);
    tallymerindex->alphasize
      = (unsigned int) decodesingleinteger(tallymerindex->mertable +
                                           numofbytes - rest +
                                           sizeof (unsigned long));
    tallymerindex->merbytes = MERBYTES(tallymerindex->mersize);
    if ((numofbytes - rest) % tallymerindex->merbytes != 0)
    {
      gt_error_set(err,"size of index is %lu which is not a multiple of %lu",
                   (unsigned long) (numofbytes - rest),
                   tallymerindex->merbytes);
      haserr = true;
    }
  }
  if (!haserr)
  {
    tallymerindex->numofmers
      = (unsigned long) (numofbytes - rest) / tallymerindex->merbytes;
    assert(tallymerindex->mertable != NULL);
    if (tallymerindex->numofmers == 0)
    {
      tallymerindex->lastmer = tallymerindex->mertable - 1;
    } else
    {
      tallymerindex->lastmer
        = tallymerindex->mertable + numofbytes - rest - tallymerindex->merbytes;
    }
  }
  if (haserr)
  {
    if (tallymerindex->mappedfileptr != NULL)
    {
      gt_fa_xmunmap(tallymerindex->mappedfileptr);
      tallymerindex->mappedfileptr = NULL;
    }
    FREESPACE(tallymerindex);
  }
  return haserr ? NULL : tallymerindex;
}

int mapMCTinfo(MCTinfo *mctinfo,size_t numofmers,
               const GtStr *tallymerindexname,GtError *err)
{
  size_t numofbytes;
  void *tmp;
  bool haserr = false;

  mctinfo->indexfilename = tallymerindexname;
  mctinfo->mappedmctfileptr = genericmaponlytable(tallymerindexname,
                                                  COUNTSSUFFIX,&numofbytes,err);
  if (mctinfo->mappedmctfileptr == NULL)
  {
    mctinfo->smallcounts = NULL;
    haserr = true;
  } else
  {
    mctinfo->smallcounts = (Uchar *) mctinfo->mappedmctfileptr;
    tmp = &mctinfo->smallcounts[numofmers];
    mctinfo->largecounts = (Largecount *) tmp;
    if (numofbytes < numofmers)
    {
      gt_error_set(err,"size of file \"%s.%s\" is smaller than minimum size "
                       "%lu",gt_str_get(tallymerindexname),COUNTSSUFFIX,
                       (unsigned long) numofmers);
      haserr = true;
    }
  }
  if ((numofbytes - numofmers) % sizeof (Largecount) != 0)
  {
    gt_error_set(err,"(numofbytes - numofmers) = %lu must be a multiple of %lu",
           (unsigned long) (numofbytes - numofmers),
           (unsigned long) sizeof (Largecount));
    haserr = true;
  }
  mctinfo->numoflargecounts = (unsigned long) (numofbytes - numofmers)/
                              (unsigned long) sizeof (Largecount);
  if (haserr && mctinfo->mappedmctfileptr != NULL)
  {
    gt_fa_xmunmap(mctinfo->mappedmctfileptr);
    mctinfo->mappedmctfileptr = NULL;
  }
  return haserr ? -1 : 0;
}

void tallymerindex_show(const Tallymerindex *tallymerindex)
{
  printf("# indexfilename = %s\n",gt_str_get(tallymerindex->indexfilename));
  printf("# alphasize = %u\n",tallymerindex->alphasize);
  printf("# mersize = %lu\n",tallymerindex->mersize);
  printf("# numofmers = %lu\n",tallymerindex->numofmers);
  printf("# merbytes = %lu\n",tallymerindex->merbytes);
}

void tallymerindex_delete(Tallymerindex **tallymerindexptr)
{
  Tallymerindex *tallymerindex = *tallymerindexptr;

  gt_fa_xmunmap(tallymerindex->mappedfileptr);
  tallymerindex->mappedfileptr = NULL;
  FREESPACE(tallymerindex);
  *tallymerindexptr = NULL;
}
