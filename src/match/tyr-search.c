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
#include "core/unused_api.h"
#include "divmodmul.h"
#include "opensfxfile.h"
#include "tyr-search.h"
#include "spacedef.h"

struct Tallymercountinfo
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

Tallymercountinfo *tallymercountinfo_new(size_t numofmers,
                                         const GtStr *tallymerindexname,
                                         GtError *err)
{
  size_t numofbytes;
  void *tmp;
  bool haserr = false;
  Tallymercountinfo *tallymercountinfo;

  ALLOCASSIGNSPACE(tallymercountinfo,NULL,Tallymercountinfo,1);
  tallymercountinfo->indexfilename = tallymerindexname;
  tallymercountinfo->mappedmctfileptr
    = genericmaponlytable(tallymerindexname,COUNTSSUFFIX,&numofbytes,err);
  if (tallymercountinfo->mappedmctfileptr == NULL)
  {
    tallymercountinfo->smallcounts = NULL;
    haserr = true;
  } else
  {
    tallymercountinfo->smallcounts
      = (Uchar *) tallymercountinfo->mappedmctfileptr;
    tmp = &tallymercountinfo->smallcounts[numofmers];
    tallymercountinfo->largecounts = (Largecount *) tmp;
    if (numofbytes < numofmers)
    {
      gt_error_set(err,"size of file \"%s.%s\" is smaller than minimum size "
                       "%lu",gt_str_get(tallymerindexname),COUNTSSUFFIX,
                       (unsigned long) numofmers);
      haserr = true;
    }
  }
  if (!haserr && (numofbytes - numofmers) % sizeof (Largecount) != 0)
  {
    gt_error_set(err,"(numofbytes - numofmers) = %lu must be a multiple of %lu",
           (unsigned long) (numofbytes - numofmers),
           (unsigned long) sizeof (Largecount));
    haserr = true;
  }
  if (!haserr)
  {
    tallymercountinfo->numoflargecounts
      = (unsigned long) (numofbytes - numofmers)/
        (unsigned long) sizeof (Largecount);
  }
  if (haserr)
  {
    if (tallymercountinfo->mappedmctfileptr != NULL)
    {
      gt_fa_xmunmap(tallymercountinfo->mappedmctfileptr);
      tallymercountinfo->mappedmctfileptr = NULL;
    }
    FREESPACE(tallymercountinfo);
  }
  return haserr ? NULL : tallymercountinfo;
}

void tallymercountinfo_delete(Tallymercountinfo **tallymercountinfoptr)
{
  Tallymercountinfo *tallymercountinfo = *tallymercountinfoptr;

  gt_fa_xmunmap(tallymercountinfo->mappedmctfileptr);
  tallymercountinfo->mappedmctfileptr = NULL;
  FREESPACE(tallymercountinfo);
  *tallymercountinfoptr = NULL;
}

int tallymersearch(const GtStr *tallymerindexname,
                   GT_UNUSED const GtStrArray *queryfilenames,
                   GT_UNUSED unsigned int showmode,
                   GT_UNUSED unsigned int strand,
                   bool verbose,
                   GtError *err)
{
  Tallymerindex *tallymerindex;
  bool haserr = false;

  tallymerindex = tallymerindex_new(tallymerindexname,err);
  if (tallymerindex == NULL)
  {
    haserr = true;
  }
  if (verbose)
  {
    tallymerindex_show(tallymerindex);
  }
  if (tallymerindex != NULL)
  {
    tallymerindex_delete(&tallymerindex);
  }
  return haserr ? -1 : 0;
}
