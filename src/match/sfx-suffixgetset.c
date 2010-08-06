/*
  Copyright (c) 2010 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2010 Center for Bioinformatics, University of Hamburg

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

#include "core/ma_api.h"
#include "core/fa.h"
#include "core/assert_api.h"
#include "core/unused_api.h"
#include "suffixptr.h"
#include "sfx-suffixgetset.h"

struct Suffixsortspace
{
  Suffixptr *sortspace;
  unsigned long sortspaceoffset,
                bucketleftidx;
  bool freesortspace,
       sortspacestart,
       unmapsortspace;
};

Suffixsortspace *suffixsortspace_new(unsigned long numofentries)
{
  Suffixsortspace *suffixsortspace;

  suffixsortspace = gt_malloc(sizeof(*suffixsortspace));
  if (numofentries == 0)
  {
    suffixsortspace->sortspace = NULL;
    suffixsortspace->freesortspace = false;
    suffixsortspace->sortspacestart = false;
  } else
  {
    suffixsortspace->sortspace = gt_malloc(sizeof(*suffixsortspace->sortspace) *
                                           numofentries);
    /*
    fprintf(stderr,"sortspace of size %lu at %lu\n",
            numofentries,(unsigned long) suffixsortspace->sortspace);
    */
    suffixsortspace->sortspacestart = true;
    suffixsortspace->freesortspace = true;
  }
  suffixsortspace->sortspaceoffset = 0;
  suffixsortspace->bucketleftidx = 0;
  suffixsortspace->unmapsortspace = false;
  return suffixsortspace;
}

Suffixsortspace *suffixsortspace_new_fromfile(int filedesc,
                                              const char *filename,
                                              unsigned long numofentries)
{
  Suffixsortspace *suffixsortspace;
  suffixsortspace = gt_malloc(sizeof(*suffixsortspace));
  suffixsortspace->sortspace
    = gt_fa_mmap_generic_fd(filedesc,filename,
                            (size_t) numofentries * sizeof (Suffixptr),
                            (size_t) 0,false,false,NULL);
  suffixsortspace->sortspaceoffset = 0;
  suffixsortspace->freesortspace = false;
  suffixsortspace->sortspacestart = true;
  suffixsortspace->bucketleftidx = 0;
  suffixsortspace->unmapsortspace = true;
  return suffixsortspace;
}

void suffixsortspace_delete(Suffixsortspace *suffixsortspace)
{
  if (suffixsortspace != NULL)
  {
    if (suffixsortspace->freesortspace)
    {
      gt_free(suffixsortspace->sortspace);
    }
    if (suffixsortspace->unmapsortspace)
    {
      gt_fa_xmunmap(suffixsortspace->sortspace);
    }
    gt_free(suffixsortspace);
  }
}

static void suffixptrassert(const Suffixsortspace *sssp,
                            const Suffixptr *subbucket,
                            unsigned long subbucketleft,
                            unsigned long idx)
{
  gt_assert(sssp != NULL);
  gt_assert(sssp->sortspace != NULL);
  gt_assert(sssp->sortspaceoffset <= sssp->bucketleftidx + subbucketleft + idx);
  gt_assert(subbucket != NULL);
  if (subbucket + idx != sssp->sortspace +
                         sssp->bucketleftidx + subbucketleft + idx)
  {
    fprintf(stderr,"idx=%lu,subbucket=%lu,sssp->sortspace=%lu,"
           "bucketleftidx=%lu,subbucketleft=%lu,sortspaceoffset=%lu\n",idx,
           (unsigned long) subbucket,
           (unsigned long) sssp->sortspace,
           sssp->bucketleftidx,
           subbucketleft,
           sssp->sortspaceoffset);
    exit(EXIT_FAILURE);
  }
  gt_assert(subbucket + idx == sssp->sortspace +
                               sssp->bucketleftidx + subbucketleft + idx);
}

unsigned long suffixptrget(const Suffixsortspace *sssp,
                           const Suffixptr *subbucket,
                           unsigned long subbucketleft,
                           unsigned long idx)
{
  suffixptrassert(sssp,subbucket,subbucketleft,idx);
  return SUFFIXPTRGET(sssp->sortspace,sssp->bucketleftidx +
                                      subbucketleft + idx -
                                      sssp->sortspaceoffset);
}

void suffixptrset(Suffixsortspace *sssp,
                  Suffixptr *subbucket,
                  unsigned long subbucketleft,
                  unsigned long idx,
                  unsigned long value)
{
  suffixptrassert(sssp,subbucket,subbucketleft,idx);
  SUFFIXPTRSET(sssp->sortspace,sssp->bucketleftidx + subbucketleft + idx -
                               sssp->sortspaceoffset,value);
}

void suffixptrset2(const Suffixsortspace *sssp,
                   unsigned long idx,
                   unsigned long value)
{
  SUFFIXPTRSET(sssp->sortspace,idx - sssp->sortspaceoffset,value);
}

unsigned long suffixptrget3(const Suffixsortspace *sssp,unsigned long idx)
{
  return SUFFIXPTRGET(sssp->sortspace,idx);
}

void suffixptrset3(Suffixsortspace *sssp,
                   unsigned long idx,
                   unsigned long value)
{
  SUFFIXPTRSET(sssp->sortspace,idx,value);
}

unsigned long gt_suffixsortspace_bucketleftidx_get(const Suffixsortspace *sssp)
{
  return sssp->bucketleftidx;
}

unsigned long gt_suffixsortspace_offset_get(const Suffixsortspace *sssp)
{
  return sssp->sortspaceoffset;
}

void gt_suffixsortspace_bucketleftidx_set(Suffixsortspace *sssp,
                                          unsigned long value)
{
  sssp->bucketleftidx = value;
}

void gt_suffixsortspace_offset_set(Suffixsortspace *sssp,
                                   unsigned long offset)
{
  sssp->sortspaceoffset = offset;
}

void gt_suffixsortspace_sortspace_delete(Suffixsortspace *sssp)
{
  gt_free(sssp->sortspace);
  sssp->sortspace = NULL;
}

Suffixptr *gt_suffixsortspace_sortspace_get(const Suffixsortspace *sssp)
{
  return sssp->sortspace;
}
