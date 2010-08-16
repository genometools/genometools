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
#include "sfx-suffixgetset.h"

#ifdef  SUFFIXPTRNEWVERSION

#define SUFFIXPTRGET(TAB,IDX)     TAB[IDX].value
#define SUFFIXPTRSET(TAB,IDX,VAL) TAB[IDX].value = VAL

typedef struct
{
  unsigned long value;
} Suffixptr;

#else

#define SUFFIXPTRGET(TAB,IDX)     TAB[IDX]
#define SUFFIXPTRSET(TAB,IDX,VAL) TAB[IDX] = VAL

typedef unsigned long Suffixptr;

#endif

struct Suffixsortspace
{
  Suffixptr *sortspace;
  unsigned long offset,
                bucketleftidx;
  bool unmapsortspace;
};

Suffixsortspace *suffixsortspace_new(unsigned long numofentries)
{
  Suffixsortspace *suffixsortspace;

  gt_assert(numofentries > 0);
  suffixsortspace = gt_malloc(sizeof(*suffixsortspace));
  suffixsortspace->sortspace = gt_malloc(sizeof(*suffixsortspace->sortspace) *
                                         numofentries);
  /*
  fprintf(stderr,"sortspace of size %lu at %lu\n",
          numofentries,(unsigned long) suffixsortspace->sortspace);
  */
  suffixsortspace->offset = 0;
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
  suffixsortspace->offset = 0;
  suffixsortspace->bucketleftidx = 0;
  suffixsortspace->unmapsortspace = true;
  return suffixsortspace;
}

void suffixsortspace_delete(Suffixsortspace *suffixsortspace)
{
  if (suffixsortspace != NULL)
  {
    if (suffixsortspace->unmapsortspace)
    {
      gt_fa_xmunmap(suffixsortspace->sortspace);
    } else
    {
      gt_free(suffixsortspace->sortspace);
    }
    gt_free(suffixsortspace);
  }
}

/*
static void suffixptrassert(const Suffixsortspace *sssp,
                            const Suffixptr *subbucket,
                            unsigned long subbucketleft,
                            unsigned long idx)
{
  gt_assert(sssp != NULL);
  gt_assert(sssp->sortspace != NULL);
  gt_assert(sssp->offset <= sssp->bucketleftidx + subbucketleft + idx);
  gt_assert(subbucket != NULL);
  if (subbucket + idx != sssp->sortspace +
                         sssp->bucketleftidx + subbucketleft + idx)
  {
    fprintf(stderr,"idx=%lu,subbucket=%lu,sssp->sortspace=%lu,"
           "bucketleftidx=%lu,subbucketleft=%lu,offset=%lu\n",idx,
           (unsigned long) subbucket,
           (unsigned long) sssp->sortspace,
           sssp->bucketleftidx,
           subbucketleft,
           sssp->offset);
    exit(EXIT_FAILURE);
  }
  gt_assert(subbucket + idx == sssp->sortspace +
                               sssp->bucketleftidx + subbucketleft + idx);
}
*/

unsigned long suffixptrget(const Suffixsortspace *sssp,
                           unsigned long subbucketleft,
                           unsigned long idx)
{
  /*suffixptrassert(sssp,subbucket,subbucketleft,idx);*/
  return SUFFIXPTRGET(sssp->sortspace,sssp->bucketleftidx +
                                      subbucketleft + idx -
                                      sssp->offset);
}

void suffixptrset(Suffixsortspace *sssp,
                  unsigned long subbucketleft,
                  unsigned long idx,
                  unsigned long value)
{
  /*suffixptrassert(sssp,subbucket,subbucketleft,idx);*/
  SUFFIXPTRSET(sssp->sortspace,sssp->bucketleftidx + subbucketleft + idx -
                               sssp->offset,value);
}

void suffixptrset2(const Suffixsortspace *sssp,
                   unsigned long idx,
                   unsigned long value)
{
  SUFFIXPTRSET(sssp->sortspace,idx - sssp->offset,value);
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

void gt_suffixsortspace_bucketleftidx_set(Suffixsortspace *sssp,
                                          unsigned long value)
{
  sssp->bucketleftidx = value;
}

unsigned long gt_suffixsortspace_offset_get(const Suffixsortspace *sssp)
{
  return sssp->offset;
}

void gt_suffixsortspace_offset_set(Suffixsortspace *sssp,
                                   unsigned long offset)
{
  sssp->offset = offset;
}

void gt_suffixsortspace_sortspace_delete(Suffixsortspace *sssp)
{
  gt_free(sssp->sortspace);
  sssp->sortspace = NULL;
}

unsigned long *gt_suffixsortspace_ulong_get(const Suffixsortspace *sssp)
{
  return (unsigned long *) sssp->sortspace; /* XXX remove type cast */
}
