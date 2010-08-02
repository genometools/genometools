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

#ifndef SFX_SUFFIXGETSET_H
#define SFX_SUFFIXGETSET_H

#include "core/assert_api.h"
#include "core/unused_api.h"
#include "suffixptr.h"

typedef struct
{
  Suffixptr *sortspace;
  unsigned long sortspaceoffset,
                bucketleftidx;
  bool freesortspace;
} Suffixsortspace;

typedef void (*Dc_processunsortedrange)(void *,
                                        Suffixptr *,
                                        unsigned long,
                                        unsigned long,
                                        unsigned long);

/*@unused@*/ static inline Suffixsortspace
              *suffixsortspace_new(unsigned long numofentries)
{
  Suffixsortspace *suffixsortspace;

  suffixsortspace = gt_malloc(sizeof(*suffixsortspace));
  if (numofentries == 0)
  {
    suffixsortspace->sortspace = NULL;
    suffixsortspace->freesortspace = false;
  } else
  {
    suffixsortspace->sortspace = gt_malloc(sizeof(*suffixsortspace->sortspace) *
                                           numofentries);
    suffixsortspace->freesortspace = true;
  }
  suffixsortspace->sortspaceoffset = 0;
  suffixsortspace->bucketleftidx = 0;
  return suffixsortspace;
}

/*@unused@*/ static inline void
             suffixsortspace_delete(Suffixsortspace *suffixsortspace)
{
  if (suffixsortspace != NULL)
  {
    if (suffixsortspace->freesortspace)
    {
      gt_free(suffixsortspace->sortspace);
    }
    gt_free(suffixsortspace);
  }
}

/*@unused@*/ static inline void suffixptrassert(const Suffixsortspace *sssp,
                                                unsigned long subbucketleft,
                                                unsigned long idx)
{
  gt_assert(sssp != NULL);
  gt_assert(sssp->sortspaceoffset <= sssp->bucketleftidx + subbucketleft + idx);
  /*
  fprintf(stderr,"idx=%lu,bucketleftidx=%lu,subbucketleft=%lu,"
                 "sortspaceoffset=%lu\n",
          idx,sssp->bucketleftidx,subbucketleft,sssp->sortspaceoffset);
  */
}

/*@unused@*/ static inline unsigned long suffixptrget(
                                         const Suffixsortspace *sssp,
                                         const Suffixptr *subbucket,
                                         unsigned long subbucketleft,
                                         unsigned long idx)
{
  suffixptrassert(sssp,subbucketleft,idx);
  return SUFFIXPTRGET(subbucket,idx);
}

/*@unused@*/ static inline void suffixptrset(Suffixsortspace *sssp,
                                             Suffixptr *subbucket,
                                             unsigned long subbucketleft,
                                             unsigned long idx,
                                             unsigned long value)
{
  suffixptrassert(sssp,subbucketleft,idx);
  SUFFIXPTRSET(subbucket,idx,value);
}

/*@unused@*/ static inline void suffixptrset2(const Suffixsortspace *sssp,
                                              unsigned long idx,
                                              unsigned long value)
{
  SUFFIXPTRSET(sssp->sortspace,idx - sssp->sortspaceoffset,value);
}

/*@unused@*/ static inline unsigned long suffixptrget3(
                                         const Suffixsortspace *sssp,
                                         unsigned long idx)
{
  return SUFFIXPTRGET(sssp->sortspace,idx);
}

/*@unused@*/ static inline void suffixptrset3(const Suffixsortspace *sssp,
                                              unsigned long idx,
                                              unsigned long value)
{
  SUFFIXPTRSET(sssp->sortspace,idx,value);
}

#endif
