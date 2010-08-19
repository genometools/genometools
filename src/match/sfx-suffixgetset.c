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

#include <errno.h>
#include <string.h>
#include "core/ma_api.h"
#include "core/fa.h"
#include "core/assert_api.h"
#include "sfx-suffixgetset.h"

#define SUFFIXPTRGET(TAB,IDX)     TAB[IDX]
#define SUFFIXPTRSET(TAB,IDX,VAL) TAB[IDX] = VAL

struct GtSuffixsortspace
{
  unsigned long *ulongtab;
  unsigned long offset,
                bucketleftidx;
  bool unmapsortspace;
};

GtSuffixsortspace *gt_suffixsortspace_new(unsigned long numofentries)
{
  GtSuffixsortspace *suffixsortspace;

  gt_assert(numofentries > 0);
  suffixsortspace = gt_malloc(sizeof(*suffixsortspace));
  suffixsortspace->ulongtab
    = gt_malloc(sizeof(*suffixsortspace->ulongtab) * numofentries);
  suffixsortspace->offset = 0;
  suffixsortspace->bucketleftidx = 0;
  suffixsortspace->unmapsortspace = false;
  return suffixsortspace;
}

GtSuffixsortspace *gt_suffixsortspace_new_fromfile(int filedesc,
                                                   const char *filename,
                                                   unsigned long numofentries)
{
  GtSuffixsortspace *suffixsortspace;

  suffixsortspace = gt_malloc(sizeof(*suffixsortspace));
  suffixsortspace->ulongtab
    = gt_fa_mmap_generic_fd(filedesc,filename,
                            (size_t) numofentries *
                            sizeof (*suffixsortspace->ulongtab),
                            (size_t) 0,false,false,NULL);
  suffixsortspace->offset = 0;
  suffixsortspace->bucketleftidx = 0;
  suffixsortspace->unmapsortspace = true;
  return suffixsortspace;
}

void gt_suffixsortspace_delete(GtSuffixsortspace *suffixsortspace)
{
  if (suffixsortspace != NULL)
  {
    if (suffixsortspace->unmapsortspace)
    {
      gt_fa_xmunmap(suffixsortspace->ulongtab);
    } else
    {
      gt_free(suffixsortspace->ulongtab);
    }
    gt_free(suffixsortspace);
  }
}

/*
static void suffixptrassert(const GtSuffixsortspace *sssp,
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

unsigned long gt_suffixsortspace_getdirect(const GtSuffixsortspace *sssp,
                                           unsigned long idx)
{
  return SUFFIXPTRGET(sssp->ulongtab,idx);
}

void gt_suffixsortspace_setdirect(GtSuffixsortspace *sssp,
                        unsigned long idx,
                        unsigned long value)
{
  SUFFIXPTRSET(sssp->ulongtab,idx,value);
}

unsigned long suffixptrget(const GtSuffixsortspace *sssp,
                           unsigned long subbucketleft,
                           unsigned long idx)
{
  /*suffixptrassert(sssp,subbucket,subbucketleft,idx);*/
  return SUFFIXPTRGET(sssp->ulongtab,sssp->bucketleftidx +
                                     subbucketleft + idx -
                                     sssp->offset);
}

void suffixptrset(GtSuffixsortspace *sssp,
                  unsigned long subbucketleft,
                  unsigned long idx,
                  unsigned long value)
{
  /*suffixptrassert(sssp,subbucket,subbucketleft,idx);*/
  SUFFIXPTRSET(sssp->ulongtab,sssp->bucketleftidx + subbucketleft + idx -
                              sssp->offset,value);
}

void gt_suffixsortspace_setdirectwithoffset(const GtSuffixsortspace *sssp,
                                            unsigned long idx,
                                            unsigned long value)
{
  SUFFIXPTRSET(sssp->ulongtab,idx - sssp->offset,value);
}

unsigned long gt_suffixsortspace_bucketleftidx_get(const GtSuffixsortspace
                                                   *sssp)
{
  return sssp->bucketleftidx;
}

void gt_suffixsortspace_bucketleftidx_set(GtSuffixsortspace *sssp,
                                          unsigned long value)
{
  sssp->bucketleftidx = value;
}

unsigned long gt_suffixsortspace_offset_get(const GtSuffixsortspace *sssp)
{
  return sssp->offset;
}

void gt_suffixsortspace_offset_set(GtSuffixsortspace *sssp,
                                   unsigned long offset)
{
  sssp->offset = offset;
}

void gt_suffixsortspace_sortspace_delete(GtSuffixsortspace *sssp)
{
  gt_free(sssp->ulongtab);
  sssp->ulongtab = NULL;
}

unsigned long *gt_suffixsortspace_ulong_get(const GtSuffixsortspace *sssp)
{
  return (unsigned long *) sssp->ulongtab; /* XXX constrain the type cast */
}

int gt_suffixsortspace_to_file (FILE *outfpsuftab,
                                const GtSuffixsortspace *suffixsortspace,
                                unsigned long numberofsuffixes,
                                GtError *err)
{
  bool haserr = false;

/*
  unsigned long idx;

  for (idx = 0; !haserr && idx < numberofsuffixes; idx++)
  {
    unsigned long value = gt_suffixsortspace_getdirect(suffixsortspace,idx);
    if (fwrite(&value,
               sizeof (value),
               (size_t) 1,
               outfpsuftab)
               != (size_t) 1)
    {
      gt_error_set(err,"cannot write one  item of size %u: errormsg=\"%s\"",
                   (unsigned int) sizeof (value),
                   strerror(errno));
      haserr = true;
    }
  }
*/
  if (fwrite(suffixsortspace->ulongtab,
             sizeof (*suffixsortspace->ulongtab),
             (size_t) numberofsuffixes,
             outfpsuftab)
             != (size_t) numberofsuffixes)
  {
    gt_error_set(err,"cannot write %lu items of size %u: errormsg=\"%s\"",
                 numberofsuffixes,
                 (unsigned int) sizeof (*suffixsortspace->ulongtab),
                 strerror(errno));
    haserr = true;
  }
  return haserr ? -1 : 0;
}
