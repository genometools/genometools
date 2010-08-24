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
#include <limits.h>
#include "core/assert_api.h"
#include "core/bitpackarray.h"
#include "core/fa.h"
#include "core/log_api.h"
#include "core/ma_api.h"
#include "core/mathsupport.h"
#include "sfx-suffixgetset.h"

struct GtSuffixsortspace
{
  unsigned long(*getdirect)(const GtSuffixsortspace *,unsigned long);
  void (*setdirect)(const GtSuffixsortspace *,unsigned long,unsigned long);
  BitPackArray *bitpackarray;
  bool unmapsortspace;
  unsigned long maxindex,
                maxvalue,
                offset,
                bucketleftidx,
                *ulongtab;
};

static unsigned long getdirect_bitpackarray(const GtSuffixsortspace *sssp,
                                            unsigned long idx)
{
  gt_assert(sssp->bitpackarray != NULL);
  return (unsigned long)
#ifdef _LP64
                         bitpackarray_get_uint64
#else
                         bitpackarray_get_uint32
#endif
                                                 (sssp->bitpackarray,
                                                  (BitOffset) idx);
}

static void setdirect_bitpackarray(const GtSuffixsortspace *sssp,
                                   unsigned long idx,
                                   unsigned long value)
{
  gt_assert(sssp->bitpackarray != NULL);
#ifdef _LP64
  bitpackarray_store_uint64
#else
  bitpackarray_store_uint32
#endif
                            (sssp->bitpackarray,(BitOffset) idx,
#ifdef _LP64
                             (uint64_t)
#else
                             (uint32_t)
#endif
                            value);
}

static unsigned long getdirect_ulong(const GtSuffixsortspace *sssp,
                                     unsigned long idx)
{
  gt_assert(sssp->ulongtab != NULL);
  return (unsigned long) sssp->ulongtab[idx];
}

static void setdirect_ulong(const GtSuffixsortspace *sssp,
                            unsigned long idx,
                            unsigned long value)
{
  gt_assert(sssp->ulongtab != NULL);
  sssp->ulongtab[idx] = value;
}

GtSuffixsortspace *gt_suffixsortspace_new(unsigned long numofentries,
                                          unsigned long maxvalue,
                                          bool suftabasulongarray)
{
  GtSuffixsortspace *suffixsortspace;

  gt_assert(numofentries > 0);
  suffixsortspace = gt_malloc(sizeof(*suffixsortspace));
  suffixsortspace->maxindex = numofentries-1;
  suffixsortspace->maxvalue = maxvalue;
#define GT_SUFTABASULONGARRAY
#ifdef GT_SUFTABASULONGARRAY
  suftabasulongarray = true;
#endif
  if (suftabasulongarray)
  {
    gt_log_log("suftab as array: maxvalue=%lu,numofentries=%lu\n",
               maxvalue,numofentries);
    suffixsortspace->bitpackarray = NULL;
    suffixsortspace->ulongtab
      = gt_malloc(sizeof(*suffixsortspace->ulongtab) * numofentries);
    suffixsortspace->getdirect = getdirect_ulong;
    suffixsortspace->setdirect = setdirect_ulong;
  } else
  {
    unsigned int bitspervalue
      = gt_determinebitspervalue((uint64_t) maxvalue);

    gt_log_log("suftab as bitpackarray: maxvalue=%lu,numofentries=%lu,"
               "bitspervalue=%u\n",maxvalue,numofentries,bitspervalue);

    suffixsortspace->ulongtab = NULL;
    suffixsortspace->bitpackarray
      = bitpackarray_new(bitspervalue,(BitOffset) numofentries,true);
    suffixsortspace->getdirect = getdirect_bitpackarray;
    suffixsortspace->setdirect = setdirect_bitpackarray;
  }
  suffixsortspace->offset = 0;
  suffixsortspace->bucketleftidx = 0;
  suffixsortspace->unmapsortspace = false;
  return suffixsortspace;
}

GtSuffixsortspace *gt_suffixsortspace_new_fromfile(int filedesc,
                                                   const char *filename,
                                                   unsigned long numofentries,
                                                   unsigned long maxvalue)
{
  GtSuffixsortspace *suffixsortspace;

  suffixsortspace = gt_malloc(sizeof(*suffixsortspace));
  suffixsortspace->bitpackarray = NULL;
  suffixsortspace->ulongtab
    = gt_fa_mmap_generic_fd(filedesc,filename,
                            (size_t) numofentries *
                            sizeof (*suffixsortspace->ulongtab),
                            (size_t) 0,false,false,NULL);
  suffixsortspace->offset = 0;
  suffixsortspace->bucketleftidx = 0;
  suffixsortspace->unmapsortspace = true;
  suffixsortspace->maxindex = numofentries - 1;
  suffixsortspace->maxvalue = maxvalue;
  suffixsortspace->getdirect = getdirect_ulong;
  suffixsortspace->setdirect = setdirect_ulong;
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
      bitpackarray_delete(suffixsortspace->bitpackarray);
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
  gt_assert(idx <= sssp->maxindex);
  /*
  printf("idx=%lu\n",idx);
  */
#ifdef GT_SUFTABASULONGARRAY
  return sssp->ulongtab[idx];
#else
  return sssp->getdirect(sssp,idx);
#endif
}

void gt_suffixsortspace_setdirect(GtSuffixsortspace *sssp,
                                  unsigned long idx,
                                  unsigned long value)
{
  gt_assert(idx <= sssp->maxindex);
  gt_assert(value <= sssp->maxvalue);
  /*
  printf("idx=%lu,value=%lu\n",idx,value);
  */
#ifdef GT_SUFTABASULONGARRAY
  sssp->ulongtab[idx] = value;
#else
  sssp->setdirect(sssp,idx,value);
#endif
}

unsigned long gt_suffixsortspace_get(const GtSuffixsortspace *sssp,
                                     unsigned long subbucketleft,
                                     unsigned long idx)
{
  /*suffixptrassert(sssp,subbucket,subbucketleft,idx);*/
  return gt_suffixsortspace_getdirect(sssp, sssp->bucketleftidx
                                              + subbucketleft
                                              + idx
                                              - sssp->offset);
}

void gt_suffixsortspace_set(GtSuffixsortspace *sssp,
                            unsigned long subbucketleft,
                            unsigned long idx,
                            unsigned long value)
{
  /*suffixptrassert(sssp,subbucket,subbucketleft,idx);*/
  gt_suffixsortspace_setdirect(sssp, sssp->bucketleftidx
                                       + subbucketleft
                                       + idx
                                       - sssp->offset,value);
}

void gt_suffixsortspace_setdirectwithoffset(GtSuffixsortspace *sssp,
                                            unsigned long idx,
                                            unsigned long value)
{
  gt_suffixsortspace_setdirect(sssp,idx - sssp->offset,value);
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
  bitpackarray_delete(sssp->bitpackarray);
  sssp->bitpackarray = NULL;
  gt_free(sssp->ulongtab);
  sssp->ulongtab = NULL;
}

unsigned long *gt_suffixsortspace_ulong_get(const GtSuffixsortspace *sssp)
{
  gt_assert(sssp->ulongtab != NULL);
  return (unsigned long *) sssp->ulongtab; /* XXX constrain the type cast */
}

int gt_suffixsortspace_to_file (FILE *outfpsuftab,
                                const GtSuffixsortspace *suffixsortspace,
                                unsigned long numberofsuffixes,
                                GtError *err)
{
  bool haserr = false;

  if (suffixsortspace->ulongtab != NULL)
  {
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
  } else
  {
    unsigned long idx;

    gt_assert(suffixsortspace->ulongtab == NULL);
    for (idx = 0; !haserr && idx < numberofsuffixes; idx++)
    {
      unsigned long value = getdirect_bitpackarray(suffixsortspace,idx);
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
  }
  return haserr ? -1 : 0;
}
