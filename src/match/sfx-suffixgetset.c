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
#include "core/defined-types.h"
#include "core/fa.h"
#include "core/log_api.h"
#include "core/ma_api.h"
#include "core/mathsupport.h"
#include "core/safearith.h"
#include "core/unused_api.h"
#include "sfx-suffixgetset.h"

struct GtSuffixsortspace
{
  unsigned long(*getdirect)(const GtSuffixsortspace *,unsigned long);
  void (*setdirect)(const GtSuffixsortspace *,unsigned long,unsigned long);
  BitPackArray *bitpackarray;
  bool unmapsortspace;
  Definedunsignedlong longestidx;
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

size_t gt_suffixsortspace_requiredspace(unsigned long numofentries,
                                        unsigned long maxvalue,
                                        bool suftabcompressedbytes)
{
  size_t requiredspace = sizeof (GtSuffixsortspace);

  if (suftabcompressedbytes)
  {
    unsigned int bitspervalue
      = gt_determinebitspervalue((uint64_t) maxvalue);
    requiredspace += sizeofbitarray(bitspervalue,(BitOffset) numofentries);
  } else
  {
    requiredspace += numofentries * sizeof (unsigned long);
  }
  return requiredspace;
}

static void gt_suffixsortspace_overflow_abort(GT_UNUSED const char *f,
                                              GT_UNUSED int l,
                                              void *data)
{
  fprintf(stderr, "error: overflow detected while calculating size of "
                  "suffix sorting space: %lu * %lu bytes is too large for "
                  "the current platform, please recompile GenomeTools with "
                  "support for a larger address space to prevent this (e.g. "
                  "64 bit instead of 32 bit) or use the `-parts' option.\n",
                  (unsigned long) sizeof (unsigned long),
                  *(unsigned long*) data);
  exit(GT_EXIT_PROGRAMMING_ERROR);
}

GtSuffixsortspace *gt_suffixsortspace_new(unsigned long numofentries,
                                          unsigned long maxvalue,
                                          bool suftabcompressedbytes)
{
  GtSuffixsortspace *suffixsortspace;

  gt_assert(numofentries > 0);
  suffixsortspace = gt_malloc(sizeof (*suffixsortspace));
  suffixsortspace->maxindex = numofentries-1;
  suffixsortspace->maxvalue = maxvalue;
  suffixsortspace->longestidx.defined = false;
  suffixsortspace->longestidx.valueunsignedlong = 0;
#define GT_SUFTABASULONGARRAY
#ifdef GT_SUFTABASULONGARRAY
  suftabcompressedbytes = false;
#endif
  if (suftabcompressedbytes)
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
  } else
  {
    size_t sufspacesize;
    gt_log_log("suftab as array: maxvalue=%lu,numofentries=%lu",
               maxvalue,numofentries);
    suffixsortspace->bitpackarray = NULL;
    sufspacesize = (size_t) gt_safe_mult_ulong_check((unsigned long)
                                            sizeof (*suffixsortspace->ulongtab),
                                            numofentries,
                                            gt_suffixsortspace_overflow_abort,
                                            &numofentries);
    suffixsortspace->ulongtab = gt_malloc(sufspacesize);
    gt_log_log("sizeof (ulongtab)=%lu bytes",
               sizeof (*suffixsortspace->ulongtab) * numofentries);
    suffixsortspace->getdirect = getdirect_ulong;
    suffixsortspace->setdirect = setdirect_ulong;
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

  suffixsortspace = gt_malloc(sizeof (*suffixsortspace));
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
  suffixsortspace->longestidx.defined = false;
  suffixsortspace->longestidx.valueunsignedlong = 0;
  suffixsortspace->getdirect = getdirect_ulong;
  suffixsortspace->setdirect = setdirect_ulong;
  return suffixsortspace;
}

void gt_suffixsortspace_delete(GtSuffixsortspace *suffixsortspace,
                               bool checklongestdefined)
{
  if (suffixsortspace != NULL)
  {
    if (checklongestdefined)
    {
      gt_assert(suffixsortspace->longestidx.defined);
    }
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
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }
  gt_assert(subbucket + idx == sssp->sortspace +
                               sssp->bucketleftidx + subbucketleft + idx);
}
*/

unsigned long gt_suffixsortspace_getdirect(const GtSuffixsortspace *sssp,
                                           unsigned long idx)
{
  if (idx > sssp->maxindex)
  {
    fprintf(stderr,"line %d: idx=%lu > %lu = maxindex\n",__LINE__,
                                                    idx,sssp->maxindex);
    exit(EXIT_FAILURE); /* XXX remove later */
  }
  gt_assert(idx <= sssp->maxindex);
#ifdef GT_SUFTABASULONGARRAY
  /*
  printf("getdirect(%lu)=%lu\n",idx,sssp->ulongtab[idx]);
  */
  return sssp->ulongtab[idx];
#else
  return sssp->getdirect(sssp,idx);
#endif
}

void gt_suffixsortspace_setdirect(GtSuffixsortspace *sssp,
                                  unsigned long idx,
                                  unsigned long value)
{
  if (idx > sssp->maxindex)
  {
    fprintf(stderr,"line %d: idx=%lu > %lu = maxindex\n",__LINE__,
                                                    idx,sssp->maxindex);
    exit(EXIT_FAILURE); /* XXX remove later */
  }
  gt_assert(idx <= sssp->maxindex);
  gt_assert(value <= sssp->maxvalue);
  /*printf("idx=%lu,value=%lu\n",idx,value);*/
  if (value == 0)
  {
    sssp->longestidx.defined = true;
    sssp->longestidx.valueunsignedlong = idx + sssp->offset;
  }
#ifdef GT_SUFTABASULONGARRAY
  sssp->ulongtab[idx] = value;
#else
  sssp->setdirect(sssp,idx,value);
#endif
}

void gt_suffixsortspace_showrange(const GtSuffixsortspace *sssp,
                                  unsigned long subbucketleft,
                                  unsigned long width)
{
  unsigned long idx;

  printf("%lu,%lu=",sssp->bucketleftidx+subbucketleft-sssp->offset,
                    sssp->bucketleftidx+subbucketleft+width-1-sssp->offset);
  for (idx=sssp->bucketleftidx+subbucketleft-sssp->offset;
       idx<sssp->bucketleftidx+subbucketleft+width-sssp->offset;
       idx++)
  {
    printf(" %lu",sssp->ulongtab[idx]);
  }
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

unsigned long gt_suffixsortspace_longest(const GtSuffixsortspace *sssp)
{
  gt_assert(sssp->longestidx.defined);
  return sssp->longestidx.valueunsignedlong;
}

int gt_suffixsortspace_to_file (FILE *outfpsuftab,
                                const GtSuffixsortspace *sssp,
                                unsigned long numberofsuffixes,
                                GtError *err)
{
  bool haserr = false;

  if (sssp->ulongtab != NULL)
  {
    if (fwrite(sssp->ulongtab,
               sizeof (*sssp->ulongtab),
               (size_t) numberofsuffixes,
               outfpsuftab)
               != (size_t) numberofsuffixes)
    {
      gt_error_set(err,"cannot write %lu items of size %u: errormsg=\"%s\"",
                   numberofsuffixes,
                   (unsigned int) sizeof (*sssp->ulongtab),
                   strerror(errno));
      haserr = true;
    }
  } else
  {
    unsigned long idx;

    gt_assert(sssp->ulongtab == NULL);
    for (idx = 0; !haserr && idx < numberofsuffixes; idx++)
    {
      unsigned long value = getdirect_bitpackarray(sssp,idx);
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
