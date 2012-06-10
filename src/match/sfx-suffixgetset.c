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
#include "core/defined-types.h"
#include "core/fa.h"
#include "core/log_api.h"
#include "core/logger_api.h"
#include "core/ma_api.h"
#include "core/mathsupport.h"
#include "core/safearith.h"
#include "core/unused_api.h"
#include "core/xansi_api.h"
#include "sfx-suffixgetset.h"

struct GtSuffixsortspace
{
  bool unmapsortspace, currentexport;
  Definedunsignedlong longestidx;
  uint32_t *uinttab;
  size_t basesize;
  GtSuffixsortspace_exportptr exportptr;
  unsigned long maxindex,
                maxvalue,
                partoffset,
                bucketleftidx,
                *ulongtab;
};

static bool gt_decide_to_use_uint(bool useuint,unsigned long maxvalue)
{
  if (useuint && maxvalue <= (unsigned long) UINT_MAX)
  {
    return true;
  }
  return false;
}

uint64_t gt_suffixsortspace_requiredspace(unsigned long numofentries,
                                          unsigned long maxvalue,
                                          bool useuint)
{
  uint64_t requiredspace = (uint64_t) sizeof (GtSuffixsortspace);

  if (gt_decide_to_use_uint(useuint,maxvalue))
  {
    gt_assert(maxvalue <= (unsigned long) UINT_MAX);
    requiredspace += (uint64_t) numofentries * (uint64_t) sizeof (uint32_t);
  } else
  {
    requiredspace += (uint64_t) numofentries *
                     (uint64_t) sizeof (unsigned long);
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
                                          bool useuint,
                                          GT_UNUSED GtLogger *logger)
{
  GtSuffixsortspace *suffixsortspace;
  unsigned long sufspacesize;

  gt_assert(numofentries > 0);
  suffixsortspace = gt_malloc(sizeof (*suffixsortspace));
  suffixsortspace->maxindex = numofentries-1;
  suffixsortspace->maxvalue = maxvalue;
  suffixsortspace->longestidx.defined = false;
  suffixsortspace->longestidx.valueunsignedlong = 0;
  suffixsortspace->exportptr.ulongtabsectionptr = NULL;
  suffixsortspace->exportptr.uinttabsectionptr = NULL;
  suffixsortspace->currentexport = false;
#ifdef _LP64
  gt_logger_log(logger,"suftab uses %dbit values: "
                         "maxvalue=%lu,numofentries=%lu",
                         gt_decide_to_use_uint(useuint,maxvalue) ? 32 : 64,
                         maxvalue,numofentries);
#endif
  suffixsortspace->basesize = gt_decide_to_use_uint(useuint,maxvalue)
                                ? sizeof (*suffixsortspace->uinttab)
                                : sizeof (*suffixsortspace->ulongtab);
  sufspacesize
    = gt_safe_mult_ulong_check((unsigned long) suffixsortspace->basesize,
                               numofentries,
                               gt_suffixsortspace_overflow_abort,
                               &numofentries);
  gt_log_log("sizeof (suftab)=%lu bytes",sufspacesize);
  if (gt_decide_to_use_uint(useuint,maxvalue))
  {
    suffixsortspace->ulongtab = NULL;
    suffixsortspace->uinttab = gt_malloc((size_t) sufspacesize);
  } else
  {
    suffixsortspace->uinttab = NULL;
    suffixsortspace->ulongtab = gt_malloc((size_t) sufspacesize);
  }
  suffixsortspace->partoffset = 0;
  suffixsortspace->bucketleftidx = 0;
  suffixsortspace->unmapsortspace = false;
  return suffixsortspace;
}

GtSuffixsortspace *gt_suffixsortspace_new_fromfile(int filedesc,
                                                   const char *filename,
                                                   unsigned long numofentries,
                                                   unsigned long maxvalue,
                                                   bool useuint)
{
  GtSuffixsortspace *suffixsortspace;
  void *ptr;

  suffixsortspace = gt_malloc(sizeof (*suffixsortspace));
  suffixsortspace->basesize = gt_decide_to_use_uint(useuint,maxvalue)
                                ? sizeof (*suffixsortspace->uinttab)
                                : sizeof (*suffixsortspace->ulongtab);
  ptr = gt_fa_mmap_generic_fd(filedesc,filename,
                              (size_t) numofentries * suffixsortspace->basesize,
                              (size_t) 0,false,false,NULL);
  if (gt_decide_to_use_uint(useuint,maxvalue))
  {
    suffixsortspace->uinttab = ptr;
    suffixsortspace->ulongtab = NULL;
  } else
  {
    suffixsortspace->ulongtab = ptr;
    suffixsortspace->uinttab = NULL;
  }
  suffixsortspace->partoffset = 0;
  suffixsortspace->bucketleftidx = 0;
  suffixsortspace->unmapsortspace = true;
  suffixsortspace->maxindex = numofentries - 1;
  suffixsortspace->maxvalue = maxvalue;
  suffixsortspace->longestidx.defined = false;
  suffixsortspace->longestidx.valueunsignedlong = 0;
  suffixsortspace->exportptr.ulongtabsectionptr = NULL;
  suffixsortspace->exportptr.uinttabsectionptr = NULL;
  suffixsortspace->currentexport = false;
  return suffixsortspace;
}

void gt_suffixsortspace_delete(GtSuffixsortspace *suffixsortspace,
                               GT_UNUSED bool checklongestdefined)
{
  if (suffixsortspace != NULL)
  {
    gt_assert(!checklongestdefined || suffixsortspace->longestidx.defined);
    if (suffixsortspace->unmapsortspace)
    {
      gt_fa_xmunmap(suffixsortspace->ulongtab);
      gt_fa_xmunmap(suffixsortspace->uinttab);
    } else
    {
      gt_free(suffixsortspace->uinttab);
      gt_free(suffixsortspace->ulongtab);
    }
    gt_free(suffixsortspace);
  }
}

void gt_suffixsortspace_nooffsets(GT_UNUSED const GtSuffixsortspace *sssp)
{
  gt_assert(sssp->partoffset == 0);
  gt_assert(sssp->bucketleftidx == 0);
}

unsigned long gt_suffixsortspace_getdirect(const GtSuffixsortspace *sssp,
                                           unsigned long idx)
{
  gt_assert(idx <= sssp->maxindex);
  if (sssp->ulongtab != NULL)
  {
    return sssp->ulongtab[idx];
  }
  gt_assert(sssp->uinttab != NULL);
  return (unsigned long) sssp->uinttab[idx];
}

void gt_suffixsortspace_updatelongest(GtSuffixsortspace *sssp,unsigned long idx)
{
  sssp->longestidx.defined = true;
  sssp->longestidx.valueunsignedlong = sssp->bucketleftidx + idx;
}

void gt_suffixsortspace_setdirect(GtSuffixsortspace *sssp,
                                  unsigned long idx,
                                  unsigned long value)
{
  gt_assert(idx <= sssp->maxindex && value <= sssp->maxvalue);
  if (value == 0)
  {
    sssp->longestidx.defined = true;
    sssp->longestidx.valueunsignedlong = sssp->partoffset + idx;
  }
  if (sssp->ulongtab != NULL)
  {
    sssp->ulongtab[idx] = value;
  } else
  {
    gt_assert (sssp->uinttab != NULL);
    gt_assert (value <= (unsigned long) UINT_MAX);
    sssp->uinttab[idx] = (uint32_t) value;
  }
}

void gt_suffixsortspace_showrange(const GtSuffixsortspace *sssp,
                                  unsigned long subbucketleft,
                                  unsigned long width)
{
  unsigned long idx;

  printf("%lu,%lu=",sssp->bucketleftidx+subbucketleft-sssp->partoffset,
                    sssp->bucketleftidx+subbucketleft+width-1-sssp->partoffset);
  for (idx=sssp->bucketleftidx+subbucketleft-sssp->partoffset;
       idx<sssp->bucketleftidx+subbucketleft+width-sssp->partoffset;
       idx++)
  {
    printf(" %lu", gt_suffixsortspace_getdirect(sssp,idx));
  }
}

void gt_suffixsortspace_checkorder(const GtSuffixsortspace *sssp,
                                   unsigned long subbucketleft,
                                   unsigned long width)
{
  unsigned long idx, currentpos;
  GT_UNUSED unsigned long prevpos;

  gt_assert(width > 0);
  prevpos = gt_suffixsortspace_getdirect(sssp,
                                         sssp->bucketleftidx+subbucketleft -
                                         sssp->partoffset);
  for (idx=sssp->bucketleftidx+subbucketleft - sssp->partoffset + 1;
       idx<sssp->bucketleftidx+subbucketleft + width - sssp->partoffset;
       idx++)
  {
    currentpos = gt_suffixsortspace_getdirect(sssp,idx);
    gt_assert(prevpos > currentpos);
    prevpos = currentpos;
  }
}

GtSuffixsortspace_exportptr *gt_suffixsortspace_exportptr(
                                  unsigned long subbucketleft,
                                  GtSuffixsortspace *sssp)
{
  if (sssp->ulongtab != NULL)
  {
    sssp->exportptr.ulongtabsectionptr = sssp->ulongtab + sssp->bucketleftidx
                                                        + subbucketleft
                                                        - sssp->partoffset;
    sssp->exportptr.uinttabsectionptr = NULL;
  } else
  {
    gt_assert(sssp->uinttab != NULL);
    sssp->exportptr.uinttabsectionptr = sssp->uinttab + sssp->bucketleftidx
                                                      + subbucketleft
                                                      - sssp->partoffset;
    sssp->exportptr.ulongtabsectionptr = NULL;
  }
  sssp->currentexport = true;
  return &sssp->exportptr;
}

void gt_suffixsortspace_export_done(GtSuffixsortspace *sssp)
{
  sssp->exportptr.ulongtabsectionptr = NULL;
  sssp->exportptr.uinttabsectionptr = NULL;
  sssp->currentexport = false;
}

unsigned long gt_suffixsortspace_get(const GtSuffixsortspace *sssp,
                                     unsigned long subbucketleft,
                                     unsigned long idx)
{
  return gt_suffixsortspace_getdirect(sssp, sssp->bucketleftidx
                                              + subbucketleft
                                              + idx
                                              - sssp->partoffset);
}

const unsigned long *gt_suffixsortspace_getptr_ulong(
                                               const GtSuffixsortspace *sssp,
                                               unsigned long subbucketleft)
{
  if (sssp->ulongtab != NULL)
  {
    return sssp->ulongtab + sssp->bucketleftidx + subbucketleft
                          - sssp->partoffset;
  } else
  {
    return NULL;
  }
}

const uint32_t *gt_suffixsortspace_getptr_uint32(const GtSuffixsortspace *sssp,
                                                 unsigned long subbucketleft)
{
  if (sssp->uinttab != NULL)
  {
    return sssp->uinttab + sssp->bucketleftidx + subbucketleft
                         - sssp->partoffset;
  } else
  {
    return NULL;
  }
}

void gt_suffixsortspace_set(GtSuffixsortspace *sssp,
                            unsigned long subbucketleft,
                            unsigned long idx,
                            unsigned long value)
{
  gt_suffixsortspace_setdirect(sssp, sssp->bucketleftidx
                                       + subbucketleft
                                       + idx
                                       - sssp->partoffset,value);
}

unsigned long gt_suffixsortspace_bucketleftidx_get(const GtSuffixsortspace
                                                   *sssp)
{
  gt_assert(sssp != NULL);
  return sssp->bucketleftidx;
}

void gt_suffixsortspace_bucketleftidx_set(GtSuffixsortspace *sssp,
                                          unsigned long value)
{
  gt_assert(sssp->bucketleftidx == value || !sssp->currentexport);
  sssp->bucketleftidx = value;
}

void gt_suffixsortspace_partoffset_set(GtSuffixsortspace *sssp,
                                       unsigned long partoffset)
{
  gt_assert(sssp->partoffset == partoffset || !sssp->currentexport);
  sssp->partoffset = partoffset;
}

void gt_suffixsortspace_sortspace_delete(GtSuffixsortspace *sssp)
{
  gt_free(sssp->ulongtab);
  sssp->ulongtab = NULL;
  gt_free(sssp->uinttab);
  sssp->uinttab = NULL;
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
                                GT_UNUSED GtError *err)
{
  bool haserr = false;
  size_t basesize = sssp->ulongtab != NULL ? sizeof (*sssp->ulongtab)
                                           : sizeof (*sssp->uinttab);

  gt_xfwrite(sssp->ulongtab != NULL ? (void *) sssp->ulongtab
                                    : (void *) sssp->uinttab,
             basesize,
             (size_t) numberofsuffixes,
             outfpsuftab);
  return haserr ? -1 : 0;
}
