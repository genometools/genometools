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
#include "core/intbits.h"
#include "core/arraydef.h"
#include "sfx-suffixgetset.h"

struct GtSuffixsortspace
{
  bool unmapsortspace, currentexport;
  Definedunsignedlong longestidx;
  uint32_t *uinttab;
  size_t basesize;
  GtSuffixsortspace_exportptr exportptr;
  GtUword maxindex,
                maxvalue,
                partoffset,
                bucketleftidx,
                *ulongtab;
};

static bool gt_decide_to_use_uint(bool useuint,GtUword maxvalue)
{
  if (useuint && maxvalue <= (GtUword) UINT_MAX)
  {
    return true;
  }
  return false;
}

uint64_t gt_suffixsortspace_requiredspace(GtUword numofentries,
                                          GtUword maxvalue,
                                          bool useuint)
{
  uint64_t requiredspace = (uint64_t) sizeof (GtSuffixsortspace);

  if (gt_decide_to_use_uint(useuint,maxvalue))
  {
    gt_assert(maxvalue <= (GtUword) UINT_MAX);
    requiredspace += (uint64_t) numofentries * (uint64_t) sizeof (uint32_t);
  } else
  {
    requiredspace += (uint64_t) numofentries *
                     (uint64_t) sizeof (GtUword);
  }
  return requiredspace;
}

static void gt_suffixsortspace_overflow_abort(GT_UNUSED const char *f,
                                              GT_UNUSED int l,
                                              void *data)
{
  fprintf(stderr, "error: overflow detected while calculating size of "
                  "suffix sorting space: "GT_WU" * "GT_WU" bytes is too large "
                  "for " "the current platform, please recompile GenomeTools "
                  "with support for a larger address space to prevent this "
                  "(e.g. 64 bit instead of 32 bit) or use the `-parts' "
                  "option.\n",
                  (GtUword) sizeof (GtUword),
                  *(GtUword*) data);
  exit(GT_EXIT_PROGRAMMING_ERROR);
}

GtSuffixsortspace *gt_suffixsortspace_new(GtUword numofentries,
                                          GtUword maxvalue,
                                          bool useuint,
                                          GT_UNUSED GtLogger *logger)
{
  GtSuffixsortspace *suffixsortspace;
  GtUword sufspacesize;

  gt_assert(numofentries > 0);
  suffixsortspace = gt_malloc(sizeof (*suffixsortspace));
  suffixsortspace->maxindex = numofentries-1;
  suffixsortspace->maxvalue = maxvalue;
  suffixsortspace->longestidx.defined = false;
  suffixsortspace->longestidx.valueunsignedlong = 0;
  suffixsortspace->exportptr.ulongtabsectionptr = NULL;
  suffixsortspace->exportptr.uinttabsectionptr = NULL;
  suffixsortspace->currentexport = false;
#if defined (_LP64) || defined (_WIN64)
  gt_logger_log(logger,"suftab uses %dbit values: "
                         "maxvalue="GT_WU",numofentries="GT_WU"",
                         gt_decide_to_use_uint(useuint,maxvalue) ? 32 : 64,
                         maxvalue,numofentries);
#endif
  suffixsortspace->basesize = gt_decide_to_use_uint(useuint,maxvalue)
                                ? sizeof (*suffixsortspace->uinttab)
                                : sizeof (*suffixsortspace->ulongtab);
  sufspacesize
    = gt_safe_mult_ulong_check((GtUword) suffixsortspace->basesize,
                               numofentries,
                               gt_suffixsortspace_overflow_abort,
                               &numofentries);
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
                                                   GtUword numofentries,
                                                   GtUword maxvalue,
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

GtUword gt_suffixsortspace_getdirect(const GtSuffixsortspace *sssp,
                                           GtUword idx)
{
  gt_assert(idx <= sssp->maxindex);
  if (sssp->ulongtab != NULL)
  {
    return sssp->ulongtab[idx];
  }
  gt_assert(sssp->uinttab != NULL);
  return (GtUword) sssp->uinttab[idx];
}

void gt_suffixsortspace_updatelongest(GtSuffixsortspace *sssp,GtUword idx)
{
  sssp->longestidx.defined = true;
  sssp->longestidx.valueunsignedlong = sssp->bucketleftidx + idx;
}

void gt_suffixsortspace_setdirect(GtSuffixsortspace *sssp,
                                  GtUword idx,
                                  GtUword value)
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
    gt_assert (value <= (GtUword) UINT_MAX);
    sssp->uinttab[idx] = (uint32_t) value;
  }
}

void gt_suffixsortspace_showrange(const GtSuffixsortspace *sssp,
                                  GtUword subbucketleft,
                                  GtUword width)
{
  GtUword idx;

  printf(""GT_WU","GT_WU"=",sssp->bucketleftidx+subbucketleft-sssp->partoffset,
                    sssp->bucketleftidx+subbucketleft+width-1-sssp->partoffset);
  for (idx=sssp->bucketleftidx+subbucketleft-sssp->partoffset;
       idx<sssp->bucketleftidx+subbucketleft+width-sssp->partoffset;
       idx++)
  {
    printf(" "GT_WU"", gt_suffixsortspace_getdirect(sssp,idx));
  }
}

void gt_suffixsortspace_checkorder(const GtSuffixsortspace *sssp,
                                   GtUword subbucketleft,
                                   GtUword width)
{
  GtUword idx, currentpos;
  GT_UNUSED GtUword prevpos;

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
                                  GtUword subbucketleft,
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

GtUword gt_suffixsortspace_get(const GtSuffixsortspace *sssp,
                                     GtUword subbucketleft,
                                     GtUword idx)
{
  return gt_suffixsortspace_getdirect(sssp, sssp->bucketleftidx
                                              + subbucketleft
                                              + idx
                                              - sssp->partoffset);
}

const GtUword *gt_suffixsortspace_getptr_ulong(
                                               const GtSuffixsortspace *sssp,
                                               GtUword subbucketleft)
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
                                                 GtUword subbucketleft)
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
                            GtUword subbucketleft,
                            GtUword idx,
                            GtUword value)
{
  gt_suffixsortspace_setdirect(sssp, sssp->bucketleftidx
                                       + subbucketleft
                                       + idx
                                       - sssp->partoffset,value);
}

GtUword gt_suffixsortspace_bucketleftidx_get(const GtSuffixsortspace
                                                   *sssp)
{
  gt_assert(sssp != NULL);
  return sssp->bucketleftidx;
}

void gt_suffixsortspace_bucketleftidx_set(GtSuffixsortspace *sssp,
                                          GtUword value)
{
  gt_assert(sssp->bucketleftidx == value || !sssp->currentexport);
  sssp->bucketleftidx = value;
}

void gt_suffixsortspace_partoffset_set(GtSuffixsortspace *sssp,
                                       GtUword partoffset)
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

GtUword *gt_suffixsortspace_ulong_get(const GtSuffixsortspace *sssp)
{
  gt_assert(sssp->ulongtab != NULL);
  return (GtUword *) sssp->ulongtab; /* XXX constrain the type cast */
}

GtUword gt_suffixsortspace_longest(const GtSuffixsortspace *sssp)
{
  gt_assert(sssp->longestidx.defined);
  return sssp->longestidx.valueunsignedlong;
}

void gt_suffixsortspace_to_file (FILE *outfpsuftab,
                                const GtSuffixsortspace *sssp,
                                GtUword numberofsuffixes)
{
  size_t basesize = sssp->ulongtab != NULL ? sizeof (*sssp->ulongtab)
                                           : sizeof (*sssp->uinttab);

  gt_xfwrite(sssp->ulongtab != NULL ? (void *) sssp->ulongtab
                                    : (void *) sssp->uinttab,
             basesize,
             (size_t) numberofsuffixes,
             outfpsuftab);
}

void gt_suffixsortspace_compressed_to_file (const GtSuffixsortspace *sssp,
                                            GtBitbuffer *bb,
                                            GtUword numberofsuffixes)
{
  if (sssp->ulongtab != NULL)
  {
    gt_bitbuffer_next_ulongtab(bb,sssp->ulongtab,numberofsuffixes);
  } else
  {
    gt_assert (sssp->uinttab != NULL);
    gt_bitbuffer_next_uint32tab(bb,sssp->uinttab,numberofsuffixes);
  }
}
