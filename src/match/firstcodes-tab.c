/*
  Copyright (c) 2011 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2011 Center for Bioinformatics, University of Hamburg

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

#include "core/assert_api.h"
#include "core/divmodmul.h"
#include "core/ma.h"
#include "core/disc_distri.h"
#include "core/log.h"
#include "core/spacecalc.h"
#include "core/hashmap-generic.h"
#include "firstcodes-tab.h"

DECLARE_HASHMAP(unsigned long, ul, uint32_t, u32, static, inline)
DEFINE_HASHMAP(unsigned long, ul, uint32_t, u32, gt_ht_ul_elem_hash,
               gt_ht_ul_elem_cmp, NULL_DESTRUCTOR, NULL_DESTRUCTOR, static,
               inline)

void gt_firstcodes_countocc_new(GtFirstcodestab *fct,
                                unsigned long numofsequences)
{
  fct->countocc = gt_calloc((size_t) (numofsequences+1),
                            sizeof (*fct->countocc));
  fct->countocc_small = gt_calloc((size_t) (numofsequences+1),
                                  sizeof (*fct->countocc_small));
  fct->countocc_exceptions = ul_u32_gt_hashmap_new();
  fct->usesample = false;
  gt_assert(fct->countocc_exceptions != NULL);
}

void gt_firstcodes_countocc_resize(GtFirstcodestab *fct,
                                   unsigned long numofdifferentcodes)
{
  fct->countocc = gt_realloc(fct->countocc,sizeof (*fct->countocc) *
                                           (numofdifferentcodes+1));
  fct->countocc_small = gt_realloc(fct->countocc_small,
                                   sizeof (*fct->countocc_small) *
                                           (numofdifferentcodes+1));
}

typedef struct
{
  unsigned long smallcount, smallsum,
                largecount, largesum,
                hugecount, hugesum;
} GtCountdistri_info;

static void gt_firstcodes_evaluate_distvalue(unsigned long key,
                                             unsigned long long value,
                                             void *data)
{
  GtCountdistri_info *cdi = (GtCountdistri_info *) data;

  gt_assert(value > 0);
  if (key <= UINT8_MAX)
  {
    cdi->smallcount++;
    cdi->smallsum += value;
  } else
  {
    if (key <= UINT16_MAX)
    {
      cdi->largecount++;
      cdi->largesum += value;
    } else
    {
      cdi->hugecount++;
      cdi->hugesum += value;
    }
  }
}

static void gt_firstcodes_evaluate_countdistri(const GtDiscDistri *countdistri)
{
  GtCountdistri_info cdi;
  unsigned long sum;
  size_t spacenow, spacedirectstore, spaceopt, spacewithhash;;

  cdi.smallcount = 0;
  cdi.smallsum = 0;
  cdi.largecount = 0;
  cdi.largesum = 0;
  cdi.hugecount = 0;
  cdi.hugesum = 0;
  gt_disc_distri_foreach(countdistri,gt_firstcodes_evaluate_distvalue,&cdi);
  sum = cdi.smallsum + cdi.largesum + cdi.hugesum;
  gt_log_log("small=%lu,%lu (%.2f)",cdi.smallcount,cdi.smallsum,
          (double) cdi.smallsum/sum);
  gt_log_log("large=%lu,%lu (%.2f)",cdi.largecount,cdi.largesum,
          (double) cdi.largesum/sum);
  gt_log_log("huge=%lu,%lu (%.2f)",cdi.hugecount,cdi.hugesum,
          (double) cdi.largesum/sum);
  spacenow = sizeof (uint32_t) * sum;
  spaceopt = sizeof (uint8_t) * sum;
  spacedirectstore = sizeof (uint32_t) * cdi.largesum;
  spacewithhash = (2 * sizeof (void *)) * cdi.largecount;
  gt_log_log("spacenow=%.2f, spaceopt (direct)=%2.f spaceopt (hash) =%.2f",
          GT_MEGABYTES(spacenow),
          GT_MEGABYTES(spaceopt+spacedirectstore),
          GT_MEGABYTES(spaceopt+spacewithhash));
}

#define GT_PARTIALSUM_COUNT_GET(IDX)             fct->countocc[IDX]
#define GT_PARTIALSUM_LEFTBORDER_SET(IDX,VALUE)  fct->countocc[IDX] = VALUE

unsigned long gt_firstcodes_partialsums(GtFirstcodestab *fct)
{
  unsigned long idx, partsum, maxbucketsize, largevalues = 0, samplecount = 0;
  uint32_t currentcount;
  GtDiscDistri *countdistri = gt_disc_distri_new();
  unsigned long bitmask;
  const unsigned long maxvalue = UINT8_MAX;

  gt_assert(fct->differentcodes < UINT32_MAX);
  for (idx = 0; idx < fct->differentcodes; idx++)
  {
    currentcount = GT_PARTIALSUM_COUNT_GET(idx);
    gt_assert(currentcount > 0);
    if (currentcount <= GT_FIRSTCODES_MAXSMALL)
    {
      fct->countocc_small[idx] = (uint8_t) currentcount;
    } else
    {
      uint32_t *valueptr;

      gt_assert(fct->countocc_exceptions != NULL);
      gt_assert(ul_u32_gt_hashmap_get(fct->countocc_exceptions,idx) == NULL);
      ul_u32_gt_hashmap_add(fct->countocc_exceptions, idx, currentcount);
      valueptr = ul_u32_gt_hashmap_get(fct->countocc_exceptions, idx);
      gt_assert(valueptr != NULL);
      gt_assert(*valueptr == currentcount);
      fct->countocc_small[idx] = 0;
      largevalues++;
    }
  }
  fct->overflow_index = 0;
  currentcount = GT_PARTIALSUM_COUNT_GET(0);
  partsum = (unsigned long) currentcount;
  maxbucketsize = (unsigned long) currentcount;
  gt_disc_distri_add(countdistri,(unsigned long) currentcount);
  fct->sampleshift = 9U;
  while (true)
  {
    fct->sampledistance = 1UL << fct->sampleshift;
    if (fct->sampledistance < fct->differentcodes)
    {
      break;
    }
    fct->sampleshift--;
  }
  bitmask = fct->sampledistance - 1;
  fct->numofsamples = 1UL + 1UL + fct->differentcodes/fct->sampledistance;
  fct->countocc_samples = gt_malloc(sizeof(*fct->countocc_samples) *
                                    fct->numofsamples);
  gt_assert(samplecount < fct->numofsamples);
  fct->countocc_samples[samplecount++] = partsum;
  for (idx = 1UL; idx < fct->differentcodes; idx++)
  {
    currentcount = GT_PARTIALSUM_COUNT_GET(idx);
    gt_disc_distri_add(countdistri,(unsigned long) currentcount);
    if (maxbucketsize < (unsigned long) currentcount)
    {
      maxbucketsize = (unsigned long) currentcount;
    }
    partsum += (unsigned long) currentcount;
    if ((idx & bitmask) == 0)
    {
      gt_assert(samplecount < fct->numofsamples);
      fct->countocc_samples[samplecount++] = partsum;
    }
    if (partsum <= maxvalue)
    {
      GT_PARTIALSUM_LEFTBORDER_SET(idx,(uint32_t) partsum);
    } else
    {
      gt_assert(idx > 0);
      partsum -= (unsigned long) currentcount;
      fct->overflow_index = idx;
      break;
    }
  }
  if (fct->overflow_index == 0)
  {
    gt_assert(partsum <= maxvalue);
    GT_PARTIALSUM_LEFTBORDER_SET(fct->differentcodes,(uint32_t) partsum);
  } else
  {
    unsigned long overflowcells = fct->differentcodes - fct->overflow_index + 1;
    unsigned long *overflow_leftborder_ptr;

    fct->overflow_leftborder
      = gt_malloc(sizeof (*fct->overflow_leftborder) * overflowcells);
    overflow_leftborder_ptr = fct->overflow_leftborder - fct->overflow_index;
    currentcount = GT_PARTIALSUM_COUNT_GET(fct->overflow_index);
    gt_disc_distri_add(countdistri,(unsigned long) currentcount);
    partsum += (unsigned long) currentcount;
    overflow_leftborder_ptr[fct->overflow_index] = partsum;
    for (idx = fct->overflow_index+1; idx < fct->differentcodes; idx++)
    {
      currentcount = GT_PARTIALSUM_COUNT_GET(idx);
      gt_disc_distri_add(countdistri,(unsigned long) currentcount);
      if (maxbucketsize < (unsigned long) currentcount)
      {
        maxbucketsize = (unsigned long) currentcount;
      }
      partsum += currentcount;
      if ((idx & bitmask) == 0)
      {
        gt_assert(samplecount < fct->numofsamples);
        fct->countocc_samples[samplecount++] = partsum;
      }
      overflow_leftborder_ptr[idx] = partsum;
    }
    overflow_leftborder_ptr[fct->differentcodes] = partsum;
  }
  gt_assert(idx > 0);
  if (partsum > fct->countocc_samples[samplecount-1])
  {
    gt_assert(samplecount < fct->numofsamples);
    fct->countocc_samples[samplecount] = partsum;
    fct->numofsamples = samplecount;
  } else
  {
    gt_assert(partsum == fct->countocc_samples[samplecount-1]);
    fct->numofsamples = samplecount - 1;
  }
  gt_firstcodes_evaluate_countdistri(countdistri);
  gt_disc_distri_delete(countdistri);
  return maxbucketsize;
}

unsigned long gt_firstcodes_get_leftborder(const GtFirstcodestab *fct,
                                           unsigned long idx)
{
  if (fct->usesample)
  {
    gt_assert(idx <= fct->numofsamples);
    return fct->countocc_samples[idx];
  } else
  {
    gt_assert(idx <= fct->differentcodes);
    if (fct->overflow_index == 0 || idx < fct->overflow_index)
    {
      return (unsigned long) fct->countocc[idx];
    } else
    {
      gt_assert(idx >= fct->overflow_index);
      return fct->overflow_leftborder[idx - fct->overflow_index];
    }
  }
}

unsigned long gt_firstcodes_numofallcodes(const GtFirstcodestab *fct)
{
  if (fct->usesample)
  {
    return fct->numofsamples;
  }
  return fct->differentcodes;
}

unsigned long gt_firstcodes_findfirstlarger(const GtFirstcodestab *fct,
                                            unsigned long suftaboffset)
{
  unsigned long left = 0, right, mid, midval, found;

  if (fct->usesample)
  {
    right = found = fct->numofsamples;
  } else
  {
    right = found = fct->differentcodes;
  }
  while (left+1 < right)
  {
    mid = GT_DIV2(left+right);
    midval = gt_firstcodes_get_leftborder(fct,mid);
    if (suftaboffset == midval)
    {
      return mid;
    }
    if (suftaboffset < midval)
    {
      found = mid;
      right = mid - 1;
    } else
    {
      left = mid + 1;
    }
  }
  gt_assert(suftaboffset <= gt_firstcodes_get_leftborder(fct,found));
  return found;
}

unsigned long gt_firstcodes_sample2full(const GtFirstcodestab *fct,
                                        unsigned long idx)
{
  gt_assert(idx <= fct->numofsamples);
  if (idx < fct->numofsamples)
  {
    return idx << fct->sampleshift;
  }
  return fct->differentcodes - 1;
}

unsigned long gt_firstcodes_idx2code(const GtFirstcodestab *fct,
                                     unsigned long idx)
{
  gt_assert(idx <= fct->differentcodes);
  if (idx == fct->differentcodes)
  {
    return fct->allfirstcodes[idx-1];
  }
  return fct->allfirstcodes[idx];
}

void gt_firstcodes_countocc_delete(GtFirstcodestab *fct)
{
  gt_free(fct->countocc);
  fct->countocc = NULL;
  gt_free(fct->countocc_small);
  fct->countocc_small = NULL;
  gt_hashtable_delete(fct->countocc_exceptions);
  gt_free(fct->countocc_samples);
}

void gt_firstcodes_countocc_setnull(GtFirstcodestab *fct)
{
  fct->countocc = NULL;
}

void **gt_firstcodes_countocc_address(GtFirstcodestab *fct)
{
  return (void **) &fct->countocc;
}

void gt_firstcodes_countocc_remap(GtFirstcodestab *fct,uint32_t *ptr)
{
  fct->countocc = ptr;
}
