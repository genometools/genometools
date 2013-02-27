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

#include "core/xansi_api.h"
#include "core/fa.h"
#include "core/assert_api.h"
#include "core/divmodmul.h"
#include "core/ma.h"
#include "core/mathsupport.h"
#ifdef SKDEBUG
#include "core/disc_distri_api.h"
#endif
#include "core/log.h"
#include "core/spacecalc.h"
#include "core/hashmap-generic.h"
#include "core/arraydef.h"
#include "randomcodes-tab.h"
#include "firstcodes-psbuf.h"

void gt_randomcodes_countocc_new(GtFirstcodesspacelog *fcsl,
                                       GtRandomcodestab *rct,
                                       unsigned long numofsequences)
{
  rct->countocc_small = gt_calloc((size_t) (numofsequences+1),
                                  sizeof (*rct->countocc_small));
  GT_FCI_ADDWORKSPACE(fcsl,"countocc_small",
                      sizeof (*rct->countocc_small) *
                      (numofsequences+1));
  rct->countocc_exceptions = ul_u32_gt_hashmap_new();
  gt_assert(rct->countocc_exceptions != NULL);
  rct->outfilenameleftborder = NULL;
  rct->leftborder_samples = NULL;
#ifdef _LP64
  rct->modvaluebits = 32U; /* XXX remove the following later */
  if (rct->modvaluebits == 32U)
  {
    rct->modvaluemask = UINT32_MAX;
  } else
  {
    rct->modvaluemask = (uint32_t) ((1UL << rct->modvaluebits) - 1);
  }
  GT_INITARRAY(&rct->bitchangepoints,GtUlong);
#endif
  rct->differentcodes = numofsequences;
}

void gt_randomcodes_countocc_resize(GtFirstcodesspacelog *fcsl,
                                          GtRandomcodestab *rct,
                                          unsigned long numofdifferentcodes)
{
  rct->countocc_small = gt_realloc(rct->countocc_small,
                                   sizeof (*rct->countocc_small) *
                                           (numofdifferentcodes+1));
  GT_FCI_SUBTRACTADDWORKSPACE(fcsl,"countocc_small",
                              sizeof (*rct->countocc_small) *
                              (numofdifferentcodes+1));
}

#ifdef SKDEBUG

typedef struct
{
  unsigned long smallcount, smallsum,
                largecount, largesum,
                hugecount, hugesum;
} GtCountdistri_info;

static void gt_randomcodes_evaluate_distvalue(unsigned long key,
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

static void gt_randomcodes_evaluate_countdistri(const GtDiscDistri *countdistri)
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
  gt_disc_distri_foreach(countdistri,gt_randomcodes_evaluate_distvalue,&cdi);
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
#endif

#ifdef SKDEBUG
static void checkcodesorder(const unsigned long *tab,unsigned long len,
                            bool allowequal)
{
  unsigned long idx;

  for (idx=1UL; idx < len; idx++)
  {
    gt_assert(tab[idx-1] < tab[idx] || (allowequal && tab[idx-1] == tab[idx]));
  }
}
#endif

unsigned long gt_randomcodes_remdups(unsigned long *allrandomcodes,
    unsigned int codesize, unsigned long numofcodes, GtLogger *logger)
{
  unsigned long numofdifferentcodes = 0,
                shift = (unsigned long)GT_MULT2(GT_UNITSIN2BITENC - codesize);
  if (numofcodes != 0)
  {
    unsigned long *storeptr, *readptr;

    for (storeptr = allrandomcodes, readptr = allrandomcodes+1;
         readptr < allrandomcodes + numofcodes;
         readptr++)
    {
      if ((*storeptr ^ *readptr) << shift)
      {
        storeptr++;
        *storeptr = *readptr;
      }
    }
    numofdifferentcodes = (unsigned long) (storeptr - allrandomcodes + 1);
    if (numofdifferentcodes < numofcodes)
    {
#ifdef SKDEBUG
      checkcodesorder(allrandomcodes,numofdifferentcodes,false);
#endif
    }
  }
  gt_logger_log(logger,"number of different bucket codes=%lu (%.2f%%) "
                       "of %lu sampled codes",
                numofdifferentcodes,
                100.00 * (double) numofdifferentcodes/numofcodes,
                numofcodes);
  return numofdifferentcodes;
}

static uint32_t gt_randomcodes_countocc_get(const GtRandomcodestab *rct,
                                           unsigned long idx)
{
  if (rct->countocc_small[idx] != GT_RANDOMCODES_COUNTOCC_OVERFLOW)
  {
    return (uint32_t) rct->countocc_small[idx];
  } else
  {
    uint32_t *valueptr = ul_u32_gt_hashmap_get(rct->countocc_exceptions,idx);

    gt_assert(valueptr != NULL);
    return *valueptr + (uint32_t) GT_RANDOMCODES_MAXSMALL;
  }
}

#ifdef _LP64
#define GT_RCT_PARTIALSUM_LEFTBORDER_SET(BUF,VALUE)\
        if ((BUF)->nextfree == (BUF)->allocated)\
        {\
          gt_leftborderbuffer_flush(BUF);\
        }\
        (BUF)->spaceuint32_t[(BUF)->nextfree++]\
          = (uint32_t) ((VALUE) & rct->modvaluemask)
#else
#define GT_RCT_PARTIALSUM_LEFTBORDER_SET(BUF,VALUE)\
        if ((BUF)->nextfree == (BUF)->allocated)\
        {\
          gt_leftborderbuffer_flush(BUF);\
        }\
        (BUF)->spaceuint32_t[(BUF)->nextfree++] = (uint32_t) (VALUE)
#endif

#define GT_RANDOMCODES_ADD_SAMPLE(PARTSUM)\
        gt_assert(samplecount < rct->numofsamples);\
        rct->leftborder_samples[samplecount++] = PARTSUM

unsigned long gt_randomcodes_partialsums(GtFirstcodesspacelog *fcsl,
                                        GtRandomcodestab *rct,
                                        GT_UNUSED unsigned long
                                                            expectedlastpartsum)
{
  unsigned long idx, partsum, maxbucketsize, bitmask, samplecount = 0,
                spacewithhashmap = 0, spacewithouthashmap = 0;
  uint32_t currentcount;
  GtLeftborderOutbuffer *leftborderbuffer_all = NULL;
#ifdef _LP64
  const unsigned int btp = gt_determinebitspervalue(expectedlastpartsum);
  unsigned long exceedvalue = 1UL << rct->modvaluebits;
#endif
#ifdef SKDEBUG
  GtDiscDistri *countdistri = gt_disc_distri_new();
#endif

  gt_assert(rct->differentcodes < UINT32_MAX);
  gt_log_log("hashmap_addcount=%lu (%.2f%%)",rct->hashmap_addcount,
                  100.0 * (double) rct->hashmap_addcount/
                                   rct->differentcodes);
  gt_log_log("hashmap_incrementcount=%lu (%.2f%%)",
                  rct->hashmap_incrementcount,
                  100.0 * (double) rct->hashmap_incrementcount/
                                   rct->all_incrementcount);
  gt_log_log("hashmap_getcount=%lu (%.2f%%)",
                  rct->hashmap_getcount,
                  100.0 * (double) rct->hashmap_getcount/
                                   rct->all_incrementcount);

#ifdef _LP64
  if (btp <= rct->modvaluebits)
  {
    rct->bitchangepoints.allocatedGtUlong = 0;
    rct->bitchangepoints.spaceGtUlong = NULL;
  } else
  {
    rct->bitchangepoints.allocatedGtUlong = 1UL << (btp - rct->modvaluebits);
    gt_log_log("lastpartsum=%lu, bitchangepoints.allocated=%lu",
              expectedlastpartsum,rct->bitchangepoints.allocatedGtUlong);
    rct->bitchangepoints.spaceGtUlong
      = gt_malloc(sizeof (*rct->bitchangepoints.spaceGtUlong)
                  * rct->bitchangepoints.allocatedGtUlong);
  }
  rct->bitchangepoints.nextfreeGtUlong = 0;
#endif
  currentcount = gt_randomcodes_countocc_get(rct,0);
  partsum = (unsigned long) currentcount;
  maxbucketsize = (unsigned long) currentcount;
#ifdef SKDEBUG
  gt_disc_distri_add(countdistri,(unsigned long) currentcount);
#endif
  rct->sampleshift = 9U;
  while (true)
  {
    rct->sampledistance = 1UL << rct->sampleshift;
    if (rct->sampledistance < rct->differentcodes)
    {
      break;
    }
    rct->sampleshift--;
  }
  bitmask = rct->sampledistance - 1;
  rct->numofsamples = 1UL + 1UL + rct->differentcodes/rct->sampledistance;
  rct->leftborder_samples = gt_malloc(sizeof (*rct->leftborder_samples) *
                                      rct->numofsamples);
  GT_FCI_ADDWORKSPACE(fcsl,"leftborder_samples",
                      sizeof (*rct->leftborder_samples) * rct->numofsamples);
  GT_RANDOMCODES_ADD_SAMPLE(partsum);
  leftborderbuffer_all = gt_leftborderbuffer_new("leftborder",fcsl);
  GT_RCT_PARTIALSUM_LEFTBORDER_SET(leftborderbuffer_all,partsum);
  for (idx = 1UL; idx < rct->differentcodes; idx++)
  {
    currentcount = gt_randomcodes_countocc_get(rct,idx);
#ifdef _LP64
    gt_assert(currentcount <= rct->modvaluemask);
#endif
#ifdef SKDEBUG
    gt_disc_distri_add(countdistri,(unsigned long) currentcount);
#endif
    if (maxbucketsize < (unsigned long) currentcount)
    {
      maxbucketsize = (unsigned long) currentcount;
    }
    partsum += currentcount;
#ifdef _LP64
    if (rct->bitchangepoints.allocatedGtUlong > 0 && partsum >= exceedvalue)
    {
      gt_assert(idx > 0 && rct->bitchangepoints.nextfreeGtUlong <
                           rct->bitchangepoints.allocatedGtUlong);
      gt_assert(rct->bitchangepoints.spaceGtUlong != NULL);
      rct->bitchangepoints.spaceGtUlong
        [rct->bitchangepoints.nextfreeGtUlong++] = idx-1;
      exceedvalue
        = ((exceedvalue >> rct->modvaluebits) + 1) << rct->modvaluebits;
    }
#endif
    if ((idx & bitmask) == 0)
    {
      GT_RANDOMCODES_ADD_SAMPLE(partsum);
    }
    GT_RCT_PARTIALSUM_LEFTBORDER_SET(leftborderbuffer_all,partsum);
  }
  GT_RCT_PARTIALSUM_LEFTBORDER_SET(leftborderbuffer_all,partsum);
  rct->outfilenameleftborder
      = gt_leftborderbuffer_delete(leftborderbuffer_all,fcsl,
                                   gt_randomcodes_leftborder_entries(rct));
  if (partsum > rct->leftborder_samples[samplecount-1])
  {
    GT_RANDOMCODES_ADD_SAMPLE(partsum);
  } else
  {
    gt_assert(partsum == rct->leftborder_samples[samplecount-1]);
  }
  gt_assert(expectedlastpartsum == partsum);
  rct->numofsamples = samplecount-1;
#ifdef SKDEBUG
  gt_randomcodes_evaluate_countdistri(countdistri);
  gt_disc_distri_delete(countdistri);
#endif
  gt_assert (rct->countocc_small != NULL);
  gt_free(rct->countocc_small);
  GT_FCI_SUBTRACTWORKSPACE(fcsl,"countocc_small");
  rct->countocc_small = NULL;
  if (rct->hashmap_addcount > 0 && gt_ma_bookkeeping_enabled())
  {
    spacewithhashmap = gt_ma_get_space_current() + gt_fa_get_space_current();
  }
  gt_hashtable_delete(rct->countocc_exceptions);
  if (rct->hashmap_addcount > 0 && gt_ma_bookkeeping_enabled())
  {
    unsigned long hashmapspace;

    spacewithouthashmap = gt_ma_get_space_current() + gt_fa_get_space_current();
    gt_assert(spacewithouthashmap < spacewithhashmap);
    hashmapspace = spacewithhashmap - spacewithouthashmap;
    gt_log_log("space for hashmap=%.2f (%lu bytes per entry)",
               GT_MEGABYTES(hashmapspace),hashmapspace/rct->hashmap_addcount);
  }
  rct->countocc_exceptions = NULL;
  return maxbucketsize;
}

unsigned long gt_randomcodes_get_sample(const GtRandomcodestab *rct,
                                       unsigned long idx)
{
  gt_assert(idx <= rct->numofsamples);
  return rct->leftborder_samples[idx];
}

unsigned long gt_randomcodes_get_leftborder(const GtRandomcodestab *rct,
                                           unsigned long idx)
{
#ifdef _LP64
  GT_CHANGEPOINT_GET_RCT(changepoint);

  return (unsigned long) rct->leftborder[idx]
                         + (changepoint << rct->modvaluebits);
#else
  return (unsigned long) rct->leftborder[idx];
#endif
}

unsigned long gt_randomcodes_leftborder_entries(const GtRandomcodestab *rct)
{
  return rct->differentcodes + 1;
}

unsigned long gt_randomcodes_numofsamples(const GtRandomcodestab *rct)
{
  return rct->numofsamples;
}

unsigned long gt_randomcodes_findfirstsamplelarger(const GtRandomcodestab *rct,
                                                  unsigned long suftaboffset)
{
  unsigned long left = 0, right, mid, midval, found;

  right = found = rct->numofsamples;
  while (left+1 < right)
  {
    mid = GT_DIV2(left+right);
    midval = gt_randomcodes_get_sample(rct,mid);
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
  gt_assert(suftaboffset <= gt_randomcodes_get_sample(rct,found));
  return found;
}

unsigned long gt_randomcodes_sample2full(const GtRandomcodestab *rct,
                                        unsigned long idx)
{
  gt_assert(idx <= rct->numofsamples);
  if (idx < rct->numofsamples)
  {
    return idx << rct->sampleshift;
  }
  return rct->differentcodes - 1;
}

void gt_randomcodes_samples_delete(GtFirstcodesspacelog *fcsl,
                                  GtRandomcodestab *rct)
{
  if (rct->leftborder_samples != NULL)
  {
    gt_free(rct->leftborder_samples);
    GT_FCI_SUBTRACTWORKSPACE(fcsl,"leftborder_samples");
    rct->leftborder_samples = NULL;
  }
}

void gt_randomcodes_countocc_delete(GtFirstcodesspacelog *fcsl,
                                   GtRandomcodestab *rct)
{
  if (rct->countocc_small != NULL)
  {
    GT_FCI_SUBTRACTWORKSPACE(fcsl,"countocc_small");
    gt_free(rct->countocc_small);
    rct->countocc_small = NULL;
  }
  gt_hashtable_delete(rct->countocc_exceptions);
  rct->countocc_exceptions = NULL;
}

void gt_randomcodes_tab_delete(GtFirstcodesspacelog *fcsl,GtRandomcodestab *rct)
{
  gt_randomcodes_samples_delete(fcsl,rct);
  gt_str_delete(rct->outfilenameleftborder);
  rct->outfilenameleftborder = NULL;
#ifdef _LP64
  GT_FREEARRAY(&rct->bitchangepoints,GtUlong);
#endif
}

void gt_randomcodes_countocc_setnull(GtRandomcodestab *rct)
{
  rct->leftborder = NULL;
  rct->countocc_small = NULL;
  rct->leftborder_samples = NULL;
  rct->countocc_exceptions = NULL;
  rct->differentcodes = 0;
  rct->lastincremented_idx = 0;
  rct->lastincremented_valueptr = NULL;
  rct->hashmap_addcount = 0;
  rct->hashmap_incrementcount = 0;
  rct->all_incrementcount = 0;
  rct->hashmap_getcount = 0;
  rct->outfilenameleftborder = NULL;
#ifdef _LP64
  GT_INITARRAY(&rct->bitchangepoints,GtUlong);
#endif
}

uint32_t **gt_randomcodes_leftborder_address(GtRandomcodestab *rct)
{
  return &rct->leftborder;
}

void gt_randomcodes_leftborder_remap(GtRandomcodestab *rct,uint32_t *ptr)
{
  rct->leftborder = ptr;
}

const GtStr *gt_randomcodes_outfilenameleftborder(const GtRandomcodestab *rct)
{
  return rct->outfilenameleftborder;
}
