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

#ifndef RANDOMCODES_TAB_H
#define RANDOMCODES_TAB_H

#include <inttypes.h>
#include "core/unused_api.h"
#include "core/str_api.h"
#include "core/hashmap-generic.h"
#include "core/logger_api.h"
#include "core/arraydef.h"
#include "firstcodes-tab.h"
#include "firstcodes-spacelog.h"

#define GT_RANDOMCODES_MAXSMALL (UINT8_MAX - 1)
#define GT_RANDOMCODES_COUNTOCC_OVERFLOW UINT8_MAX

typedef struct
{
  unsigned long differentcodes,
                numofsamples,
                sampledistance,
                hashmap_addcount,
                hashmap_incrementcount,
                hashmap_getcount,
                all_incrementcount;
  unsigned int sampleshift;
  uint32_t *leftborder;
  uint8_t *countocc_small;
  GtHashtable *countocc_exceptions;
  unsigned long *leftborder_samples;
  GtStr *outfilenameleftborder;
  unsigned long lastincremented_idx;
  uint32_t *lastincremented_valueptr;
#ifdef _LP64
  uint32_t modvaluemask;
  unsigned int modvaluebits;
  GtArrayGtUlong bitchangepoints;
#endif
} GtRandomcodestab;

GT_UNUSED
static inline void gt_randomcodes_countocc_increment(GtRandomcodestab *rct,
                                                    unsigned long idx)
{
  rct->all_incrementcount++;
  if (rct->countocc_small[idx] != GT_RANDOMCODES_COUNTOCC_OVERFLOW)
  {
    if (rct->countocc_small[idx] < GT_RANDOMCODES_MAXSMALL)
    {
      rct->countocc_small[idx]++;
    } else
    {
      gt_assert (rct->countocc_small[idx] == GT_RANDOMCODES_MAXSMALL);
      rct->countocc_small[idx] = GT_RANDOMCODES_COUNTOCC_OVERFLOW;
      rct->lastincremented_valueptr
        = ul_u32_gt_hashmap_add_and_return_storage(rct->countocc_exceptions,
            idx, (uint32_t) 1);
      rct->lastincremented_idx = idx;
      rct->hashmap_addcount++;
    }
  } else
  {
    /* there is already an overflow for this index */
    if (rct->lastincremented_valueptr != NULL &&
        rct->lastincremented_idx == idx)
    {
      /* last index is identical to current index. */
      gt_assert(*rct->lastincremented_valueptr < UINT32_MAX);
      (*rct->lastincremented_valueptr)++;
    } else
    {
      uint32_t *valueptr
        = ul_u32_gt_hashmap_get(rct->countocc_exceptions,idx);

      rct->hashmap_getcount++;
      gt_assert(valueptr != NULL && *valueptr < UINT32_MAX);
      (*valueptr)++;
      rct->lastincremented_idx = idx;
      rct->lastincremented_valueptr = valueptr;
    }
    rct->hashmap_incrementcount++;
  }
}

#ifdef _LP64
#define GT_CHANGEPOINT_GET_RCT(CP)\
        unsigned long CP;\
        for (CP = 0; CP < rct->bitchangepoints.nextfreeGtUlong &&\
                     idx > rct->bitchangepoints.spaceGtUlong[CP]; CP++)\
            /* Nothing */ ;
#endif

GT_UNUSED
static inline unsigned long gt_randomcodes_insertionindex(GtRandomcodestab *rct,
                                                         unsigned long idx)
{
#ifdef _LP64
  GT_CHANGEPOINT_GET_RCT(changepoint);
  gt_assert(idx < rct->differentcodes);
  if (rct->leftborder[idx] > 0)
  {
    return (unsigned long) --rct->leftborder[idx]
                           + (changepoint << rct->modvaluebits);
  } else
  {
    gt_assert(changepoint > 0);
    changepoint--;
    rct->bitchangepoints.spaceGtUlong[changepoint]++;
    rct->leftborder[idx] = rct->modvaluemask;
    return (unsigned long)
           rct->leftborder[idx] + (changepoint << rct->modvaluebits);
  }
#else
  gt_assert(idx < rct->differentcodes && rct->leftborder[idx] > 0);
  return (unsigned long) --rct->leftborder[idx];
#endif
}

unsigned long gt_randomcodes_partialsums(GtFirstcodesspacelog *fcsl,
                                        GtRandomcodestab *rct,
                                        unsigned long expectedlastpartsum);

unsigned long gt_randomcodes_get_leftborder(const GtRandomcodestab *rct,
                                           unsigned long idx);

unsigned long gt_randomcodes_numofsamples(const GtRandomcodestab *rct);

unsigned long gt_randomcodes_findfirstsamplelarger(const GtRandomcodestab *rct,
                                                  unsigned long suftaboffset);

void gt_randomcodes_samples_delete(GtFirstcodesspacelog *fcsl,
                                  GtRandomcodestab *rct);

void gt_randomcodes_countocc_delete(GtFirstcodesspacelog *fcsl,
                                   GtRandomcodestab *rct);

void gt_randomcodes_tab_delete(GtFirstcodesspacelog *fcsl,
                              GtRandomcodestab *rct);

void gt_randomcodes_countocc_new(GtFirstcodesspacelog *fcsl,
                                       GtRandomcodestab *rct,
                                       unsigned long numofsequences);

void gt_randomcodes_countocc_resize(GtFirstcodesspacelog *fcsl,
                                          GtRandomcodestab *rct,
                                          unsigned long numofdifferentcodes);

void gt_randomcodes_countocc_setnull(GtRandomcodestab *rct);

uint32_t **gt_randomcodes_leftborder_address(GtRandomcodestab *rct);

void gt_randomcodes_leftborder_remap(GtRandomcodestab *rct,uint32_t *ptr);

const GtStr *gt_randomcodes_outfilenameleftborder(const GtRandomcodestab *rct);

unsigned long gt_randomcodes_sample2full(const GtRandomcodestab *rct,
                                        unsigned long idx);

unsigned long gt_randomcodes_leftborder_entries(const GtRandomcodestab *rct);

unsigned long gt_randomcodes_leftborder_entries(const GtRandomcodestab *rct);

unsigned long gt_randomcodes_get_sample(const GtRandomcodestab *rct,
                                       unsigned long idx);

unsigned long gt_randomcodes_remdups(unsigned long *allrandomcodes,
    unsigned int codesize, unsigned long numofcodes, GtLogger *logger);

#endif
