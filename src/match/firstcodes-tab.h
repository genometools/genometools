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

#ifndef FIRSTCODES_TAB_H
#define FIRSTCODES_TAB_H

#include <inttypes.h>
#include "core/unused_api.h"
#include "core/str_api.h"
#include "core/hashmap-generic.h"
#include "core/logger_api.h"
#include "core/arraydef.h"
#include "marksubstring.h"
#include "firstcodes-spacelog.h"

#define GT_FIRSTCODES_MAXSMALL UINT8_MAX

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
  unsigned int modvaluebits;
  uint32_t modvaluemask;
  GtArrayGtUlong bitchangepoints;
} GtFirstcodestab;

DECLARE_HASHMAP(unsigned long, ul, uint32_t, u32, static, inline)
DEFINE_HASHMAP(unsigned long, ul, uint32_t, u32, gt_ht_ul_elem_hash,
               gt_ht_ul_elem_cmp, NULL_DESTRUCTOR, NULL_DESTRUCTOR, static,
               inline)

GT_UNUSED
static inline void gt_firstcodes_countocc_increment(GtFirstcodestab *fct,
                                                    unsigned long idx,
                                                    bool firstincrement)
{
  if (firstincrement)
  {
    fct->countocc_small[idx] = (uint8_t) 1;
  } else
  {
    fct->all_incrementcount++;
    if (fct->countocc_small[idx] > 0)
    {
      if (fct->countocc_small[idx] < GT_FIRSTCODES_MAXSMALL)
      {
        fct->countocc_small[idx]++;
      } else
      {
        if (fct->countocc_small[idx] == GT_FIRSTCODES_MAXSMALL)
        {
          fct->countocc_small[idx] = 0;
          ul_u32_gt_hashmap_add(fct->countocc_exceptions, idx, (uint32_t) 1);
          fct->hashmap_addcount++;
        }
      }
    } else
    {
      /* there is already an overflow for this index */
      if (fct->lastincremented_valueptr != NULL &&
          fct->lastincremented_idx == idx)
      {
        /* last index is identucal to current index. */
        gt_assert(*fct->lastincremented_valueptr < UINT32_MAX);
        (*fct->lastincremented_valueptr)++;
      } else
      {
        uint32_t *valueptr
          = ul_u32_gt_hashmap_get(fct->countocc_exceptions,idx);

        fct->hashmap_getcount++;
        gt_assert(valueptr != NULL);
        gt_assert(*valueptr < UINT32_MAX);
        (*valueptr)++;
        fct->lastincremented_idx = idx;
        fct->lastincremented_valueptr = valueptr;
      }
      fct->hashmap_incrementcount++;
    }
  }
}

#define GT_CHANGEPOINT_GET(CP)\
        unsigned long CP;\
        for (CP = 0; CP < fct->bitchangepoints.nextfreeGtUlong &&\
                     idx > fct->bitchangepoints.spaceGtUlong[CP]; CP++)\
            /* Nothing */ ;

GT_UNUSED
static inline unsigned long gt_firstcodes_insertionindex(GtFirstcodestab *fct,
                                                         unsigned long idx)
{
  GT_CHANGEPOINT_GET(changepoint);
  gt_assert(idx < fct->differentcodes);
  if (fct->leftborder[idx] > 0)
  {
    return (unsigned long) --fct->leftborder[idx]
                           + (changepoint << fct->modvaluebits);
  } else
  {
    gt_assert(changepoint > 0);
    changepoint--;
    fct->bitchangepoints.spaceGtUlong[changepoint]++;
    fct->leftborder[idx] = fct->modvaluemask;
    return (unsigned long)
           fct->leftborder[idx] + (changepoint << fct->modvaluebits);
  }
}

unsigned long gt_firstcodes_partialsums(GtFirstcodesspacelog *fcsl,
                                        GtFirstcodestab *fct,
                                        unsigned long expectedlastpartsum);

unsigned long gt_firstcodes_get_leftborder(const GtFirstcodestab *fct,
                                           unsigned long idx);

unsigned long gt_firstcodes_numofsamples(const GtFirstcodestab *fct);

unsigned long gt_firstcodes_findfirstsamplelarger(const GtFirstcodestab *fct,
                                                  unsigned long suftaboffset);

void gt_firstcodes_samples_delete(GtFirstcodesspacelog *fcsl,
                                  GtFirstcodestab *fct);

void gt_firstcodes_countocc_delete(GtFirstcodesspacelog *fcsl,
                                   GtFirstcodestab *fct);

void gt_firstcodes_tab_delete(GtFirstcodesspacelog *fcsl,
                              GtFirstcodestab *fct);

void gt_firstcodes_countocc_setnull(GtFirstcodestab *fct);

uint32_t **gt_firstcodes_leftborder_address(GtFirstcodestab *fct);

void gt_firstcodes_leftborder_remap(GtFirstcodestab *fct,uint32_t *ptr);

const GtStr *gt_firstcodes_outfilenameleftborder(const GtFirstcodestab *fct);

unsigned long gt_firstcodes_sample2full(const GtFirstcodestab *fct,
                                        unsigned long idx);

unsigned long gt_firstcodes_leftborder_entries(const GtFirstcodestab *fct);

unsigned long gt_firstcodes_leftborder_entries(const GtFirstcodestab *fct);

unsigned long gt_firstcodes_get_sample(const GtFirstcodestab *fct,
                                       unsigned long idx);

unsigned long gt_firstcodes_remdups(unsigned long *allfirstcodes,
                                    GtFirstcodesspacelog *fcsl,
                                    GtFirstcodestab *fct,
                                    unsigned long numofsequences,
                                    Gtmarksubstring *markprefix,
                                    Gtmarksubstring *marksuffix,
                                    GtLogger *logger);

#endif
