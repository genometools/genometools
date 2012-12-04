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
#include "firstcodes-cache.h"

typedef uint8_t GtCountAFCtype;
#define GT_FIRSTCODES_MAXSMALL UINT8_MAX

typedef struct
{
  unsigned long differentcodes,
                numofsamples,
                sampledistance,
                hashmap_addcount,
                hashmap_getcount;
  unsigned int sampleshift;
  uint32_t *leftborder;
  GtCountAFCtype *countocc_small;
  GtHashtable *countocc_exceptions;
  unsigned long *leftborder_samples;
  GtStr *outfilenameleftborder;
  unsigned long differencemask, /* for extracting the difference */
                countmax;
  unsigned int shiftforcounts;
#ifdef _LP64
  GtArrayGtUlong bitchangepoints;
#endif
} GtFirstcodestab;

DECLARE_HASHMAP(unsigned long, ul, uint32_t, u32, static, inline)
DEFINE_HASHMAP(unsigned long, ul, uint32_t, u32, gt_ht_ul_elem_hash,
               gt_ht_ul_elem_cmp, NULL_DESTRUCTOR, NULL_DESTRUCTOR, static,
               inline)

#ifdef _LP64
#define GT_CHANGEPOINT_GET(CP)\
        unsigned long CP;\
        for (CP = 0; CP < fct->bitchangepoints.nextfreeGtUlong &&\
                     idx > fct->bitchangepoints.spaceGtUlong[CP]; CP++)\
            /* Nothing */ ;
#endif

#define GT_MODVALUEBITS 32U
#define GT_MODVALUEMASK UINT32_MAX

GT_UNUSED
static inline unsigned long gt_firstcodes_insertionindex(GtFirstcodestab *fct,
                                                         unsigned long idx)
{
#ifdef _LP64
  GT_CHANGEPOINT_GET(changepoint);
  gt_assert(idx < fct->differentcodes);
  if (fct->leftborder[idx] > 0)
  {
    return (unsigned long) --fct->leftborder[idx]
                           + (changepoint << GT_MODVALUEBITS);
  } else
  {
    gt_assert(changepoint > 0);
    changepoint--;
    fct->bitchangepoints.spaceGtUlong[changepoint]++;
    fct->leftborder[idx] = GT_MODVALUEMASK;
    return (unsigned long)
           fct->leftborder[idx] + (changepoint << GT_MODVALUEBITS);
  }
#else
  gt_assert(idx < fct->differentcodes && fct->leftborder[idx] > 0);
  return (unsigned long) --fct->leftborder[idx];
#endif
}

unsigned long gt_firstcodes_partialsums(GtFirstcodesspacelog *fcsl,
                                        GtFirstcodestab *fct,
                                        const unsigned long *differences,
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
                                    GtArrayGtIndexwithcode **binsearchcache,
                                    unsigned int addbscache_depth,
                                    bool withdistbits,
                                    GtLogger *logger);

unsigned long gt_firstcodes_accumulatecounts_merge(
                                        GtFirstcodestab *tab,
                                        unsigned long *differences,
                                        unsigned long differentcodes,
                                        const unsigned long *querystream_fst,
                                        const unsigned long *querystream_lst,
                                        unsigned long subjectindex,
                                        unsigned long subjectcode);

#endif
