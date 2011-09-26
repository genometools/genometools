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
#include "core/hashmap-generic.h"

#define GT_FIRSTCODES_MAXSMALL UINT8_MAX

typedef struct
{
  unsigned long differentcodes,
                *allfirstcodes,
                overflow_index,
                numofsamples,
                sampledistance;
  uint32_t *countocc;
  uint8_t *countocc_small;
  GtHashtable *countocc_exceptions;
  unsigned long *overflow_leftborder;
  unsigned long *countocc_samples;
  bool usesample;
} GtFirstcodestab;

void gt_firstcodes_countocc_new(GtFirstcodestab *fct,
                                unsigned long numofsequences);

void gt_firstcodes_countocc_resize(GtFirstcodestab *fct,
                                   unsigned long numofdifferentcodes);

GT_UNUSED
static inline void gt_firstcodes_countocc_increment(GtFirstcodestab *fct,
                                                    unsigned long idx)
{
  fct->countocc[idx]++;
}

GT_UNUSED
static inline unsigned long gt_firstcodes_insertionindex(GtFirstcodestab *fct,
                                                         unsigned long idx)
{
  gt_assert(idx < fct->differentcodes);
  if (fct->overflow_index == 0 || idx < fct->overflow_index)
  {
    return (unsigned long) --fct->countocc[idx];
  } else
  {
    return --fct->overflow_leftborder[idx - fct->overflow_index];
  }
}

unsigned long gt_firstcodes_partialsums(GtFirstcodestab *fct);

unsigned long gt_firstcodes_get_leftborder(const GtFirstcodestab *fct,
                                           unsigned long idx);

unsigned long gt_firstcodes_numofallcodes(const GtFirstcodestab *fct);

unsigned long gt_firstcodes_findfirstlarger(const GtFirstcodestab *fct,
                                            unsigned long suftaboffset);

unsigned long gt_firstcodes_idx2code(const GtFirstcodestab *fct,
                                     unsigned long idx);

void gt_firstcodes_countocc_delete(GtFirstcodestab *fct);

void gt_firstcodes_countocc_setnull(GtFirstcodestab *fct);

void **gt_firstcodes_countocc_address(GtFirstcodestab *fct);

void gt_firstcodes_countocc_remap(GtFirstcodestab *fct,uint32_t *ptr);

unsigned long gt_firstcodes_sample2full(const GtFirstcodestab *fct,
                                        unsigned long idx);

#endif
