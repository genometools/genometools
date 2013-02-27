/*
  Copyright (c) 2008 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2008 Center for Bioinformatics, University of Hamburg

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

#include "annotationsketch/text_width_calculator_rep.h"
#include "annotationsketch/text_width_calculator.h"
#include "core/assert_api.h"
#include "core/class_alloc.h"
#include "core/ma.h"
#include "core/thread_api.h"
#include "core/unused_api.h"

struct GtTextWidthCalculatorMembers {
  unsigned int reference_count;
  GtRWLock *lock;
};

struct GtTextWidthCalculatorClass {
  size_t size;
  TextWidthCalculatorGetTextWidth get_text_width;
  TextWidthCalculatorFreeFunc free;
};

const GtTextWidthCalculatorClass* gt_text_width_calculator_class_new(
                                 size_t size,
                                 TextWidthCalculatorGetTextWidth get_text_width,
                                 TextWidthCalculatorFreeFunc free)
{
  GtTextWidthCalculatorClass *c_class = gt_class_alloc(sizeof *c_class);
  c_class->size = size;
  c_class->get_text_width = get_text_width;
  c_class->free = free;
  return c_class;
}

GtTextWidthCalculator* gt_text_width_calculator_create(
                                         const GtTextWidthCalculatorClass *twcc)
{
  GtTextWidthCalculator *twc;
  gt_assert(twcc && twcc->size);
  twc = gt_calloc(1, twcc->size);
  twc->c_class = twcc;
  twc->pvt = gt_calloc(1, sizeof (GtTextWidthCalculatorMembers));
  twc->pvt->lock = gt_rwlock_new();
  return twc;
}

GtTextWidthCalculator* gt_text_width_calculator_ref(GtTextWidthCalculator *twc)
{
  gt_assert(twc);
  gt_rwlock_wrlock(twc->pvt->lock);
  twc->pvt->reference_count++;
  gt_rwlock_unlock(twc->pvt->lock);
  return twc;
}

void gt_text_width_calculator_delete(GtTextWidthCalculator *twc)
{
  if (!twc) return;
  gt_rwlock_wrlock(twc->pvt->lock);
  if (twc->pvt->reference_count) {
    twc->pvt->reference_count--;
    gt_rwlock_unlock(twc->pvt->lock);
    return;
  }
  gt_assert(twc->c_class);
  if (twc->c_class->free)
    twc->c_class->free(twc);
  gt_rwlock_unlock(twc->pvt->lock);
  gt_rwlock_delete(twc->pvt->lock);
  gt_free(twc->pvt);
  gt_free(twc);
}

double gt_text_width_calculator_get_text_width(GtTextWidthCalculator *twc,
                                               const char* text,
                                               GtError *err)
{
  double width;
  gt_assert(twc && text);
  gt_rwlock_rdlock(twc->pvt->lock);
  gt_assert(twc->c_class);
  width = twc->c_class->get_text_width(twc, text, err);
  gt_rwlock_unlock(twc->pvt->lock);
  return width;
}

void* gt_text_width_calculator_cast(GT_UNUSED
                                    const GtTextWidthCalculatorClass *twcc,
                                    GtTextWidthCalculator *twc)
{
  gt_assert(twcc && twc&& twc->c_class == twcc);
  return twc;
}

void* gt_text_width_calculator_try_cast(const GtTextWidthCalculatorClass *twcc,
                                        GtTextWidthCalculator *twc)
{
  gt_assert(twcc && twc);
  if (twc->c_class == twcc)
    return twc;
  return NULL;
}
