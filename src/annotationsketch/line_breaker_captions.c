/*
  Copyright (c) 2008-2009 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2008-2009 Center for Bioinformatics, University of Hamburg

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

#include <math.h>
#include "core/class_alloc_lock.h"
#include "core/hashmap.h"
#include "core/ma.h"
#include "core/mathsupport.h"
#include "core/str.h"
#include "annotationsketch/coords.h"
#include "annotationsketch/default_formats.h"
#include "annotationsketch/drawing_range.h"
#include "annotationsketch/line_breaker_captions.h"
#include "annotationsketch/line_breaker_rep.h"

struct GtLineBreakerCaptions {
  const GtLineBreaker parent_instance;
  GtLayout *layout;
  unsigned long width;
  double margins;
  GtHashmap *linepositions;
};

#define gt_line_breaker_captions_cast(LB)\
        gt_line_breaker_cast(gt_line_breaker_captions_class(), LB)

int calculate_drawing_range(GtLineBreakerCaptions *lbc, GtDrawingRange *rng,
                            GtBlock* block, GtError *err)
{
  double textwidth = 0.0;
  GtDrawingRange drange;
  gt_assert(block && lbc);
  drange = gt_coords_calc_generic_range(gt_block_get_range(block),
                                        gt_layout_get_range(lbc->layout));
  drange.start *= lbc->width-2*lbc->margins;
  drange.end *= lbc->width-2*lbc->margins;
  if (gt_block_get_caption(block))
  {
    textwidth = gt_text_width_calculator_get_text_width(
                                      gt_layout_get_twc(lbc->layout),
                                      gt_str_get(gt_block_get_caption(block)),
                                      err);
    if (gt_double_smaller_double(textwidth, 0))
      return -1;
    if (gt_double_smaller_double(gt_drawing_range_length(drange), textwidth))
      drange.end = drange.start + textwidth;
  }
  rng->start = drange.start;
  rng->end = drange.end;
  return 0;
}

int gt_line_breaker_captions_is_line_occupied(GtLineBreaker* lb, bool *result,
                                              GtLine *line, GtBlock *block,
                                              GtError *err)
{
  GtDrawingRange dr;
  GtLineBreakerCaptions *lbcap;
  int had_err = 0;
  double *num;
  gt_assert(lb && block && line);
  lbcap = gt_line_breaker_captions_cast(lb);
  had_err = calculate_drawing_range(lbcap, &dr, block, err);
  if (!had_err) {
    if (!(num = gt_hashmap_get(lbcap->linepositions, line)))
      *result = false;
    else
      *result = (dr.start <= *num);
  }
  return had_err;
}

int gt_line_breaker_captions_register_block(GtLineBreaker *lb,
                                            GtLine *line,
                                            GtBlock *block,
                                            GtError *err)
{
  GtDrawingRange dr;
  GtLineBreakerCaptions *lbcap;
  int had_err = 0;
  double *num;
  gt_assert(lb && block && line);
  lbcap = gt_line_breaker_captions_cast(lb);
  if (!(num = gt_hashmap_get(lbcap->linepositions, line)))
  {
    num = gt_calloc(1, sizeof (double));
    gt_hashmap_add(lbcap->linepositions, line, num);
  }
  had_err = calculate_drawing_range(lbcap, &dr, block, err);
  if (!had_err)
    *num = floor(dr.end);
  return had_err;
}

void gt_line_breaker_captions_delete(GtLineBreaker *lb)
{
  GtLineBreakerCaptions *lbcap;
  if (!lb) return;
  lbcap = gt_line_breaker_captions_cast(lb);
  gt_hashmap_delete(lbcap->linepositions);
}

const GtLineBreakerClass* gt_line_breaker_captions_class(void)
{
  static const GtLineBreakerClass *lbc = NULL;
  gt_class_alloc_lock_enter();
  if (!lbc)
  {
    lbc = gt_line_breaker_class_new(sizeof (GtLineBreakerCaptions),
                                   gt_line_breaker_captions_is_line_occupied,
                                   gt_line_breaker_captions_register_block,
                                   gt_line_breaker_captions_delete);
  }
  gt_class_alloc_lock_leave();
  return lbc;
}

GtLineBreaker* gt_line_breaker_captions_new(GtLayout *layout,
                                            unsigned long width,
                                            GtStyle *style)
{
  GtLineBreakerCaptions *lbcap;
  GtLineBreaker *lb;
  gt_assert(layout);
  lb = gt_line_breaker_create(gt_line_breaker_captions_class());
  lbcap = gt_line_breaker_captions_cast(lb);
  lbcap->layout = layout;
  lbcap->width = width;
  if (!gt_style_get_num(style, "format", "margins", &lbcap->margins,
                        NULL, NULL)) {
    lbcap->margins = MARGINS_DEFAULT;
  }
  lbcap->linepositions = gt_hashmap_new(GT_HASH_DIRECT, NULL, gt_free_func);
  return lb;
}
