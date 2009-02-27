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
#include "core/hashmap.h"
#include "core/ma.h"
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

static GtDrawingRange calculate_drawing_range(GtLineBreakerCaptions *lbc,
                                              GtBlock* block)
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
                                      gt_str_get(gt_block_get_caption(block)));
  if (textwidth > gt_drawing_range_length(drange))
    drange.end = drange.start + textwidth;
  }
  return drange;
}

bool gt_line_breaker_captions_is_line_occupied(GtLineBreaker* lb, GtLine *line,
                                               GtBlock *block)
{
  GtDrawingRange dr;
  GtLineBreakerCaptions *lbcap;
  double *num;
  gt_assert(lb && block && line);
  lbcap = gt_line_breaker_captions_cast(lb);
  dr = calculate_drawing_range(lbcap, block);
  if (!(num = gt_hashmap_get(lbcap->linepositions, line)))
    return false;
  else
    return (dr.start <= *num);
}

void gt_line_breaker_captions_register_block(GtLineBreaker *lb,
                                             GtLine *line,
                                             GtBlock *block)
{
  GtDrawingRange dr;
  GtLineBreakerCaptions *lbcap;
  double *num;
  gt_assert(lb && block && line);
  lbcap = gt_line_breaker_captions_cast(lb);
  if (!(num = gt_hashmap_get(lbcap->linepositions, line)))
  {
    num = gt_calloc(1, sizeof (double));
    gt_hashmap_add(lbcap->linepositions, line, num);
  }
  dr = calculate_drawing_range(lbcap, block);
  *num = floor(dr.end);
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
  if (!lbc)
  {
    lbc = gt_line_breaker_class_new(sizeof (GtLineBreakerCaptions),
                                   gt_line_breaker_captions_is_line_occupied,
                                   gt_line_breaker_captions_register_block,
                                   gt_line_breaker_captions_delete);
  }
  return lbc;
}

GtLineBreaker* gt_line_breaker_captions_new(GtLayout *layout,
                                            unsigned long width,
                                            GtStyle *style)
{
  gt_assert(layout);
  GtLineBreakerCaptions *lbcap;
  GtLineBreaker *lb;
  lb = gt_line_breaker_create(gt_line_breaker_captions_class());
  lbcap = gt_line_breaker_captions_cast(lb);
  lbcap->layout = layout;
  lbcap->width = width;
  if (!gt_style_get_num(style, "format", "margins", &lbcap->margins, NULL))
    lbcap->margins = MARGINS_DEFAULT;
  lbcap->linepositions = gt_hashmap_new(HASH_DIRECT, NULL, gt_free_func);
  return lb;
}
