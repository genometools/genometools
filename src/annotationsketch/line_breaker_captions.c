/*
  Copyright (c) 2008 Sascha Steinbiss <ssteinbiss@stud.zbh.uni-hamburg.de>
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

#include "core/hashmap.h"
#include "core/ma.h"
#include "core/str.h"
#include "annotationsketch/line_breaker_captions.h"
#include "annotationsketch/line_breaker_rep.h"

struct GtLineBreakerCaptions {
  const GtLineBreaker parent_instance;
  GtCanvas *canvas;
  Hashmap *linepositions;
};

#define gt_line_breaker_captions_cast(LB)\
        gt_line_breaker_cast(gt_line_breaker_captions_class(), LB)

static GtDrawingRange calculate_drawing_range(GtLineBreakerCaptions *lcb,
                                            GtBlock* block)
{
  double textwidth = 0.0;
  GtDrawingRange drange;
  assert(block && lcb);
  drange = gt_canvas_convert_coords(lcb->canvas, gt_block_get_range(block));
  if (gt_block_get_caption(block))
  {
    textwidth = gt_canvas_get_text_width(lcb->canvas,
                                      gt_str_get(gt_block_get_caption(block)));
  if (textwidth > gt_drawing_range_length(drange))
    drange.end = drange.start + textwidth;
  }
  return drange;
}

bool gt_line_breaker_captions_is_gt_line_occupied(GtLineBreaker* lb,
                                            GtLine *line,
                                            GtBlock *block)
{
  GtDrawingRange dr;
  GtLineBreakerCaptions *lbcap;
  unsigned long *num;
  assert(lb && block && line);
  lbcap = gt_line_breaker_captions_cast(lb);
  dr = calculate_drawing_range(lbcap, block);
  if (!(num = hashmap_get(lbcap->linepositions, line)))
    return false;
  else
    return (dr.start < *num);
}

void gt_line_breaker_captions_register_block(GtLineBreaker *lb,
                                          GtLine *line,
                                          GtBlock *block)
{
  GtDrawingRange dr;
  GtLineBreakerCaptions *lbcap;
  unsigned long *num;
  assert(lb && block && line);
  lbcap = gt_line_breaker_captions_cast(lb);
  if (!(num = hashmap_get(lbcap->linepositions, line)))
  {
    num = gt_malloc(sizeof (unsigned long));
    hashmap_add(lbcap->linepositions, line, num);
  }
  dr = calculate_drawing_range(lbcap, block);
  *num = dr.end + 1;
}

void gt_line_breaker_captions_delete(GtLineBreaker *lb)
{
  GtLineBreakerCaptions *lbcap;
  if (!lb) return;
  lbcap = gt_line_breaker_captions_cast(lb);
  hashmap_delete(lbcap->linepositions);
}

const GtLineBreakerClass* gt_line_breaker_captions_class(void)
{
  static const GtLineBreakerClass gt_line_breaker_class =
    { sizeof (GtLineBreakerCaptions),
      gt_line_breaker_captions_is_gt_line_occupied,
      gt_line_breaker_captions_register_block,
      gt_line_breaker_captions_delete };
  return &gt_line_breaker_class;
}

GtLineBreaker* gt_line_breaker_captions_new(GtCanvas *canvas)
{
  assert(canvas);
  GtLineBreakerCaptions *lbcap;
  GtLineBreaker *lb;
  lb = gt_line_breaker_create(gt_line_breaker_captions_class());
  lbcap = gt_line_breaker_captions_cast(lb);
  lbcap->canvas = canvas;
  lbcap->linepositions = hashmap_new(HASH_DIRECT, NULL, gt_free_func);
  return lb;
}
