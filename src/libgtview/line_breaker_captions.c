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

#include "libgtcore/hashmap.h"
#include "libgtcore/ma.h"
#include "libgtcore/str.h"
#include "libgtview/line_breaker_captions.h"
#include "libgtview/line_breaker_rep.h"

struct LineBreakerCaptions {
  const LineBreaker parent_instance;
  Canvas *canvas;
  Hashmap *linepositions;
};

#define line_breaker_captions_cast(LB)\
        line_breaker_cast(line_breaker_captions_class(), LB)

static DrawingRange calculate_drawing_range(LineBreakerCaptions *lcb,
                                            Block* block)
{
  double textwidth = 0.0;
  DrawingRange drange;
  assert(block && lcb);
  drange = canvas_convert_coords(lcb->canvas, block_get_range(block));
  if (block_get_caption(block))
  {
    textwidth = canvas_get_text_width(lcb->canvas,
                                      str_get(block_get_caption(block)));
  if (textwidth > drawing_range_length(drange))
    drange.end = drange.start + textwidth;
  }
  return drange;
}

bool line_breaker_captions_is_line_occupied(LineBreaker* lb,
                                            Line *line,
                                            Block *block)
{
  DrawingRange dr;
  LineBreakerCaptions *lbcap;
  unsigned long *num;
  assert(lb && block && line);
  lbcap = line_breaker_captions_cast(lb);
  dr = calculate_drawing_range(lbcap, block);
  if (!(num = hashmap_get(lbcap->linepositions, line)))
    return false;
  else
    return (dr.start < *num);
}

void line_breaker_captions_register_block(LineBreaker *lb,
                                          Line *line,
                                          Block *block)
{
  DrawingRange dr;
  LineBreakerCaptions *lbcap;
  unsigned long *num;
  assert(lb && block && line);
  lbcap = line_breaker_captions_cast(lb);
  if (!(num = hashmap_get(lbcap->linepositions, line)))
  {
    num = ma_malloc(sizeof (unsigned long));
    hashmap_add(lbcap->linepositions, line, num);
  }
  dr = calculate_drawing_range(lbcap, block);
  *num = dr.end + 1;
}

void line_breaker_captions_delete(LineBreaker *lb)
{
  LineBreakerCaptions *lbcap;
  if (!lb) return;
  lbcap = line_breaker_captions_cast(lb);
  hashmap_delete(lbcap->linepositions);
}

const LineBreakerClass* line_breaker_captions_class(void)
{
  static const LineBreakerClass line_breaker_class =
    { sizeof (LineBreakerCaptions),
      line_breaker_captions_is_line_occupied,
      line_breaker_captions_register_block,
      line_breaker_captions_delete };
  return &line_breaker_class;
}

LineBreaker* line_breaker_captions_new(Canvas *canvas)
{
  assert(canvas);
  LineBreakerCaptions *lbcap;
  LineBreaker *lb;
  lb = line_breaker_create(line_breaker_captions_class());
  lbcap = line_breaker_captions_cast(lb);
  lbcap->canvas = canvas;
  lbcap->linepositions = hashmap_new(HASH_DIRECT, NULL, ma_free_func);
  return lb;
}
