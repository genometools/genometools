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

#include <assert.h>
#include "core/ma.h"
#include "core/unused_api.h"
#include "annotationsketch/graphics_rep.h"

GtGraphics* gt_graphics_create(const GtGraphicsClass *gc)
{
  GtGraphics *g;
  assert(gc && gc->size);
  g = gt_calloc(1, gc->size);
  g->c_class = gc;
  return g;
}

void gt_graphics_delete(GtGraphics *g)
{
  if (!g) return;
  assert(g->c_class);
  if (g->c_class->free)
    g->c_class->free(g);
  gt_free(g);
}

void* gt_graphics_cast(GT_UNUSED const GtGraphicsClass *gc,
                    GtGraphics *g)
{
  assert(gc && g && g->c_class == gc);
  return g;
}

void gt_graphics_draw_text(GtGraphics *g, double x, double y, const char* txt)
{
  assert(g && g->c_class && txt);
  g->c_class->draw_text(g, x, y, txt);
}

void gt_graphics_draw_text_centered(GtGraphics *g, double x, double y,
                                    const char *t)
{
  assert(g && g->c_class && t);
  g->c_class->draw_text_centered(g, x, y, t);
}

void gt_graphics_draw_text_right(GtGraphics *g, double x, double y,
                                 const char *txt)
{
  assert(g && g->c_class && txt);
  g->c_class->draw_text_right(g, x, y, txt);
}

void gt_graphics_draw_colored_text(GtGraphics *g, double x, double y,
                                   GtColor col, const char *txt)
{
  assert(g && g->c_class && txt);
  g->c_class->draw_colored_text(g, x, y, col, txt);
}

double gt_graphics_get_text_height(GtGraphics *g)
{
  assert(g && g->c_class);
  return g->c_class->get_text_height(g);
}

double gt_graphics_get_text_width(GtGraphics *g, const char *txt)
{
  assert(g && g->c_class && txt);
  return g->c_class->get_text_width(g, txt);
}

void gt_graphics_set_font(GtGraphics *g, const char *family,
                       FontSlant slant, FontWeight weight)
{
  assert(g && g->c_class && family);
  g->c_class->set_font(g, family, slant, weight);
}

double gt_graphics_get_image_height(GtGraphics *g)
{
  assert(g && g->c_class);
  return g->c_class->get_image_height(g);
}

double gt_graphics_get_image_width(GtGraphics *g)
{
  assert(g && g->c_class);
  return g->c_class->get_image_width(g);
}

void gt_graphics_set_margins(GtGraphics *g, double margin_x, double margin_y)
{
  assert(g && g->c_class);
  g->c_class->set_margins(g, margin_x, margin_y);
}

void gt_graphics_draw_horizontal_line(GtGraphics *g, double x, double y,
                                   double width)
{
  assert(g && g->c_class);
  g->c_class->draw_horizontal_line(g, x, y, width);
}

void gt_graphics_draw_vertical_line(GtGraphics *g, double x, double y,
                                 GtColor color, double length)
{
  assert(g && g->c_class);
  g->c_class->draw_vertical_line(g, x, y, color, length);
}

void gt_graphics_draw_box(GtGraphics *g, double x, double y, double width,
                       double height, GtColor fill_color,
                       ArrowStatus arrow_status, double arrow_width,
                       double stroke_width, GtColor stroke_color,
                       bool dashed)
{
  assert(g && g->c_class);
  g->c_class->draw_box(g, x, y, width, height, fill_color, arrow_status,
                       arrow_width, stroke_width, stroke_color, dashed);
}

void gt_graphics_draw_dashes(GtGraphics *g, double x, double y, double width,
                          double height, ArrowStatus arrow_status,
                          double arrow_width, double stroke_width,
                          GtColor stroke_color)
{
  assert(g && g->c_class);
  g->c_class->draw_dashes(g, x, y, width, height, arrow_status, arrow_width,
                          stroke_width, stroke_color);
}

void gt_graphics_draw_caret(GtGraphics *g, double x, double y, double width,
                         double height, ArrowStatus arrow_status,
                         double arrow_width,  double stroke_width,
                         GtColor stroke_color)
{
  assert(g && g->c_class);
  g->c_class->draw_caret(g, x, y, width, height, arrow_status, arrow_width,
                         stroke_width, stroke_color);
}

void gt_graphics_draw_rectangle(GtGraphics *g, double x, double y,
                             bool filled, GtColor fill_color, bool outlined,
                             GtColor outline_color, double outline_width,
                             double width)
{
  assert(g && g->c_class);
  g->c_class->draw_rectangle(g, x, y, filled, fill_color, outlined,
                             outline_color, outline_width, width);
}

void gt_graphics_draw_arrowhead(GtGraphics *g, double x, double y,
                                GtColor col, ArrowStatus arrow_status)
{
  assert(g && g->c_class);
  g->c_class->draw_arrowhead(g, x, y, col, arrow_status);
}

int gt_graphics_save_to_file(const GtGraphics *g, const char *filename,
                             GtError *err)
{
  gt_error_check(err);
  assert(g && g->c_class);
  return g->c_class->save_to_file(g, filename, err);
}

void gt_graphics_save_to_stream(const GtGraphics *g, GtStr *stream)
{
  assert(g && g->c_class);
  return g->c_class->save_to_stream(g, stream);
}
