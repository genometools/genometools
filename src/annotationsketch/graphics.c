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

#include "core/assert_api.h"
#include "core/class_alloc.h"
#include "core/ma.h"
#include "core/thread_api.h"
#include "core/unused_api.h"
#include "annotationsketch/graphics_rep.h"

struct GtGraphicsMembers {
  unsigned int reference_count;
  GtRWLock *lock;
};

struct GtGraphicsClass {
  size_t size;
  GtGraphicsDrawTextFunc draw_text,
                         draw_text_clip,
                         draw_text_centered,
                         draw_text_right;
  GtGraphicsDrawColoredTextFunc draw_colored_text;
  GtGraphicsGetSingleExtentFunc get_text_height;
  GtGraphicsGetTextWidthFunc get_text_width;
  GtGraphicsSetColorFunc set_background_color;
  GtGraphicsSetFontFunc set_font;
  GtGraphicsGetSingleExtentFunc get_image_width,
                                get_image_height,
                                get_xmargins,
                                get_ymargins;
  GtGraphicsSetMarginsFunc set_margins;
  GtGraphicsDrawLineFunc draw_horizontal_line,
                         draw_vertical_line;
  GtGraphicsDrawLineToFunc draw_line;
  GtGraphicsDrawBoxFunc draw_box;
  GtGraphicsDrawSimpleFunc draw_dashes,
                           draw_caret;
  GtGraphicsDrawRectFunc draw_rectangle;
  GtGraphicsDrawArrowheadFunc draw_arrowhead;
  GtGraphicsDrawCurveDataFunc draw_curve;
  GtGraphicsSaveToFileFunc save_to_file;
  GtGraphicsSaveToStreamFunc save_to_stream;
  GtGraphicsFreeFunc free;
};

const GtGraphicsClass* gt_graphics_class_new(size_t size,
                                         GtGraphicsDrawTextFunc draw_text,
                                         GtGraphicsDrawTextFunc draw_text_clip,
                                         GtGraphicsDrawTextFunc
                                                     draw_text_centered,
                                         GtGraphicsDrawTextFunc draw_text_right,
                                         GtGraphicsDrawColoredTextFunc
                                                     draw_colored_text,
                                         GtGraphicsGetSingleExtentFunc
                                                     get_text_height,
                                         GtGraphicsGetTextWidthFunc
                                                     get_text_width,
                                         GtGraphicsSetColorFunc
                                                     set_background_color,
                                         GtGraphicsSetFontFunc
                                                     set_font,
                                         GtGraphicsGetSingleExtentFunc
                                                     get_image_width,
                                         GtGraphicsGetSingleExtentFunc
                                                     get_image_height,
                                         GtGraphicsGetSingleExtentFunc
                                                     get_xmargins,
                                         GtGraphicsGetSingleExtentFunc
                                                     get_ymargins,
                                         GtGraphicsSetMarginsFunc
                                                     set_margins,
                                         GtGraphicsDrawLineToFunc
                                                     draw_line,
                                         GtGraphicsDrawLineFunc
                                                     draw_horizontal_line,
                                         GtGraphicsDrawLineFunc
                                                     draw_vertical_line,
                                         GtGraphicsDrawBoxFunc draw_box,
                                         GtGraphicsDrawSimpleFunc draw_dashes,
                                         GtGraphicsDrawSimpleFunc draw_caret,
                                         GtGraphicsDrawRectFunc draw_rectangle,
                                         GtGraphicsDrawArrowheadFunc
                                                     draw_arrowhead,
                                         GtGraphicsDrawCurveDataFunc
                                                     draw_curve,
                                         GtGraphicsSaveToFileFunc save_to_file,
                                         GtGraphicsSaveToStreamFunc
                                                     save_to_stream,
                                         GtGraphicsFreeFunc free)
{
  GtGraphicsClass *c_class = gt_class_alloc(sizeof *c_class);
  c_class->size = size;
  c_class->draw_text = draw_text;
  c_class->draw_text_clip = draw_text_clip;
  c_class->draw_text_right = draw_text_right;
  c_class->draw_text_centered = draw_text_centered;
  c_class->draw_colored_text = draw_colored_text;
  c_class->get_text_height = get_text_height;
  c_class->get_text_width = get_text_width;
  c_class->set_background_color = set_background_color;
  c_class->set_font = set_font;
  c_class->get_image_width = get_image_width;
  c_class->get_image_height = get_image_height;
  c_class->get_xmargins = get_xmargins;
  c_class->get_ymargins = get_ymargins;
  c_class->set_margins = set_margins;
  c_class->draw_line = draw_line;
  c_class->draw_horizontal_line = draw_horizontal_line;
  c_class->draw_vertical_line = draw_vertical_line;
  c_class->draw_box = draw_box;
  c_class->draw_dashes = draw_dashes;
  c_class->draw_caret = draw_caret;
  c_class->draw_rectangle = draw_rectangle;
  c_class->draw_arrowhead = draw_arrowhead;
  c_class->draw_curve = draw_curve;
  c_class->save_to_file = save_to_file;
  c_class->save_to_stream = save_to_stream;
  c_class->free = free;
  return c_class;
}

GtGraphics* gt_graphics_ref(GtGraphics *g)
{
  gt_assert(g && g->pvt);
  gt_rwlock_wrlock(g->pvt->lock);
  g->pvt->reference_count++;
  gt_rwlock_unlock(g->pvt->lock);
  return g;
}

GtGraphics* gt_graphics_create(const GtGraphicsClass *gc)
{
  GtGraphics *g;
  gt_assert(gc && gc->size);
  g = gt_calloc(1, gc->size);
  g->c_class = gc;
  g->pvt = gt_calloc(1, sizeof (GtGraphicsMembers));
  g->pvt->lock = gt_rwlock_new();
  return g;
}

void gt_graphics_delete(GtGraphics *g)
{
  if (!g) return;
  gt_rwlock_wrlock(g->pvt->lock);
  if (g->pvt->reference_count)
  {
    g->pvt->reference_count--;
    gt_rwlock_unlock(g->pvt->lock);
    return;
  }
  gt_assert(g->c_class);
  if (g->c_class->free)
    g->c_class->free(g);
  gt_rwlock_unlock(g->pvt->lock);
  gt_rwlock_delete(g->pvt->lock);
  gt_free(g->pvt);
  gt_free(g);
}

void* gt_graphics_cast(GT_UNUSED const GtGraphicsClass *gc,
                       GtGraphics *g)
{
  gt_assert(gc && g);
  gt_assert(g->c_class == gc);
  return g;
}

void gt_graphics_draw_text(GtGraphics *g, double x, double y, const char* txt)
{
  gt_assert(g && g->c_class && txt);
  gt_rwlock_wrlock(g->pvt->lock);
  g->c_class->draw_text(g, x, y, txt);
  gt_rwlock_unlock(g->pvt->lock);
}

void gt_graphics_draw_text_clip(GtGraphics *g, double x, double y,
                                const char* txt)
{
  gt_assert(g && g->c_class && txt);
  gt_rwlock_wrlock(g->pvt->lock);
  g->c_class->draw_text_clip(g, x, y, txt);
  gt_rwlock_unlock(g->pvt->lock);
}

void gt_graphics_draw_text_centered(GtGraphics *g, double x, double y,
                                    const char *t)
{
  gt_assert(g && g->c_class && t);
  gt_rwlock_wrlock(g->pvt->lock);
  g->c_class->draw_text_centered(g, x, y, t);
  gt_rwlock_unlock(g->pvt->lock);
}

void gt_graphics_draw_text_right(GtGraphics *g, double x, double y,
                                 const char *txt)
{
  gt_assert(g && g->c_class && txt);
  gt_rwlock_wrlock(g->pvt->lock);
  g->c_class->draw_text_right(g, x, y, txt);
  gt_rwlock_unlock(g->pvt->lock);
}

void gt_graphics_draw_colored_text(GtGraphics *g, double x, double y,
                                   GtColor col, const char *txt)
{
  gt_assert(g && g->c_class && txt);
  gt_rwlock_wrlock(g->pvt->lock);
  g->c_class->draw_colored_text(g, x, y, col, txt);
  gt_rwlock_unlock(g->pvt->lock);
}

double gt_graphics_get_text_height(GtGraphics *g)
{
  double ret;
  gt_assert(g && g->c_class);
  gt_rwlock_rdlock(g->pvt->lock);
  ret = g->c_class->get_text_height(g);
  gt_rwlock_unlock(g->pvt->lock);
  return ret;
}

double gt_graphics_get_text_width(GtGraphics *g, const char *text)
{
  double ret;
  gt_assert(g && g->c_class && text);
  gt_rwlock_rdlock(g->pvt->lock);
  ret = g->c_class->get_text_width(g, text);
  gt_rwlock_unlock(g->pvt->lock);
  return ret;
}

int gt_graphics_set_background_color(GtGraphics *g, GtColor color)
{
  int ret;
  gt_assert(g && g->c_class);
  gt_rwlock_wrlock(g->pvt->lock);
  ret = g->c_class->set_background_color(g, color);
  gt_rwlock_unlock(g->pvt->lock);
  return ret;
}

void gt_graphics_set_font(GtGraphics *g, const char *family,
                       FontSlant slant, FontWeight weight, double size)
{
  gt_assert(g && g->c_class && family);
  gt_rwlock_wrlock(g->pvt->lock);
  g->c_class->set_font(g, family, slant, weight, size);
  gt_rwlock_unlock(g->pvt->lock);
}

double gt_graphics_get_image_height(GtGraphics *g)
{
  double ret;
  gt_assert(g && g->c_class);
  gt_rwlock_rdlock(g->pvt->lock);
  ret = g->c_class->get_image_height(g);
  gt_rwlock_unlock(g->pvt->lock);
  return ret;
}

double gt_graphics_get_image_width(GtGraphics *g)
{
  double ret;
  gt_assert(g && g->c_class);
  gt_rwlock_rdlock(g->pvt->lock);
  ret = g->c_class->get_image_width(g);
  gt_rwlock_unlock(g->pvt->lock);
  return ret;
}

double gt_graphics_get_xmargins(GtGraphics *g)
{
  double ret;
  gt_assert(g && g->c_class);
  gt_rwlock_rdlock(g->pvt->lock);
  ret = g->c_class->get_xmargins(g);
  gt_rwlock_unlock(g->pvt->lock);
  return ret;
}

double gt_graphics_get_ymargins(GtGraphics *g)
{
  double ret;
  gt_assert(g && g->c_class);
  gt_rwlock_rdlock(g->pvt->lock);
  ret = g->c_class->get_ymargins(g);
  gt_rwlock_unlock(g->pvt->lock);
  return ret;
}

void gt_graphics_set_margins(GtGraphics *g, double margin_x, double margin_y)
{
  gt_assert(g && g->c_class);
  gt_rwlock_wrlock(g->pvt->lock);
  g->c_class->set_margins(g, margin_x, margin_y);
  gt_rwlock_unlock(g->pvt->lock);
}

void gt_graphics_draw_line(GtGraphics *g, double x, double y,
                           double xto, double yto, GtColor color,
                           double stroke_width)
{
  gt_assert(g && g->c_class);
  gt_rwlock_wrlock(g->pvt->lock);
  g->c_class->draw_line(g, x, y, xto, yto, color, stroke_width);
  gt_rwlock_unlock(g->pvt->lock);
}

void gt_graphics_draw_horizontal_line(GtGraphics *g, double x, double y,
                                      GtColor color, double width,
                                      double stroke_width)
{
  gt_assert(g && g->c_class);
  gt_rwlock_wrlock(g->pvt->lock);
  g->c_class->draw_horizontal_line(g, x, y, color, width, stroke_width);
  gt_rwlock_unlock(g->pvt->lock);
}

void gt_graphics_draw_vertical_line(GtGraphics *g, double x, double y,
                                    GtColor color, double length,
                                    double stroke_width)
{
  gt_assert(g && g->c_class);
  gt_rwlock_wrlock(g->pvt->lock);
  g->c_class->draw_vertical_line(g, x, y, color, length, stroke_width);
  gt_rwlock_unlock(g->pvt->lock);
}

void gt_graphics_draw_box(GtGraphics *g, double x, double y, double width,
                          double height, GtColor fill_color,
                          ArrowStatus arrow_status, double arrow_width,
                          double stroke_width, GtColor stroke_color,
                          bool dashed)
{
  gt_assert(g && g->c_class);
  gt_rwlock_wrlock(g->pvt->lock);
  g->c_class->draw_box(g, x, y, width, height, fill_color, arrow_status,
                       arrow_width, stroke_width, stroke_color, dashed);
  gt_rwlock_unlock(g->pvt->lock);
}

void gt_graphics_draw_dashes(GtGraphics *g, double x, double y, double width,
                             double height, ArrowStatus arrow_status,
                             double arrow_width, double stroke_width,
                             GtColor stroke_color)
{
  gt_assert(g && g->c_class);
  gt_rwlock_wrlock(g->pvt->lock);
  g->c_class->draw_dashes(g, x, y, width, height, arrow_status, arrow_width,
                          stroke_width, stroke_color);
  gt_rwlock_unlock(g->pvt->lock);
}

void gt_graphics_draw_caret(GtGraphics *g, double x, double y, double width,
                            double height, ArrowStatus arrow_status,
                            double arrow_width,  double stroke_width,
                            GtColor stroke_color)
{
  gt_assert(g && g->c_class);
  gt_rwlock_wrlock(g->pvt->lock);
  g->c_class->draw_caret(g, x, y, width, height, arrow_status, arrow_width,
                         stroke_width, stroke_color);
  gt_rwlock_unlock(g->pvt->lock);
}

void gt_graphics_draw_rectangle(GtGraphics *g, double x, double y,
                                bool filled, GtColor fill_color, bool stroked,
                                GtColor stroke_color, double stroke_width,
                                double width, double height)
{
  gt_assert(g && g->c_class);
  gt_rwlock_wrlock(g->pvt->lock);
  g->c_class->draw_rectangle(g, x, y, filled, fill_color, stroked,
                             stroke_color, stroke_width, width, height);
  gt_rwlock_unlock(g->pvt->lock);
}

void gt_graphics_draw_arrowhead(GtGraphics *g, double x, double y,
                                GtColor col, ArrowStatus arrow_status)
{
  gt_assert(g && g->c_class);
  gt_rwlock_wrlock(g->pvt->lock);
  g->c_class->draw_arrowhead(g, x, y, col, arrow_status);
  gt_rwlock_unlock(g->pvt->lock);
}

void gt_graphics_draw_curve_data(GtGraphics *g, double x, double y,
                                 GtColor color,
                                 double data[], unsigned long ndata,
                                 GtRange valrange, unsigned long height)
{
  gt_assert(g && g->c_class);
  gt_rwlock_wrlock(g->pvt->lock);
  g->c_class->draw_curve(g, x, y, color, data, ndata, valrange, height);
  gt_rwlock_unlock(g->pvt->lock);
}

int gt_graphics_save_to_file(const GtGraphics *g, const char *filename,
                             GtError *err)
{
  int ret;
  gt_error_check(err);
  gt_assert(g && g->c_class);
  gt_rwlock_wrlock(g->pvt->lock);
  ret = g->c_class->save_to_file(g, filename, err);
  gt_rwlock_unlock(g->pvt->lock);
  return ret;
}

void gt_graphics_save_to_stream(const GtGraphics *g, GtStr *stream)
{
  gt_assert(g && g->c_class);
  gt_rwlock_wrlock(g->pvt->lock);
  g->c_class->save_to_stream(g, stream);
  gt_rwlock_unlock(g->pvt->lock);
}

/* the following are functions for the Ruby wrappers as workarounds for
   64bit issues */

typedef struct {
  double x,
         y;
  const char *txt;
} GraphicsDrawTextFuncParams;

void gt_graphics_draw_text_p(GtGraphics *g, GraphicsDrawTextFuncParams *params)
{
  gt_assert(g && g->c_class && params && params->txt);
  gt_rwlock_wrlock(g->pvt->lock);
  g->c_class->draw_text(g, params->x, params->y, params->txt);
  gt_rwlock_unlock(g->pvt->lock);
}

void gt_graphics_draw_text_clip_p(GtGraphics *g,
                                  GraphicsDrawTextFuncParams *params)
{
  gt_assert(g && g->c_class && params && params->txt);
  gt_rwlock_wrlock(g->pvt->lock);
  g->c_class->draw_text_clip(g, params->x, params->y, params->txt);
  gt_rwlock_unlock(g->pvt->lock);
}

void gt_graphics_draw_text_centered_p(GtGraphics *g,
                                      GraphicsDrawTextFuncParams *params)
{
  gt_assert(g && g->c_class && params && params->txt);
  gt_rwlock_wrlock(g->pvt->lock);
  g->c_class->draw_text_centered(g, params->x, params->y, params->txt);
  gt_rwlock_unlock(g->pvt->lock);
}

void gt_graphics_draw_text_right_p(GtGraphics *g,
                                   GraphicsDrawTextFuncParams *params)
{
  gt_assert(g && g->c_class && params && params->txt);
  gt_rwlock_wrlock(g->pvt->lock);
  g->c_class->draw_text_right(g, params->x, params->y, params->txt);
  gt_rwlock_unlock(g->pvt->lock);
}

typedef struct {
  double x,
         y;
  GtColor col;
  const char *txt;
} GraphicsDrawColoredTextFuncParams;

void gt_graphics_draw_colored_text_p(GtGraphics *g,
                                     GraphicsDrawColoredTextFuncParams *params)
{
  gt_assert(g && g->c_class && params && params->txt);
  gt_rwlock_wrlock(g->pvt->lock);
  g->c_class->draw_colored_text(g, params->x, params->y, params->col,
                                params->txt);
  gt_rwlock_unlock(g->pvt->lock);
}

typedef struct {
  double margin_x,
         margin_y;
} GraphicsSetMarginsFuncParams;

void gt_graphics_set_margins_p(GtGraphics *g,
                               GraphicsSetMarginsFuncParams *params)
{
  gt_assert(g && g->c_class && params);
  gt_rwlock_wrlock(g->pvt->lock);
  g->c_class->set_margins(g, params->margin_x, params->margin_y);
  gt_rwlock_unlock(g->pvt->lock);
}

typedef struct {
  double x,
         y,
         xto,
         yto;
  GtColor color;
  double stroke_width;
} GraphicsDrawLineToFuncParams;

void gt_graphics_draw_line_p(GtGraphics *g,
                             GraphicsDrawLineToFuncParams *params)
{
  gt_assert(g && g->c_class && params);
  gt_rwlock_wrlock(g->pvt->lock);
  g->c_class->draw_line(g, params->x, params->y, params->xto, params->yto,
                        params->color, params->stroke_width);
  gt_rwlock_unlock(g->pvt->lock);
}

typedef struct {
  double x,
         y;
  GtColor color;
  double len,
         stroke_width;
} GraphicsDrawLineFuncParams;

void gt_graphics_draw_horizontal_line_p(GtGraphics *g,
                                        GraphicsDrawLineFuncParams *params)
{
  gt_assert(g && g->c_class && params);
  gt_rwlock_wrlock(g->pvt->lock);
  g->c_class->draw_horizontal_line(g, params->x, params->y, params->color,
                                   params->len, params->stroke_width);
  gt_rwlock_unlock(g->pvt->lock);
}

void gt_graphics_draw_vertical_line_p(GtGraphics *g,
                                      GraphicsDrawLineFuncParams *params)
{
  gt_assert(g && g->c_class && params);
  gt_rwlock_wrlock(g->pvt->lock);
  g->c_class->draw_vertical_line(g, params->x, params->y, params->color,
                                 params->len, params->stroke_width);
  gt_rwlock_unlock(g->pvt->lock);
}

typedef struct {
  double x,
         y,
         width,
         height;
  GtColor fill_color;
  double arrow_width,
         stroke_width;
  GtColor stroke_color;
  int dashed;
  int arrow_status;
} GraphicsDrawBoxFuncParams;

void gt_graphics_draw_box_p(GtGraphics *g, GraphicsDrawBoxFuncParams *params)
{
  gt_assert(g && g->c_class && params);
  gt_rwlock_wrlock(g->pvt->lock);
  g->c_class->draw_box(g, params->x, params->y, params->width, params->height,
                       params->fill_color, params->arrow_status,
                       params->arrow_width, params->stroke_width,
                       params->stroke_color, params->dashed);
  gt_rwlock_unlock(g->pvt->lock);
}

typedef struct {
  double x,
         y,
         width,
         height;
  double arrow_width,
         stroke_width;
  GtColor stroke_color;
  ArrowStatus arrow_status;
} GraphicsDrawSimpleFuncParams;

void gt_graphics_draw_dashes_p(GtGraphics *g,
                               GraphicsDrawSimpleFuncParams *params)
{
  gt_assert(g && g->c_class && params);
  gt_rwlock_wrlock(g->pvt->lock);
  g->c_class->draw_dashes(g, params->x, params->y, params->width,
                          params->height, params->arrow_status,
                          params->arrow_width, params->stroke_width,
                          params->stroke_color);
  gt_rwlock_unlock(g->pvt->lock);
}

void gt_graphics_draw_caret_p(GtGraphics *g,
                              GraphicsDrawSimpleFuncParams *params)
{
  gt_assert(g && g->c_class && params);
  gt_rwlock_wrlock(g->pvt->lock);
  g->c_class->draw_caret(g, params->x, params->y, params->width,
                         params->height, params->arrow_status,
                         params->arrow_width, params->stroke_width,
                         params->stroke_color);
  gt_rwlock_unlock(g->pvt->lock);
}

typedef struct {
  double x,
         y;
  GtColor fill_color;
  GtColor outline_color;
  double outline_width,
         width,
         height;
  bool filled,
       outlined;
} GraphicsDrawRectFuncParams;

void gt_graphics_draw_rectangle_p(GtGraphics *g,
                                  GraphicsDrawRectFuncParams *params)
{
  gt_assert(g && g->c_class && params);
  gt_rwlock_wrlock(g->pvt->lock);
  g->c_class->draw_rectangle(g, params->x, params->y, params->filled,
                             params->fill_color, params->outlined,
                             params->outline_color, params->outline_width,
                             params->width, params->height);
  gt_rwlock_unlock(g->pvt->lock);
}

typedef struct {
  double x,
         y;
  GtColor col;
  ArrowStatus arrow_status;
} GraphicsDrawArrowheadFuncParams;

void gt_graphics_draw_arrowhead_p(GtGraphics *g,
                                  GraphicsDrawArrowheadFuncParams *params)
{
  gt_assert(g && g->c_class && params);
  gt_rwlock_wrlock(g->pvt->lock);
  g->c_class->draw_arrowhead(g, params->x, params->y, params->col,
                             params->arrow_status);
  gt_rwlock_unlock(g->pvt->lock);
}

typedef struct {
  double x,
         y;
  GtColor color;
} GraphicsDrawCurveFuncParams;

void gt_graphics_draw_curve_data_p(GtGraphics *g,
                                   GraphicsDrawCurveFuncParams *params,
                                   double *data,
                                   unsigned long rngstart,
                                   unsigned long rngend,
                                   unsigned long ndata,
                                   unsigned long height)
{
  GtRange rng;
  gt_assert(g && g->c_class && params);
  gt_rwlock_wrlock(g->pvt->lock);
  rng.start = rngstart;
  rng.end = rngend;
  g->c_class->draw_curve(g, params->x, params->y, params->color, data,
                         ndata, rng, height);
  gt_rwlock_unlock(g->pvt->lock);
}
