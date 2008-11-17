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

#include <math.h>
#include <string.h>
#include "core/bittab.h"
#include "core/ensure.h"
#include "core/ma.h"
#include "core/minmax.h"
#include "core/range.h"
#include "core/unused_api.h"
#include "annotationsketch/canvas_members.h"
#include "annotationsketch/canvas_rep.h"
#include "annotationsketch/canvas_cairo.h"
#include "annotationsketch/canvas_cairo_context.h"
#include "annotationsketch/default_formats.h"
#include "annotationsketch/graphics_cairo.h"
#include "annotationsketch/style.h"

struct GtCanvasCairoContext {
  const GtCanvas parent_instance;
  cairo_t *context;
};

#define canvas_cairo_context_cast(C)\
        gt_canvas_cast(gt_canvas_cairo_context_class(), C)

const GtCanvasClass* gt_canvas_cairo_context_class(void)
{
  static const GtCanvasClass *canvas_class = NULL;
  if (!canvas_class) {
    canvas_class = gt_canvas_class_new(sizeof (GtCanvasCairoContext),
                                       gt_canvas_cairo_visit_layout_pre,
                                       gt_canvas_cairo_visit_layout_post,
                                       gt_canvas_cairo_visit_track_pre,
                                       gt_canvas_cairo_visit_track_post,
                                       gt_canvas_cairo_visit_line_pre,
                                       gt_canvas_cairo_visit_line_post,
                                       gt_canvas_cairo_visit_block,
                                       gt_canvas_cairo_visit_element,
                                       gt_canvas_cairo_draw_ruler,
                                       NULL);
  }
  return canvas_class;
}

GtCanvas* gt_canvas_cairo_context_new(GtStyle *sty, cairo_t *context,
                                      unsigned long width,
                                      unsigned long height,
                                      GtImageInfo *ii)
{
  GtCanvas *canvas;
  GtCanvasCairoContext *ccc;
  double margins = 10.0;
  gt_assert(sty && width > 0 && height > 0);
  canvas = gt_canvas_create(gt_canvas_cairo_context_class());
  canvas->pvt->y += HEADER_SPACE;
  canvas->pvt->g = gt_graphics_cairo_new_from_context(context,
                                                      width,
                                                      height);
  gt_style_get_num(sty, "format", "margins", &margins, NULL);
  gt_graphics_set_margins(canvas->pvt->g, margins, 0);
  canvas->pvt->margins = margins;
  if (ii)
    gt_image_info_set_height(ii, height);
  canvas->pvt->sty = sty;
  canvas->pvt->ii = ii;
  canvas->pvt->width = width;
  canvas->pvt->height = height;
  canvas->pvt->bt = NULL;
  /* 0.5 displacement to eliminate fuzzy horizontal lines */
  canvas->pvt->y = 0.5 + HEADER_SPACE;
  ccc = canvas_cairo_context_cast(canvas);
  ccc->context = context;
  return canvas;
}
