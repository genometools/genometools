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

#include <math.h>
#include <string.h>
#include "core/bittab.h"
#include "core/ensure.h"
#include "core/ma.h"
#include "core/minmax.h"
#include "core/range.h"
#include "core/unused_api.h"
#include "annotationsketch/canvas_rep.h"
#include "annotationsketch/canvas_cairo_context.h"
#include "annotationsketch/graphics_cairo.h"
#include "annotationsketch/style.h"

#define MARGINS_DEFAULT           10
#define HEADER_SPACE              70

struct GT_CanvasCairoContext {
  const GT_Canvas parent_instance;
  cairo_t *context;
};

#define canvas_cairo_context_cast(C)\
        gt_canvas_cast(canvas_cairo_context_class(), C)

int canvas_cairo_context_visit_diagram_pre(GT_Canvas *canvas, GT_Diagram *dia)
{
  double margins;

  assert(canvas && dia);

  if (gt_style_get_num(canvas->sty, "format", "margins", &margins, NULL))
    canvas->margins = margins;
  else
    canvas->margins = MARGINS_DEFAULT;

  if (!gt_style_get_bool(canvas->sty, "format", "show_track_captions",
                         &canvas->show_track_captions, NULL))
    canvas->show_track_captions = true;

  canvas->viewrange = gt_diagram_get_range(dia);
  if (canvas->g)
  {
    gt_graphics_delete(canvas->g);
    canvas->g = NULL;
  }
  canvas->g = gt_graphics_cairo_new(GT_GRAPHICS_PNG, canvas->width, 1);

  /* calculate scaling factor */
  canvas->factor = ((double) canvas->width
                     -(2*canvas->margins))
                    / gt_range_length(canvas->viewrange);
  return 0;
}

int canvas_cairo_context_visit_diagram_post(GT_Canvas *canvas, GT_Diagram *dia)
{
  int had_err = 0;

  assert(canvas && dia);

  /* set initial image-specific values */
  canvas->y += HEADER_SPACE;
  canvas->height = gt_canvas_calculate_height(canvas, dia);
  if (canvas->ii)
    gt_image_info_set_height(canvas->ii, canvas->height);
  if (canvas->g)
  {
    gt_graphics_delete(canvas->g);
    canvas->g = NULL;
  }
  canvas->g = gt_graphics_cairo_new_from_context(
                                 ((GT_CanvasCairoContext*) canvas)->context,
                                 canvas->width, canvas->height);
  gt_graphics_set_margins(canvas->g, canvas->margins, 0);

  /* Add ruler/scale to the image */
  gt_canvas_draw_ruler(canvas);

  return had_err;
}

const GT_CanvasClass* canvas_cairo_context_class(void)
{
  static const GT_CanvasClass canvas_class =
    { sizeof (GT_CanvasCairoContext),
      canvas_cairo_context_visit_diagram_pre,
      canvas_cairo_context_visit_diagram_post,
      NULL };
  return &canvas_class;
}

GT_Canvas* canvas_cairo_context_new(GT_Style *sty, cairo_t *context,
                                    unsigned long width, GT_ImageInfo *ii)
{
  GT_Canvas *canvas;
  GT_CanvasCairoContext *ccc;
  assert(sty && width > 0);
  canvas = gt_canvas_create(canvas_cairo_context_class());
  canvas->sty = sty;
  canvas->ii = ii;
  canvas->width = width;
  canvas->bt = NULL;
  canvas->y = 0.5; /* 0.5 displacement to eliminate fuzzy horizontal lines */
  ccc = canvas_cairo_context_cast(canvas);
  ccc->context = context;
  return canvas;
}
