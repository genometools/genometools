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
#include "core/unused.h"
#include "annotationsketch/canvas_rep.h"
#include "annotationsketch/canvas_cairo_context.h"
#include "annotationsketch/graphics_cairo.h"

#define MARGINS_DEFAULT           10
#define HEADER_SPACE              70

struct CanvasCairoContext {
  const Canvas parent_instance;
  cairo_t *context;
};

#define canvas_cairo_context_cast(C)\
        canvas_cast(canvas_cairo_context_class(), C)

int canvas_cairo_context_visit_diagram_pre(Canvas *canvas, Diagram *dia)
{
  double margins;

  assert(canvas && dia);

  if (style_get_num(canvas->sty, "format", "margins", &margins, NULL))
    canvas->margins = margins;
  else
    canvas->margins = MARGINS_DEFAULT;

  if (!style_get_bool(canvas->sty, "format", "show_track_captions",
                       &canvas->show_track_captions, NULL))
    canvas->show_track_captions = true;

  canvas->viewrange = diagram_get_range(dia);
  if (canvas->g)
  {
    graphics_delete(canvas->g);
    canvas->g = NULL;
  }
  canvas->g = graphics_cairo_new(GRAPHICS_PNG, canvas->width, 1);

  /* calculate scaling factor */
  canvas->factor = ((double) canvas->width
                     -(2*canvas->margins))
                    / range_length(canvas->viewrange);
  return 0;
}

int canvas_cairo_context_visit_diagram_post(Canvas *canvas, Diagram *dia)
{
  int had_err = 0;

  assert(canvas && dia);

  /* set initial image-specific values */
  canvas->y += HEADER_SPACE;
  canvas->height = calculate_height(canvas, dia);
  if (canvas->ii)
    image_info_set_height(canvas->ii, canvas->height);
  if (canvas->g)
  {
    graphics_delete(canvas->g);
    canvas->g = NULL;
  }
  canvas->g = graphics_cairo_new_from_context(
                                 ((CanvasCairoContext*) canvas)->context,
                                 canvas->width, canvas->height);
  graphics_set_margins(canvas->g, canvas->margins, 0);

  /* Add ruler/scale to the image */
  draw_ruler(canvas);

  return had_err;
}

const CanvasClass* canvas_cairo_context_class(void)
{
  static const CanvasClass canvas_class =
    { sizeof (CanvasCairoContext),
      canvas_cairo_context_visit_diagram_pre,
      canvas_cairo_context_visit_diagram_post,
      NULL };
  return &canvas_class;
}

Canvas* canvas_cairo_context_new(Style *sty, cairo_t *context,
                                 unsigned long width, ImageInfo *ii)
{
  Canvas *canvas;
  CanvasCairoContext *ccc;
  assert(sty && width > 0);
  canvas = canvas_create(canvas_cairo_context_class());
  canvas->sty = sty;
  canvas->ii = ii;
  canvas->width = width;
  canvas->bt = NULL;
  canvas->y = 0.5; /* 0.5 displacement to eliminate fuzzy horizontal lines */
  ccc = canvas_cairo_context_cast(canvas);
  ccc->context = context;
  return canvas;
}
