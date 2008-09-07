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
#include "annotationsketch/canvas.h"
#include "annotationsketch/canvas_cairo_file.h"
#include "annotationsketch/canvas_rep.h"
#include "annotationsketch/graphics_cairo.h"

#define MARGINS_DEFAULT           10
#define HEADER_SPACE              70

struct CanvasCairoFile {
  const Canvas parent_instance;
  GraphicsOutType type;
};

#define canvas_cairo_file_cast(C)\
        canvas_cast(canvas_cairo_file_class(), C)

int canvas_cairo_file_visit_diagram_pre(Canvas *canvas, Diagram *dia)
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
  canvas->g = graphics_cairo_new(((CanvasCairoFile*) canvas)->type,
                                  canvas->width, 1);

  /* calculate scaling factor */
  canvas->factor = ((double) canvas->width
                     -(2*canvas->margins))
                    / range_length(canvas->viewrange);
  return 0;
}

int canvas_cairo_file_visit_diagram_post(Canvas *canvas, Diagram *dia)
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
  canvas->g = graphics_cairo_new(((CanvasCairoFile*) canvas)->type,
                                 canvas->width, canvas->height);
  graphics_set_margins(canvas->g, canvas->margins, 0);

  /* Add ruler/scale to the image */
  draw_ruler(canvas);

  return had_err;
}

int canvas_cairo_file_to_file(CanvasCairoFile *canvas, const char *filename,
                              Error *err)
{
  int had_err = 0;
  Canvas *c = (Canvas*) canvas;
  assert(canvas && filename && err);
  /* write out result file */
  if (c->g)
    had_err = graphics_save_to_file(c->g, filename, err);
  else
  {
    error_set(err, "No graphics has been created yet!");
    had_err = -1;
  }

  return had_err;
}

int canvas_cairo_file_to_stream(CanvasCairoFile *canvas, Str *stream)
{
  int had_err = 0;
  Canvas *c = (Canvas*) canvas;
  assert(canvas && stream);

  /* write out result file */
  if (c->g)
    graphics_save_to_stream(c->g, stream);

  return had_err;
}

const CanvasClass* canvas_cairo_file_class(void)
{
  static const CanvasClass canvas_class =
    { sizeof (CanvasCairoFile),
      canvas_cairo_file_visit_diagram_pre,
      canvas_cairo_file_visit_diagram_post,
      NULL };
  return &canvas_class;
}

Canvas* canvas_cairo_file_new(Style *sty, GraphicsOutType type,
                              unsigned long width, ImageInfo *ii)
{
  Canvas *canvas;
  CanvasCairoFile *ccf;
  assert(sty && width > 0);
  canvas = canvas_create(canvas_cairo_file_class());
  canvas->sty = sty;
  canvas->ii = ii;
  canvas->width = width;
  canvas->bt = NULL;
  canvas->y = 0.5; /* 0.5 displacement to eliminate fuzzy horizontal lines */
  ccf = canvas_cairo_file_cast(canvas);
  ccf->type = type;
  return canvas;
}
