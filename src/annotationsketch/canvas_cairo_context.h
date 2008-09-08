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

#ifndef CANVAS_CAIRO_CONTEXT_H
#define CANVAS_CAIRO_CONTEXT_H

#include <cairo.h>
#include "annotationsketch/canvas.h"

/* Implements the Canvas interface.
   This Canvas uses the GraphicsCairo class.  */
typedef struct GT_CanvasCairoContext GT_CanvasCairoContext;

const GT_CanvasClass* gt_canvas_cairo_context_class(void);
/* Create a new Canvas object tied to the cairo_t <context> and <width>
   using the style given in <style>. The optional <image_info> is filled when
   the created Canvas object is used to render a Diagram object. */
GT_Canvas* canvas_cairo_context_new(GT_Style *style, cairo_t *context,
                                    unsigned long width,
                                    GT_ImageInfo *image_info);
#endif
