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

#ifndef CANVAS_CAIRO_FILE_API_H
#define CANVAS_CAIRO_FILE_API_H

#include "annotationsketch/canvas.h"
#include "annotationsketch/graphics.h"
#include "annotationsketch/image_info.h"

/* Implements the Canvas interface.
   This Canvas uses the GraphicsCairo class.  */
typedef struct GT_CanvasCairoFile GT_CanvasCairoFile;

/* Create a new Canvas object with given <output_type> and <width> using the
   configuration given in <style>. The optional <image_info> is filled when
   the created Canvas object is used to render a Diagram object. */
GT_Canvas* gt_canvas_cairo_file_new(GT_Style *style,
                                    GraphicsOutType output_type,
                                    unsigned long width,
                                    GT_ImageInfo *image_info);
/* Write rendered <canvas> to file with name <filename>. */
int     gt_canvas_cairo_file_to_file(GT_CanvasCairoFile *canvas,
                                     const char *filename, GT_Error *err);
/* Append rendered <canvas> to given <stream>. */
int     gt_canvas_cairo_file_to_stream(GT_CanvasCairoFile *canvas,
                                       GT_Str *stream);

#endif
