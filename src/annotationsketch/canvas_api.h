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

#ifndef CANVAS_API_H
#define CANVAS_API_H

#include "annotationsketch/graphics.h"
#include "annotationsketch/image_info.h"
#include "annotationsketch/style.h"

typedef struct GT_Canvas GT_Canvas;

/* Create a new GT_Canvas object with given <output_type> and <width> using the
   configuration given in <style>. The optional <image_info> is filled when
   the created GT_Canvas object is used to render a GT_Diagram object. */
GT_Canvas*    gt_canvas_new(GT_Style *style, GraphicsOutType output_type,
                            unsigned long width, GT_ImageInfo *image_info);
/* Returns the height of the given <canvas>. */
unsigned long gt_canvas_get_height(GT_Canvas *canvas);
/* Write rendered <canvas> to file with name <filename>. */
int           gt_canvas_to_file(GT_Canvas *canvas, const char *filename,
                                Error*);
/* Append rendered <canvas> to given <stream>. */
int           gt_canvas_to_stream(GT_Canvas *canvas, Str *stream);
/* Delete the given <canvas>. */
void          gt_canvas_delete(GT_Canvas *canvas);

#endif
