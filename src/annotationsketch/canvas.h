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

#ifndef CANVAS_H
#define CANVAS_H

typedef struct GT_Canvas GT_Canvas;

#include "annotationsketch/block.h"
#include "annotationsketch/diagram.h"
#include "annotationsketch/drawing_range.h"
#include "annotationsketch/element.h"
#include "annotationsketch/graphics.h"
#include "annotationsketch/line.h"
#include "annotationsketch/image_info.h"
#include "annotationsketch/style.h"
#include "annotationsketch/track.h"
#include "core/range.h"

/* Create a new GT_Canvas object with given <output_type> and <width> using the
   configuration given in <style>. The optional <image_info> is filled when
   the created GT_Canvas object is used to render a GT_Diagram object. */
GT_Canvas*    gt_canvas_new(GT_Style *style, GraphicsOutType output_type,
                            unsigned long width, ImageInfo *image_info);
/* Returns a pixel-based range for a nucleotide-based range
   using the scaling factor defined for the given <canvas> */
DrawingRange  gt_canvas_convert_coords(GT_Canvas* canvas, Range);
/* Returns the height of the given <canvas>. */
unsigned long gt_canvas_get_height(GT_Canvas *canvas);
/* Returns rendered width in pixels of the given text. */
double        gt_canvas_get_text_width(GT_Canvas*, const char *text);
/* Callback function for GT_Diagram rendering. */
int           gt_canvas_visit_gt_diagram_pre(GT_Canvas*, GT_Diagram*);
/* Callback function for GT_Diagram rendering. */
int           gt_canvas_visit_gt_diagram_post(GT_Canvas*, GT_Diagram*);
/* Callback function for GT_Diagram rendering. */
int           gt_canvas_visit_track_pre(GT_Canvas*, Track*);
/* Callback function for GT_Diagram rendering. */
int           gt_canvas_visit_track_post(GT_Canvas*, Track*);
/* Callback function for GT_Diagram rendering. */
int           gt_canvas_visit_line_pre(GT_Canvas*, Line*);
/* Callback function for GT_Diagram rendering. */
int           gt_canvas_visit_line_post(GT_Canvas*, Line*);
/* Callback function for GT_Diagram rendering. */
int           gt_canvas_visit_block(GT_Canvas*, GT_Block*);
/* Callback function for GT_Diagram rendering. */
int           gt_canvas_visit_element(GT_Canvas*, Element*);
/* Write rendered <canvas> to file with name <filename>. */
int           gt_canvas_to_file(GT_Canvas *canvas, const char *filename,
                                Error*);
/* Append rendered <canvas> to given <stream>. */
int           gt_canvas_to_stream(GT_Canvas *canvas, Str *stream);
/* Delete the given <canvas>. */
void          gt_canvas_delete(GT_Canvas *canvas);

#endif
