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

typedef struct Canvas Canvas;

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

/* Create a new Canvas object with given <output_type> and <width> using the
   configuration given in <style>. The optional <image_info> is filled when
   the created Canvas object is used to render a Diagram object. */
Canvas*       canvas_new(Style *style, GraphicsOutType output_type,
                         unsigned long width, ImageInfo *image_info);
/* Returns a pixel-based range for a nucleotide-based range
   using the scaling factor defined for the given <canvas> */
DrawingRange  canvas_convert_coords(Canvas* canvas, Range);
/* Returns the height of the given <canvas>. */
unsigned long canvas_get_height(Canvas *canvas);
/* Returns rendered width in pixels of the given text. */
double        canvas_get_text_width(Canvas*, const char *text);
/* Callback function for Diagram rendering. */
int           canvas_visit_diagram_pre(Canvas*, Diagram*);
/* Callback function for Diagram rendering. */
int           canvas_visit_diagram_post(Canvas*, Diagram*);
/* Callback function for Diagram rendering. */
int           canvas_visit_track_pre(Canvas*, Track*);
/* Callback function for Diagram rendering. */
int           canvas_visit_track_post(Canvas*, Track*);
/* Callback function for Diagram rendering. */
int           canvas_visit_line_pre(Canvas*, Line*);
/* Callback function for Diagram rendering. */
int           canvas_visit_line_post(Canvas*, Line*);
/* Callback function for Diagram rendering. */
int           canvas_visit_block(Canvas*, GT_Block*);
/* Callback function for Diagram rendering. */
int           canvas_visit_element(Canvas*, Element*);
/* Write rendered <canvas> to file with name <filename>. */
int           canvas_to_file(Canvas *canvas, const char *filename, Error*);
/* Append rendered <canvas> to given <stream>. */
int           canvas_to_stream(Canvas *canvas, Str *stream);
/* Delete the given <canvas>. */
void          canvas_delete(Canvas *canvas);

#endif
