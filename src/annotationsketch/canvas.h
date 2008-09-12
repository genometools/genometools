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

typedef struct GtCanvasClass GtCanvasClass;

#include "annotationsketch/canvas_api.h"
#include "annotationsketch/block.h"
#include "annotationsketch/diagram.h"
#include "annotationsketch/drawing_range.h"
#include "annotationsketch/element.h"
#include "annotationsketch/track.h"

unsigned long gt_canvas_calculate_height(GtCanvas *canvas, GT_Diagram *dia);
void          gt_canvas_draw_ruler(GtCanvas *canvas);

/* Returns a pixel-based range for a nucleotide-based range
   using the scaling factor defined for the given <canvas> */
GT_DrawingRange  gt_canvas_convert_coords(GtCanvas *canvas, GT_Range);
/* Returns rendered width in pixels of the given text. */
double        gt_canvas_get_text_width(GtCanvas*, const char *text);
/* Callback function for Diagram rendering. */
int           gt_canvas_visit_diagram_pre(GtCanvas*, GT_Diagram*);
/* Callback function for Diagram rendering. */
int           gt_canvas_visit_diagram_post(GtCanvas*, GT_Diagram*);
/* Callback function for Diagram rendering. */
int           gt_canvas_visit_track_pre(GtCanvas*, GtTrack*);
/* Callback function for Diagram rendering. */
int           gt_canvas_visit_track_post(GtCanvas*, GtTrack*);
/* Callback function for Diagram rendering. */
int           gt_canvas_visit_line_pre(GtCanvas*, GtLine*);
/* Callback function for Diagram rendering. */
int           gt_canvas_visit_line_post(GtCanvas*, GtLine*);
/* Callback function for Diagram rendering. */
int           gt_canvas_visit_block(GtCanvas*, GtBlock*);
/* Callback function for Diagram rendering. */
int           gt_canvas_visit_element(GtCanvas*, GtElement*);

#endif
