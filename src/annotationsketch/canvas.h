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

#ifndef CANVAS_H
#define CANVAS_H

#include <stdio.h>
#include "core/error_api.h"
#include "annotationsketch/canvas_api.h"
#include "annotationsketch/block.h"
#include "annotationsketch/custom_track.h"
#include "annotationsketch/diagram.h"
#include "annotationsketch/drawing_range.h"
#include "annotationsketch/element.h"
#include "annotationsketch/layout.h"
#include "annotationsketch/line.h"
#include "annotationsketch/style.h"
#include "annotationsketch/track.h"

typedef struct GtCanvasClass GtCanvasClass;

void*           gt_canvas_cast(const GtCanvasClass *cc, GtCanvas *c);
void*           gt_canvas_try_cast(const GtCanvasClass *cc, GtCanvas *c);

GtUword   gt_canvas_calculate_height(GtCanvas*, GtDiagram*);
int             gt_canvas_draw_ruler(GtCanvas*, GtRange, GtError *err);

void            gt_format_ruler_label(char *txt, GtWord pos,
                                      const char *unitstr, size_t buflen);
GtStyle*        gt_canvas_get_style(GtCanvas *canvas);

/* Callback function for rendering. */
int             gt_canvas_visit_layout_pre(GtCanvas*, GtLayout*, GtError*);
/* Callback function for rendering. */
int             gt_canvas_visit_layout_post(GtCanvas*, GtLayout*, GtError*);
/* Callback function for rendering. */
int             gt_canvas_visit_track_pre(GtCanvas*, GtTrack*, GtError*);
/* Callback function for rendering. */
int             gt_canvas_visit_track_post(GtCanvas*, GtTrack*, GtError*);
/* Callback function for rendering. */
int             gt_canvas_visit_line_pre(GtCanvas*, GtLine*, GtError*);
/* Callback function for rendering. */
int             gt_canvas_visit_line_post(GtCanvas*, GtLine*, GtError*);
/* Callback function for rendering. */
int             gt_canvas_visit_block(GtCanvas*, GtBlock*, GtError*);
/* Callback function for rendering. */
int             gt_canvas_visit_element(GtCanvas*, GtElement*, GtError*);
/* Callback function for rendering. */
int             gt_canvas_visit_custom_track(GtCanvas*, GtCustomTrack*,
                                             GtError*);

#endif
