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

#ifndef CANVAS_CAIRO_H
#define CANVAS_CAIRO_H

#include "annotationsketch/block.h"
#include "annotationsketch/canvas.h"
#include "annotationsketch/element.h"
#include "annotationsketch/layout.h"
#include "annotationsketch/track.h"
#include "core/error.h"

int  gt_canvas_cairo_visit_layout_pre(GtCanvas*, GtLayout*, GtError*);
int  gt_canvas_cairo_visit_layout_post(GtCanvas*, GtLayout*, GtError*);
int  gt_canvas_cairo_visit_track_pre(GtCanvas*, GtTrack*, GtError*);
int  gt_canvas_cairo_visit_track_post(GtCanvas*, GtTrack*, GtError*);
int  gt_canvas_cairo_visit_line_pre(GtCanvas*, GtLine*, GtError*);
int  gt_canvas_cairo_visit_line_post(GtCanvas*, GtLine*, GtError*);
int  gt_canvas_cairo_visit_block(GtCanvas*, GtBlock*, GtError*);
int  gt_canvas_cairo_visit_element(GtCanvas*, GtElement*, GtError*);
int  gt_canvas_cairo_visit_custom_track(GtCanvas*, GtCustomTrack*, GtError*);
/* Renders a ruler with dynamic scale labeling and optional grid. */
int  gt_canvas_cairo_draw_ruler(GtCanvas*, GtRange, GtError*);

#endif
