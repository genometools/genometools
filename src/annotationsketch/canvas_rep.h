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

#ifndef CANVAS_REP_H
#define CANVAS_REP_H

#include <stdio.h>
#include "core/bittab.h"
#include "core/range.h"
#include "annotationsketch/canvas.h"
#include "annotationsketch/custom_track.h"
#include "annotationsketch/graphics.h"
#include "annotationsketch/layout.h"

typedef int  (*GtCanvasVisitLayoutFunc)(GtCanvas*, GtLayout*, GtError*);
typedef int  (*GtCanvasVisitTrackFunc)(GtCanvas*, GtTrack*, GtError*);
typedef int  (*GtCanvasVisitLineFunc)(GtCanvas*, GtLine*, GtError*);
typedef int  (*GtCanvasVisitBlockFunc)(GtCanvas*, GtBlock*, GtError*);
typedef int  (*GtCanvasVisitElementFunc)(GtCanvas*, GtElement*, GtError*);
typedef int  (*GtCanvasVisitCustomTrackFunc)(GtCanvas*, GtCustomTrack*,
                                             GtError*);
typedef int  (*GtCanvasDrawRulerFunc)(GtCanvas*, GtRange, GtError*);
typedef void (*GtCanvasFreeFunc)(GtCanvas*);

typedef struct GtCanvasMembers GtCanvasMembers;

struct GtCanvas {
  const GtCanvasClass *c_class;
  GtCanvasMembers *pvt;
};
const GtCanvasClass* gt_canvas_class_new(size_t size,
                                         GtCanvasVisitLayoutFunc la_visit_pre,
                                         GtCanvasVisitLayoutFunc la_visit_post,
                                         GtCanvasVisitTrackFunc t_visit_pre,
                                         GtCanvasVisitTrackFunc t_visit_post,
                                         GtCanvasVisitLineFunc l_visit_pre,
                                         GtCanvasVisitLineFunc l_visit_post,
                                         GtCanvasVisitBlockFunc b_visit,
                                         GtCanvasVisitElementFunc e_visit,
                                         GtCanvasVisitCustomTrackFunc ct_visit,
                                         GtCanvasDrawRulerFunc draw_ruler_func,
                                         GtCanvasFreeFunc free);
GtCanvas* gt_canvas_create(const GtCanvasClass*);
void*     gt_canvas_cast(const GtCanvasClass *cc, GtCanvas *c);
void*     gt_canvas_try_cast(const GtCanvasClass *cc, GtCanvas *c);

#endif
