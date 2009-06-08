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

#include "annotationsketch/custom_track_rep.h"
#include "annotationsketch/custom_track.h"
#include "core/assert_api.h"
#include "core/class_alloc.h"
#include "core/ma.h"
#include "core/unused_api.h"

struct GtCustomTrackMembers {
  unsigned int reference_count;
};

struct GtCustomTrackClass {
  size_t size;
  GtCustomTrackRenderFunc render;
  GtCustomTrackGetHeightFunc get_height;
  GtCustomTrackGetTitleFunc get_title;
  GtCustomTrackFreeFunc free;
};

const GtCustomTrackClass* gt_custom_track_class_new(size_t size,
                                          GtCustomTrackRenderFunc render,
                                          GtCustomTrackGetHeightFunc get_height,
                                          GtCustomTrackGetTitleFunc get_title,
                                          GtCustomTrackFreeFunc free)
{
  GtCustomTrackClass *c_class = gt_class_alloc(sizeof *c_class);
  c_class->size = size;
  c_class->render = render;
  c_class->get_height = get_height;
  c_class->get_title = get_title;
  c_class->free = free;
  return c_class;
}

GtCustomTrack* gt_custom_track_create(const GtCustomTrackClass *ctc)
{
  GtCustomTrack *ct;
  gt_assert(ctc && ctc->size);
  ct = gt_calloc(1, ctc->size);
  ct->c_class = ctc;
  ct->pvt = gt_calloc(1, sizeof (GtCustomTrackMembers));
  return ct;
}

GtCustomTrack* gt_custom_track_ref(GtCustomTrack *ctrack)
{
  gt_assert(ctrack);
  ctrack->pvt->reference_count++;
  return ctrack;
}

void gt_custom_track_delete(GtCustomTrack *ctrack)
{
  if (!ctrack) return;
  if (ctrack->pvt->reference_count) {
    ctrack->pvt->reference_count--;
    return;
  }
  gt_assert(ctrack->c_class);
  if (ctrack->c_class->free)
    ctrack->c_class->free(ctrack);
  gt_free(ctrack->pvt);
  gt_free(ctrack);
}

int gt_custom_track_render(GtCustomTrack *ctrack, GtGraphics *graphics,
                           unsigned int start_ypos, GtRange viewrange,
                           GtStyle *style, GtError *err)
{
  gt_assert(ctrack && ctrack->c_class && graphics && style);
  return ctrack->c_class->render(ctrack, graphics, start_ypos, viewrange, style,
                                 err);
}

int gt_custom_track_sketch(GtCustomTrack *ctrack, GtCanvas *canvas,
                           GtError *err)
{
  gt_assert(ctrack && canvas);
  return gt_canvas_visit_custom_track(canvas, ctrack, err);
}

unsigned long gt_custom_track_get_height(GtCustomTrack *ctrack)
{
  gt_assert(ctrack && ctrack->c_class);
  return ctrack->c_class->get_height(ctrack);
}

const char* gt_custom_track_get_title(GtCustomTrack *ctrack)
{
  gt_assert(ctrack && ctrack->c_class);
  return ctrack->c_class->get_title(ctrack);
}

void* gt_custom_track_cast(GT_UNUSED const GtCustomTrackClass *ctc,
                           GtCustomTrack *ct)
{
  gt_assert(ctc && ct && ct->c_class == ctc);
  return ct;
}
