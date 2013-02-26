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

#include <math.h>
#include "core/class_alloc_lock.h"
#include "core/unused_api.h"
#include "annotationsketch/coords.h"
#include "annotationsketch/custom_track_example.h"
#include "annotationsketch/custom_track_rep.h"

struct GtCustomTrackExample {
  const GtCustomTrack parent_instance;
  unsigned long height;
  GtStr *title;
};

#define gt_custom_track_example_cast(ct)\
        gt_custom_track_cast(gt_custom_track_example_class(), ct)

int gt_custom_track_example_sketch(GtCustomTrack *ct, GtGraphics *graphics,
                                   unsigned int start_ypos,
                                   GtRange viewrange,
                                   GtStyle *style, GT_UNUSED GtError *err)
{
  int had_err = 0;
  GtCustomTrackExample *cte;
  double margins;
  GtColor color;
  char buffer[BUFSIZ];
  cte = gt_custom_track_example_cast(ct);

  color.red = color.blue = color.green = 0.9;
  color.alpha = 0.3;

  margins = gt_graphics_get_xmargins(graphics);
  gt_graphics_draw_rectangle(graphics,
                             margins,
                             start_ypos,
                             true,
                             color,
                             true,
                             color,
                             1.0,
                             gt_graphics_get_image_width(graphics)
                               - 2*margins,
                             cte->height);

  sprintf(buffer, "Style %p", style);
  gt_graphics_draw_text_centered(graphics,
                            margins+((gt_graphics_get_image_width(graphics)
                              - 2*margins) /2),
                             start_ypos + gt_graphics_get_text_height(graphics),
                             ((const char*) buffer));

  sprintf(buffer, "Range %lu - %lu", viewrange.start, viewrange.end);
  gt_graphics_draw_text_centered(graphics,
                           margins+((gt_graphics_get_image_width(graphics)
                             - 2*margins) /2),
                           start_ypos + 2*gt_graphics_get_text_height(graphics),
                           ((const char*) buffer));

  sprintf(buffer, "Start position y=%u", start_ypos);
  gt_graphics_draw_text_centered(graphics,
                           margins+((gt_graphics_get_image_width(graphics)
                             - 2*margins) /2),
                           start_ypos + 3*gt_graphics_get_text_height(graphics),
                           ((const char*) buffer));

  return had_err;
}

unsigned long gt_custom_track_example_get_height(GtCustomTrack *ct)
{
  GtCustomTrackExample *cte;
  cte = gt_custom_track_example_cast(ct);
  return cte->height;
}

const char* gt_custom_track_example_get_title(GtCustomTrack *ct)
{
  GtCustomTrackExample *cte;
  cte = gt_custom_track_example_cast(ct);
  return gt_str_get(cte->title);
}

void gt_custom_track_example_delete(GtCustomTrack *ct)
{
  GtCustomTrackExample *cte;
  if (!ct) return;
  cte = gt_custom_track_example_cast(ct);
  gt_str_delete(cte->title);
}

const GtCustomTrackClass* gt_custom_track_example_class(void)
{
  static const GtCustomTrackClass *ctc = NULL;
  gt_class_alloc_lock_enter();
  if (!ctc)
  {
    ctc = gt_custom_track_class_new(sizeof (GtCustomTrackExample),
                                    gt_custom_track_example_sketch,
                                    gt_custom_track_example_get_height,
                                    gt_custom_track_example_get_title,
                                    gt_custom_track_example_delete);
  }
  gt_class_alloc_lock_leave();
  return ctc;
}

GtCustomTrack* gt_custom_track_example_new()
{
  GtCustomTrackExample *cte;
  GtCustomTrack *ct;
  ct = gt_custom_track_create(gt_custom_track_example_class());
  cte = gt_custom_track_example_cast(ct);
  cte->title = gt_str_new_cstr("Example custom track for AnnotationSketch");
  cte->height = 40;
  return ct;
}
