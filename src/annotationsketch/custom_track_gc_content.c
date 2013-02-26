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
#include "annotationsketch/coords.h"
#include "annotationsketch/default_formats.h"
#include "annotationsketch/custom_track_gc_content.h"
#include "annotationsketch/custom_track_rep.h"
#include "core/class_alloc_lock.h"
#include "core/log.h"
#include "core/ma.h"
#include "core/mathsupport.h"
#include "core/minmax.h"
#include "core/unused_api.h"

struct GtCustomTrackGcContent {
  const GtCustomTrack parent_instance;
  unsigned long windowsize,
                height;
  double avg;
  bool show_scale;
  GtStr *title;
  const char *seq;
  unsigned long seqlen;
};

#define gt_custom_track_gc_content_cast(ct)\
        gt_custom_track_cast(gt_custom_track_gc_content_class(), ct)

static double get_val_for_pos(GtCustomTrackGcContent *ctgc, unsigned long pos)
{
  unsigned long i,
                gc_count = 0,
                bases = 0;
  for (i=0;i<ctgc->windowsize;i++)
  {
    if (pos + i > ctgc->seqlen)
      return -1;
    if (ctgc->seq[pos+i] == 'g' || ctgc->seq[pos+i] == 'c'
         || ctgc->seq[pos+i] == 'G' || ctgc->seq[pos+i] == 'C')
    {
      gc_count++;
    }
    bases++;
  }
  return ((double) gc_count)/((double) MIN(ctgc->windowsize, bases));
}

int gt_custom_track_gc_content_sketch(GtCustomTrack *ct, GtGraphics *graphics,
                                      unsigned int start_ypos,
                                      GtRange viewrange,
                                      GtStyle *style, GT_UNUSED GtError *err)
{
  int had_err = 0;
  GtCustomTrackGcContent *ctgc;
  double iter, iter_step, GT_UNUSED value, *data;
  unsigned long n;
  GtRange value_range = {0, 1};
  GtColor color, grey, black;
  gt_assert(ct && graphics && viewrange.start <= viewrange.end);
  gt_assert(gt_double_smaller_double(0, gt_graphics_get_image_width(graphics)
                                        -2*gt_graphics_get_xmargins(graphics)));

  ctgc = gt_custom_track_gc_content_cast(ct);

  (void) gt_style_get_color(style, "GC_content", "stroke", &color, NULL, NULL);
  grey.red = grey.blue = grey.green = 0.8;
  grey.alpha = 0.9;
  black.red = black.blue = black.green = 0.0;
  black.alpha = 0.9;
  value = get_val_for_pos(ctgc, viewrange.start);

  iter_step = (double) gt_range_length(&viewrange)
                / ((double) gt_graphics_get_image_width(graphics)
                   - 2*gt_graphics_get_xmargins(graphics));

  gt_log_log("len=%lu, iter_step = %f, width = %f, margins = %f\n",
                      gt_range_length(&viewrange),
                      iter_step, gt_graphics_get_image_width(graphics),
                      gt_graphics_get_xmargins(graphics));

  data = gt_calloc(ceil(gt_range_length(&viewrange)/iter_step)+1,
                   sizeof (double));
  n = 0;
  for (iter=viewrange.start+1;
       gt_double_smaller_double(iter, viewrange.end-ctgc->windowsize);
       iter+=iter_step)
  {
    if (floor(iter) >= ctgc->seqlen) break;
    data[n++] = get_val_for_pos(ctgc, floor(iter));
  }

  gt_log_log("i=%lu, widthval = %f\n", n,
                 (gt_graphics_get_image_width(graphics)
                    - 2*gt_graphics_get_xmargins(graphics)));
  if (ctgc->show_scale)
  {
    gt_graphics_draw_horizontal_line(graphics,
                                     gt_graphics_get_xmargins(graphics) + 1,
                                     start_ypos+1,
                                     black,
                                     2,
                                     1.0);
    gt_graphics_draw_horizontal_line(graphics,
                                     gt_graphics_get_xmargins(graphics) + 1,
                                     start_ypos + ctgc->height,
                                     black,
                                     2,
                                     1.0);
    gt_graphics_draw_text(graphics,
                          gt_graphics_get_xmargins(graphics)+ 5,
                          start_ypos
                            + gt_graphics_get_text_height(graphics)/2 - 1,
                          "100\%");
    gt_graphics_draw_text(graphics,
                          gt_graphics_get_xmargins(graphics)+ 5,
                          start_ypos
                            + ctgc->height
                            + gt_graphics_get_text_height(graphics)/2 - 1,
                          "0\%");
  }
  gt_graphics_draw_horizontal_line(graphics,
                                   gt_graphics_get_xmargins(graphics),
                                   start_ypos + (1.0-ctgc->avg) * ctgc->height,
                                   grey,
                                   (gt_graphics_get_image_width(graphics)
                                     - 2*gt_graphics_get_xmargins(graphics)),
                                   1.0);
  if (ctgc->show_scale)
  {
    gt_graphics_draw_vertical_line(graphics,
                                   gt_graphics_get_xmargins(graphics),
                                   start_ypos,
                                   black,
                                   ctgc->height,
                                   1.0);
  }
  gt_graphics_draw_curve_data(graphics,
                              gt_graphics_get_xmargins(graphics),
                              start_ypos,
                              color,
                              data,
                              n,
                              value_range,
                              ctgc->height);
  gt_free(data);
  return had_err;
}

unsigned long gt_custom_track_gc_content_get_height(GtCustomTrack *ct)
{
  GtCustomTrackGcContent *ctgc;
  ctgc = gt_custom_track_gc_content_cast(ct);
  return ctgc->height;
}

const char* gt_custom_track_gc_content_get_title(GtCustomTrack *ct)
{
  GtCustomTrackGcContent *ctgc;
  ctgc = gt_custom_track_gc_content_cast(ct);
  return gt_str_get(ctgc->title);
}

void gt_custom_track_gc_content_delete(GtCustomTrack *ct)
{
  GtCustomTrackGcContent *ctgc;
  if (!ct) return;
  ctgc = gt_custom_track_gc_content_cast(ct);
  gt_str_delete(ctgc->title);
}

const GtCustomTrackClass* gt_custom_track_gc_content_class(void)
{
  static const GtCustomTrackClass *ctc = NULL;
  gt_class_alloc_lock_enter();
  if (!ctc)
  {
    ctc = gt_custom_track_class_new(sizeof (GtCustomTrackGcContent),
                                    gt_custom_track_gc_content_sketch,
                                    gt_custom_track_gc_content_get_height,
                                    gt_custom_track_gc_content_get_title,
                                    gt_custom_track_gc_content_delete);
  }
  gt_class_alloc_lock_leave();
  return ctc;
}

GtCustomTrack* gt_custom_track_gc_content_new(const char *seq,
                                              unsigned long seqlen,
                                              unsigned long windowsize,
                                              unsigned long height,
                                              double avg,
                                              bool show_scale)
{
  GtCustomTrackGcContent *ctgc;
  GtCustomTrack *ct;
  char buf[BUFSIZ];
  ct = gt_custom_track_create(gt_custom_track_gc_content_class());
  ctgc = gt_custom_track_gc_content_cast(ct);
  ctgc->windowsize = windowsize;
  ctgc->height = height;
  ctgc->seq = seq;
  ctgc->seqlen = seqlen;
  ctgc->avg = avg;
  ctgc->show_scale = show_scale;
  ctgc->title = gt_str_new_cstr("GC content (window size ");
  gt_str_append_ulong(ctgc->title, ctgc->windowsize);
  if (gt_double_smaller_double(0, avg))
  {
    (void) snprintf(buf, BUFSIZ, ", average: %.1f%%", avg*100);
    gt_str_append_cstr(ctgc->title, buf);
  }
  gt_str_append_cstr(ctgc->title, ")");
  return ct;
}
