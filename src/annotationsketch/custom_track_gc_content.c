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
#include "core/ma.h"
#include "core/minmax.h"
#include "core/unused_api.h"
#include "annotationsketch/coords.h"
#include "annotationsketch/custom_track_gc_content.h"
#include "annotationsketch/custom_track_rep.h"

struct GtCustomTrackGcContent {
  const GtCustomTrack parent_instance;
  unsigned long windowsize,
                height;
  double avg;
  GtStr *title;
  GtSeq *seq;
};

#define gt_custom_track_gc_content_cast(ct)\
        gt_custom_track_cast(gt_custom_track_gc_content_class(), ct)

static inline double get_val_for_pos(GtCustomTrackGcContent *ctgc,
                                     unsigned long pos)
{
  const char *seq = gt_seq_get_orig(ctgc->seq);
  unsigned long i,
                gc_count = 0,
                bases = 0;
  for (i=0;i<ctgc->windowsize;i++)
  {
    if (pos + i > gt_seq_length(ctgc->seq))
      break;
    if (seq[pos+i] == 'g' || seq[pos+i] == 'c'
         || seq[pos+i] == 'G' || seq[pos+i] == 'C')
    {
      gc_count++;
    }
    bases++;
  }
  return (double) gc_count/MIN((double) ctgc->windowsize, bases);
}

int gt_custom_track_gc_content_sketch(GtCustomTrack *ct, GtGraphics *graphics,
                                      unsigned int start_ypos,
                                      GtRange viewrange,
                                      GT_UNUSED GtStyle *style)
{
  GtCustomTrackGcContent *ctgc;
  double iter,
         iter_step,
         drawpos,
         value,
         *data;
  unsigned long i;
  GtColor color;
  ctgc = gt_custom_track_gc_content_cast(ct);

  color.red = 1;
  color.green = color.blue = 0;
  color.alpha = .8;
  iter_step = gt_range_length(&viewrange)/gt_graphics_get_image_width(graphics);

  value = get_val_for_pos(ctgc, viewrange.start);
  drawpos = gt_graphics_get_xmargins(graphics)
                + gt_coords_convert_point(viewrange, viewrange.start)
                  * (gt_graphics_get_image_width(graphics)
                      - 2*gt_graphics_get_xmargins(graphics));

  data = gt_calloc(ceil(gt_range_length(&viewrange)/iter_step)+1,
                   sizeof (double));
  i = 0;
  for (iter=viewrange.start+1; iter<viewrange.end; iter+=iter_step)
  {
    data[i++] = get_val_for_pos(ctgc, floor(iter));
  }
  gt_graphics_draw_horizontal_line(graphics, gt_graphics_get_xmargins(graphics),
                                   start_ypos + (1.0-ctgc->avg) * ctgc->height,
                                   (gt_graphics_get_image_width(graphics)
                                     - 2*gt_graphics_get_xmargins(graphics)));
  gt_graphics_draw_curve_data(graphics, gt_graphics_get_xmargins(graphics),
                              start_ypos, color, data, i, ctgc->height);

  return 0;
}

unsigned long gt_custom_track_gc_content_get_height(GtCustomTrack *ct)
{
  GtCustomTrackGcContent *ctgc;
  ctgc = gt_custom_track_gc_content_cast(ct);
  return ctgc->height;
}

GtStr* gt_custom_track_gc_content_get_title(GtCustomTrack *ct)
{
  GtCustomTrackGcContent *ctgc;
  ctgc = gt_custom_track_gc_content_cast(ct);
  return ctgc->title;
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
  if (!ctc)
  {
    ctc = gt_custom_track_class_new(sizeof (GtCustomTrackGcContent),
                                    gt_custom_track_gc_content_sketch,
                                    gt_custom_track_gc_content_get_height,
                                    gt_custom_track_gc_content_get_title,
                                    gt_custom_track_gc_content_delete);
  }
  return ctc;
}

GtCustomTrack* gt_custom_track_gc_content_new(GtSeq *seq,
                                              unsigned long windowsize,
                                              unsigned long height,
                                              double avg)
{
  GtCustomTrackGcContent *ctgc;
  GtCustomTrack *ct;
  char buf[BUFSIZ];
  ct = gt_custom_track_create(gt_custom_track_gc_content_class());
  ctgc = gt_custom_track_gc_content_cast(ct);
  ctgc->windowsize = windowsize;
  ctgc->height = height;
  ctgc->seq = seq;
  ctgc->avg = avg;
  /* ctgc->title = gt_str_new_cstr(gt_seq_get_description(seq)); */
  ctgc->title = gt_str_new_cstr("GC content");
  if (avg > 0)
  {
    snprintf(buf, BUFSIZ, " (average: %.1f%%)", avg*100);
    gt_str_append_cstr(ctgc->title, buf);
  }
  return ct;
}
