/*
  Copyright (c) 2008-2009 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2008-2009 Center for Bioinformatics, University of Hamburg

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

#include "annotationsketch/block.h"
#include "annotationsketch/canvas.h"
#include "annotationsketch/cliptype.h"
#include "annotationsketch/default_formats.h"
#include "annotationsketch/diagram.h"
#include "annotationsketch/layout.h"
#include "annotationsketch/line_breaker_captions.h"
#include "annotationsketch/style.h"
#include "annotationsketch/text_width_calculator_cairo.h"
#include "annotationsketch/track.h"
#include "core/basename_api.h"
#include "core/cstr.h"
#include "core/ensure.h"
#include "core/hashmap.h"
#include "core/log.h"
#include "core/ma.h"
#include "core/minmax.h"
#include "core/msort.h"
#include "core/str.h"
#include "core/undef.h"
#include "core/unused_api.h"
#include "extended/feature_node.h"
#include "extended/genome_node.h"

/* used to separate a filename from the type in a track name */
#define FILENAME_TYPE_SEPARATOR  '|'

typedef struct {
  GtTextWidthCalculator *twc;
  GtLayout *layout;
} GtLayoutTraverseInfo;

typedef struct {
  GtCanvas *canvas;
  GtLayout *layout;
} GtRenderTraverseInfo;

typedef struct {
  unsigned long height;
  GtStyle *style;
} GtTracklineInfo;

struct GtLayout {
  GtStyle *style;
  GtTextWidthCalculator *twc;
  bool own_twc;
  GtArray *custom_tracks;
  GtHashmap *tracks;
  GtRange viewrange;
  unsigned long nof_tracks;
  unsigned int width;
};

static int blocklist_block_compare(const void *item1, const void *item2)
{
  gt_assert(item1 && item2);
  return gt_block_compare(*(GtBlock**) item1, *(GtBlock**) item2);
}

static int add_tracklines(GT_UNUSED void *key, void *value,
                          void *data, GT_UNUSED GtError *err)
{
  GtTracklineInfo *add = (GtTracklineInfo*) data;
  add->height += gt_track_get_height((GtTrack*) value, add->style);
  return 0;
}

static GtStr* track_key_new(const char *filename, const char *type)
{
  GtStr *track_key;
  track_key = gt_str_new_cstr(filename);
  gt_str_append_char(track_key, FILENAME_TYPE_SEPARATOR);
  gt_str_append_cstr(track_key, type);
  return track_key;
}

static int layout_tracks(void *key, void *value, void *data,
                         GT_UNUSED GtError *err)
{
  unsigned long i, max;
  GtTrack *track;
  GtLayoutTraverseInfo *lti = (GtLayoutTraverseInfo*) data;
  GtArray *list = (GtArray*) value;
  char *filename;
  GtStr *gt_track_key;
  const char *type = key;
  GtBlock *block;
  bool split;
  double tmp;
  gt_assert(type && list);

  /* to get a deterministic layout, we sort the GtBlocks for each type */
  gt_array_sort_stable(list, blocklist_block_compare);
  /* we take the basename of the filename to have nicer output in the
     generated graphic. this might lead to ``collapsed'' tracks, if two files
     with different paths have the same basename. */
  block = *(GtBlock**) gt_array_get(list, 0);
  filename = gt_basename(gt_genome_node_get_filename(
                                     (GtGenomeNode*)
                                       gt_block_get_top_level_feature(block)));
  gt_track_key = track_key_new(filename, type);
  gt_free(filename);

  if (!gt_style_get_bool(lti->layout->style, "format", "split_lines", &split,
                         NULL))
    split = true;
  if (split)
    if (!gt_style_get_bool(lti->layout->style, type, "split_lines", &split,
                           NULL))
      split = true;
  if (gt_style_get_num(lti->layout->style, type, "max_num_lines", &tmp, NULL))
    max = tmp;
  else
    max = 50;

  track = gt_track_new(gt_track_key, max, split,
                       gt_line_breaker_captions_new(lti->layout,
                                                    lti->layout->width,
                                                    lti->layout->style));
  lti->layout->nof_tracks++;
  for (i = 0; i < gt_array_size(list); i++) {
    block = *(GtBlock**) gt_array_get(list, i);
    gt_track_insert_block(track, block);
  }
  gt_hashmap_add(lti->layout->tracks, gt_cstr_dup(gt_str_get(gt_track_key)),
                 track);
  gt_str_delete(gt_track_key);
  return 0;
}

static int render_tracks(GT_UNUSED void *key, void *value, void *data,
                         GtError *err)
{
  GtRenderTraverseInfo *rti = (GtRenderTraverseInfo*) data;
  GtTrack *track = (GtTrack*) value;
  int had_err = 0;
  gt_assert(rti && track);
  had_err = gt_track_sketch(track, rti->canvas, err);
  return had_err;
}

static int render_custom_tracks(GT_UNUSED void *key, void *value, void *data,
                                GtError *err)
{
  GtRenderTraverseInfo *rti = (GtRenderTraverseInfo*) data;
  GtCustomTrack *ctrack = (GtCustomTrack*) value;
  int had_err = 0;
  gt_assert(rti && ctrack);
  had_err = gt_custom_track_sketch(ctrack, rti->canvas, err);
  return had_err;
}

static int check_width(unsigned int width,
                       GtStyle *style,
                       GtError *err)
{
  int had_err = 0;
  double margins = MARGINS_DEFAULT;
  gt_style_get_num(style, "format", "margins", &margins, NULL);
  if (width - 2*margins < 0)
  {
    gt_error_set(err, "layout width must at least be twice the x-margin size "
                      "(2*%.1f=%.1f) but was %u",
                      margins,
                      2*margins,
                      width);
    had_err = -1;
  }
  return had_err;
}

GtLayout* gt_layout_new(GtDiagram *diagram,
                        unsigned int width,
                        GtStyle *style,
                        GtError *err)
{
  GtLayout *layout;
  GtTextWidthCalculator *twc;
  gt_assert(diagram && width > 0 && style && err);
  if (check_width(width, style, err) < 0)
    return NULL;
  twc = gt_text_width_calculator_cairo_new(NULL, style);
  layout = gt_layout_new_with_twc(diagram, width, style, twc, err);
  layout->own_twc = true;
  return layout;
}

GtLayout* gt_layout_new_with_twc(GtDiagram *diagram,
                                 unsigned int width,
                                 GtStyle *style,
                                 GtTextWidthCalculator *twc,
                                 GtError *err)
{
  GtLayout *layout;
  GtLayoutTraverseInfo lti;
  gt_assert(diagram && style && twc && err);
  if (check_width(width, style, err) < 0)
    return NULL;
  layout = gt_calloc(1, sizeof (GtLayout));
  layout->twc = twc;
  layout->style = style;
  layout->width = width;
  layout->viewrange = gt_diagram_get_range(diagram);
  layout->nof_tracks = 0;
  lti.layout = layout;
  lti.twc = twc;
  layout->own_twc = false;
  gt_diagram_build(diagram);
  layout->custom_tracks = gt_diagram_get_custom_tracks(diagram);
  /* XXX: use other container type here! */
  layout->tracks = gt_hashmap_new(HASH_STRING, gt_free_func,
                                  (GtFree) gt_track_delete);
  (void) gt_hashmap_foreach(gt_diagram_get_blocks(diagram),
                            layout_tracks,
                            &lti,
                            NULL);
  return layout;
}

void gt_layout_delete(GtLayout *layout)
{
  if (!layout) return;
  if (layout->twc && layout->own_twc)
    gt_text_width_calculator_delete(layout->twc);
  gt_hashmap_delete(layout->tracks);
  gt_free(layout);
}

int gt_layout_sketch(GtLayout *layout, GtCanvas *canvas, GtError *err)
{
  int had_err = 0;
  unsigned long i;
  GtRenderTraverseInfo rti;
  gt_assert(layout && canvas);
  rti.layout = layout;
  rti.canvas = canvas;
  had_err = gt_canvas_visit_layout_pre(canvas, layout, err);
  if (had_err) return had_err;
  had_err = gt_hashmap_foreach_in_key_order(layout->tracks, render_tracks,
                                            &rti, err);
  if (had_err) return had_err;
  had_err = gt_canvas_visit_layout_post(canvas, layout, err);
  if (had_err) return had_err;
  for (i=0;i<gt_array_size(layout->custom_tracks);i++)
  {
    GtCustomTrack *ct = *(GtCustomTrack**) gt_array_get(layout->custom_tracks,
                                                        i);
    had_err = render_custom_tracks(NULL, ct, &rti, err);
  }
  return had_err;
}

int gt_layout_get_number_of_tracks(const GtLayout *layout)
{
  gt_assert(layout);
  return layout->nof_tracks;
}

GtRange gt_layout_get_range(const GtLayout *layout)
{
  gt_assert(layout);
  return layout->viewrange;
}

unsigned int gt_layout_get_width(const GtLayout *layout)
{
  gt_assert(layout);
  return layout->width;
}

GtTextWidthCalculator* gt_layout_get_twc(const GtLayout *layout)
{
  gt_assert(layout);
  return layout->twc;
}

unsigned long gt_layout_get_height(const GtLayout *layout)
{
  GtTracklineInfo lines;
  double tmp, head_track_space = HEAD_TRACK_SPACE_DEFAULT;
  bool show_track_captions;
  unsigned long height,
                line_height,
                i;
  gt_assert(layout);

  /* get dynamic heights from tracks */
  lines.style = layout->style; lines.height = 0;
  gt_hashmap_foreach(layout->tracks, add_tracklines,
                     &lines, NULL);
  height = lines.height;

  /* obtain line height and spacer from style */
  if (gt_style_get_num(layout->style, "format", "bar_height", &tmp, NULL))
    line_height = tmp;
  else
    line_height = BAR_HEIGHT_DEFAULT;
  if (gt_style_get_num(layout->style, "format", "bar_vspace", &tmp, NULL))
    line_height += tmp;
  else
    line_height += BAR_VSPACE_DEFAULT;

  if (!(gt_style_get_bool(layout->style, "format","show_track_captions",
                          &show_track_captions, NULL)))
    show_track_captions = true;

  /* add custom track space allotment */
  if (show_track_captions)
  {
    double theight = TOY_TEXT_HEIGHT,
           captionspace = CAPTION_BAR_SPACE_DEFAULT;
    gt_style_get_num(layout->style, "format", "track_caption_font_size",
                     &theight, NULL);
    gt_style_get_num(layout->style, "format", "track_caption_space",
                     &captionspace, NULL);
    height += gt_array_size(layout->custom_tracks)
                  * (theight + captionspace);
  }

  for (i=0;i<gt_array_size(layout->custom_tracks);i++)
  {
    GtCustomTrack *ct = *(GtCustomTrack**) gt_array_get(layout->custom_tracks,
                                                        i);
    height += gt_custom_track_get_height(ct);
    gt_style_get_num(layout->style, "format", "track_vspace", &tmp, NULL);
    height += tmp;

  }

  /* add header space and footer */
  gt_style_get_num(layout->style, "format", "ruler_space", &head_track_space,
                   NULL);
  height += HEADER_SPACE + head_track_space + FOOTER_SPACE;
  return height;
}
