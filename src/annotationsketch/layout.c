/*
  Copyright (c) 2008-2010 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2008-2010 Center for Bioinformatics, University of Hamburg

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
#include "core/cstr_api.h"
#include "core/ensure.h"
#include "core/hashmap.h"
#include "core/log.h"
#include "core/ma.h"
#include "core/mathsupport.h"
#include "core/minmax.h"
#include "core/msort.h"
#include "core/str.h"
#include "core/thread_api.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "extended/feature_node.h"
#include "extended/genome_node.h"

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
  bool own_twc,
       layout_done;
  GtArray *custom_tracks;
  GtHashmap *tracks,
            *blocks;
  GtRange viewrange;
  unsigned long nof_tracks;
  unsigned int width;
  GtRWLock *lock;
  GtTrackOrderingFunc track_ordering_func;
  GtBlockOrderingFunc block_ordering_func;
  void *t_cmp_data,
       *b_cmp_data;
};

static int blocklist_block_compare(const void *item1, const void *item2,
                                   void *data)
{
  gt_assert(item1 && item2);
  GtLayout *l = (GtLayout*) data;
  return l->block_ordering_func(*(GtBlock**) item1, *(GtBlock**) item2,
                                l->b_cmp_data);
}

static int add_tracklines(GT_UNUSED void *key, void *value,
                          void *data, GtError *err)
{
  GtTracklineInfo *add = (GtTracklineInfo*) data;
  double height;
  if (gt_track_get_height((GtTrack*) value, &height, add->style, err) < 0) {
    return 1;
  }
  add->height += height;
  return 0;
}

static int layout_tracks(void *key, void *value, void *data,
                         GtError *err)
{
  unsigned long i,
                max = 50;
  GtTrack *track = NULL;
  GtLayoutTraverseInfo *lti = (GtLayoutTraverseInfo*) data;
  GtArray *list = (GtArray*) value;
  GtStr *gt_track_key;
  GtBlock *block;
  int had_err = 0;
  bool split = true;
  double tmp = 50;
  gt_assert(list);

  /* to get a deterministic layout, we sort the GtBlocks for each type */
  if (lti->layout->block_ordering_func) {
    gt_array_sort_stable_with_data(list, blocklist_block_compare,
                                   lti->layout);
  }

  /* XXX: get first block for track property lookups, this should be reworked
     to allow arbitrary track keys! */
  block = *(GtBlock**) gt_array_get(list, 0);
  gt_track_key = gt_str_new_cstr((char*) key);

  /* obtain default settings*/
  if (gt_style_get_bool(lti->layout->style, "format", "split_lines", &split,
                         NULL, err) == GT_STYLE_QUERY_ERROR) {
    had_err = 1;
  }
  if (!had_err) {
    if (gt_style_get_num(lti->layout->style,
                         "format", "max_num_lines",
                         &tmp, NULL, err) == GT_STYLE_QUERY_ERROR) {
      had_err = 1;
    }
  }
  /* obtain track-specific settings, should be changed to query arbitrary
     track keys! */
  if (!had_err) {
    if (gt_style_get_bool(lti->layout->style, gt_block_get_type(block),
                          "split_lines",  &split, NULL,
                          err) == GT_STYLE_QUERY_ERROR) {
      had_err = 1;
    }
  }
  if (!had_err) {
    if (gt_style_get_num(lti->layout->style, gt_block_get_type(block),
                         "max_num_lines", &tmp, NULL,
                         err) == GT_STYLE_QUERY_ERROR) {
      had_err = 1;
    }
  }

  if (!had_err) {
    max = (unsigned long) tmp;
    track = gt_track_new(gt_track_key, max, split,
                         gt_line_breaker_captions_new(lti->layout,
                                                      lti->layout->width,
                                                      lti->layout->style));
    lti->layout->nof_tracks++;
    for (i = 0; !had_err && i < gt_array_size(list); i++) {
      block = *(GtBlock**) gt_array_get(list, i);
      had_err = gt_track_insert_block(track, block, err);
    }
  }
  if (!had_err) {
    gt_hashmap_add(lti->layout->tracks, gt_cstr_dup(gt_str_get(gt_track_key)),
                   track);
  }
  else
  {
    gt_track_delete(track);
  }

  gt_str_delete(gt_track_key);
  return had_err;
}

static int layout_all_tracks(GtLayout *layout, GtError *err)
{
  int had_err = 0;
  GtLayoutTraverseInfo lti;
  gt_assert(layout && layout->blocks);
  gt_error_check(err);

  if (!layout->layout_done) {
    lti.layout = layout;
    lti.twc = layout->twc;
    had_err = gt_hashmap_foreach(layout->blocks, layout_tracks, &lti, err);
    layout->layout_done = true;
  }
  return had_err;
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

static int check_width(unsigned int width, GtStyle *style, GtError *err)
{
  int had_err = 0;
  double margins = MARGINS_DEFAULT;
  if (gt_style_get_num(style,
                       "format", "margins",
                       &margins, NULL, err) == GT_STYLE_QUERY_ERROR) {
    had_err = -1;
  }
  if (!had_err && gt_double_smaller_double(width - 2*margins, 0))
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
  twc = gt_text_width_calculator_cairo_new(NULL, style, err);
  if (!twc)
    return NULL;
  layout = gt_layout_new_with_twc(diagram, width, style, twc, err);
  if (layout)
    layout->own_twc = true;
  else
    gt_text_width_calculator_delete(twc);
  return layout;
}

GtLayout* gt_layout_new_with_twc(GtDiagram *diagram,
                                 unsigned int width,
                                 GtStyle *style,
                                 GtTextWidthCalculator *twc,
                                 GtError *err)
{
  GtLayout *layout;
  GtHashmap *blocks;
  gt_assert(diagram);
  gt_assert(style);
  gt_assert(twc);
  gt_assert(err);
  if (check_width(width, style, err) < 0)
    return NULL;
  layout = gt_calloc(1, sizeof (GtLayout));
  layout->twc = twc;
  layout->style = style;
  layout->width = width;
  layout->blocks = NULL;
  layout->viewrange = gt_diagram_get_range(diagram);
  layout->nof_tracks = 0;
  layout->track_ordering_func = NULL;
  layout->block_ordering_func = gt_block_compare;
  layout->t_cmp_data = NULL;
  layout->b_cmp_data = NULL;
  layout->lock = gt_rwlock_new();
  layout->own_twc = false;
  layout->layout_done = false;
  layout->custom_tracks = gt_array_ref(gt_diagram_get_custom_tracks(diagram));
  /* XXX: use other container type here! */
  layout->tracks = gt_hashmap_new(GT_HASH_STRING, gt_free_func,
                                  (GtFree) gt_track_delete);
  blocks = gt_diagram_get_blocks(diagram, err);
  if (!blocks) {
    gt_array_delete(layout->custom_tracks);
    gt_hashmap_delete(layout->tracks);
    gt_free(layout);
    return NULL;
  } else
    layout->blocks = gt_hashmap_ref(blocks);
  return layout;
}

void gt_layout_delete(GtLayout *layout)
{
  if (!layout) return;
  gt_rwlock_wrlock(layout->lock);
  if (layout->twc && layout->own_twc)
    gt_text_width_calculator_delete(layout->twc);
  gt_hashmap_delete(layout->tracks);
  gt_array_delete(layout->custom_tracks);
  if (layout->blocks)
    gt_hashmap_delete(layout->blocks);
  gt_rwlock_unlock(layout->lock);
  gt_rwlock_delete(layout->lock);
  gt_free(layout);
}

static int track_cmp_wrapper(const void *t1, const void *t2, void *data)
{
  const char *s1 = (const char*) t1, *s2 = (const char*) t2;
  GtLayoutTraverseInfo *lti = (GtLayoutTraverseInfo*) data;
  gt_assert(t1 && t2 && lti && lti->layout);
  return lti->layout->track_ordering_func(s1, s2, lti->layout->t_cmp_data);
}

int gt_layout_sketch(GtLayout *layout, GtCanvas *target_canvas, GtError *err)
{
  int had_err = 0;
  unsigned long i;
  GtRenderTraverseInfo rti;
  gt_assert(layout && target_canvas);

  had_err = layout_all_tracks(layout, err);
  if (!had_err) {
    rti.layout = layout;
    rti.canvas = target_canvas;
    had_err = gt_canvas_visit_layout_pre(target_canvas, layout, err);
    if (had_err) return had_err;

    if (layout->track_ordering_func == NULL) {
      had_err = gt_hashmap_foreach_in_key_order(layout->tracks, render_tracks,
                                                &rti, err);
    } else {
      had_err = gt_hashmap_foreach_ordered(layout->tracks, render_tracks,
                                           &rti, (GtCompare) track_cmp_wrapper,
                                           err);
    }
    if (had_err) return had_err;
    had_err = gt_canvas_visit_layout_post(target_canvas, layout, err);
    if (had_err) return had_err;
    for (i=0;i<gt_array_size(layout->custom_tracks);i++)
    {
      GtCustomTrack *ct = *(GtCustomTrack**) gt_array_get(layout->custom_tracks,
                                                          i);
      had_err = render_custom_tracks(NULL, ct, &rti, err);
    }
  }
  return had_err ? -1 : 0;
}

void gt_layout_set_track_ordering_func(GtLayout *layout,
                                       GtTrackOrderingFunc track_ordering_func,
                                       void *data)
{
  gt_assert(layout && track_ordering_func);
  layout->track_ordering_func = track_ordering_func;
  layout->t_cmp_data = data;
}

void gt_layout_unset_track_ordering_func(GtLayout *layout)
{
  gt_assert(layout);
  layout->track_ordering_func = NULL;
  layout->t_cmp_data = NULL;
}

void gt_layout_set_block_ordering_func(GtLayout *layout,
                                       GtBlockOrderingFunc block_ordering_func,
                                       void *data)
{
  gt_assert(layout && block_ordering_func);
  layout->block_ordering_func = block_ordering_func;
  layout->b_cmp_data = data;
}

void gt_layout_unset_block_ordering_func(GtLayout *layout)
{
  gt_assert(layout);
  layout->block_ordering_func = NULL;
  layout->b_cmp_data = NULL;
}

GtRange gt_layout_get_range(const GtLayout *layout)
{
  gt_assert(layout);
  return layout->viewrange;
}

GtTextWidthCalculator* gt_layout_get_twc(const GtLayout *layout)
{
  gt_assert(layout);
  return layout->twc;
}

int gt_layout_get_height(GtLayout *layout, unsigned long *result,
                         GtError *err)
{
  int had_err = 0;
  GtTracklineInfo lines;
  double tmp, head_track_space = HEAD_TRACK_SPACE_DEFAULT;
  bool show_track_captions = true;
  unsigned long height,
                line_height,
                i;
  gt_assert(layout);

  had_err = layout_all_tracks(layout, err);
  if (!had_err) {
    /* get dynamic heights from tracks */
    lines.style = layout->style;
    lines.height = 0;
    if (gt_hashmap_foreach(layout->tracks, add_tracklines, &lines, err) < 0) {
      return -1;
    }
    height = lines.height;

    /* obtain line height and spacer from style */
    tmp = BAR_HEIGHT_DEFAULT;
    if (gt_style_get_num(layout->style,
                         "format", "bar_height",
                         &tmp, NULL, err) == GT_STYLE_QUERY_ERROR) {
      return -1;
    }
    line_height = tmp;

    tmp = BAR_VSPACE_DEFAULT;
    if (gt_style_get_num(layout->style,
                         "format", "bar_vspace",
                         &tmp, NULL, err) == GT_STYLE_QUERY_ERROR) {
      return -1;
    }
    line_height += tmp;

    if (gt_style_get_bool(layout->style,
                          "format","show_track_captions",
                          &show_track_captions,
                          NULL, err) == GT_STYLE_QUERY_ERROR) {
      return -1;
    }

    /* add custom track space allotment */
    if (show_track_captions)
    {
      double theight = TEXT_SIZE_DEFAULT,
             captionspace = CAPTION_BAR_SPACE_DEFAULT;
      if (gt_style_get_num(layout->style,
                           "format", "track_caption_font_size",
                           &theight, NULL, err) == GT_STYLE_QUERY_ERROR) {
        return -1;
      }
      if (gt_style_get_num(layout->style,
                           "format", "track_caption_space",
                           &captionspace, NULL, err) == GT_STYLE_QUERY_ERROR) {
        return -1;
      }
      height += gt_array_size(layout->custom_tracks)
                    * (theight + captionspace);
    }

    for (i=0;i<gt_array_size(layout->custom_tracks);i++)
    {
      GtCustomTrack *ct = *(GtCustomTrack**) gt_array_get(layout->custom_tracks,
                                                          i);
      height += gt_custom_track_get_height(ct);
      if (gt_style_get_num(layout->style, "format", "track_vspace", &tmp,
                           NULL, err) == GT_STYLE_QUERY_ERROR) {
        return -1;
      }
      height += tmp;
    }

    /* add header space and footer */
    if (gt_style_get_num(layout->style, "format", "ruler_space",
                         &head_track_space, NULL,
                         err) == GT_STYLE_QUERY_ERROR) {
      return -1;
    }
    height += HEADER_SPACE + head_track_space + FOOTER_SPACE;

    *result = height;
    return 0;
  }
  return had_err;
}
