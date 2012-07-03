/*
  Copyright (c) 2007      Christin Schaerfer <schaerfer@zbh.uni-hamburg.de>
  Copyright (c) 2008-2010 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2007-2010 Center for Bioinformatics, University of Hamburg

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

#include "annotationsketch/default_formats.h"
#include "annotationsketch/line.h"
#include "annotationsketch/line_breaker_bases.h"
#include "annotationsketch/style.h"
#include "annotationsketch/track.h"
#include "core/ensure.h"
#include "core/hashtable.h"
#include "core/ma.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "extended/gff3_defines.h"

struct GtTrack {
  GtStr *title;
  unsigned long max_num_lines, discarded_blocks;
  GtLineBreaker *lb;
  bool split;
  unsigned long y_index;
  GtArray *lines;
};

GtTrack* gt_track_new(GtStr *title, unsigned long max_num_lines,
                      bool split_lines, GtLineBreaker *lb)
{
  GtTrack *track;
  gt_assert(title && lb);
  track = gt_calloc(1, sizeof (GtTrack));
  gt_assert(track);
  track->title = gt_str_ref(title);
  track->lines = gt_array_new(sizeof (GtLine*));
  track->max_num_lines = max_num_lines;
  track->split = split_lines;
  track->lb = lb;
  track->y_index = GT_UNDEF_ULONG;
  return track;
}

static int get_next_free_line(GtTrack *track, GtLine **result, GtBlock *block,
                              GtError *err)
{
  unsigned long i;
  GtLine* line;
  int had_err = 0;
  bool is_occupied;
  gt_assert(track);

  /* find unoccupied line -- may need optimisation */
  for (i = 0; i < gt_array_size(track->lines); i++) {
    line = *(GtLine**) gt_array_get(track->lines, i);
    had_err = gt_line_breaker_line_is_occupied(track->lb, &is_occupied, line,
                                               block, err);
    if (had_err)
      break;
    if (!is_occupied) {
      *result = line;
      return 0;
    }
  }
  /* all lines are occupied, we need o create a new one */
  if (!had_err) {
    /* if line limit is hit, do not create any more lines! */
    if (track->max_num_lines != GT_UNDEF_ULONG
          && gt_array_size(track->lines) == track->max_num_lines)
    {
      track->discarded_blocks++;
      *result = NULL;
    }
    /* make sure there is only one line if 'split_lines' is set to false */
    if (!track->split)
    {
      if (gt_array_size(track->lines) < 1)
      {
        line = gt_line_new();
        gt_array_add(track->lines, line);
      }
      else
        line = *(GtLine**) gt_array_get(track->lines, 0);
      gt_assert(gt_array_size(track->lines) == 1);
    }
    else
    {
      line = gt_line_new();
      gt_array_add(track->lines, line);
    }
    gt_assert(line);
  }
  *result = line;
  return had_err;
}

unsigned long gt_track_get_y_index(const GtTrack *track)
{
  gt_assert(track);
  return track->y_index;
}

void gt_track_set_y_index(GtTrack *track, unsigned long y_index)
{
  gt_assert(track);
  track->y_index = y_index;
}

unsigned long gt_track_get_number_of_discarded_blocks(GtTrack *track)
{
  gt_assert(track);
  return track->discarded_blocks;
}

int gt_track_insert_block(GtTrack *track, GtBlock *block, GtError *err)
{
  GtLine *line = NULL;
  int had_err = 0;

  gt_assert(track && block);
  had_err = get_next_free_line(track, &line, block, err);
  if (!had_err)
  {
    if (line) {
      block = gt_block_ref(block);
      gt_line_insert_block(line, block);
      had_err = gt_line_breaker_register_block(track->lb, line, block, err);
    }
  }
  return had_err;
}

GtStr* gt_track_get_title(const GtTrack *track)
{
  gt_assert(track && track->title);
  return track->title;
}

static unsigned long gt_track_get_number_of_lines(const GtTrack *track)
{
  gt_assert(track);
  return gt_array_size(track->lines);
}

int gt_track_sketch(GtTrack* track, GtCanvas *canvas, GtError *err)
{
  int i = 0, had_err = 0;
  gt_assert(track && canvas);
  had_err = gt_canvas_visit_track_pre(canvas, track, err);
  if (!had_err) {
    for (i = 0; i < gt_array_size(track->lines); i++) {
      had_err = gt_line_sketch(*(GtLine**) gt_array_get(track->lines, i),
                               canvas,
                               err);
      if (had_err)
        break;
    }
  }
  if (!had_err)
    had_err = gt_canvas_visit_track_post(canvas, track, err);
  return had_err;
}

int gt_track_get_height(const GtTrack *track, double *height,
                        const GtStyle *sty, GtError *err)
{
  unsigned long i;
  double track_height = 0,
         bheight = TEXT_SIZE_DEFAULT,
         theight = TEXT_SIZE_DEFAULT,
         tcaptionspace = CAPTION_BAR_SPACE_DEFAULT,
         bcaptionspace = CAPTION_BAR_SPACE_DEFAULT,
         tmp = TRACK_VSPACE_DEFAULT;
  bool show_track_captions = true,
       show_block_captions = true;
  gt_assert(track && sty);
  if (gt_style_get_num(sty, "format", "block_caption_font_size", &bheight,
                       NULL, err) == GT_STYLE_QUERY_ERROR) {
    return -1;
  }
  if (gt_style_get_num(sty, "format", "track_caption_font_size", &theight,
                       NULL, err) == GT_STYLE_QUERY_ERROR) {
    return -1;
  }
  if (gt_style_get_num(sty, "format", "track_caption_space", &tcaptionspace,
                       NULL, err) == GT_STYLE_QUERY_ERROR) {
    return -1;
  }
  if (gt_style_get_num(sty, "format", "block_caption_space", &bcaptionspace,
                       NULL, err) == GT_STYLE_QUERY_ERROR) {
    return -1;
  }
  if (gt_style_get_num(sty, "format", "block_caption_font_size", &bheight,
                       NULL, err) == GT_STYLE_QUERY_ERROR) {
    return -1;
  }
  if (gt_style_get_num(sty, "format", "track_caption_font_size", &theight,
                       NULL, err) == GT_STYLE_QUERY_ERROR) {
    return -1;
  }
  for (i = 0; i < gt_array_size(track->lines); i++)
  {
    double itmp = BAR_VSPACE_DEFAULT,
           tmp = 0.0;
    int rval;
    GtLine *line;

    line = *(GtLine**) gt_array_get(track->lines, i);
    rval = gt_line_get_height(line, &tmp, sty, err);
    if (rval < 0)
      return -1;
    track_height += tmp;
    show_block_captions = true;

    if (gt_style_get_bool(sty,
                          "format","show_block_captions",
                          &show_block_captions,
                          NULL, err) == GT_STYLE_QUERY_ERROR) {
      return -1;
    }

    /* add caption space if necessary */
    if (gt_line_has_captions(line) && show_block_captions)
    {
      track_height += bheight + bcaptionspace;
    }
    /* add vertical spacer */
    if (gt_style_get_num(sty,
                         "format", "bar_vspace",
                         &itmp, NULL, err) == GT_STYLE_QUERY_ERROR) {
      return -1;
    }
    track_height += itmp;
  }

  /* determine display of track captions */
  if (gt_style_get_bool(sty,
                        "format","show_track_captions",
                        &show_track_captions,
                        NULL, err) == GT_STYLE_QUERY_ERROR) {
    return -1;
  }

  /* add track caption height and spacer */
  if (show_track_captions)
    track_height += theight + tcaptionspace;
  if (gt_style_get_num(sty,
                       "format", "track_vspace",
                       &tmp, NULL, err) == GT_STYLE_QUERY_ERROR) {
    return -1;
  }
  track_height += tmp;

  *height = track_height;
  return 0;
}

int gt_track_unit_test(GtError *err)
{
  int had_err = 0;
  GtBlock *b[4];
  GtRange r[4];
  GtTrack *track;
  GtGenomeNode *parent[4], *gn[4];
  GtStr *title;
  double height, tmp;
  GtStyle *sty;
  unsigned long i;
  GtLineBreaker *lb;
  double t_rest = 0,
         l_rest = 0;
  gt_error_check(err);

  title = gt_str_new_cstr("test");

  r[0].start=100UL;  r[0].end=1000UL;
  r[1].start=1001UL; r[1].end=1500UL;
  r[2].start=700UL;  r[2].end=1200UL;
  r[3].start=10UL;   r[3].end=200UL;

  for (i=0;i<4;i++)
  {
    parent[i] = gt_feature_node_new(title, gt_ft_gene, r[i].start, r[i].end,
                                    GT_STRAND_FORWARD);
    gn[i] = gt_feature_node_new(title, gt_ft_exon, r[i].start, r[i].end,
                                GT_STRAND_FORWARD);

    gt_feature_node_add_child((GtFeatureNode*) parent[i],
                              (GtFeatureNode*) gn[i]);

    gt_feature_node_add_attribute((GtFeatureNode*) parent[i], GT_GFF_NAME,
                                  "parent");
    gt_feature_node_add_attribute((GtFeatureNode*) gn[i], GT_GFF_NAME, "child");
  }

  for (i=0;i<4;i++)
  {
    b[i] = gt_block_new();
    gt_block_set_range(b[i], r[i]);
    gt_block_insert_element(b[i], (GtFeatureNode*) parent[i]);
    gt_block_insert_element(b[i], (GtFeatureNode*) gn[i]);
  }

  lb = gt_line_breaker_bases_new();

  sty = gt_style_new(err);

  if (gt_style_get_num(sty, "format", "track_caption_font_size", &tmp,
                   NULL, err) == GT_STYLE_QUERY_NOT_SET) {
    tmp = TEXT_SIZE_DEFAULT;
  }
  t_rest += tmp;
  if (gt_style_get_num(sty, "format", "track_caption_space", &tmp,
                       NULL, err) == GT_STYLE_QUERY_NOT_SET) {
    tmp = CAPTION_BAR_SPACE_DEFAULT;
  }
  t_rest += tmp;
  if (gt_style_get_num(sty, "format", "track_vspace", &tmp,
                       NULL, err) == GT_STYLE_QUERY_NOT_SET) {
    tmp = TRACK_VSPACE_DEFAULT;
  }
  t_rest += tmp;
  if (gt_style_get_num(sty, "format", "bar_vspace", &l_rest,
                       NULL, err) == GT_STYLE_QUERY_NOT_SET) {
    l_rest = BAR_VSPACE_DEFAULT;
  }

  track = gt_track_new(title, GT_UNDEF_ULONG, true, lb);
  gt_ensure(had_err, track);
  gt_ensure(had_err, gt_track_get_title(track) == title);

  gt_ensure(had_err, gt_track_get_number_of_lines(track) == 0);
  gt_ensure(had_err, gt_track_get_height(track, &height, sty, err) == 0);
  gt_ensure(had_err, height == t_rest);
  gt_ensure(had_err, !gt_error_is_set(err));

  gt_ensure(had_err, gt_track_insert_block(track, b[0], err) == 0);
  gt_ensure(had_err, !gt_error_is_set(err));
  gt_ensure(had_err, gt_track_get_number_of_lines(track) == 1);
  gt_ensure(had_err, gt_track_get_height(track, &height, sty, err) == 0);
  gt_ensure(had_err, height == t_rest + l_rest + BAR_HEIGHT_DEFAULT);
  gt_ensure(had_err, !gt_error_is_set(err));

  gt_ensure(had_err, gt_track_insert_block(track, b[1], err) == 0);
  gt_ensure(had_err, !gt_error_is_set(err));
  gt_ensure(had_err, gt_track_get_number_of_lines(track) == 1);
  gt_ensure(had_err, gt_track_get_height(track, &height, sty, err) == 0);
  gt_ensure(had_err, height == t_rest + l_rest + BAR_HEIGHT_DEFAULT);
  gt_ensure(had_err, !gt_error_is_set(err));

  gt_ensure(had_err, gt_track_insert_block(track, b[2], err) == 0);
  gt_ensure(had_err, !gt_error_is_set(err));
  gt_ensure(had_err, gt_track_get_number_of_lines(track) == 2);
  gt_ensure(had_err, gt_track_insert_block(track, b[3], err) == 0);
  gt_ensure(had_err, !gt_error_is_set(err));
  gt_ensure(had_err, gt_track_get_number_of_lines(track) == 2);
  gt_ensure(had_err, gt_track_get_height(track, &height, sty, err) == 0);
  gt_ensure(had_err, height == t_rest + 2*(l_rest + BAR_HEIGHT_DEFAULT));
  gt_ensure(had_err, !gt_error_is_set(err));

  gt_style_set_num(sty, "exon", "bar_height", 42);
  gt_ensure(had_err, gt_track_get_height(track, &height, sty, err) == 0);
  gt_ensure(had_err, height == t_rest + 2*(l_rest+42));
  gt_ensure(had_err, !gt_error_is_set(err));
  gt_style_set_num(sty, "gene", "bar_height", 23);
  gt_ensure(had_err, gt_track_get_height(track, &height, sty, err) == 0);
  gt_ensure(had_err, height == t_rest + 2*(l_rest+42));
  gt_ensure(had_err, !gt_error_is_set(err));
  gt_style_unset(sty, "exon", "bar_height");
  gt_ensure(had_err, gt_track_get_height(track, &height, sty, err) == 0);
  gt_ensure(had_err, height == t_rest + 2*(l_rest+23));
  gt_ensure(had_err, !gt_error_is_set(err));
  gt_style_unset(sty, "gene", "bar_height");
  gt_style_set_num(sty, "format", "bar_height", 99);
  gt_ensure(had_err, gt_track_get_height(track, &height, sty, err) == 0);
  gt_ensure(had_err, height == t_rest + 2*(l_rest+99));
  gt_ensure(had_err, !gt_error_is_set(err));

  gt_ensure(had_err, gt_track_get_number_of_discarded_blocks(track) == 0);

  gt_track_delete(track);
  gt_str_delete(title);
  gt_style_delete(sty);
  for (i=0;i<4;i++)
  {
    gt_block_delete(b[i]);
    gt_genome_node_delete(parent[i]);
  }
  return had_err;
}

void gt_track_delete(GtTrack *track)
{
  unsigned long i;
  if (!track) return;
  if (track->lb)
    gt_line_breaker_delete(track->lb);
  for (i = 0; i < gt_array_size(track->lines); i++)
    gt_line_delete(*(GtLine**) gt_array_get(track->lines, i));
  gt_array_delete(track->lines);
  gt_str_delete(track->title);
  gt_free(track);
}
