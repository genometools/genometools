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

#include <math.h>
#include <string.h>
#include "core/bittab.h"
#include "core/ensure.h"
#include "core/log.h"
#include "core/ma.h"
#include "core/mathsupport.h"
#include "core/minmax.h"
#include "core/range.h"
#include "core/unused_api.h"
#include "core/undef_api.h"
#include "annotationsketch/default_formats.h"
#include "annotationsketch/canvas_members.h"
#include "annotationsketch/canvas_rep.h"
#include "annotationsketch/cliptype.h"
#include "annotationsketch/coords.h"
#include "annotationsketch/graphics_cairo_api.h"
#include "annotationsketch/style.h"

/* this is a Cairo-based implementation of the GtCanvas interface hooks */

int gt_canvas_cairo_visit_layout_pre(GtCanvas *canvas,
                                     GtLayout *layout,
                                     GtError *err)
{
  double head_track_space = HEAD_TRACK_SPACE_DEFAULT;
  /* get displayed range for internal use */
  canvas->pvt->viewrange = gt_layout_get_range(layout);
  if (gt_canvas_draw_ruler(canvas, canvas->pvt->viewrange, err)) {
    return -1;
  }
  if (gt_style_get_num(canvas->pvt->sty, "format", "ruler_space",
                       &head_track_space, NULL,err) == GT_STYLE_QUERY_ERROR) {
    return -1;
  }
  canvas->pvt->y += HEADER_SPACE + head_track_space;
  return 0;
}

int gt_canvas_cairo_visit_layout_post(GT_UNUSED GtCanvas *canvas,
                                      GT_UNUSED GtLayout *layout,
                                      GT_UNUSED GtError *err)
{
  return 0;
}

int gt_canvas_cairo_visit_track_pre(GtCanvas *canvas, GtTrack *track,
                                    GtError *err)
{
  int had_err = 0;
  GtUword exceeded;
  bool show_track_captions = true;

  GtColor color;

  gt_assert(canvas && track);

  if (gt_style_get_bool(canvas->pvt->sty,
                        "format", "show_track_captions",
                        &show_track_captions,
                        NULL, err) == GT_STYLE_QUERY_ERROR) {
    return -1;
  }
  if (gt_style_get_color(canvas->pvt->sty,
                         "format", "track_title_color",
                         &color, NULL, err) == GT_STYLE_QUERY_ERROR) {
    return -1;
  }

  /* debug */
  gt_log_log("processing track %s", gt_str_get(gt_track_get_title(track)));

  if (show_track_captions)
  {
    double theight      = gt_graphics_get_text_height(canvas->pvt->g),
           captionspace = CAPTION_BAR_SPACE_DEFAULT;
    if (gt_style_get_num(canvas->pvt->sty,
                          "format", "track_caption_font_size",
                          &theight, NULL, err) == GT_STYLE_QUERY_ERROR) {
      return -1;
    }
    if (gt_style_get_num(canvas->pvt->sty,
                         "format", "track_caption_space",
                         &captionspace, NULL, err ) == GT_STYLE_QUERY_ERROR) {
      return -1;
    }
    gt_graphics_set_font(canvas->pvt->g,
                         "Sans",
                         SLANT_NORMAL,
                         WEIGHT_NORMAL,
                         theight);
    canvas->pvt->y += theight;
    /* draw track title */
    gt_graphics_draw_colored_text(canvas->pvt->g,
                                  canvas->pvt->margins,
                                  canvas->pvt->y,
                                  color,
                                  gt_str_get(gt_track_get_title(track)));
    /* draw 'line maximum exceeded' message */
    if ((exceeded = gt_track_get_number_of_discarded_blocks(track)) > 0)
    {
      char buf[BUFSIZ];
      const char *msg;
      double width;
      GtColor red;
      red.red   = LINE_EXCEEDED_MSG_R;
      red.green = LINE_EXCEEDED_MSG_G;
      red.blue  = LINE_EXCEEDED_MSG_B;
      red.alpha = 1.0;
      if (exceeded == 1) {
        msg = "(1 block not shown due to exceeded line limit)";
        strncpy(buf, msg, BUFSIZ);
      } else {
        msg = "("GT_WU" blocks not shown due to exceeded line limit)";
        /*@ignore@*/
        snprintf(buf, BUFSIZ, msg, exceeded);
        /*@end@*/
      }
      width = gt_graphics_get_text_width(canvas->pvt->g,
                                         gt_str_get(gt_track_get_title(track)));
      gt_graphics_draw_colored_text(canvas->pvt->g,
                                    canvas->pvt->margins+width+10.0,
                                    canvas->pvt->y,
                                    red,
                                    buf);
    }
    canvas->pvt->y += captionspace;
  }
  return had_err;
}

int gt_canvas_cairo_visit_track_post(GtCanvas *canvas, GT_UNUSED GtTrack *track,
                                     GtError *err)
{
  double vspace = TRACK_VSPACE_DEFAULT;
  gt_assert(canvas && track);
  /* put track spacer after track */
  if (gt_style_get_num(canvas->pvt->sty, "format", "track_vspace", &vspace,
                       NULL, err) == GT_STYLE_QUERY_ERROR) {
    return -1;
  }
  canvas->pvt->y += vspace;
  return 0;
}

int gt_canvas_cairo_visit_line_pre(GtCanvas *canvas, GtLine *line,
                                   GtError *err)
{
  int had_err = 0;
  double lheight, bar_vspace = BAR_VSPACE_DEFAULT,
         theight = gt_graphics_get_text_height(canvas->pvt->g),
         captionspace = CAPTION_BAR_SPACE_DEFAULT;
  bool show_block_captions = true;
  gt_assert(canvas && line);
  if (gt_style_get_num(canvas->pvt->sty, "format", "block_caption_space",
                       &captionspace, NULL, err) == GT_STYLE_QUERY_ERROR) {
    return -1;
  }
  if (gt_style_get_bool(canvas->pvt->sty,
                        "format", "show_block_captions",
                        &show_block_captions,
                        NULL, err) == GT_STYLE_QUERY_ERROR) {
    return -1;
  }
  had_err = gt_line_get_height(line, &lheight, canvas->pvt->sty, err);
  if (!had_err) {
    if (gt_style_get_num(canvas->pvt->sty,
                         "format", "bar_vspace",
                         &bar_vspace, NULL, err) == GT_STYLE_QUERY_ERROR) {
      return -1;
    }
    if (gt_style_get_num(canvas->pvt->sty,
                         "format", "block_caption_font_size",
                         &theight, NULL, err) == GT_STYLE_QUERY_ERROR) {
      return -1;
    }
    if (gt_line_has_captions(line) && show_block_captions)
      canvas->pvt->y +=  theight + captionspace;
    canvas->pvt->bt = gt_bittab_new(canvas->pvt->width);
    gt_bittab_unset(canvas->pvt->bt);
    canvas->pvt->y += lheight/2;
  }
  return had_err;
}

int gt_canvas_cairo_visit_line_post(GtCanvas *canvas, GtLine *line,
                                    GtError *err)
{
  int had_err = 0;
  double lheight,
         bar_vspace = BAR_VSPACE_DEFAULT;
  gt_assert(canvas && line);
  had_err = gt_line_get_height(line, &lheight, canvas->pvt->sty, err);
  if (!had_err) {
    if (gt_style_get_num(canvas->pvt->sty,
                         "format", "bar_vspace",
                         &bar_vspace, NULL, err) == GT_STYLE_QUERY_ERROR) {
      return -1;
    }
    canvas->pvt->y += bar_vspace + lheight/2;
    if (canvas->pvt->bt)
      gt_bittab_delete(canvas->pvt->bt);
    canvas->pvt->bt = NULL;
  }
  return had_err;
}

int gt_canvas_cairo_visit_block(GtCanvas *canvas, GtBlock *block,
                                GtError *err)
{
  int had_err = 0, arrow_status = ARROW_NONE;
  GtRange block_range;
  GtDrawingRange draw_range;
  GtColor grey, fillcolor, strokecolor;
  double block_start,
         block_width,
         bar_height = BAR_HEIGHT_DEFAULT,
         min_len_block = MIN_LEN_BLOCK_DEFAULT,
         arrow_width = ARROW_WIDTH_DEFAULT,
         stroke_width = STROKE_WIDTH_DEFAULT;
  const char *caption, *btype;
  GtStrand strand;
  double maxbarheight;
  btype = gt_block_get_type(block);

  gt_assert(canvas && block);

  grey.red = grey.green = grey.blue = DEFAULT_GREY_TONE;
  grey.alpha = 0.5;
  strand = gt_block_get_strand(block);
  block_range = gt_block_get_range(block);
  if (gt_block_get_max_height(block, &maxbarheight, canvas->pvt->sty, err)) {
    return -1;
  }
  if (gt_style_get_num(canvas->pvt->sty, "format", "bar_height", &bar_height,
                       NULL, err) == GT_STYLE_QUERY_ERROR) {
    return -1;
  }
  if (gt_style_get_num(canvas->pvt->sty, "format", "min_len_block",
                       &min_len_block,NULL, err) == GT_STYLE_QUERY_ERROR) {
    return -1;
  }
  if (gt_style_get_num(canvas->pvt->sty, "format", "arrow_width", &arrow_width,
                       NULL, err) == GT_STYLE_QUERY_ERROR) {
    return -1;
  }
  if (gt_style_get_num(canvas->pvt->sty, "format", "stroke_width",
                       &stroke_width, NULL, err) == GT_STYLE_QUERY_ERROR) {
    return -1;
  }
  if (gt_style_get_color_with_track(canvas->pvt->sty,
                                 btype, "stroke", &strokecolor,
                                 gt_block_get_top_level_feature(block),
                                 gt_track_get_title(canvas->pvt->current_track),
                                 err) == GT_STYLE_QUERY_ERROR) {
    return -1;
  }

  if (strand == GT_STRAND_REVERSE || strand == GT_STRAND_BOTH)
    arrow_status = ARROW_LEFT;
  if (strand == GT_STRAND_FORWARD || strand == GT_STRAND_BOTH)
    arrow_status = (arrow_status == ARROW_LEFT ? ARROW_BOTH : ARROW_RIGHT);

  /* calculate scaled drawing coordinates for this canvas */
  draw_range = gt_coords_calc_generic_range(block_range,
                                            canvas->pvt->viewrange);
  draw_range.start *= canvas->pvt->width-2*canvas->pvt->margins;
  draw_range.end *= canvas->pvt->width-2*canvas->pvt->margins;
  block_start = draw_range.start + canvas->pvt->margins;
  block_width = draw_range.end - draw_range.start;

  /* draw block caption */
  if (gt_block_caption_is_visible(block))
  {
    caption = gt_str_get(gt_block_get_caption(block));
    if (caption)
    {
      double theight = gt_graphics_get_text_height(canvas->pvt->g),
             captionspace = CAPTION_BAR_SPACE_DEFAULT;
      if (gt_style_get_num(canvas->pvt->sty,
                           "format", "block_caption_space",
                           &captionspace, NULL, err) == GT_STYLE_QUERY_ERROR) {
        return -1;
      }
      if (gt_style_get_num(canvas->pvt->sty,
                          "format", "block_caption_font_size",
                          &theight, NULL, err) == GT_STYLE_QUERY_ERROR) {
        return -1;
      }
      gt_graphics_set_font(canvas->pvt->g,
                           "Sans",
                           SLANT_NORMAL,
                           WEIGHT_NORMAL,
                           theight);
      gt_graphics_draw_text_clip(canvas->pvt->g,
                                 block_start,
                                 canvas->pvt->y - bar_height/2
                                                - captionspace,
                                 caption);
    }
  }

  /* do not draw further details in very small blocks */
  if (!gt_block_has_only_one_fullsize_element(block)
       && gt_double_smaller_double(block_width, min_len_block))
  {
    /* optimise drawing for very narrow blocks,
       if <= 1px/pt width, draw as simple lines */
    if (canvas->pvt->bt && draw_range.end-draw_range.start <= 1.1)
    {
      if (!gt_bittab_bit_is_set(canvas->pvt->bt,
                                (GtUword) draw_range.start)) {
        gt_graphics_draw_vertical_line(canvas->pvt->g,
                                       block_start,
                                       canvas->pvt->y - bar_height/2,
                                       strokecolor,
                                       bar_height,
                                       stroke_width);
        gt_bittab_set_bit(canvas->pvt->bt, (GtUword) draw_range.start);
      }
    } else {
      if (gt_style_get_color_with_track(canvas->pvt->sty,
                             btype, "fill",
                             &fillcolor,
                             gt_block_get_top_level_feature(block),
                             gt_track_get_title(canvas->pvt->current_track),
                             err) == GT_STYLE_QUERY_ERROR) {
        return -1;
      }

      gt_graphics_draw_box(canvas->pvt->g,
                           block_start,
                           canvas->pvt->y - bar_height/2,
                           block_width,
                           bar_height,
                           fillcolor,
                           arrow_status,
                           arrow_width,
                           stroke_width,
                           strokecolor,
                           true);

      /* draw arrowheads at clipped margins */
      if (draw_range.clip == CLIPPED_LEFT || draw_range.clip == CLIPPED_BOTH)
          gt_graphics_draw_arrowhead(canvas->pvt->g,
                                     canvas->pvt->margins - 10,
                                     canvas->pvt->y - 4,
                                     grey,
                                     ARROW_LEFT);
      if (draw_range.clip == CLIPPED_RIGHT || draw_range.clip == CLIPPED_BOTH)
          gt_graphics_draw_arrowhead(canvas->pvt->g,
                                     canvas->pvt->width-canvas->pvt->margins
                                       + 10,
                                     canvas->pvt->y - 4,
                                     grey,
                                     ARROW_RIGHT);
    }

    /* register coordinates in GtImageInfo object if available */
    if (canvas->pvt->ii)
    {
      GtRecMap *rm = gt_rec_map_new(block_start,
                                    canvas->pvt->y - bar_height/2,
                                    block_start + block_width,
                                    canvas->pvt->y + bar_height/2,
                                    (GtFeatureNode*) /* XXX */
                                      gt_block_get_top_level_feature(block));
      gt_image_info_add_rec_map(canvas->pvt->ii, rm);
      gt_rec_map_set_omitted_children(rm, true);
    }
    /* signal break */
    return 1;
  }

  /* get stroke color */
  if (gt_style_get_color(canvas->pvt->sty,
                         "format", "default_stroke_color",
                         &strokecolor, NULL, err) == GT_STYLE_QUERY_ERROR) {
    return -1;
  }

  /* draw parent block boundaries */
  gt_graphics_draw_dashes(canvas->pvt->g,
                          block_start,
                          canvas->pvt->y - bar_height/2,
                          block_width,
                          bar_height,
                          ARROW_NONE,
                          arrow_width,
                          stroke_width,
                          strokecolor);
  return had_err;
}

int gt_canvas_cairo_visit_element(GtCanvas *canvas, GtElement *elem,
                                  GtError *err)
{
  int had_err = 0, arrow_status = ARROW_NONE;
  GtRange elem_range = gt_element_get_range(elem);
  GtDrawingRange draw_range;
  double elem_start = GT_UNDEF_DOUBLE,
         elem_width = GT_UNDEF_DOUBLE,
         stroke_width = STROKE_WIDTH_DEFAULT,
         bar_height = BAR_HEIGHT_DEFAULT,
         arrow_width = ARROW_WIDTH_DEFAULT;
  GtColor elem_color, grey, fill_color;
  const char *type;
  GtStyleQueryStatus rval;
  GtStr *style;
  GtStrand strand = gt_element_get_strand(elem);

  gt_assert(canvas && elem);

  /* This shouldn't happen. */
  if (!gt_range_overlap(&elem_range, &canvas->pvt->viewrange))
    return -1;

  type = gt_element_get_type(elem);
  grey.red = grey.green = grey.blue = .85;
  grey.alpha = 0.5;

  /* get default or image-wide bar height */
  if (gt_style_get_num(canvas->pvt->sty, "format", "bar_height", &bar_height,
                       NULL, err) == GT_STYLE_QUERY_ERROR) {
    return -1;
  }
  /* try to get type-specific bar height */
  if (gt_style_get_num_with_track(canvas->pvt->sty, type, "bar_height",
                       &bar_height,
                       gt_element_get_node_ref(elem),
                       gt_track_get_title(canvas->pvt->current_track),
                       err) == GT_STYLE_QUERY_ERROR) {
    return -1;
  }
  /* get default or image-wide arrow width */
  if (gt_style_get_num(canvas->pvt->sty, "format", "arrow_width", &arrow_width,
                        NULL, err)== GT_STYLE_QUERY_ERROR) {
    return -1;
  }
  /* try to get type-specific arrow width */
  if (gt_style_get_num_with_track(canvas->pvt->sty, type, "arrow_width",
                       &arrow_width,
                       gt_element_get_node_ref(elem),
                       gt_track_get_title(canvas->pvt->current_track),
                       err) == GT_STYLE_QUERY_ERROR) {
    return -1;
  }

  if ((strand == GT_STRAND_REVERSE || strand == GT_STRAND_BOTH)
         /*&& delem == gt_dlist_first(elems)*/)
    arrow_status = ARROW_LEFT;
  if ((strand == GT_STRAND_FORWARD || strand == GT_STRAND_BOTH)
         /*&& gt_dlistelem_next(delem) == NULL*/)
    arrow_status = (arrow_status == ARROW_LEFT ? ARROW_BOTH : ARROW_RIGHT);

  gt_log_log("processing element from "GT_WU" to "GT_WU", strand %d\n",
             elem_range.start, elem_range.end, (int) strand);

  draw_range = gt_coords_calc_generic_range(elem_range, canvas->pvt->viewrange);
  draw_range.start *= (canvas->pvt->width-2*canvas->pvt->margins);
  draw_range.end *= (canvas->pvt->width-2*canvas->pvt->margins);
  elem_start = draw_range.start + canvas->pvt->margins;
  elem_width = draw_range.end - draw_range.start;
  gt_assert(elem_start != GT_UNDEF_DOUBLE && elem_width != GT_UNDEF_DOUBLE);

  if (gt_element_is_marked(elem)) {
    if (gt_style_get_color_with_track(canvas->pvt->sty, type, "stroke_marked",
                           &elem_color, gt_element_get_node_ref(elem),
                           gt_track_get_title(canvas->pvt->current_track),
                           err) == GT_STYLE_QUERY_ERROR) {
      return -1;
    }
    if (gt_style_get_num_with_track(canvas->pvt->sty, "format",
                          "stroke_marked_width",
                          &stroke_width, gt_element_get_node_ref(elem),
                          gt_track_get_title(canvas->pvt->current_track),
                          err) == GT_STYLE_QUERY_ERROR) {
    return -1;
    }
  }
  else {
    if (gt_style_get_color_with_track(canvas->pvt->sty, type, "stroke",
                              &elem_color,
                              gt_element_get_node_ref(elem),
                              gt_track_get_title(canvas->pvt->current_track),
                              err) == GT_STYLE_QUERY_ERROR) {
      return -1;
    }
    if (gt_style_get_num_with_track(canvas->pvt->sty, "format", "stroke_width",
                         &stroke_width,
                         gt_element_get_node_ref(elem),
                         gt_track_get_title(canvas->pvt->current_track),
                         err) == GT_STYLE_QUERY_ERROR) {
      return -1;
    }
    if (gt_style_get_num_with_track(canvas->pvt->sty, type, "stroke_width",
                         &stroke_width,
                         gt_element_get_node_ref(elem),
                         gt_track_get_title(canvas->pvt->current_track),
                         err) == GT_STYLE_QUERY_ERROR) {
      return -1;
    }
  }
  if (gt_style_get_color_with_track(canvas->pvt->sty, type, "fill",
                                 &fill_color,
                                 gt_element_get_node_ref(elem),
                                 gt_track_get_title(canvas->pvt->current_track),
                                 err) == GT_STYLE_QUERY_ERROR) {
    return -1;
  }

  if (canvas->pvt->bt &&
          gt_double_smaller_double(draw_range.end-draw_range.start, 1.1))
  {
    if ((GtUword) draw_range.start > gt_bittab_size(canvas->pvt->bt))
      return had_err;
    if (gt_bittab_bit_is_set(canvas->pvt->bt, (GtUword) draw_range.start))
      return had_err;
    gt_graphics_draw_vertical_line(canvas->pvt->g,
                                   elem_start,
                                   canvas->pvt->y - bar_height/2,
                                   elem_color,
                                   bar_height,
                                   stroke_width);
    gt_bittab_set_bit(canvas->pvt->bt, (GtUword) draw_range.start);
  }

  /* register coordinates in GtImageInfo object if available */
  if (canvas->pvt->ii)
  {
    GtRecMap *rm = gt_rec_map_new(elem_start, canvas->pvt->y - bar_height/2,
                                  elem_start+elem_width,
                                  canvas->pvt->y+bar_height/2,
                                  (GtFeatureNode*)
                                    gt_element_get_node_ref(elem));
    gt_image_info_add_rec_map(canvas->pvt->ii, rm);
    }

  if (canvas->pvt->bt && draw_range.end-draw_range.start <= 1.1)
  {
    return had_err;
  }

  gt_log_log("drawing element from %f to %f, arrow status: %d",
             draw_range.start, draw_range.end, arrow_status);

  /* draw each element according to style set in the style */
  style = gt_str_new();
  rval = gt_style_get_str_with_track(canvas->pvt->sty, type, "style", style,
                          gt_element_get_node_ref(elem),
                          gt_track_get_title(canvas->pvt->current_track),
                          err);
  switch (rval) {
    case GT_STYLE_QUERY_NOT_SET:
      gt_str_set(style, "box");    /* default style */
      break;
    case GT_STYLE_QUERY_ERROR:
      gt_str_delete(style);
      gt_assert(gt_error_is_set(err));
      return -1;
    default:
      break;
  }

  if (strcmp(gt_str_get(style), "box") == 0)
  {
    gt_graphics_draw_box(canvas->pvt->g,
                         elem_start,
                         canvas->pvt->y - bar_height/2,
                         elem_width,
                         bar_height,
                         fill_color,
                         arrow_status,
                         arrow_width,
                         stroke_width,
                         elem_color,
                         false);
  }
  else if (strcmp(gt_str_get(style), "rectangle") == 0)
  {
    gt_graphics_draw_box(canvas->pvt->g,
                         elem_start,
                         canvas->pvt->y - bar_height/2,
                         elem_width,
                         bar_height,
                         fill_color,
                         ARROW_NONE,
                         arrow_width,
                         stroke_width,
                         elem_color,
                         false);
  }
  else if (strcmp(gt_str_get(style), "caret") == 0)
  {
    gt_graphics_draw_caret(canvas->pvt->g,
                           elem_start,
                           canvas->pvt->y - bar_height/2,
                           elem_width,
                           bar_height,
                           ARROW_NONE,
                           arrow_width,
                           stroke_width,
                           elem_color);
  }
  else if (strcmp(gt_str_get(style), "dashes") == 0)
  {
    gt_graphics_draw_dashes(canvas->pvt->g,
                            elem_start,
                            canvas->pvt->y - bar_height/2,
                            elem_width,
                            bar_height,
                            arrow_status,
                            arrow_width,
                            stroke_width,
                            elem_color);
  }
  else if (strcmp(gt_str_get(style), "line") == 0)
  {
    gt_graphics_draw_horizontal_line(canvas->pvt->g,
                                     elem_start,
                                     canvas->pvt->y - bar_height/2,
                                     elem_color,
                                     elem_width,
                                     1.0);
  }
  else
  {
     gt_graphics_draw_box(canvas->pvt->g,
                          elem_start,
                          canvas->pvt->y - bar_height/2,
                          elem_width,
                          bar_height,
                          fill_color,
                          arrow_status,
                          arrow_width,
                          stroke_width,
                          elem_color,
                          false);
  }
  gt_str_delete(style);

  /* draw arrowheads at clipped margins */
  if (draw_range.clip == CLIPPED_LEFT || draw_range.clip == CLIPPED_BOTH)
      gt_graphics_draw_arrowhead(canvas->pvt->g,
                                 canvas->pvt->margins - 10,
                                 canvas->pvt->y - 4,
                                 grey,
                                 ARROW_LEFT);
  if (draw_range.clip == CLIPPED_RIGHT || draw_range.clip == CLIPPED_BOTH)
      gt_graphics_draw_arrowhead(canvas->pvt->g,
                                 canvas->pvt->width-canvas->pvt->margins + 10,
                                 canvas->pvt->y - 4,
                                 grey,
                                 ARROW_RIGHT);
  return had_err;
}

int gt_canvas_cairo_visit_custom_track(GtCanvas *canvas,
                                       GtCustomTrack *ct,
                                       GtError *err)
{
  bool show_track_captions = true;
  double space;
  int had_err = 0;
  GtColor color;
  gt_assert(canvas && ct);

  if (gt_style_get_bool(canvas->pvt->sty,
                        "format", "show_track_captions",
                        &show_track_captions,
                        NULL, err) == GT_STYLE_QUERY_ERROR) {
    return -1;
  }
  if (gt_style_get_color(canvas->pvt->sty, "format", "track_title_color",
                         &color, NULL, err) == GT_STYLE_QUERY_ERROR) {
    return -1;
  }

  if (show_track_captions)
  {
    double theight = gt_graphics_get_text_height(canvas->pvt->g),
           captionspace = CAPTION_BAR_SPACE_DEFAULT;
    if (gt_style_get_num(canvas->pvt->sty,
                         "format", "track_caption_font_size",
                         &theight, NULL, err) == GT_STYLE_QUERY_ERROR) {
      return -1;
    }
    /* draw track title */
    gt_graphics_set_font(canvas->pvt->g,
                           "Sans",
                           SLANT_NORMAL,
                           WEIGHT_NORMAL,
                           theight);
    gt_graphics_draw_colored_text(canvas->pvt->g,
                                  canvas->pvt->margins,
                                  canvas->pvt->y,
                                  color,
                                  gt_custom_track_get_title(ct));
    if (gt_style_get_num(canvas->pvt->sty,
                         "format", "track_caption_space",
                         &captionspace, NULL, err) == GT_STYLE_QUERY_ERROR) {
      return -1;
    }
    canvas->pvt->y += theight + captionspace;
  }

  /* call rendering function */
  had_err = gt_custom_track_render(ct,
                                   canvas->pvt->g,
                                   canvas->pvt->y,
                                   canvas->pvt->viewrange,
                                   canvas->pvt->sty,
                                   err);
  canvas->pvt->y += gt_custom_track_get_height(ct);

  /* put spacers after track */
  space = TRACK_VSPACE_DEFAULT;
  if (gt_style_get_num(canvas->pvt->sty, "format", "track_vspace", &space,
                       NULL, err) == GT_STYLE_QUERY_ERROR) {
    return -1;
  }
  canvas->pvt->y += space;

  return had_err;
}

/* Renders a ruler with dynamic scale labeling and optional grid. */
int gt_canvas_cairo_draw_ruler(GtCanvas *canvas, GtRange viewrange,
                                GtError *err)
{
  double step, minorstep, vmajor, vminor, offset,
         theight = gt_graphics_get_text_height(canvas->pvt->g);
  GtWord base_length, tick;
  GtColor rulercol, gridcol;
  GtStr *left_str, *right_str, *unit;
  char str[BUFSIZ];
  GtStyleQueryStatus rval;
  bool showgrid = true;
  gt_assert(canvas);

  if (gt_style_get_bool(canvas->pvt->sty, "format", "show_grid", &showgrid,
                        NULL, err) == GT_STYLE_QUERY_ERROR) {
    return -1;
  }
  if (gt_style_get_num(canvas->pvt->sty, "format", "ruler_font_size",
                       &theight, NULL, err) == GT_STYLE_QUERY_ERROR) {
    return -1;
  }

  /* get offset value from style, default: 0 */
  offset = 0;
  if (gt_style_get_num(canvas->pvt->sty, "format", "ruler_offset",
                       &offset, NULL, err) == GT_STYLE_QUERY_ERROR) {
    return -1;
  }

  /* get unit value from style, default: empty */
  unit = gt_str_new();
  if (gt_style_get_str(canvas->pvt->sty,
                       "format", "unit",
                       unit, NULL, err) == GT_STYLE_QUERY_ERROR) {
    gt_str_delete(unit);
    return -1;
  }

  /* get additional description texts from style */
  left_str = gt_str_new();
  rval = gt_style_get_str(canvas->pvt->sty, "format",
                          "ruler_left_text", left_str, NULL, err);
  switch (rval) {
    case GT_STYLE_QUERY_NOT_SET:
      gt_str_append_cstr(left_str, FIVE_PRIME_STRING);
      break;
    case GT_STYLE_QUERY_ERROR:
      gt_str_delete(unit);
      gt_str_delete(left_str);
      return -1;
      break;  /* shouldn't reach this */
    default:
      break;
  }
  right_str = gt_str_new();
  rval = gt_style_get_str(canvas->pvt->sty, "format",
                        "ruler_right_text", right_str, NULL, err);
  switch (rval) {
    case GT_STYLE_QUERY_NOT_SET:
      gt_str_append_cstr(right_str, THREE_PRIME_STRING);
      break;
    case GT_STYLE_QUERY_ERROR:
      gt_str_delete(unit);
      gt_str_delete(left_str);
      gt_str_delete(right_str);
      return -1;
      break;  /* shouldn't reach this */
    default:
      break;
  }

  /* reset font to default */
  gt_graphics_set_font(canvas->pvt->g,
                       "Sans",
                       SLANT_NORMAL,
                       WEIGHT_NORMAL,
                       theight);

  rulercol.red = rulercol.green = rulercol.blue = RULER_GREY;
  rulercol.alpha = 1.0;
  gridcol.red = gridcol.green = gridcol.blue = GRID_GREY;
  gridcol.alpha = 1.0;

  /* determine range and step of the scale */
  base_length = gt_range_length(&viewrange);

  /* determine tick steps */
  step = pow(10,ceil(log10(base_length))-1);
  minorstep = step/10.0;

  /* calculate starting positions */
  vminor = (double) (floor(viewrange.start / minorstep))*minorstep;
  vmajor = (double) (floor(viewrange.start / step))*step;

  /* draw major ticks */
  for (tick = vmajor; tick <= viewrange.end; tick += step)
  {
    double drawtick = (gt_coords_convert_point(viewrange, tick)
                       * (canvas->pvt->width-2*canvas->pvt->margins))
                       + canvas->pvt->margins;
    if (tick < viewrange.start) continue;
    gt_graphics_draw_vertical_line(canvas->pvt->g,
                                   drawtick,
                                   canvas->pvt->y + 30,
                                   rulercol,
                                   10,
                                   1.0);
    gt_format_ruler_label(str, tick + (GtWord)offset, gt_str_get(unit), BUFSIZ);
    gt_graphics_draw_text_centered(canvas->pvt->g,
                                   drawtick,
                                   canvas->pvt->y + 20,
                                   str);
  }
  /* draw minor ticks */
  if (minorstep >= 1)
  {
    for (tick = vminor; tick <= viewrange.end; tick += minorstep)
    {
      double drawtick;
      if (tick < viewrange.start)
        continue;
      drawtick = (gt_coords_convert_point(viewrange, tick)
                    * (canvas->pvt->width-2*canvas->pvt->margins))
                  + canvas->pvt->margins;
      if (showgrid)
      {
        gt_graphics_draw_vertical_line(canvas->pvt->g,
                                       drawtick,
                                       canvas->pvt->y + 40,
                                       gridcol,
                                       canvas->pvt->height - 40 - 15,
                                       1.0);
      }
      gt_graphics_draw_vertical_line(canvas->pvt->g,
                                     drawtick,
                                     canvas->pvt->y + 35,
                                     rulercol,
                                     5,
                                     1.0);
    }
  }
  /* draw ruler line */
  gt_graphics_draw_horizontal_line(canvas->pvt->g,
                                   canvas->pvt->margins,
                                   canvas->pvt->y + 40,
                                   rulercol,
                                   canvas->pvt->width - 2
                                     * canvas->pvt->margins,
                                   1.25);

    gt_graphics_draw_text_right(canvas->pvt->g,
                              canvas->pvt->margins - 10,
                              canvas->pvt->y + 39 + (theight/2),
                              gt_str_get(left_str));
  gt_graphics_draw_text(canvas->pvt->g,
                        canvas->pvt->width - canvas->pvt->margins + 10,
                        canvas->pvt->y + 39 + (theight/2),
                        gt_str_get(right_str));

  gt_str_delete(unit);
  gt_str_delete(left_str);
  gt_str_delete(right_str);

  return 0;
}
