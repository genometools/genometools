/*
  Copyright (c) 2008 Sascha Steinbiss <ssteinbiss@stud.zbh.uni-hamburg.de>
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
#include <string.h>
#include "core/bittab.h"
#include "core/ensure.h"
#include "core/ma.h"
#include "core/minmax.h"
#include "core/range.h"
#include "core/unused.h"
#include "annotationsketch/block.h"
#include "annotationsketch/canvas.h"
#include "annotationsketch/graphics.h"
#include "annotationsketch/line.h"

#define MARGINS_DEFAULT           10
#define BAR_HEIGHT_DEFAULT        15
#define BAR_VSPACE_DEFAULT        10
#define TOY_TEXT_HEIGHT          8.0
#define TRACK_VSPACE_DEFAULT      20
#define CAPTION_BAR_SPACE_DEFAULT  7
#define MIN_LEN_BLOCK_DEFAULT     30
#define ARROW_WIDTH_DEFAULT        6
#define STROKE_WIDTH_DEFAULT     0.6
#define FONT_SIZE_DEFAULT         10

#define HEADER_SPACE              70
#define FOOTER_SPACE              20

struct GT_Canvas {
  GT_Range viewrange;
  double factor, y, margins;
  unsigned long width, height;
  GT_Style *sty;
  bool show_track_captions;
  Bittab *bt;
  Graphics *g;
  GraphicsOutType type;
  GT_ImageInfo *ii;
};

typedef enum
{
  CLIPPED_RIGHT,
  CLIPPED_LEFT,
  CLIPPED_NONE,
  CLIPPED_BOTH
} ClipType;

/* Calculate the final height of the image to be created. */
static unsigned long calculate_height(GT_Canvas *canvas, GT_Diagram *dia)
{
  TracklineInfo lines;
  double tmp;
  unsigned long height;
  unsigned long line_height;
  assert(dia && canvas);

  /* get line information for height calculation */
  lines.total_captionlines = lines.total_lines = 0;
  gt_diagram_get_lineinfo(dia, &lines);

  /* obtain line height and spacer from style */
  if (gt_style_get_num(canvas->sty, "format", "bar_height", &tmp, NULL))
    line_height = tmp;
  else
    line_height = BAR_HEIGHT_DEFAULT;
  if (gt_style_get_num(canvas->sty, "format", "bar_vspace", &tmp, NULL))
    line_height += tmp;
  else
    line_height += BAR_VSPACE_DEFAULT;

  /* get total height of all lines */
  height  = lines.total_lines * line_height;
  height += lines.total_captionlines * (TOY_TEXT_HEIGHT
                                          + CAPTION_BAR_SPACE_DEFAULT);
  /* add track caption height and spacer */
  if (canvas->show_track_captions)
  {
    if (gt_style_get_num(canvas->sty, "format", "track_vspace", &tmp, NULL))
      height += gt_diagram_get_number_of_tracks(dia)
                  * (TOY_TEXT_HEIGHT
                      + CAPTION_BAR_SPACE_DEFAULT
                      + tmp);
    else
      height += gt_diagram_get_number_of_tracks(dia)
                  * (TOY_TEXT_HEIGHT
                      + CAPTION_BAR_SPACE_DEFAULT
                      + TRACK_VSPACE_DEFAULT);
  }

  /* add header space and footer */
  height += HEADER_SPACE + FOOTER_SPACE;
  if (gt_style_get_verbose(canvas->sty))
    fprintf(stderr, "calculated height: %lu\n", height);
  return height;
}

double gt_canvas_get_text_width(GT_Canvas *canvas, const char *text)
{
  assert(canvas);
  if (!text) return 0.0;
  return graphics_get_text_width(canvas->g, text);
}

static double convert_point(GT_Canvas *canvas, long pos)
{
  return (double) ((canvas->factor *
                      MAX(0,(pos-(long) canvas->viewrange.start)))
                      + canvas->margins);
}

/* Converts base range <node_range> into a pixel range.
   If the range exceeds visibility boundaries, clipping info is set. */
DrawingRange gt_canvas_convert_coords(GT_Canvas *canvas, GT_Range node_range)
{
  DrawingRange converted_range;
  converted_range.clip = CLIPPED_NONE;
  node_range.end++;
  /* scale coordinates to target image width */
  /* first, check if left side has to be clipped */
  if ((long) node_range.start < (long) canvas->viewrange.start )
  {
    converted_range.clip = CLIPPED_LEFT;
    converted_range.start = MAX(0.0, canvas->margins - 5);
  }
  else
  {
    converted_range.start = convert_point(canvas, node_range.start);
  }
  /* then, check right side. */
  if ((long) node_range.end > (long) canvas->viewrange.end+1)
  {
    converted_range.clip = (converted_range.clip == CLIPPED_LEFT ?
                                                      CLIPPED_BOTH :
                                                      CLIPPED_RIGHT);
    converted_range.end = (double) canvas->width - canvas->margins + 5;
  }
  else
  {
    converted_range.end = convert_point(canvas, node_range.end);
  }
  return converted_range;
}

/* Formats a given position number for short display in the ruler. */
static void format_ruler_label(char *txt, unsigned long pos, size_t buflen)
{
  assert(txt);
  double fpos;
  if (pos > 1000000000)
  {
    fpos = (double) pos / 1000000000;
    (void) snprintf(txt, buflen, "%.2fG", fpos);
  }
  else if (pos > 1000000)
  {
    fpos = (double) pos / 1000000;
    (void) snprintf(txt, buflen, "%.2fM", fpos);
  }
  else if (pos > 1000)
  {
    fpos = (double) pos / 1000;
    (void) snprintf(txt, buflen, "%.2fK", fpos);
  } else
    (void) snprintf(txt, buflen, "%li", pos);
}

/* Renders a ruler with dynamic scale labeling and optional grid. */
static void draw_ruler(GT_Canvas *canvas)
{
  double step, minorstep, vmajor, vminor, margins;
  long base_length, tick;
  GT_Color rulercol, gridcol;
  char str[BUFSIZ];
  bool showgrid;

  assert(canvas);

  margins = canvas->margins;

  if (!(gt_style_get_bool(canvas->sty, "format","show_grid", &showgrid, NULL)))
    showgrid = true;

  rulercol.red = rulercol.green = rulercol.blue = .2;
  gridcol.red = gridcol.green = gridcol.blue = .93;

  /* determine range and step of the scale */
  base_length = gt_range_length(canvas->viewrange);

  /* determine tick steps */
  step = pow(10,ceil(log10(base_length))-1);
  minorstep = step/10.0;

  /* calculate starting positions */
  vminor = (double) (floor(canvas->viewrange.start / minorstep))*minorstep;
  vmajor = (double) (floor(canvas->viewrange.start / step))*step;

  /* draw major ticks */
  for (tick = vmajor; tick <= canvas->viewrange.end; tick += step)
  {
    if (tick < canvas->viewrange.start) continue;
    graphics_draw_vertical_line(canvas->g,
                                convert_point(canvas, tick),
                                30,
                                rulercol,
                                10);
    format_ruler_label(str, tick, BUFSIZ);
    graphics_draw_text_centered(canvas->g,
                                convert_point(canvas, tick),
                                20,
                                str);
  }
  /* draw minor ticks */
  if (minorstep >= 1)
  {
    for (tick = vminor; tick <= canvas->viewrange.end; tick += minorstep)
    {
      if (tick < canvas->viewrange.start) continue;
      if (showgrid)
      {
        graphics_draw_vertical_line(canvas->g,
                                    convert_point(canvas, tick),
                                    40,
                                    gridcol,
                                    canvas->height-40-15);
      }
      graphics_draw_vertical_line(canvas->g,
                                  convert_point(canvas, tick),
                                  35,
                                  rulercol,
                                  5);
    }
  }
  /* draw ruler line */
  graphics_draw_horizontal_line(canvas->g,
                                canvas->margins,
                                40,
                                canvas->width-2*margins);
  /* put 3' and 5' captions at the ends */
  graphics_draw_text_centered(canvas->g,
                              canvas->margins-10,
                              45-(TOY_TEXT_HEIGHT/2),
                              "5'");
  graphics_draw_text_centered(canvas->g,
                              canvas->width-canvas->margins+10,
                              45-(TOY_TEXT_HEIGHT/2),
                              "3'");
}

GT_Canvas* gt_canvas_new(GT_Style *sty, GraphicsOutType type,
                   unsigned long width, GT_ImageInfo *ii)
{
  assert(sty && width > 0);
  GT_Canvas *canvas;
  canvas = ma_calloc(1, sizeof (GT_Canvas));
  canvas->sty = sty;
  canvas->ii = ii;
  canvas->width = width;
  canvas->bt = NULL;
  canvas->type = type;
  canvas->y = 0.5; /* 0.5 displacement to eliminate fuzzy horizontal lines */
  return canvas;
}

unsigned long gt_canvas_get_height(GT_Canvas *canvas)
{
  assert(canvas);
  return canvas->height;
}

int gt_canvas_visit_gt_diagram_pre(GT_Canvas *canvas, GT_Diagram *dia)
{
  double margins;

  assert(canvas && dia);

  if (gt_style_get_num(canvas->sty, "format", "margins", &margins, NULL))
    canvas->margins = margins;
  else
    canvas->margins = MARGINS_DEFAULT;

  if (!gt_style_get_bool(canvas->sty, "format", "show_track_captions",
                       &canvas->show_track_captions, NULL))
    canvas->show_track_captions = true;

  canvas->viewrange = gt_diagram_get_range(dia);
  if (canvas->g)
  {
    graphics_delete(canvas->g);
    canvas->g = NULL;
  }
  canvas->g = graphics_new(canvas->type, canvas->width, 1);

  /* calculate scaling factor */
  canvas->factor = ((double) canvas->width
                     -(2*canvas->margins))
                    / gt_range_length(canvas->viewrange);
  return 0;
}

int gt_canvas_visit_gt_diagram_post(GT_Canvas *canvas, GT_Diagram *dia)
{
  int had_err = 0;

  assert(canvas && dia);

  /* set initial image-specific values */
  canvas->y += HEADER_SPACE;
  canvas->height = calculate_height(canvas, dia);
  if (canvas->ii)
    gt_image_info_set_height(canvas->ii, canvas->height);
  if (canvas->g)
  {
    graphics_delete(canvas->g);
    canvas->g = NULL;
  }
  canvas->g = graphics_new(canvas->type, canvas->width, canvas->height);
  graphics_set_margins(canvas->g, canvas->margins, 0);

  /* Add ruler/scale to the image */
  draw_ruler(canvas);

  return had_err;
}

int gt_canvas_visit_track_pre(GT_Canvas *canvas, Track *track)
{
  int had_err = 0;
  unsigned long exceeded;
  GT_Color color;

  assert(canvas && track);

  gt_style_get_color(canvas->sty, "format", "track_title_color", &color, NULL);

  /* debug */
  if (gt_style_get_verbose(canvas->sty))
    fprintf(stderr, "processing track %s\n", str_get(track_get_title(track)));

  if (canvas->show_track_captions)
  {
    /* draw track title */
    graphics_draw_colored_text(canvas->g,
                               canvas->margins,
                               canvas->y,
                               color,
                               str_get(track_get_title(track)));

    /* draw 'line maximum exceeded' message */
    if ((exceeded = track_get_number_of_discarded_blocks(track)) > 0)
    {
      char buf[BUFSIZ];
      const char *msg;
      double width;
      GT_Color red;
      red.red   = 0.7;
      red.green = red.blue  = 0.4;
      if (exceeded == 1)
        msg = "(1 block not shown due to exceeded line limit)";
      else
      {
        msg = "(%lu blocks not shown due to exceeded line limit)";
        snprintf(buf, BUFSIZ, msg, exceeded);
      }
      width = graphics_get_text_width(canvas->g,
                                      str_get(track_get_title(track)));
      graphics_draw_colored_text(canvas->g,
                                 canvas->margins+width+10.0,
                                 canvas->y,
                                 red,
                                 buf);
    }
    canvas->y += TOY_TEXT_HEIGHT + CAPTION_BAR_SPACE_DEFAULT;
  }
  return had_err;
}

int gt_canvas_visit_track_post(GT_Canvas *canvas, UNUSED Track *track)
{
  double vspace;
  assert(canvas && track);
  /* put track spacer after track */
  if (gt_style_get_num(canvas->sty, "format", "track_vspace", &vspace, NULL))
    canvas->y += vspace;
  else
    canvas->y += TRACK_VSPACE_DEFAULT;
  return 0;
}

int gt_canvas_visit_line_pre(GT_Canvas *canvas, Line *line)
{
  int had_err = 0;
  assert(canvas && line);
  canvas->bt = bittab_new(canvas->width);
  if (line_has_captions(line))
    canvas->y += TOY_TEXT_HEIGHT + CAPTION_BAR_SPACE_DEFAULT;
  return had_err;
}

int gt_canvas_visit_line_post(GT_Canvas *canvas, UNUSED Line *line)
{
  int had_err = 0;
  double tmp;
  assert(canvas && line);
  if (gt_style_get_num(canvas->sty, "format", "bar_height", &tmp, NULL))
    canvas->y += tmp;
  else
    canvas->y += BAR_HEIGHT_DEFAULT;
  if (gt_style_get_num(canvas->sty, "format", "bar_vspace", &tmp, NULL))
    canvas->y += tmp;
  else
    canvas->y += BAR_VSPACE_DEFAULT;
  bittab_delete(canvas->bt);
  canvas->bt = NULL;
  return had_err;
}

int gt_canvas_visit_block(GT_Canvas *canvas, GT_Block *block)
{
  int had_err = 0, arrow_status = ARROW_NONE;
  GT_Range block_range;
  DrawingRange draw_range;
  GT_Color grey, fillcolor, strokecolor;
  double bar_height, min_len_block, arrow_width, stroke_width;
  const char* caption;
  GT_Strand strand;

  assert(canvas && block);

  grey.red = grey.green = grey.blue = .85;
  strand = gt_block_get_strand(block);
  block_range = gt_block_get_range(block);
  if (!gt_style_get_num(canvas->sty, "format", "bar_height", &bar_height, NULL))
    bar_height = BAR_HEIGHT_DEFAULT;
  if (!gt_style_get_num(canvas->sty, "format", "min_len_block", &min_len_block,
                     NULL))
    min_len_block = MIN_LEN_BLOCK_DEFAULT;
  if (!gt_style_get_num(canvas->sty, "format", "arrow_width", &arrow_width,
                        NULL)) {
    arrow_width = ARROW_WIDTH_DEFAULT;
  }
  if (!gt_style_get_num(canvas->sty, "format", "stroke_width", &stroke_width,
                     NULL))
    stroke_width = STROKE_WIDTH_DEFAULT;

  if (strand == GT_STRAND_REVERSE || strand == GT_STRAND_BOTH)
    arrow_status = ARROW_LEFT;
  if (strand == GT_STRAND_FORWARD || strand == GT_STRAND_BOTH)
    arrow_status = (arrow_status == ARROW_LEFT ? ARROW_BOTH : ARROW_RIGHT);

  /* draw block caption */
  draw_range = gt_canvas_convert_coords(canvas, block_range);
  if (gt_block_caption_is_visible(block)) {
    caption = str_get(gt_block_get_caption(block));
    if (caption)
    {
      graphics_draw_text(canvas->g,
                         MAX(canvas->margins, draw_range.start),
                         canvas->y -CAPTION_BAR_SPACE_DEFAULT,
                         caption);
    }
  }

  /* do not draw further details in very small blocks */
  if (!gt_block_has_only_one_fullsize_element(block)
       && draw_range.end-draw_range.start < min_len_block)
  {
    GT_GenomeFeatureType *btype = gt_block_get_type(block);
    gt_style_get_color(canvas->sty, gt_genome_feature_type_get_cstr(btype),
                           "fill", &fillcolor,
                           gt_block_get_top_level_feature(block));
    gt_style_get_color(canvas->sty, gt_genome_feature_type_get_cstr(btype),
                           "stroke", &strokecolor,
                           gt_block_get_top_level_feature(block));
    graphics_draw_box(canvas->g,
                      draw_range.start,
                      canvas->y,
                      draw_range.end-draw_range.start+1,
                      bar_height,
                      fillcolor,
                      arrow_status,
                      arrow_width,
                      1,
                      strokecolor,
                      true);
    /* draw arrowheads at clipped margins */
    if (draw_range.clip == CLIPPED_LEFT || draw_range.clip == CLIPPED_BOTH)
        graphics_draw_arrowhead(canvas->g,
                                canvas->margins-10,
                                canvas->y+((bar_height-8)/2),
                                grey,
                                ARROW_LEFT);
    if (draw_range.clip == CLIPPED_RIGHT || draw_range.clip == CLIPPED_BOTH)
        graphics_draw_arrowhead(canvas->g,
                                canvas->width-canvas->margins+10,
                                canvas->y+((bar_height-8)/2),
                                grey,
                                ARROW_RIGHT);
    /* register coordinates in GT_ImageInfo object if available */
    if (canvas->ii)
    {
      GT_RecMap *rm = gt_recmap_new(draw_range.start, canvas->y,
                                    draw_range.end, canvas->y+bar_height,
                                    (GT_GenomeFeature*) /* XXX */
                                    gt_block_get_top_level_feature(block));
      gt_image_info_add_recmap(canvas->ii, rm);
      rm->has_omitted_children = true;
    }
    return -1;
  }

  gt_style_get_color(canvas->sty, "format", "default_stroke_color",
                     &strokecolor, NULL);

  /* draw parent block boundaries */
  graphics_draw_dashes(canvas->g,
                       draw_range.start,
                       canvas->y,
                       draw_range.end - draw_range.start,
                       bar_height,
                       ARROW_NONE,
                       arrow_width,
                       stroke_width,
                       strokecolor);
  return had_err;
}

int gt_canvas_visit_element(GT_Canvas *canvas, Element *elem)
{
  int had_err = 0, arrow_status = ARROW_NONE;
  GT_Range elem_range = element_get_range(elem);
  DrawingRange draw_range;
  double elem_start, elem_width, stroke_width, bar_height, arrow_width;
  GT_Color elem_color, grey, fill_color;
  const char *type;
  GT_Str *style;
  GT_Strand strand = element_get_strand(elem);

  assert(canvas && elem);

  /* This shouldn't happen. */
  if (!gt_range_overlap(elem_range, canvas->viewrange))
    return -1;

  type = gt_genome_feature_type_get_cstr(element_get_type(elem));
  grey.red = grey.green = grey.blue = .85;
  if (!gt_style_get_num(canvas->sty, "format", "bar_height", &bar_height, NULL))
    bar_height = BAR_HEIGHT_DEFAULT;
  if (!gt_style_get_num(canvas->sty, "format", "arrow_width", &arrow_width,
                        NULL)) {
    arrow_width = ARROW_WIDTH_DEFAULT;
  }

  if ((strand == GT_STRAND_REVERSE || strand == GT_STRAND_BOTH)
         /*&& delem == dlist_first(elems)*/)
    arrow_status = ARROW_LEFT;
  if ((strand == GT_STRAND_FORWARD || strand == GT_STRAND_BOTH)
         /*&& dlistelem_next(delem) == NULL*/)
    arrow_status = (arrow_status == ARROW_LEFT ? ARROW_BOTH : ARROW_RIGHT);

  if (gt_style_get_verbose(canvas->sty))
    fprintf(stderr, "processing element from %lu to %lu, strand %d\n",
            elem_range.start,
            elem_range.end,
            (int) strand);

  draw_range = gt_canvas_convert_coords(canvas, elem_range);
  elem_start = draw_range.start;
  elem_width = draw_range.end - draw_range.start;

  if (element_is_marked(elem)) {
    gt_style_get_color(canvas->sty, type, "stroke_marked", &elem_color,
                    element_get_node_ref(elem));
    if (!gt_style_get_num(canvas->sty, "format", "stroke_marked_width",
                       &stroke_width, element_get_node_ref(elem)))
    stroke_width = STROKE_WIDTH_DEFAULT;
  }
  else {
    gt_style_get_color(canvas->sty, type, "stroke", &elem_color,
                    element_get_node_ref(elem));
    if (!gt_style_get_num(canvas->sty, "format", "stroke_width", &stroke_width,
                       element_get_node_ref(elem)))
    stroke_width = STROKE_WIDTH_DEFAULT;
  }
  gt_style_get_color(canvas->sty, type, "fill", &fill_color,
                  element_get_node_ref(elem));

  if (draw_range.end-draw_range.start <= 1.1)
  {
    if (bittab_bit_is_set(canvas->bt, (unsigned long) draw_range.start))
      return had_err;
    graphics_draw_vertical_line(canvas->g,
                                draw_range.start,
                                canvas->y,
                                elem_color,
                                bar_height);
    bittab_set_bit(canvas->bt, (unsigned long) draw_range.start);
  }

  /* register coordinates in GT_ImageInfo object if available */
  if (canvas->ii)
  {
    GT_RecMap *rm = gt_recmap_new(elem_start, canvas->y,
                                  elem_start+elem_width, canvas->y+bar_height,
                                  (GT_GenomeFeature*) /* XXX */
                                  element_get_node_ref(elem));
    gt_image_info_add_recmap(canvas->ii, rm);
  }

  if (draw_range.end-draw_range.start <= 1.1)
  {
    return had_err;
  }

  if (gt_style_get_verbose(canvas->sty))
    fprintf(stderr, "drawing element from %f to %f, arrow status: %d\n",
            draw_range.start,
            draw_range.end,
            arrow_status);

  /* draw each element according to style set in the style */
  style = str_new();
  if (!gt_style_get_str(canvas->sty, type, "style", style,
                     element_get_node_ref(elem)))
    str_set(style, "box");

  if (strcmp(str_get(style), "box")==0)
  {
    graphics_draw_box(canvas->g,
                      elem_start,
                      canvas->y,
                      elem_width,
                      bar_height,
                      fill_color,
                      arrow_status,
                      arrow_width,
                      stroke_width,
                      elem_color,
                      false);
  }
  else if (strcmp(str_get(style), "caret")==0)
  {
    graphics_draw_caret(canvas->g,
                        elem_start,
                        canvas->y,
                        elem_width,
                        bar_height,
                        ARROW_NONE,
                        arrow_width,
                        stroke_width,
                        elem_color);
  }
  else if (strcmp(str_get(style), "dashes")==0)
  {
    graphics_draw_dashes(canvas->g,
                         elem_start,
                         canvas->y,
                         elem_width,
                         bar_height,
                         arrow_status,
                         arrow_width,
                         stroke_width,
                         elem_color);
  }
  else if (strcmp(str_get(style), "line")==0)
  {
    graphics_draw_horizontal_line(canvas->g,
                                  elem_start,
                                  canvas->y,
                                  elem_width);
  }
  else
  {
     graphics_draw_box(canvas->g,
                       elem_start,
                       canvas->y,
                       elem_width,
                       bar_height,
                       fill_color,
                       arrow_status,
                       arrow_width,
                       stroke_width,
                       elem_color,
                       false);
  }
  str_delete(style);

  /* draw arrowheads at clipped margins */
  if (draw_range.clip == CLIPPED_LEFT || draw_range.clip == CLIPPED_BOTH)
      graphics_draw_arrowhead(canvas->g,
                              canvas->margins-10,
                              canvas->y+((bar_height-8)/2),
                              grey,
                              ARROW_LEFT);
  if (draw_range.clip == CLIPPED_RIGHT || draw_range.clip == CLIPPED_BOTH)
      graphics_draw_arrowhead(canvas->g,
                              canvas->width-canvas->margins+10,
                              canvas->y+((bar_height-8)/2),
                              grey,
                              ARROW_RIGHT);
  return had_err;
}

int gt_canvas_to_file(GT_Canvas *canvas, const char *filename, GT_Error *err)
{
  int had_err = 0;
  assert(canvas && filename && err);

  /* write out result file */
  if (canvas->g)
    had_err = graphics_save_to_file(canvas->g, filename, err);
  else
  {
    gt_error_set(err, "No graphics has been created yet!");
    had_err = -1;
  }

  return had_err;
}

int gt_canvas_to_stream(GT_Canvas *canvas, GT_Str *stream)
{
  int had_err = 0;
  assert(canvas && stream);

  /* write out result file */
  if (canvas->g)
    graphics_save_to_stream(canvas->g, stream);

  return had_err;
}

void gt_canvas_delete(GT_Canvas *canvas)
{
  if (!canvas) return;
  if (canvas->g)
    graphics_delete(canvas->g);
  if (canvas->bt)
    bittab_delete(canvas->bt);
  ma_free(canvas);
  canvas = NULL;
}
