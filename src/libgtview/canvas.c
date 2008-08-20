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
#include "libgtcore/bittab.h"
#include "libgtcore/ensure.h"
#include "libgtcore/ma.h"
#include "libgtcore/minmax.h"
#include "libgtcore/range.h"
#include "libgtcore/unused.h"
#include "libgtview/block.h"
#include "libgtview/canvas.h"
#include "libgtview/graphics.h"
#include "libgtview/line.h"

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

struct Canvas {
  Range viewrange;
  double factor, y;
  unsigned long width, height, margins;
  Config *cfg;
  bool show_track_captions;
  Bittab *bt;
  Graphics *g;
  GraphicsOutType type;
  ImageInfo *ii;
};

typedef enum
{
  CLIPPED_RIGHT,
  CLIPPED_LEFT,
  CLIPPED_NONE,
  CLIPPED_BOTH
} ClipType;

/* Calculate the final height of the image to be created. */
static unsigned long calculate_height(Canvas *canvas, Diagram *dia)
{
  TracklineInfo lines;
  double tmp;
  unsigned long height;
  unsigned long line_height;
  assert(dia && canvas);

  /* get line information for height calculation */
  lines.total_captionlines = lines.total_lines = 0;
  diagram_get_lineinfo(dia, &lines);

  /* obtain line height and spacer from configuration settings */
  if (config_get_num(canvas->cfg, "format", "bar_height", &tmp))
    line_height = tmp;
  else
    line_height = BAR_HEIGHT_DEFAULT;
  if (config_get_num(canvas->cfg, "format", "bar_vspace", &tmp))
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
    if (config_get_num(canvas->cfg, "format", "track_vspace", &tmp))
      height += diagram_get_number_of_tracks(dia)
                  * (TOY_TEXT_HEIGHT
                      + CAPTION_BAR_SPACE_DEFAULT
                      + tmp);
    else
      height += diagram_get_number_of_tracks(dia)
                  * (TOY_TEXT_HEIGHT
                      + CAPTION_BAR_SPACE_DEFAULT
                      + TRACK_VSPACE_DEFAULT);
  }

  /* add header space and footer */
  height += HEADER_SPACE + FOOTER_SPACE;
  if (config_get_verbose(canvas->cfg))
    fprintf(stderr, "calculated height: %lu\n", height);
  return height;
}

double canvas_get_text_width(Canvas *canvas, const char *text)
{
  assert(canvas);
  if (!text) return 0.0;
  return graphics_get_text_width(canvas->g, text);
}

static double convert_point(Canvas *canvas, long pos)
{
  return (double) ((canvas->factor *
                      MAX(0,(pos-(long) canvas->viewrange.start)))
                      + canvas->margins);
}

/* Converts base range <node_range> into a pixel range.
   If the range exceeds visibility boundaries, clipping info is set. */
DrawingRange canvas_convert_coords(Canvas *canvas, Range node_range)
{
  DrawingRange converted_range;
  converted_range.clip = CLIPPED_NONE;
  node_range.end++;
  /* scale coordinates to target image width */
  /* first, check if left side has to be clipped */
  if ((long) node_range.start < (long) canvas->viewrange.start )
  {
    converted_range.clip = CLIPPED_LEFT;
    converted_range.start = MAX(0, canvas->margins-5);
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
    converted_range.end = (double) canvas->width - canvas->margins+5;
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
static void draw_ruler(Canvas *canvas)
{
  double step, minorstep, vmajor, vminor, margins;
  long base_length, tick;
  Color rulercol, gridcol;
  char str[BUFSIZ];
  bool showgrid;

  assert(canvas);

  margins = canvas->margins;

  if (!(config_get_bool(canvas->cfg, "format","show_grid", &showgrid)))
    showgrid = true;

  rulercol.red = rulercol.green = rulercol.blue = .2;
  gridcol.red = gridcol.green = gridcol.blue = .93;

  /* determine range and step of the scale */
  base_length = range_length(canvas->viewrange);

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

Canvas* canvas_new(Config *cfg, GraphicsOutType type,
                   unsigned long width, ImageInfo *ii)
{
  assert(cfg && width > 0);
  Canvas *canvas;
  canvas = ma_calloc(1, sizeof (Canvas));
  canvas->cfg = cfg;
  canvas->ii = ii;
  canvas->width = width;
  canvas->bt = NULL;
  canvas->type = type;
  canvas->y = 0.5; /* 0.5 displacement to eliminate fuzzy horizontal lines */
  return canvas;
}

unsigned long canvas_get_height(Canvas *canvas)
{
  assert(canvas);
  return canvas->height;
}

int canvas_visit_diagram_pre(Canvas *canvas, Diagram *dia)
{
  double margins;

  assert(canvas && dia);

  if (config_get_num(canvas->cfg, "format", "margins", &margins))
    canvas->margins = margins;
  else
    canvas->margins = MARGINS_DEFAULT;

  if (!config_get_bool(canvas->cfg, "format", "show_track_captions",
                       &canvas->show_track_captions))
    canvas->show_track_captions = true;

  canvas->viewrange = diagram_get_range(dia);
  if (canvas->g)
  {
    graphics_delete(canvas->g);
    canvas->g = NULL;
  }
  canvas->g = graphics_new(canvas->type, canvas->width, 1);

  /* calculate scaling factor */
  canvas->factor = ((double) canvas->width
                     -(2*canvas->margins))
                    / range_length(canvas->viewrange);
  return 0;
}

int canvas_visit_diagram_post(Canvas *canvas, Diagram *dia)
{
  int had_err = 0;

  assert(canvas && dia);

  /* set initial image-specific values */
  canvas->y += HEADER_SPACE;
  canvas->height = calculate_height(canvas, dia);
  if (canvas->ii)
    image_info_set_height(canvas->ii, canvas->height);
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

int canvas_visit_track_pre(Canvas *canvas, Track *track)
{
  int had_err = 0;
  unsigned long exceeded;
  Color color;

  assert(canvas && track);

  config_get_color(canvas->cfg, "format", "track_title_color", &color);

  /* debug */
  if (config_get_verbose(canvas->cfg))
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
      Color red;
      red.red   = 0.7;
      red.green = red.blue  = 0.4;
      if (exceeded == 1)
        msg = "(1 block not shown due to exceeded line limit)";
      else
      {
        msg = "(%lu blocks not shown due to exceeded line limit)";
        snprintf(buf, BUFSIZ, msg, exceeded);
      }
      width = graphics_get_text_width(canvas->g, str_get(track_get_title(track)));
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

int canvas_visit_track_post(Canvas *canvas, Track *track)
{
  double vspace;
  assert(canvas && track);
  /* put track spacer after track */
  if (config_get_num(canvas->cfg, "format", "track_vspace", &vspace))
    canvas->y += vspace;
  else
    canvas->y += TRACK_VSPACE_DEFAULT;
  return 0;
}

int canvas_visit_line_pre(Canvas *canvas, Line *line)
{
  int had_err = 0;
  assert(canvas && line);
  canvas->bt = bittab_new(canvas->width);
  if (line_has_captions(line))
    canvas->y += TOY_TEXT_HEIGHT + CAPTION_BAR_SPACE_DEFAULT;
  return had_err;
}

int canvas_visit_line_post(Canvas *canvas, Line *line)
{
  int had_err = 0;
  double tmp;
  assert(canvas && line);
  if (config_get_num(canvas->cfg, "format", "bar_height", &tmp))
    canvas->y += tmp;
  else
    canvas->y += BAR_HEIGHT_DEFAULT;
  if (config_get_num(canvas->cfg, "format", "bar_vspace", &tmp))
    canvas->y += tmp;
  else
    canvas->y += BAR_VSPACE_DEFAULT;
  bittab_delete(canvas->bt);
  canvas->bt = NULL;
  return had_err;
}

int canvas_visit_block(Canvas *canvas, Block *block)
{
  int had_err = 0, arrow_status = ARROW_NONE;
  Range block_range;
  DrawingRange draw_range;
  Color grey, fillcolor, strokecolor;
  double bar_height, min_len_block, arrow_width, stroke_width;
  const char* caption;
  Strand strand;

  assert(canvas && block);

  grey.red = grey.green = grey.blue = .85;
  strand = block_get_strand(block);
  block_range = block_get_range(block);
  if (!config_get_num(canvas->cfg, "format", "bar_height", &bar_height))
    bar_height = BAR_HEIGHT_DEFAULT;
  if (!config_get_num(canvas->cfg, "format", "min_len_block", &min_len_block))
    min_len_block = MIN_LEN_BLOCK_DEFAULT;
  if (!config_get_num(canvas->cfg, "format", "arrow_width", &arrow_width))
    arrow_width = ARROW_WIDTH_DEFAULT;
  if (!config_get_num(canvas->cfg, "format", "stroke_width", &stroke_width))
    stroke_width = STROKE_WIDTH_DEFAULT;

  if (strand == STRAND_REVERSE || strand == STRAND_BOTH)
    arrow_status = ARROW_LEFT;
  if (strand == STRAND_FORWARD || strand == STRAND_BOTH)
    arrow_status = (arrow_status == ARROW_LEFT ? ARROW_BOTH : ARROW_RIGHT);

  /* draw block caption */
  draw_range = canvas_convert_coords(canvas, block_range);
  if (block_caption_is_visible(block)) {
    caption = str_get(block_get_caption(block));
    if (caption)
    {
      graphics_draw_text(canvas->g,
                         MAX(canvas->margins, draw_range.start),
                         canvas->y -CAPTION_BAR_SPACE_DEFAULT,
                         caption);
    }
  }

  /* do not draw further details in very small blocks */
  if (!block_has_only_one_fullsize_element(block)
       && draw_range.end-draw_range.start < min_len_block)
  {
    GenomeFeatureType *btype = block_get_type(block);
    config_get_color(canvas->cfg, genome_feature_type_get_cstr(btype),
                           "fill", &fillcolor);
    config_get_color(canvas->cfg, genome_feature_type_get_cstr(btype),
                           "stroke", &strokecolor);
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
    /* register coordinates in ImageInfo object if available */
    if (canvas->ii)
    {
      RecMap *rm = recmap_create(draw_range.start, canvas->y,
                                 draw_range.end, canvas->y+bar_height,
                                 block_get_top_level_feature(block));
      image_info_add_recmap(canvas->ii, rm);
      rm->has_omitted_children = true;
    }
    return -1;
  }

  config_get_color(canvas->cfg, "format", "default_stroke_color", &strokecolor);

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

int canvas_visit_element(Canvas *canvas, Element *elem)
{
  int had_err = 0, arrow_status = ARROW_NONE;
  Range elem_range = element_get_range(elem);
  DrawingRange draw_range;
  double elem_start, elem_width, stroke_width, bar_height, arrow_width;
  Color elem_color, grey, fill_color;
  const char *type;
  Str *style;
  Strand strand = element_get_strand(elem);

  assert(canvas && elem);

  /* This shouldn't happen. */
  if (!range_overlap(elem_range, canvas->viewrange))
    return -1;

  type = (char*) genome_feature_type_get_cstr(element_get_type(elem));
  grey.red = grey.green = grey.blue = .85;
  if (!config_get_num(canvas->cfg, "format", "bar_height", &bar_height))
    bar_height = BAR_HEIGHT_DEFAULT;
  if (!config_get_num(canvas->cfg, "format", "arrow_width", &arrow_width))
    arrow_width = ARROW_WIDTH_DEFAULT;

  if ((strand == STRAND_REVERSE || strand == STRAND_BOTH)
         /*&& delem == dlist_first(elems)*/)
    arrow_status = ARROW_LEFT;
  if ((strand == STRAND_FORWARD || strand == STRAND_BOTH)
         /*&& dlistelem_next(delem) == NULL*/)
    arrow_status = (arrow_status == ARROW_LEFT ? ARROW_BOTH : ARROW_RIGHT);

  if (config_get_verbose(canvas->cfg))
    fprintf(stderr, "processing element from %lu to %lu, strand %d\n",
            elem_range.start,
            elem_range.end,
            (int) strand);

  draw_range = canvas_convert_coords(canvas, elem_range);
  elem_start = draw_range.start;
  elem_width = draw_range.end - draw_range.start;

  if (element_is_marked(elem)) {
    config_get_color(canvas->cfg, type, "stroke_marked", &elem_color);
    if (!config_get_num(canvas->cfg, "format", "stroke_marked_width",
                       &stroke_width))
    stroke_width = STROKE_WIDTH_DEFAULT;
  }
  else {
    config_get_color(canvas->cfg, type, "stroke", &elem_color);
    if (!config_get_num(canvas->cfg, "format", "stroke_width", &stroke_width))
    stroke_width = STROKE_WIDTH_DEFAULT;
  }
  config_get_color(canvas->cfg, type, "fill", &fill_color);

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

  /* register coordinates in ImageInfo object if available */
  if (canvas->ii)
  {
    RecMap *rm = recmap_create(elem_start, canvas->y,
                               elem_start+elem_width, canvas->y+bar_height,
                               element_get_node_ref(elem));
    image_info_add_recmap(canvas->ii, rm);
  }

  if (draw_range.end-draw_range.start <= 1.1)
  {
    return had_err;
  }

  if (config_get_verbose(canvas->cfg))
    fprintf(stderr, "drawing element from %f to %f, arrow status: %d\n",
            draw_range.start,
            draw_range.end,
            arrow_status);

  /* draw each element according to style set in the config */
  style = str_new();
  if (!config_get_str(canvas->cfg, type, "style", style))
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

int canvas_to_file(Canvas *canvas, const char *filename, Error *err)
{
  int had_err = 0;
  assert(canvas && filename && err);

  /* write out result file */
  if (canvas->g)
    had_err = graphics_save_to_file(canvas->g, filename, err);
  else
  {
    error_set(err, "No graphics has been created yet!");
    had_err = -1;
  }

  return had_err;
}

int canvas_to_stream(Canvas *canvas, Str *stream)
{
  int had_err = 0;
  assert(canvas && stream);

  /* write out result file */
  if (canvas->g)
    graphics_save_to_stream(canvas->g, stream);

  return had_err;
}

void canvas_delete(Canvas *canvas)
{
  if (!canvas) return;
  if (canvas->g)
    graphics_delete(canvas->g);
  if (canvas->bt)
    bittab_delete(canvas->bt);
  ma_free(canvas);
  canvas = NULL;
}
