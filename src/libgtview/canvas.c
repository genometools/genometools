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

struct Canvas {
  Range viewrange;
  double factor, y;
  unsigned long width, height, margins;
  Config *cfg;
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

typedef struct
{
  double start, end;
  ClipType clip;
} DrawingRange;

/* Calculate the final height of the image to be created. */
static unsigned long calculate_height(Canvas *canvas, Diagram *dia)
{
  TracklineInfo lines;
  unsigned long height;
  unsigned long line_height;
  assert(dia && canvas);

  /* get line information for height calculation */
  lines.total_captionlines = lines.total_lines = 0;
  diagram_get_lineinfo(dia, &lines);

  /* obtain line height and spacer from configuration settings */
  line_height = ((unsigned long) config_get_num(canvas->cfg,
                                               "format",
                                               "bar_height",
                                               15)) +
                ((unsigned long) config_get_num(canvas->cfg,
                                               "format",
                                               "bar_vspace",
                                               10));

  /* get total height of all lines */
  height  = lines.total_lines * line_height;
  height += lines.total_captionlines * 15;
  /* add track caption height and spacer */
  height += diagram_get_number_of_tracks(dia)
            * ((config_get_num(canvas->cfg, "format","track_vspace", 20))+15);
  /* add header space and footer */
  height += 70 + 20;
  if (config_get_verbose(canvas->cfg))
    fprintf(stderr, "calculated height: %lu\n", height);
  return height;
}

static double convert_point(Canvas *canvas, long pos)
{
  return (double) ((canvas->factor *
                      MAX(0,(pos-(long) canvas->viewrange.start)))
                      + config_get_num(canvas->cfg, "format", "margins", 10));
}

/* Converts base range <node_range> into a pixel range.
   If the range exceeds visibility boundaries, clipping info is set. */
DrawingRange convert_coords(Canvas *canvas, Range node_range)
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

  assert(canvas);

  rulercol.red = rulercol.green = rulercol.blue = .2;
  gridcol.red = gridcol.green = gridcol.blue = .9;
  margins = config_get_num(canvas->cfg, "format", "margins", 10);

  /* determine range and step of the scale */
  base_length = range_length(canvas->viewrange);

  /* determine tick steps */
  step = pow(10,ceil(log10(base_length))-1);
  minorstep = step/4.0;

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
      if (strcmp(config_get_cstr(canvas->cfg, "format","show_grid", "no"),
                 "yes") == 0)
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
                                canvas->width-2*config_get_num(canvas->cfg,
                                                       "format",
                                                       "margins", 10));
  /* put 3' and 5' captions at the ends */
  graphics_draw_text_centered(canvas->g,
                              canvas->margins-10,
                              45-(graphics_get_text_height(canvas->g)/2),
                              "5'");
  graphics_draw_text_centered(canvas->g,
                              canvas->width-canvas->margins+10,
                              45-(graphics_get_text_height(canvas->g)/2),
                              "3'");
}

/* This function disables captions for blocks if they overlap with
   neighboring captions. */
void mark_caption_collisions(Canvas *canvas, Line *line)
{
  int i, j;
  Array *blocks;

  assert(canvas && line);

  blocks = line_get_blocks(line);
  for (i = 0; i < array_size(blocks)-1; i++) {
    Block *this_block = *(Block**) array_get(blocks, i);
    if (block_caption_is_visible(this_block)) {
      Range block_range = block_get_range(this_block);
      const char *caption;
      Range cur_range;
      caption = str_get(block_get_caption(this_block));
      if (!caption) caption = "";
      cur_range.start = MAX(canvas->margins,
                            convert_point(canvas, block_range.start));
      cur_range.end   = cur_range.start
                          + graphics_get_text_width(canvas->g, caption);
      for (j = i-1; j >= 0; j--) {
        Block *left_block = *(Block**) array_get(blocks, j);
        Range chk_range = block_get_range(left_block);
        caption = str_get(block_get_caption(left_block));
        if (!caption) caption = "";
        chk_range.start = convert_point(canvas, chk_range.start);
        chk_range.end   = chk_range.start
                            + graphics_get_text_width(canvas->g, caption);
        if (range_overlap(chk_range, cur_range))
          block_set_caption_visibility(left_block, false);
      }
      for (j = i+1; j < array_size(blocks); j++) {
        Block *right_block = *(Block**) array_get(blocks, j);
        Range chk_range = block_get_range(right_block);
        caption = str_get(block_get_caption(right_block));
        if (!caption) caption = "";
        chk_range.start = convert_point(canvas, chk_range.start);
        chk_range.end   = chk_range.start
                            + graphics_get_text_width(canvas->g, caption);
        if (range_overlap(chk_range, cur_range))
          block_set_caption_visibility(right_block, false);
      }
    }
  }
}

Canvas* canvas_new(Config *cfg, GraphicsOutType type,
                   unsigned long width, ImageInfo *ii)
{
  assert(cfg && width > 0);
  Canvas *canvas;
  canvas = ma_malloc(sizeof (Canvas));
  canvas->cfg = cfg;
  canvas->ii = ii;
  canvas->width = width;
  canvas->bt = NULL;
  canvas->type = type;
  return canvas;
}

unsigned long canvas_get_height(Canvas *canvas)
{
  assert(canvas);
  return canvas->height;
}

int canvas_visit_diagram(Canvas *canvas, Diagram *dia)
{
  int had_err = 0;

  assert(canvas && dia);

  canvas->margins = config_get_num(canvas->cfg, "format", "margins", 10);

  /* set initial image-specific values */
  canvas->y = 70;
  canvas->width = canvas->width;
  canvas->viewrange = diagram_get_range(dia);
  canvas->height = calculate_height(canvas, dia);
  if (canvas->ii)
    image_info_set_height(canvas->ii, canvas->height);

  /* calculate scaling factor */
  canvas->factor = ((double) canvas->width
                     -(2*canvas->margins))
                    / range_length(canvas->viewrange);

  if (config_get_verbose(canvas->cfg))
    fprintf(stderr, "scaling factor is %f\n", canvas->factor);

  canvas->g = graphics_new(canvas->type, canvas->width, canvas->height);
  graphics_set_margins(canvas->g, canvas->margins, 0);

  /* Add ruler/scale to the image */
  draw_ruler(canvas);

  return had_err;
}

int canvas_visit_track_pre(Canvas *canvas, Track *track)
{
  int had_err = 0;

  assert(canvas && track);

  /* debug */
  if (config_get_verbose(canvas->cfg))
    fprintf(stderr, "processing track %s\n", str_get(track_get_title(track)));

  /* draw track title */
  graphics_draw_colored_text(canvas->g,
                             config_get_num(canvas->cfg,
                                            "format", "margins", 10),
                             canvas->y,
                             config_get_color(canvas->cfg,
                                              "format",
                                              "track_title_color"),
                             str_get(track_get_title(track)));
  canvas->y += 15;

  return had_err;
}

int canvas_visit_track_post(Canvas *canvas, Track *track)
{
  assert(canvas && track);
  /* put track spacer after track, except if at last track */
  canvas->y += config_get_num(canvas->cfg, "format", "track_vspace", 20);
  return 0;
}

int canvas_visit_line_pre(Canvas *canvas, Line *line)
{
  int had_err = 0;
  assert(canvas && line);
  canvas->bt = bittab_new(canvas->width);
  if (line_has_captions(line))
  {
    mark_caption_collisions(canvas, line);
    canvas->y += 15;
  }
  return had_err;
}

int canvas_visit_line_post(Canvas *canvas, Line *line)
{
  int had_err = 0;
  assert(canvas && line);
  canvas->y += config_get_num(canvas->cfg, "format", "bar_height", 15) +
               config_get_num(canvas->cfg, "format", "bar_vspace", 10);

  bittab_delete(canvas->bt);
  canvas->bt = NULL;
  return had_err;
}

int canvas_visit_block(Canvas *canvas, Block *block)
{
  int had_err = 0, arrow_status = ARROW_NONE;
  Range block_range = block_get_range(block);
  DrawingRange draw_range;
  Color grey;
  double bar_height = config_get_num(canvas->cfg, "format", "bar_height", 15),
         min_len_block = config_get_num(canvas->cfg,
                                        "format", "min_len_block", 40);
  const char* caption;
  Strand strand = block_get_strand(block);
  grey.red = grey.green = grey.blue = .85;

  assert(canvas && block);

  if (strand == STRAND_REVERSE || strand == STRAND_BOTH)
    arrow_status = ARROW_LEFT;
  if (strand == STRAND_FORWARD || strand == STRAND_BOTH)
    arrow_status = (arrow_status == ARROW_LEFT ? ARROW_BOTH : ARROW_RIGHT);

  /* draw block caption */
  draw_range = convert_coords(canvas, block_range);
  if (block_caption_is_visible(block)) {
    caption = str_get(block_get_caption(block));
    if (caption)
    {
      graphics_draw_text(canvas->g,
                         MAX(canvas->margins, draw_range.start),
                         canvas->y-6,
                         caption);
    }
  }

  /* do not draw further details in very small blocks */
  if (!block_has_only_one_fullsize_element(block)
       && draw_range.end-draw_range.start < min_len_block)
  {
    GenomeFeatureType *btype = block_get_type(block);
    graphics_draw_box(canvas->g,
                      draw_range.start,
                      canvas->y,
                      draw_range.end-draw_range.start+1,
                      bar_height,
                      config_get_color(canvas->cfg,
                                       genome_feature_type_get_cstr(btype),
                                       "fill"),
                      arrow_status,
                      config_get_num(canvas->cfg, "format", "arrow_width", 6),
                      1,
                      config_get_color(canvas->cfg,
                                       genome_feature_type_get_cstr(btype),
                                       "stroke"),
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

  /* draw parent block boundaries */
  graphics_draw_dashes(canvas->g,
                       draw_range.start,
                       canvas->y,
                       draw_range.end - draw_range.start,
                       config_get_num(canvas->cfg, "format",
                                      "bar_height", 15),
                       ARROW_NONE,
                       config_get_num(canvas->cfg, "format",
                                      "arrow_width", 6),
                       config_get_num(canvas->cfg, "format",
                                      "stroke_width", 1),
                       config_get_color(canvas->cfg, "format",
                                        "default_stroke_color"));
  return had_err;
}

int canvas_visit_element(Canvas *canvas, Element *elem)
{
  int had_err = 0, arrow_status = ARROW_NONE;
  Range elem_range = element_get_range(elem);
  DrawingRange draw_range;
  double elem_start, elem_width, stroke_width,
         bar_height = config_get_num(canvas->cfg, "format", "bar_height", 15);
  Color elem_color, grey;
  grey.red = grey.green = grey.blue = .85;
  const char *style,
             *type = genome_feature_type_get_cstr(element_get_type(elem));
  Strand strand = element_get_strand(elem);

  assert(canvas && elem);

  /* This shouldn't happen. */
  if (!range_overlap(elem_range, canvas->viewrange))
    return -1;

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

  draw_range = convert_coords(canvas, elem_range);
  elem_start = draw_range.start;
  elem_width = draw_range.end - draw_range.start;

  if (element_is_marked(elem)) {
    elem_color = config_get_color(canvas->cfg, type, "stroke_marked");
    stroke_width = config_get_num(canvas->cfg, "format", "stroke_marked_width",
                                  1);
  }
  else {
    elem_color = config_get_color(canvas->cfg, type, "stroke");
    stroke_width = config_get_num(canvas->cfg, "format", "stroke_width", 1);
  }

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
  style = config_get_cstr(canvas->cfg, type, "style", "box");

  if (strcmp(style, "box")==0)
  {
    graphics_draw_box(canvas->g,
                      elem_start,
                      canvas->y,
                      elem_width,
                      bar_height,
                      config_get_color(canvas->cfg, type, "fill"),
                      arrow_status,
                      config_get_num(canvas->cfg, "format", "arrow_width", 6),
                      stroke_width,
                      elem_color,
                      false);
  }
  else if (strcmp(style, "caret")==0)
  {
    graphics_draw_caret(canvas->g,
                        elem_start,
                        canvas->y,
                        elem_width,
                        bar_height,
                        ARROW_NONE,
                        config_get_num(canvas->cfg, "format", "arrow_width", 6),
                        stroke_width,
                        elem_color);
  }
  else if (strcmp(style, "dashes")==0)
  {
    graphics_draw_dashes(canvas->g,
                         elem_start,
                         canvas->y,
                         elem_width,
                         bar_height,
                         arrow_status,
                         config_get_num(canvas->cfg,
                                        "format", "arrow_width", 6),
                         stroke_width,
                         elem_color);
  }
  else if (strcmp(style, "line")==0)
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
                       config_get_color(canvas->cfg, type, "fill"),
                       arrow_status,
                       config_get_num(canvas->cfg, "format", "arrow_width", 6),
                       stroke_width,
                       elem_color,
                       false);
  }

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
