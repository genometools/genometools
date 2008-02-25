/*
  Copyright (c) 2007      Sascha Steinbiss <ssteinbiss@zbh.uni-hamburg.de>
  Copyright (c)      2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007-2008 Center for Bioinformatics, University of Hamburg

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
#include "libgtcore/ma.h"
#include "libgtcore/minmax.h"
#include "libgtcore/unused.h"
#include "libgtview/render.h"
#include "libgtview/graphics.h"
#include "libgtview/element.h"
#include "libgtview/track.h"
#include "libgtview/line.h"
#include "libgtview/block.h"

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

struct Render {
  Diagram *dia;
  Config *cfg;
  Graphics *g;
  Range range;
  double y, margins, factor;
  unsigned int width, height, cur_track;
};

Render* render_new(Config *cfg)
{
  Render *r;
  assert(cfg);
  r = ma_malloc(sizeof (Render));
  r->dia = NULL;
  r->cfg = cfg;
  r->margins = config_get_num(r->cfg, "format", "margins", 10);
  return r;
}

/* Calculate the final height of the image to be created. */
static unsigned int render_calculate_height(Render *r)
{
  unsigned int lines = diagram_get_total_lines(r->dia);
  unsigned int height;
  /* obtain line height and spacer from configuration settings */
  unsigned int line_height = ((unsigned int) config_get_num(r->cfg,
                                                            "format",
                                                            "bar_height",
                                                            15
                                                            )) +
                             ((unsigned int) config_get_num(r->cfg,
                                                            "format",
                                                            "bar_vspace",
                                                            10
                                                            )) + 15;

  assert(r);

  /* get total height of all lines */
  height = lines * line_height;
  /* add track caption height and spacer */
  height += diagram_get_number_of_tracks(r->dia)
              * ((config_get_num(r->cfg, "format","track_vspace", 20))+15);
  /* add header space and footer */
  height += 70 + 20;
  if (config_get_verbose(r->cfg))
    fprintf(stderr, "calculated height: %u\n", height);
  return height;
}

double render_convert_point(Render *r, long pos)
{
  return (double) ((r->factor * MAX(0,(pos-(long) r->range.start)))
                      + r->margins);
}

/* Converts base range <node_range> into a pixel range.
   If the range exceeds visibility boundaries, clipping info is set. */
static DrawingRange render_convert_coords(Render *r, Range node_range)
{
  DrawingRange converted_range;
  converted_range.clip = CLIPPED_NONE;

  node_range.end++;

  /* scale coordinates to target image width */
  /* first, check if left side has to be clipped */
  if ((long) node_range.start < (long) r->range.start )
  {
    converted_range.clip = CLIPPED_LEFT;
    converted_range.start = MAX(0, r->margins-5);
  }
  else
  {
    converted_range.start = render_convert_point(r, node_range.start);
  }

  /* then, check right side. */
  if ((long) node_range.end > (long) r->range.end+1)
  {
    converted_range.clip = (converted_range.clip == CLIPPED_LEFT ?
                                                      CLIPPED_BOTH :
                                                      CLIPPED_RIGHT);
    converted_range.end = r->width - r->margins+5;
  }
  else
  {
    converted_range.end = render_convert_point(r, node_range.end);
  }
  return converted_range;
}

static void render_line(Render *r, Line *line)
{
  int i;
  Array *blocks;

  assert(r && line);

  blocks = line_get_blocks(line);
  /* begin drawing block */
  for (i = 0; i < array_size(blocks); i++) {
    Dlistelem *delem;
    Block *block = *(Block**) array_get(blocks, i);
    Dlist *elems = block_get_elements(block);
    Range block_range = block_get_range(block);
    DrawingRange draw_range;
    Color grey;
    double bar_height = config_get_num(r->cfg, "format", "bar_height", 15),
           min_len_block = config_get_num(r->cfg,
                                          "format", "min_len_block", 40);
    const char* caption;
    Strand strand = block_get_strand(block);
    grey.red = grey.green = grey.blue = .85;

    /* draw block caption */
    draw_range = render_convert_coords(r, block_range);
    if (block_caption_is_visible(block)) {
      caption = str_get(block_get_caption(block));
      if (!caption) caption = "";
      graphics_draw_text(r->g,
                         MAX(r->margins, draw_range.start),
                         r->y-6,
                         caption);
    }

    /* do not draw further details in very small blocks */
    if (!block_has_only_one_fullsize_element(block)
         && draw_range.end-draw_range.start < min_len_block)
    {
      int arrow_status = ARROW_NONE;
      GenomeFeatureType btype = block_get_type(block);
      Strand strand = block_get_strand(block);
      if (strand == STRAND_REVERSE || strand == STRAND_BOTH)
        arrow_status = ARROW_LEFT;
      if (strand == STRAND_FORWARD || strand == STRAND_BOTH)
        arrow_status = (arrow_status == ARROW_LEFT ? ARROW_BOTH : ARROW_RIGHT);
      graphics_draw_box(r->g,
                        draw_range.start,
                        r->y,
                        draw_range.end-draw_range.start,
                        bar_height,
                        config_get_color(r->cfg,
                                         genome_feature_type_get_cstr(btype)),
                        arrow_status,
                        config_get_num(r->cfg, "format", "arrow_width", 6),
                        1,
                        config_get_color(r->cfg, "stroke"),
                        true);
      /* draw arrowheads at clipped margins */
      if (draw_range.clip == CLIPPED_LEFT || draw_range.clip == CLIPPED_BOTH)
          graphics_draw_arrowhead(r->g,
                                  r->margins-10,
                                  r->y+((bar_height-8)/2),
                                  grey,
                                  ARROW_LEFT);
      if (draw_range.clip == CLIPPED_RIGHT || draw_range.clip == CLIPPED_BOTH)
          graphics_draw_arrowhead(r->g,
                                  r->width-r->margins+10,
                                  r->y+((bar_height-8)/2),
                                  grey,
                                  ARROW_RIGHT);
      continue;
    }

    /* draw parent block boundaries */
    graphics_draw_dashes(r->g,
                         draw_range.start,
                         r->y,
                         draw_range.end - draw_range.start,
                         config_get_num(r->cfg, "format",
                                        "bar_height", 15),
                         ARROW_NONE,
                         config_get_num(r->cfg, "format",
                                        "arrow_width", 6),
                         config_get_num(r->cfg, "format",
                                        "stroke_width", 1),
                         config_get_color(r->cfg, "stroke"));

    /* draw elements in block */
    for (delem = dlist_first(elems); delem; delem = dlistelem_next(delem)) {
      Element *elem = (Element*) dlistelem_get_data(delem);
      Range elem_range = element_get_range(elem);
      DrawingRange draw_range;
      double elem_start, elem_width, stroke_width;
      Color elem_color;
      int arrow_status = ARROW_NONE;
      const char *style,
                 *type = genome_feature_type_get_cstr(element_get_type(elem));

      /* This shouldn't happen. */
      if (!range_overlap(elem_range, diagram_get_range(r->dia)))
        continue;

      if ((strand == STRAND_REVERSE || strand == STRAND_BOTH)
             /*&& delem == dlist_first(elems)*/)
        arrow_status = ARROW_LEFT;
      if ((strand == STRAND_FORWARD || strand == STRAND_BOTH)
             /*&& dlistelem_next(delem) == NULL*/)
        arrow_status = (arrow_status == ARROW_LEFT ? ARROW_BOTH : ARROW_RIGHT);

      if (config_get_verbose(r->cfg))
        fprintf(stderr, "processing element from %lu to %lu, strand %d\n",
                elem_range.start,
                elem_range.end,
                (int) strand);

      draw_range = render_convert_coords(r, elem_range);
      elem_start = draw_range.start;
      elem_width = draw_range.end - draw_range.start;

      if (config_get_verbose(r->cfg))
        fprintf(stderr, "drawing element from %f to %f, arrow status: %d\n",
                draw_range.start,
                draw_range.end,
                arrow_status);

      /* draw each element according to style set in the config */
      style = config_get_cstr(r->cfg, "feature_styles", type, "box");

      if (element_is_marked(elem)) {
        elem_color = config_get_color(r->cfg, "stroke_marked");
        stroke_width = config_get_num(r->cfg, "format", "stroke_marked_width",
                                      1);
      }
      else {
        elem_color = config_get_color(r->cfg, "stroke");
        stroke_width = config_get_num(r->cfg, "format", "stroke_width", 1);
      }

      if (strcmp(style, "box")==0)
      {
        graphics_draw_box(r->g,
                       elem_start,
                       r->y,
                       elem_width,
                       bar_height,
                       config_get_color(r->cfg, type),
                       arrow_status,
                       config_get_num(r->cfg, "format", "arrow_width", 6),
                       stroke_width,
                       elem_color,
                       false);
      }
      else if (strcmp(style, "caret")==0)
      {
        graphics_draw_caret(r->g,
                       elem_start,
                       r->y,
                       elem_width,
                       bar_height,
                       ARROW_NONE,
                       config_get_num(r->cfg, "format", "arrow_width", 6),
                       stroke_width,
                       elem_color);
      }
      else if (strcmp(style, "dashes")==0)
      {
        graphics_draw_dashes(r->g,
                       elem_start,
                       r->y,
                       elem_width,
                       bar_height,
                       arrow_status,
                       config_get_num(r->cfg, "format", "arrow_width", 6),
                       stroke_width,
                       elem_color);
      }
      else if (strcmp(style, "line")==0)
      {
        graphics_draw_horizontal_line(r->g,
                       elem_start,
                       r->y,
                       elem_width);
      }
      else
      {
         graphics_draw_box(r->g,
                       elem_start,
                       r->y,
                       elem_width,
                       bar_height,
                       config_get_color(r->cfg, type),
                       arrow_status,
                       config_get_num(r->cfg, "format", "arrow_width", 6),
                       stroke_width,
                       elem_color,
                       false);
      }

      /* draw arrowheads at clipped margins */
      if (draw_range.clip == CLIPPED_LEFT || draw_range.clip == CLIPPED_BOTH)
          graphics_draw_arrowhead(r->g,
                                  r->margins-10,
                                  r->y+((bar_height-8)/2),
                                  grey,
                                  ARROW_LEFT);
      if (draw_range.clip == CLIPPED_RIGHT || draw_range.clip == CLIPPED_BOTH)
          graphics_draw_arrowhead(r->g,
                                  r->width-r->margins+10,
                                  r->y+((bar_height-8)/2),
                                  grey,
                                  ARROW_RIGHT);
    }
  }
  /* do not add line spacing after the last line of a track */
  if (i!=array_size(blocks)-1)
    r->y += config_get_num(r->cfg, "format", "bar_height", 15) +
               config_get_num(r->cfg, "format", "bar_vspace", 10) + 15;
}

/* This function disables captions for blocks if they overlap with
   neighboring captions. */
static void mark_caption_collisions(Render *r, Line *line)
{
  int i, j;
  Array *blocks;

  assert(r && line);

  blocks = line_get_blocks(line);
  for (i = 0; i < array_size(blocks)-1; i++) {
    Block *this_block = *(Block**) array_get(blocks, i);
    if (block_caption_is_visible(this_block)) {
      Range block_range = block_get_range(this_block);
      const char *caption;
      Range cur_range;
      caption = str_get(block_get_caption(this_block));
      if (!caption) caption = "";
      cur_range.start = MAX(r->margins,
                            render_convert_point(r, block_range.start));
      cur_range.end   = cur_range.start
                          + graphics_get_text_width(r->g, caption);
      for (j = i-1; j >= 0; j--) {
        Block *left_block = *(Block**) array_get(blocks, j);
        Range chk_range = block_get_range(left_block);
        caption = str_get(block_get_caption(left_block));
        if (!caption) caption = "";
        chk_range.start = render_convert_point(r, chk_range.start);
        chk_range.end   = chk_range.start
                            + graphics_get_text_width(r->g, caption);
        if (range_overlap(chk_range, cur_range))
          block_set_caption_visibility(left_block, false);
      }
      for (j = i+1; j < array_size(blocks); j++) {
        Block *right_block = *(Block**) array_get(blocks, j);
        Range chk_range = block_get_range(right_block);
        caption = str_get(block_get_caption(right_block));
        if (!caption) caption = "";
        chk_range.start = render_convert_point(r, chk_range.start);
        chk_range.end   = chk_range.start
                            + graphics_get_text_width(r->g, caption);
        if (range_overlap(chk_range, cur_range))
          block_set_caption_visibility(right_block, false);
      }
    }
  }
}

static int render_track(void *key, void *value, void *data, UNUSED Error *err)
{
  Render* r = (Render*) data;
  Track* track = (Track*) value;
  Array* lines = track_get_lines(track);
  int i;

  assert(value && key && data);

  if (config_get_verbose(r->cfg))
    fprintf(stderr, "processing track %s\n", (const char*) key);

  /* draw track title */
  graphics_draw_colored_text(r->g,
                             r->margins,
                             r->y-6,
                             config_get_color(r->cfg, "track_title"),
                             str_get(track_get_title(track)));
  r->y += 15;

  /* render each line */
  for (i = 0; i < array_size(lines); i++) {
    Line* line = *(Line**) array_get(lines, i);
    mark_caption_collisions(r, line);
    render_line(r, line);
  }

  /* put track spacer after track, except if at last track */
  if (r->cur_track++ != diagram_get_number_of_tracks(r->dia))
    r->y += config_get_num(r->cfg, "format", "track_vspace", 20);
  return 0;
}

/* Formats a given position number for short display in the ruler. */
static void format_ruler_label(char *txt, long pos)
{
  assert(txt);
  (void) snprintf(txt, BUFSIZ, "%li", pos);
}

/* Renders a ruler with dynamic scale labeling and optional grid. */
static void render_ruler(Render *r)
{
  double step, minorstep, vmajor, vminor;
  long base_length, tick;
  Color rulercol, gridcol;
  char str[BUFSIZ];

  assert(r);

  rulercol.red = rulercol.green = rulercol.blue = .2;
  gridcol.red = gridcol.green = gridcol.blue = .9;

  /* determine range and step of the scale */
  base_length = range_length(r->range);

  /* determine tick steps */
  step = pow(10,ceil(log10(base_length))-1);
  minorstep = step/4.0;

  /* calculate starting positions */
  vminor = (double) (floor(r->range.start / minorstep))*minorstep;
  vmajor = (double) (floor(r->range.start / step))*step;

  /* draw major ticks */
  for (tick = vmajor; tick <= r->range.end; tick += step)
  {
    if (tick < r->range.start) continue;
    graphics_draw_vertical_line(r->g,
                                render_convert_point(r, tick),
                                30,
                                rulercol,
                                10);
    format_ruler_label(str, tick);
    graphics_draw_text_centered(r->g,
                                render_convert_point(r, tick),
                                20,
                                str);
  }
  /* draw minor ticks */
  if (minorstep >= 1)
  {
    for (tick = vminor; tick <= r->range.end; tick += minorstep)
    {
      if (tick < r->range.start) continue;
      if (strcmp(config_get_cstr(r->cfg, "format","show_grid", "no"),
                 "yes") == 0)
        graphics_draw_vertical_line(r->g,
                                    render_convert_point(r, tick),
                                    40,
                                    gridcol,
                                    r->height);
      graphics_draw_vertical_line(r->g,
                                  render_convert_point(r, tick),
                                  35,
                                  rulercol,
                                  5);
    }
  }
  /* draw ruler line */
  graphics_draw_horizontal_line(r->g, r->margins, 40, r->width-2*r->margins);
  /* put 3' and 5' captions at the ends */
  graphics_draw_text_centered(r->g,
                              r->margins-10,
                              45-(graphics_get_text_height(r->g)/2),
                              "5'");
  graphics_draw_text_centered(r->g,
                              r->width-r->margins+10,
                              45-(graphics_get_text_height(r->g)/2),
                              "3'");
}

static void render_prepare_graphics(Render *r, Diagram *dia, unsigned int width)
{
  unsigned int height;

  assert(r && width > 1);

  /* set initial image-specific values */
  r->y = 70;
  r->width = width;
  r->dia = dia;
  r->range = diagram_get_range(dia);
  r->height = height = render_calculate_height(r);

  /* calculate scaling factor */
    r->factor = ((double) r->width
                 -(2*r->margins))
               / range_length(r->range);
  if (config_get_verbose(r->cfg))
     fprintf(stderr, "scaling factor is %f\n",r->factor);

  /* create new Graphics backend */
  r->g = graphics_new(width, height);
  graphics_set_margins(r->g, r->margins, 0);

  /* Add ruler/scale to the image */
  render_ruler(r);

  r->cur_track = 0;
  if (diagram_get_number_of_tracks(r->dia) > 0)
  {
    /* process (render) each track */
    Hashtable *tracks = diagram_get_tracks(r->dia);
    (void) hashtable_foreach_ao(tracks, render_track, r, NULL);
  }
  else if (config_get_verbose(r->cfg))
    fprintf(stderr, "diagram has no tracks!\n");

  if (config_get_verbose(r->cfg))
    fprintf(stderr, "actual used height: %f\n", r->y);
}

int render_to_png(Render *r, Diagram *dia, const char *filename,
                  unsigned int width, Error *err)
{
  int had_err;

  /* prepare graphics object */
  render_prepare_graphics(r, dia, width);

  /* write out result file */
  had_err = graphics_save_to_file(r->g, filename, err);
  graphics_delete(r->g);

  return had_err;
}

void render_to_png_stream(Render *r, Diagram *dia, Str *stream,
                          unsigned int width)
{
  /* prepare graphics object */
  render_prepare_graphics(r, dia, width);

  /* write out stream */
  graphics_save_to_stream(r->g, stream);
  graphics_delete(r->g);
}

void render_delete(Render *r)
{
  if (!r) return;
  ma_free(r);
}
