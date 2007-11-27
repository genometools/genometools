/*
  Copyright (c) 2007 Sascha Steinbiss <ssteinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>,
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

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

#include <cairo.h>
#include <cairo-pdf.h>
#include <cairo-ps.h>
#include "libgtcore/ma.h"
#include "libgtcore/minmax.h"
#include "libgtcore/fileutils.h"
#include "libgtview/graphics.h"

typedef enum {
  PDF,
  PNG,
  PS
} CairoGraphicsType;

struct Graphics {
  cairo_t *cr;
  cairo_surface_t *surf;
  CairoGraphicsType type;
  double margin_x, margin_y, height, width;
  const char* filename;
};

void graphics_initialize(Graphics *g, const char *filename, unsigned int width,
                         unsigned int height, Env *env)
{
  g->filename = filename;
  g->cr = cairo_create(g->surf);
  assert(cairo_status(g->cr) == CAIRO_STATUS_SUCCESS);
  g->width = width;
  g->height = height;
  g->margin_x = g->margin_y = 20;
  cairo_set_source_rgb(g->cr, 1, 1, 1);
  cairo_set_operator(g->cr, CAIRO_OPERATOR_SOURCE);
  cairo_paint(g->cr);
  cairo_set_line_join(g->cr, CAIRO_LINE_JOIN_ROUND);
  cairo_set_line_cap(g->cr, CAIRO_LINE_CAP_ROUND);
}

Graphics* graphics_new_png(const char *filename, unsigned int width,
                           unsigned int height, Env *env)
{
  Graphics *g = ma_malloc(sizeof (Graphics));
  g->type = PNG;
  g->surf = cairo_image_surface_create(CAIRO_FORMAT_RGB24, width, height);
  graphics_initialize(g, filename, width, height, env);
  return g;
}

Graphics* graphics_new_pdf(const char *filename, unsigned int width,
                           unsigned int height, Env *env)
{
  Graphics *g = ma_malloc(sizeof (Graphics));
  g->type = PDF;
  g->surf = cairo_pdf_surface_create(filename, width, height);
  graphics_initialize(g, filename, width, height, env);
  return g;
}

Graphics* graphics_new_ps(const char *filename, unsigned int width,
                          unsigned int height, Env *env)
{
  Graphics *g = ma_malloc(sizeof (Graphics));
  g->type = PS;
  g->surf = cairo_ps_surface_create(filename, width, height);
  graphics_initialize(g, filename, width, height, env);
  return g;
}

void graphics_draw_text(Graphics *g, double x, double y, const char *text)
{
  assert(g && text);
  cairo_set_source_rgb(g->cr, 0, 0, 0);
  cairo_move_to(g->cr, x, y);
  cairo_show_text(g->cr, text);
}

void graphics_draw_text_centered(Graphics *g, double x, double y,
                                 const char *text)
{
  cairo_text_extents_t ext;
  assert(g && text);
  cairo_set_source_rgb(g->cr, 0, 0, 0);
  /* get text extents */
  cairo_text_extents(g->cr, text, &ext);
  /* draw text w/ its center at the given coords */
  cairo_move_to(g->cr, x-(ext.width/2)-1, y);
  cairo_show_text(g->cr, text);
}

void graphics_draw_text_right(Graphics *g, double x, double y, const char *text)
{
  cairo_text_extents_t ext;
  assert(g && text);
  cairo_set_source_rgb(g->cr, 0, 0, 0);
  /* get text extents */
  cairo_text_extents(g->cr, text, &ext);
  /* draw text w/ its right end at the given coords */
  cairo_move_to(g->cr, x-(ext.width)-1, y);
  cairo_show_text(g->cr, text);
}

void graphics_draw_colored_text(Graphics *g, double x, double y, Color color,
                                const char *text)
{
  assert(g && text);
  cairo_set_source_rgb(g->cr, color.red, color.green, color.blue);
  cairo_move_to(g->cr, x, y);
  cairo_show_text(g->cr, text);
}

double graphics_get_text_height(Graphics *g)
{
  cairo_text_extents_t ext;
  assert(g);
  /* get text extents */
  cairo_text_extents(g->cr, "A", &ext);
  return ext.height;
}

void graphics_set_margins(Graphics *g, double margin_x, double margin_y)
{
  assert(g);
  g->margin_x = margin_x;
  g->margin_y = margin_y;
}

void graphics_draw_horizontal_line(Graphics *g, double x, double y,
                                   double width)
{
  assert(g);
  cairo_move_to(g->cr, x, y);
  cairo_rel_line_to(g->cr, width, 0);
  cairo_stroke(g->cr);
}

void graphics_draw_vertical_line(Graphics *g, double x, double y,
                                 Color color, double length)
{
  assert(g);
  cairo_move_to(g->cr, x, y);
  cairo_set_source_rgb(g->cr, color.red, color.green, color.blue);
  cairo_rel_line_to(g->cr, 0, length);
  cairo_stroke(g->cr);
}

double graphics_get_text_width(Graphics *g, const char* text)
{
  cairo_text_extents_t ext;
  assert(g);
  /* get text extents */
  cairo_text_extents(g->cr, text, &ext);
  return ext.width;
}

void graphics_draw_box(Graphics *g, double x, double y, double width,
                       double height, Color fill_color,
                       ArrowStatus arrow_status, double arrow_width,
                       double stroke_width, Color stroke_color)
{
  assert(g);
  /* save cairo context */
  cairo_save(g->cr);
  cairo_rectangle(g->cr, g->margin_x, g->margin_y,
                  g->width-2*g->margin_x, g->height-2*g->margin_y);
  cairo_clip(g->cr);
  /* construct shape of the box or arrow */
  switch (arrow_status)
  {
    case ARROW_RIGHT:
      cairo_move_to(g->cr, x, y);
      if (width - arrow_width > 0)
        cairo_rel_line_to(g->cr, width - arrow_width, 0);
      cairo_line_to(g->cr, x + width, y + height / 2);
      if (width - arrow_width > 0)
        cairo_line_to(g->cr, x + width - arrow_width, y + height);
      cairo_line_to(g->cr, x, y + height);
      cairo_close_path(g->cr);
      break;
    case ARROW_LEFT:
      cairo_move_to(g->cr, x + width, y);
      if (width - arrow_width > 0)
        cairo_rel_line_to(g->cr, -(width - arrow_width), 0);
      cairo_line_to(g->cr, x, y + height / 2);
      cairo_line_to(g->cr, x + MIN(width, arrow_width), y + height);
      if (width - arrow_width > 0)
        cairo_line_to(g->cr, x + width, y + height);
      cairo_close_path(g->cr);
      break;
    case ARROW_BOTH: /* XXX */
    case ARROW_NONE:
      cairo_rectangle(g->cr, x, y, width, height);
   }
   /* fill area */
   cairo_set_source_rgb(g->cr, fill_color.red,
                               fill_color.green,
                               fill_color.blue);
   cairo_fill_preserve(g->cr);
   /* draw outline */
   cairo_set_line_width(g->cr, stroke_width);
   cairo_set_source_rgb(g->cr, stroke_color.red,
                               stroke_color.green,
                               stroke_color.blue);
   cairo_stroke(g->cr);
   /* restore cairo context */
   cairo_restore(g->cr);
}

void graphics_draw_dashes(Graphics *g, double x, double y, double width,
                          double height, ArrowStatus arrow_status,
                          double arrow_width, double stroke_width,
                          Color stroke_color)
{
  double dashes[] = {3.0};
  assert(g);
  /* save cairo context */
  cairo_save(g->cr);
  cairo_rectangle(g->cr, g->margin_x, g->margin_y,
                  g->width-2*g->margin_x, g->height-2*g->margin_y);
  cairo_clip(g->cr);
  cairo_set_line_width(g->cr, stroke_width);
  cairo_set_source_rgb(g->cr, stroke_color.red,
                              stroke_color.green,
                              stroke_color.blue);
  /* create arrowhead path */
  if (width - arrow_width > 0)
  {
    switch (arrow_status)
    {
      case ARROW_RIGHT:
        cairo_move_to(g->cr, x + width - arrow_width, y);
        cairo_line_to(g->cr, x + width, y + height / 2);
        cairo_line_to(g->cr, x + width - arrow_width, y + height);
        /* draw arrowhead */
        cairo_stroke(g->cr);
        break;
      case ARROW_LEFT:
        cairo_move_to(g->cr, width - arrow_width, y);
        cairo_line_to(g->cr, x + width, y + height / 2);
        cairo_line_to(g->cr, x + width - arrow_width, y + height);
        /* draw arrowhead */
        cairo_stroke(g->cr);
        break;
      case ARROW_BOTH: /* XXX */
      case ARROW_NONE: break;
    }
  }
  /* draw dashes */
  cairo_set_dash(g->cr, dashes, 1, (double) 0);
  cairo_move_to(g->cr, x, y+(height/2));
  cairo_rel_line_to(g->cr, width, 0);
  cairo_stroke(g->cr);
  /* restore cairo context */
  cairo_restore(g->cr);
}

void graphics_draw_caret(Graphics *g, double x, double y, double width,
                         double height, ArrowStatus arrow_status,
                         double arrow_width, double stroke_width,
                         Color stroke_color)
{
  assert(g);
  /* save cairo context */
  cairo_save(g->cr);
  cairo_rectangle(g->cr, g->margin_x, g->margin_y,
                  g->width-2*g->margin_x, g->height-2*g->margin_y);
  cairo_clip(g->cr);
  /* create caret path */
  switch (arrow_status)
  {
    case ARROW_RIGHT:
      if (width - arrow_width > 0)
      {
        width -= arrow_width;
        cairo_move_to(g->cr, x, y + (height/2));
        cairo_line_to(g->cr, x + (width/2), y);
        cairo_line_to(g->cr, x + width, y + (height/2));
        cairo_line_to(g->cr, x + width + arrow_width, y + (height/2));
        /* add arrowhead */
        cairo_move_to(g->cr, x + width, y);
        cairo_line_to(g->cr, x + width + arrow_width, y + height / 2);
        cairo_line_to(g->cr, x + width, y + height);
      }
      break;
    case ARROW_LEFT:
      if (width - arrow_width > 0)
      {
        width -= arrow_width;
        /* draw arrowhead */
        cairo_move_to(g->cr, x + arrow_width , y);
        cairo_line_to(g->cr, x, y + (height/2));
        cairo_line_to(g->cr, x + arrow_width, y + height);
        cairo_move_to(g->cr, x, y + (height/2));
        cairo_line_to(g->cr, x + arrow_width, y + (height/2));
        cairo_line_to(g->cr, x + arrow_width + (width/2), y);
        cairo_line_to(g->cr, x + arrow_width + width, y + (height/2));
      }
      break;
    case ARROW_BOTH: /* XXX */
    case ARROW_NONE:
     cairo_move_to(g->cr, x, y + (height/2));
     cairo_line_to(g->cr, x + (width/2), y);
     cairo_line_to(g->cr, x + width, y + (height/2));
   }
   /* draw caret */
   cairo_set_line_width(g->cr, stroke_width);
   cairo_set_source_rgb(g->cr, stroke_color.red,
                               stroke_color.green,
                               stroke_color.blue);
   cairo_stroke(g->cr);
   /* restore cairo context */
   cairo_restore(g->cr);
}

void graphics_draw_rectangle(Graphics *g, double x, double y,
                             bool filled, Color fill_color, bool outlined,
                             Color outline_color, double outline_width,
                             double width, double height)
{
  assert(g);
  /* save cairo context */
  cairo_save(g->cr);
  cairo_new_path(g->cr);
  cairo_rectangle(g->cr, x, y, x + width, y + width);
  if (filled) {
    cairo_set_source_rgb(g->cr, fill_color.red,
                                fill_color.green,
                                fill_color.blue);
    cairo_fill_preserve(g->cr);
  }
  if (outlined) {
    cairo_set_line_width(g->cr, outline_width);
    cairo_set_source_rgb(g->cr, outline_color.red,
                                outline_color.green,
                                outline_color.blue);
    cairo_stroke(g->cr);
  }
  /* restore cairo context */
  cairo_restore(g->cr);
}

void graphics_draw_arrowhead(Graphics *g, double x, double y,
                             Color color, ArrowStatus arrow_status)
{
  double arrow_height = 8, arrow_width = 5;
  assert(g);
  /* save cairo context */
  cairo_save(g->cr);
  cairo_reset_clip(g->cr);
  cairo_set_source_rgb(g->cr, color.red, color.green, color.blue);
  switch (arrow_status)
  {
    case ARROW_LEFT:
      cairo_move_to(g->cr, x+arrow_width, y);
      cairo_line_to(g->cr, x, y+(arrow_height/2));
      cairo_line_to(g->cr, x+arrow_width, y+arrow_height);
      cairo_close_path(g->cr);
      /* fill area */
      cairo_fill_preserve(g->cr);
      cairo_stroke(g->cr);
      break;
    case ARROW_RIGHT:
      cairo_move_to(g->cr, x, y);
      cairo_line_to(g->cr, x+arrow_width, y+(arrow_height/2));
      cairo_line_to(g->cr, x, y+arrow_height);
      cairo_close_path(g->cr);
      /* fill area */
      cairo_fill_preserve(g->cr);
      cairo_stroke(g->cr);
      break;
    case ARROW_BOTH: /* XXX */
    case ARROW_NONE: break;
  }
  /* restore cairo context */
  cairo_restore(g->cr);
}

int graphics_save(const Graphics *g, Env *env)
{
  cairo_status_t rval = CAIRO_STATUS_SUCCESS;
  env_error_check(env);
  assert(g);

  switch (g->type)
  {
    case PNG:
      rval = cairo_surface_write_to_png(g->surf, g->filename);
      break;
    case PDF:
    case PS:
      cairo_show_page(g->cr);
      cairo_surface_flush(g->surf);
      if (!file_exists(g->filename))
        rval = CAIRO_STATUS_WRITE_ERROR;
      break;
  }
  assert(rval == CAIRO_STATUS_SUCCESS || rval == CAIRO_STATUS_WRITE_ERROR);
  if (rval == CAIRO_STATUS_WRITE_ERROR) {
    env_error_set(env, "an I/O error occurred while attempting to write image "
                  "file \"%s\"", g->filename);
    return -1;
  }
  return 0;
}

void graphics_delete(Graphics *g, Env *env)
{
  if (!g) return;
  cairo_surface_destroy(g->surf); /* reference counted */
  cairo_destroy(g->cr);
  ma_free(g);
}
