/*
  Copyright (c) 2007 Sascha Steinbiss <ssteinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>,
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <cairo.h>
#include <cairo-pdf.h>
#include <cairo-ps.h>
#include <libgtview/graphics.h>

#define EXON_ARROW_WIDTH        6

struct Graphics {
  cairo_t *cr;
  cairo_surface_t *surf;
  double margin_x, margin_y, height, width;
  const char* fn;
};

Graphics* graphics_new(unsigned int width, unsigned int height, Env *env)
{
  Graphics *g = env_ma_malloc(env, sizeof (Graphics));
  g->surf = cairo_image_surface_create(CAIRO_FORMAT_RGB24, width, height);
  g->cr = cairo_create(g->surf);
  assert(cairo_status(g->cr) == CAIRO_STATUS_SUCCESS);
  cairo_set_source_rgb(g->cr, 1, 1, 1);
  cairo_set_operator(g->cr, CAIRO_OPERATOR_SOURCE);
  cairo_paint(g->cr);
  return g;
}

void graphics_draw_exon_box(Graphics *g, double x, double y, double width,
                            double height, Strand strand)
{
  assert(g);
  cairo_set_source_rgb(g->cr, 0, 0, 1);
  switch (strand) {
    case STRAND_FORWARD:
      cairo_move_to(g->cr, x, y);
      if (width - EXON_ARROW_WIDTH > 0)
        cairo_rel_line_to(g->cr, width - EXON_ARROW_WIDTH, 0);
      cairo_line_to(g->cr, x + width, y + height / 2);
      if (width - EXON_ARROW_WIDTH > 0)
        cairo_line_to(g->cr, x + width - EXON_ARROW_WIDTH, y + height);
      cairo_line_to(g->cr, x, y + height);
      cairo_close_path(g->cr);
      break;
    case STRAND_REVERSE:
      cairo_move_to(g->cr, x + width, y);
      if (width - EXON_ARROW_WIDTH > 0)
        cairo_rel_line_to(g->cr, -(width - EXON_ARROW_WIDTH), 0);
      cairo_line_to(g->cr, x, y + height / 2);
      cairo_line_to(g->cr, x + MIN(width, EXON_ARROW_WIDTH), y + height);
      if (width - EXON_ARROW_WIDTH > 0)
        cairo_line_to(g->cr, x + width, y + height);
      cairo_close_path(g->cr);
      break;
    case STRAND_BOTH:
    case STRAND_UNKNOWN:
      cairo_rectangle(g->cr, x, y, width, height);
   }

   cairo_fill_preserve(g->cr);
   cairo_set_source_rgb(g->cr, 0, 0, 0);
   cairo_stroke(g->cr);
}

void graphics_draw_horizontal_line(Graphics *g, double x, double y,
                                   double width)
{
  assert(g);
  cairo_move_to(g->cr, x, y);
  cairo_rel_line_to(g->cr, width, 0);
  cairo_stroke(g->cr);
}

void graphics_draw_text(Graphics *g, double x, double y, const char *text)
{
  assert(g && text);
  cairo_set_source_rgb(g->cr, 0, 0, 0);
  cairo_move_to(g->cr, x, y);
  cairo_show_text(g->cr, text);
}

void graphics_save_as_png(const Graphics *g, const char *path)
{
  assert(g);
  cairo_surface_write_to_png(g->surf, path);
}

void graphics_delete(Graphics *g, Env *env)
{
  if (!g) return;
  cairo_surface_destroy(g->surf); /* reference counted */
  cairo_destroy(g->cr);
  env_ma_free(g, env);
}

/* new functions -------------------------------------------------------------*/
Graphics* graphics_new_png(const char *fname, unsigned int width,
                           unsigned int height, Env *env)
{
  Graphics *g = env_ma_malloc(env, sizeof (Graphics));
  g->surf = cairo_image_surface_create(CAIRO_FORMAT_RGB24, width, height);
  g->fn = fname;
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
  return g;
}

void graphics_draw_box(Graphics *g, double x, double y, double width,
                       double height, Color fill_color, int arrow_status,
                       double arrow_width, double stroke_width,
                       Color stroke_color)
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
                          double height, int arrow_status,
                          double arrow_width, double stroke_width,
                          Color stroke_color)
{
  double dashes[] = {5.0};
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
                         double height, int arrow_status,
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

void graphics_draw_vertical_line(Graphics *g, double x, double y,
                                 Color color, double length)
{
  assert(g);
  cairo_move_to(g->cr, x, y);
  cairo_set_source_rgb(g->cr, color.red, color.green, color.blue);
  cairo_rel_line_to(g->cr, 0, length);
  cairo_stroke(g->cr);
}

double graphics_get_text_height(Graphics *g)
{
  cairo_text_extents_t ext;
  assert(g);
  /* get text extents */
  cairo_text_extents(g->cr, "A", &ext);
	return ext.height;
}

double graphics_get_text_width(Graphics *g, const char* text)
{
  cairo_text_extents_t ext;
  assert(g);
  /* get text extents */
  cairo_text_extents(g->cr, text, &ext);
	return ext.width;
}

void graphics_draw_colored_text(Graphics *g,
                                double x,
                                double y,
                                Color color,
                                const char *text)
{
  assert(g && text);
  cairo_set_source_rgb(g->cr,
                       color.red,
                       color.green,
                       color.blue);
  cairo_move_to(g->cr, x, y);
  cairo_show_text(g->cr, text);
}

void graphics_draw_arrowhead(Graphics *g, double x, double y,
                             Color color, int arrow_status)
{
  assert(g);
  double arrow_height = 8, arrow_width = 5;

  /* save cairo context */
  cairo_save(g->cr);
  cairo_reset_clip(g->cr);
  cairo_set_source_rgb(g->cr, color.red, color.green, color.blue);
  if (arrow_status == ARROW_LEFT) {
    cairo_move_to(g->cr, x+arrow_width, y);
    cairo_line_to(g->cr, x, y+(arrow_height/2));
    cairo_line_to(g->cr, x+arrow_width, y+arrow_height);
    cairo_close_path(g->cr);
    /* fill area */
    cairo_fill_preserve(g->cr);
    cairo_stroke(g->cr);
  }
  if (arrow_status == ARROW_RIGHT) {
    cairo_move_to(g->cr, x, y);
    cairo_line_to(g->cr, x+arrow_width, y+(arrow_height/2));
    cairo_line_to(g->cr, x, y+arrow_height);
    cairo_close_path(g->cr);
    /* fill area */
    cairo_fill_preserve(g->cr);
    cairo_stroke(g->cr);
  }
  /* restore cairo context */
  cairo_restore(g->cr);
}

void graphics_set_margins(Graphics *g, double margin_x, double margin_y,
                          double width, double height)
{
  g->margin_x = margin_x;
  g->margin_y = margin_y;
}

void graphics_save(const Graphics *g)
{
  assert(g);
  cairo_surface_write_to_png(g->surf, g->fn);
}
