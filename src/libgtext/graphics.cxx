/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <agg_pixfmt_rgb.h>
extern "C" {
#include <libgtext/graphics.h>
}

#define EXON_ARROW_WIDTH        8

struct Graphics {
  unsigned char *pixbuf;
  agg::rendering_buffer *rbuf;
};

Graphics* graphics_new(unsigned int width, unsigned int height, Env *env)
{
  Graphics *g = (Graphics*) env_ma_malloc(env, sizeof (Graphics));
  g->pixbuf = (unsigned char*) env_ma_malloc(env, width * height * 4);
  g->rbuf = new agg::rendering_buffer;
  assert(g->rbuf);
  g->rbuf->attach(g->pixbuf, width, height, width * 4);
  return g;
}

void graphics_draw_exon_box(Graphics *g, double x, double y, double width,
                            double height, Strand strand)
{
  assert(g);

#if 0
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
#endif
}

void graphics_draw_horizontal_line(Graphics *g, double x, double y,
                                   double width)
{
  assert(g);
#if 0
  cairo_move_to(g->cr, x, y);
  cairo_rel_line_to(g->cr, width, 0);
  cairo_stroke(g->cr);
#endif
}

void graphics_draw_text(Graphics *g, double x, double y, const char *text)
{
  assert(g && text);
#if 0
  cairo_set_source_rgb(g->cr, 0, 0, 0);
  cairo_move_to(g->cr, x, y);
  cairo_show_text(g->cr, text);
#endif
}

void graphics_save_as_png(const Graphics *g, const char *path)
{
  assert(g);
#if 0
  cairo_surface_write_to_png(g->surf, path);
#endif
}

void graphics_delete(Graphics *g, Env *env)
{
  if (!g) return;
  env_ma_free(g->pixbuf, env);
  delete g->rbuf;
  env_ma_free(g, env);
}
