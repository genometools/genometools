/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#define AGG_RENDERING_BUFFER row_ptr_cache<int8u>
#include <agg_rendering_buffer.h>
#include <png.h>
extern "C" {
#include <libgtext/graphics.h>
}

#define EXON_ARROW_WIDTH        8

struct Graphics {
  unsigned int width,
               height;
  unsigned char *pixbuf;
  agg::rendering_buffer *rbuf;
};

Graphics* graphics_new(unsigned int width, unsigned int height, Env *env)
{
  Graphics *g = (Graphics*) env_ma_malloc(env, sizeof (Graphics));
  g->width = width;
  g->height = height;
  g->pixbuf = (unsigned char*) env_ma_malloc(env, width * height * 4);
  memset(g->pixbuf, 255, width * height * 4);
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

void graphics_save_as_png(const Graphics *g, const char *path, Env *env)
{
  FILE *outfp;
  png_structp png = NULL;
  png_infop info = NULL;
  png_bytepp rows;
  int has_err = 0;

  env_error_check(env);
  assert(g && path);

  rows = (png_byte**) g->rbuf->rows();
  outfp = env_fa_xfopen(env, path, "w");

  png = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
  if (!png) has_err = -1;
  if (!has_err) {
    info = png_create_info_struct(png);
    if (!info) has_err = -1;
  }
  if (!has_err) {
    if (setjmp(png_jmpbuf(png)))
      has_err = -1;
    else {
      png_init_io(png, outfp);
      png_set_IHDR(png, info, g->width, g->height, 8, PNG_COLOR_TYPE_RGB_ALPHA,
                   PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT,
                   PNG_FILTER_TYPE_DEFAULT);
      png_set_rows(png, info, rows);
      png_write_png(png, info, PNG_TRANSFORM_IDENTITY, NULL);
    }
  }
  if (png) png_destroy_write_struct(&png, &info);
  assert(!has_err);
  env_fa_xfclose(outfp, env);
}

void graphics_delete(Graphics *g, Env *env)
{
  if (!g) return;
  env_ma_free(g->pixbuf, env);
  delete g->rbuf;
  env_ma_free(g, env);
}
