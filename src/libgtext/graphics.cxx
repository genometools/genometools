/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#define AGG_RENDERING_BUFFER row_ptr_cache<int8u>
#include <agg_path_storage.h>
#include <agg_pixfmt_rgba.h>
#include <agg_rasterizer_scanline_aa.h>
#include <agg_renderer_base.h>
#include <agg_renderer_scanline.h>
#include <agg_scanline_p.h>
#include <agg_gsv_text.h>
#include <png.h>
extern "C" {
#include <libgtext/graphics.h>
}

#define EXON_ARROW_WIDTH  8

/* the used pixel format */
typedef agg::pixfmt_rgba32 PixelFormat;

/* the used base renderer */
typedef agg::renderer_base<PixelFormat> RendererBase;

/* the used solid renderer */
typedef agg::renderer_scanline_aa_solid<RendererBase> RendererSolid;

struct Graphics {
  unsigned int width,
               height;
  unsigned char *pixel_buffer; /* the actual pixel buffer */
  /* the AGG stuff */
  agg::rendering_buffer *rendering_buffer; /* the AGG rendering buffer (is
                                              attached to the <pixel_buffer>) */
  agg::rasterizer_scanline_aa<> *rasterizer_scanline;
  PixelFormat *pixel_format;               /* the AGG pixel format (is attached
                                              to the <rendering_buffer> */
  RendererBase *renderer_base;             /* the AGG base renderer (is attached
                                              to the <pixel_format> */
  RendererSolid *renderer_solid;           /* the AGG solid renderer (is
                                              attached to the <renderer_base> */
  agg::scanline_p8 *scanline;
};

Graphics* graphics_new(unsigned int width, unsigned int height, Env *env)
{
  Graphics *g = (Graphics*) env_ma_malloc(env, sizeof (Graphics));
  g->width = width;
  g->height = height;
  /* init the AGG stuff */
  g->pixel_buffer = (unsigned char*) env_ma_malloc(env, width * height * 4);
  g->rendering_buffer = new agg::rendering_buffer;
  assert(g->rendering_buffer);
  g->rendering_buffer->attach(g->pixel_buffer, width, height, width * 4);
  g->rasterizer_scanline = new agg::rasterizer_scanline_aa<>;
  assert(g->rasterizer_scanline);
  g->pixel_format = new PixelFormat(*g->rendering_buffer);
  assert(g->pixel_format);
  g->renderer_base = new RendererBase;
  assert(g->renderer_base);
  g->renderer_base->attach(*g->pixel_format);
  g->renderer_solid = new RendererSolid;
  assert(g->renderer_solid);
  g->renderer_solid->attach(*g->renderer_base);
  g->scanline = new agg::scanline_p8;
  assert(g->scanline);

  /* set white background */
  g->rasterizer_scanline->clip_box(0, 0, g->renderer_base->width(),
                                   g->renderer_base->height());
  g->renderer_base->clear(agg::rgba(1, 1, 1, 1));

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

static void graphics_set_color(Graphics *g)
{
  assert(g);
  g->renderer_solid->color(agg::rgba(0, 0, 0)); /* set to black */
}

static void graphics_render(Graphics *g)
{
  assert(g);
  agg::render_scanlines(*g->rasterizer_scanline, *g->scanline, *g->renderer_solid);
}

void graphics_draw_horizontal_line(Graphics *g, double x, double y,
                                   double width)
{
  agg::path_storage path_storage;
  assert(g);
  path_storage.move_to(x, y);
  path_storage.line_to(x + width, y);
  agg::conv_stroke<agg::path_storage> sp(path_storage);
  sp.line_cap(agg::round_cap);
  sp.line_join(agg::round_join);
  sp.width(2);
  g->rasterizer_scanline->add_path(sp);

  graphics_set_color(g);
  graphics_render(g);
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

  rows = (png_byte**) g->rendering_buffer->rows();
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
  delete g->scanline;
  delete g->renderer_solid;
  delete g->renderer_base;
  delete g->pixel_format;
  delete g->rasterizer_scanline;
  delete g->rendering_buffer;
  env_ma_free(g->pixel_buffer, env);
  env_ma_free(g, env);
}
