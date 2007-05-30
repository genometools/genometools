/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <cairo.h>
#include <libgtext/graphics.h>

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
