/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef GRAPHICS_H
#define GRAPHICS_H

#include <gtcore.h>

/* the graphics class */
typedef struct {
  cairo_t *cr;
  cairo_surface_t *surf;
} Graphics;

Graphics* graphics_new(unsigned int width, unsigned int height, Env*);
void      graphics_save_as_png(const Graphics*, const char *path);
void      graphics_delete(Graphics*, Env*);

#endif
