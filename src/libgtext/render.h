/*
  Copyright (c) 2007 Sascha Steinbiss <ssteinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef RENDER_H
#define RENDER_H

#include <gtcore.h>
#include <libgtext/diagram.h>
#include <libgtext/config.h>

/* the render class */
typedef struct Render Render;

Render* render_new(Diagram* dia, Config* cfg, Env* env);
void    render_to_png(Render* r, char* fn, unsigned int width, Env* env);
void    render_delete(Render* r, Env* env);

#endif
