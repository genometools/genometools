/*
  Copyright (c) 2007 Sascha Steinbiss <ssteinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef RENDER_H
#define RENDER_H

#include "libgtview/diagram.h"
#include "libgtview/config.h"

#define DEFAULT_RENDER_WIDTH  800

/* The Render class used for Diagram to Image conversion. */
typedef struct Render Render;

/* <cfg> is used to determine drawing options. */
Render* render_new(Config *cfg, Env*);
/* Render <diagram> to PNG file <filename> (relative to working directory) */
int     render_to_png(Render*, Diagram *diagram, const char *filename,
                      unsigned int width, Env*);
void    render_delete(Render*, Env*);

#endif
