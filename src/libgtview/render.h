/*
  Copyright (c) 2007 Sascha Steinbiss <ssteinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/
/**
 * \file render.h
 * \author Sascha Steinbiss <ssteinbiss@stud.zbh.uni-hamburg.de>
 */

#ifndef RENDER_H
#define RENDER_H

#include "libgtview/diagram.h"
#include "libgtview/config.h"

#define DEFAULT_RENDER_WIDTH  800

/* the Render class
   contains methods for Diagram->Image conversion */
typedef struct Render Render;

/*!
Creates a new Render object.
\param cfg Pointer to Config object. Used to determine
           drawing options.
\param env Pointer to Environment object.
\return Created Render object.
*/
Render* render_new(Config *cfg, Env *env);

/*!
Renders a Diagram to a PNG file.
\param r Render object.
\param dia Diagram that should be rendered.
\param filename Filename (relative to working directory)
          the image should be written to.
\param width Target image width (in pixels).
\param env Pointer to Environment object.
\return Error code.
*/
int     render_to_png(Render *r, Diagram *dia, const char *filename,
                      unsigned int width, Env *env);

/*!
Deletes a Render object
\param r Render object
\param env Pointer to Environment object.
*/
void    render_delete(Render *r, Env *env);

#endif
