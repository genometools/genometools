/*
  Copyright (c) 2007      Sascha Steinbiss <ssteinbiss@zbh.uni-hamburg.de>
  Copyright (c)      2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007-2008 Center for Bioinformatics, University of Hamburg

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

#ifndef RENDER_H
#define RENDER_H

#include "libgtview/diagram.h"
#include "libgtview/config.h"

#define DEFAULT_RENDER_WIDTH  800

/* The Render class used for Diagram to Image conversion. */
typedef struct Render Render;

/* <cfg> is used to determine drawing options. */
Render* render_new(Config*);
/* Render diagram to PNG file <filename> of given <width>
   (relative to working directory) */
int     render_to_png(Render*, Diagram*, const char *filename,
                      unsigned int width, Error*);
/* Render <diagram> to <stream> (in PNG format of given <width>) */
void    render_to_png_stream(Render*, Diagram *diagram, Str *stream,
                             unsigned int width);
void    render_delete(Render*);

#endif
