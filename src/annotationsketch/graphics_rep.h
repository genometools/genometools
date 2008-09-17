/*
  Copyright (c) 2008 Sascha Steinbiss <ssteinbiss@stud.zbh.uni-hamburg.de>
  Copyright (c) 2008 Center for Bioinformatics, University of Hamburg

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

#ifndef GRAPHICS_REP_H
#define GRAPHICS_REP_H

#include <stdio.h>
#include "annotationsketch/graphics.h"

struct GtGraphicsClass {
  size_t size;
  void    (*draw_text)(GtGraphics*, double, double, const char*);
  void    (*draw_text_centered)(GtGraphics*, double, double, const char*);
  void    (*draw_text_right)(GtGraphics*, double, double, const char*);
  void    (*draw_colored_text)(GtGraphics*, double, double, GtColor,
                               const char*);
  double  (*get_text_height)(GtGraphics*);
  double  (*get_text_width)(GtGraphics*, const char*);
  void    (*set_font)(GtGraphics*, const char*, FontSlant, FontWeight);
  double  (*get_image_width)(GtGraphics*);
  double  (*get_image_height)(GtGraphics*);
  void    (*set_margins)(GtGraphics*, double, double);
  void    (*draw_horizontal_line)(GtGraphics*, double, double, double);
  void    (*draw_vertical_line)(GtGraphics*, double, double, GtColor, double);
  void    (*draw_box)(GtGraphics*, double, double, double, double, GtColor,
                      ArrowStatus, double, double, GtColor, bool);
  void    (*draw_dashes)(GtGraphics*, double, double, double, double,
                         ArrowStatus, double, double, GtColor);
  void    (*draw_caret)(GtGraphics*, double, double, double, double,
                        ArrowStatus, double, double, GtColor);
  void    (*draw_rectangle)(GtGraphics*, double, double, bool, GtColor, bool,
                            GtColor, double, double);
  void    (*draw_arrowhead)(GtGraphics*, double, double, GtColor,
                            ArrowStatus);
  int     (*save_to_file)(const GtGraphics*, const char*, GtError*);
  void    (*save_to_stream)(const GtGraphics*, GtStr*);
  void    (*free)(GtGraphics*);
};

typedef struct {
    unsigned int reference_count;
} GtGraphicsPrivate;

struct GtGraphics {
  const GtGraphicsClass *c_class;
  GtGraphicsPrivate *pvt;
};

GtGraphics* gt_graphics_create(const GtGraphicsClass*);
void*       gt_graphics_cast(const GtGraphicsClass*, GtGraphics*);

#endif
