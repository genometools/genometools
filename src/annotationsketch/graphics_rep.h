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

struct GraphicsClass {
  size_t size;
  void    (*draw_text)(Graphics*, double, double, const char*);
  void    (*draw_text_centered)(Graphics*, double, double, const char*);
  void    (*draw_text_right)(Graphics*, double, double, const char*);
  void    (*draw_colored_text)(Graphics*, double, double, GT_Color,
                               const char*);
  double  (*get_text_height)(Graphics*);
  double  (*get_text_width)(Graphics*, const char*);
  void    (*set_font)(Graphics*, const char*, FontSlant, FontWeight);
  double  (*get_image_width)(Graphics*);
  double  (*get_image_height)(Graphics*);
  void    (*set_margins)(Graphics*, double, double);
  void    (*draw_horizontal_line)(Graphics*, double, double, double);
  void    (*draw_vertical_line)(Graphics*, double, double, GT_Color, double);
  void    (*draw_box)(Graphics*, double, double, double, double, GT_Color,
                      ArrowStatus, double, double, GT_Color, bool);
  void    (*draw_dashes)(Graphics*, double, double, double, double, ArrowStatus,
                         double, double, GT_Color);
  void    (*draw_caret)(Graphics*, double, double, double, double, ArrowStatus,
                        double, double, GT_Color);
  void    (*draw_rectangle)(Graphics*, double, double, bool, GT_Color, bool,
                            GT_Color, double, double);
  void    (*draw_arrowhead)(Graphics*, double, double, GT_Color, ArrowStatus);
  int     (*save_to_file)(const Graphics*, const char*, GT_Error*);
  void    (*save_to_stream)(const Graphics*, GT_Str*);
  void    (*free)(Graphics*);
};

struct Graphics {
  const GraphicsClass *c_class;
  unsigned int reference_count;
};

Graphics* graphics_create(const GraphicsClass*);
void*     graphics_cast(const GraphicsClass*, Graphics*);

#endif
