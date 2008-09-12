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

struct GT_GraphicsClass {
  size_t size;
  void    (*draw_text)(GT_Graphics*, double, double, const char*);
  void    (*draw_text_centered)(GT_Graphics*, double, double, const char*);
  void    (*draw_text_right)(GT_Graphics*, double, double, const char*);
  void    (*draw_colored_text)(GT_Graphics*, double, double, GT_Color,
                               const char*);
  double  (*get_text_height)(GT_Graphics*);
  double  (*get_text_width)(GT_Graphics*, const char*);
  void    (*set_font)(GT_Graphics*, const char*, FontSlant, FontWeight);
  double  (*get_image_width)(GT_Graphics*);
  double  (*get_image_height)(GT_Graphics*);
  void    (*set_margins)(GT_Graphics*, double, double);
  void    (*draw_horizontal_line)(GT_Graphics*, double, double, double);
  void    (*draw_vertical_line)(GT_Graphics*, double, double, GT_Color, double);
  void    (*draw_box)(GT_Graphics*, double, double, double, double, GT_Color,
                      ArrowStatus, double, double, GT_Color, bool);
  void    (*draw_dashes)(GT_Graphics*, double, double, double, double,
                         ArrowStatus, double, double, GT_Color);
  void    (*draw_caret)(GT_Graphics*, double, double, double, double,
                        ArrowStatus, double, double, GT_Color);
  void    (*draw_rectangle)(GT_Graphics*, double, double, bool, GT_Color, bool,
                            GT_Color, double, double);
  void    (*draw_arrowhead)(GT_Graphics*, double, double, GT_Color,
                            ArrowStatus);
  int     (*save_to_file)(const GT_Graphics*, const char*, GtError*);
  void    (*save_to_stream)(const GT_Graphics*, GtStr*);
  void    (*free)(GT_Graphics*);
};

struct GT_Graphics {
  const GT_GraphicsClass *c_class;
  unsigned int reference_count;
};

GT_Graphics* gt_graphics_create(const GT_GraphicsClass*);
void*     gt_graphics_cast(const GT_GraphicsClass*, GT_Graphics*);

#endif
