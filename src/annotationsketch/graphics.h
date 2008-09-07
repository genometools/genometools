/*
  Copyright (c) 2007      Sascha Steinbiss <ssteinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2007-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>,
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

#ifndef GRAPHICS_H
#define GRAPHICS_H

#include "core/str.h"

typedef enum
{
  ARROW_LEFT,
  ARROW_RIGHT,
  ARROW_BOTH,
  ARROW_NONE
} ArrowStatus;

typedef enum
{
  SLANT_NORMAL,
  SLANT_ITALIC
} FontSlant;

typedef enum
{
  WEIGHT_NORMAL,
  WEIGHT_BOLD
} FontWeight;

typedef enum {
  GRAPHICS_PDF,
  GRAPHICS_PNG,
  GRAPHICS_PS,
  GRAPHICS_SVG
} GraphicsOutType;

typedef struct Graphics Graphics;

#include "annotationsketch/color.h"
#include "core/error.h"

/* Create a new Graphics object, which is an abstraction of a drawing surface on
   which several relevant primitives can be drawn. */
Graphics* graphics_new(GraphicsOutType, unsigned int width,
                       unsigned int height);
void      graphics_draw_text(Graphics*, double x, double y, const char*);
#define   graphics_draw_text_left(g,x,y,t) \
          graphics_draw_text(g,x,y,t);
void      graphics_draw_text_centered(Graphics*, double x, double y,
                                      const char*);
void      graphics_draw_text_right(Graphics*, double x, double y, const char*);
void      graphics_draw_colored_text(Graphics*, double x, double y, GT_Color,
                                     const char*);
double    graphics_get_text_height(Graphics*);
double    graphics_get_text_width(Graphics*, const char *text);
void      graphics_set_font(Graphics *g, const char *family,
                            FontSlant slant, FontWeight weight);
double    graphics_get_image_width(Graphics*);
double    graphics_get_image_height(Graphics*);
/* Set margins (space to the image boundaries that are clear of elements)
   in the graphics.
   <margin_x> denotes the Margin to the left and right, in pixels.
   <margin_y> denotes the Margin to the top and bottom, in pixels. */
void      graphics_set_margins(Graphics*, double margin_x, double margin_y);
void      graphics_draw_horizontal_line(Graphics*, double x, double y,
                                        double width);
/* Draws a vertical line beginning at the given coordinates downwards. */
void      graphics_draw_vertical_line(Graphics*, double x, double y,
                                      GT_Color color, double length);
void      graphics_draw_box(Graphics*, double x, double y, double width,
                            double height, GT_Color fill_color,
                            ArrowStatus arrow_status, double arrow_width,
                            double stroke_width, GT_Color stroke_color,
                            bool dashed);
void      graphics_draw_dashes(Graphics*, double x, double y, double width,
                               double height, ArrowStatus arrow_status,
                               double arrow_width, double stroke_width,
                               GT_Color stroke_color);
/* Draws a caret (^) style glyph. */
void      graphics_draw_caret(Graphics*, double x, double y, double width,
                              double height, ArrowStatus arrow_status,
                              double arrow_width,  double stroke_width,
                              GT_Color stroke_color);
void      graphics_draw_rectangle(Graphics*, double x, double y,
                                  bool filled, GT_Color fill_color,
                                  bool outlined, GT_Color outline_color,
                                  double outline_width, double width);
void      graphics_draw_arrowhead(Graphics*, double x, double y, GT_Color,
                                  ArrowStatus);
/* Write out the Graphic to the given file with <filename>. */
int       graphics_save_to_file(const Graphics*, const char *filename, GT_Error*);
/* Write out the Graphic to the given <stream>. */
void      graphics_save_to_stream(const Graphics*, Str *stream);

void      graphics_delete(Graphics*);

#endif
