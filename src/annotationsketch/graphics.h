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

#include "annotationsketch/graphics_api.h"
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

typedef struct GT_GraphicsClass GT_GraphicsClass;
typedef struct GT_Graphics GT_Graphics;

#include "annotationsketch/color.h"
#include "core/str.h"
#include "core/error.h"

/* Create a new GT_Graphics object, which is an abstraction of a drawing
   surface on which several relevant primitives can be drawn. */
void      gt_graphics_draw_text(GT_Graphics*, double x, double y, const char*);
#define   gt_graphics_draw_text_left(g,x,y,t) \
          gt_graphics_draw_text(g,x,y,t);
void      gt_graphics_draw_text_centered(GT_Graphics*, double x, double y,
                                         const char*);
void      gt_graphics_draw_text_right(GT_Graphics*, double x, double y,
                                      const char*);
void      gt_graphics_draw_colored_text(GT_Graphics*, double x, double y,
                                        GT_Color, const char*);
double    gt_graphics_get_text_height(GT_Graphics*);
double    gt_graphics_get_text_width(GT_Graphics*, const char *text);
void      gt_graphics_set_font(GT_Graphics *g, const char *family,
                            FontSlant slant, FontWeight weight);
double    gt_graphics_get_image_width(GT_Graphics*);
double    gt_graphics_get_image_height(GT_Graphics*);
/* Set margins (space to the image boundaries that are clear of elements)
   in the graphics.
   <margin_x> denotes the Margin to the left and right, in pixels.
   <margin_y> denotes the Margin to the top and bottom, in pixels. */
void      gt_graphics_set_margins(GT_Graphics*, double margin_x,
                                  double margin_y);
void      gt_graphics_draw_horizontal_line(GT_Graphics*, double x, double y,
                                        double width);
/* Draws a vertical line beginning at the given coordinates downwards. */
void      gt_graphics_draw_vertical_line(GT_Graphics*, double x, double y,
                                         GT_Color color, double length);
void      gt_graphics_draw_box(GT_Graphics*, double x, double y, double width,
                            double height, GT_Color fill_color,
                            ArrowStatus arrow_status, double arrow_width,
                            double stroke_width, GT_Color stroke_color,
                            bool dashed);
void      gt_graphics_draw_dashes(GT_Graphics*, double x, double y,
                                  double width, double height,
                                  ArrowStatus arrow_status, double arrow_width,
                                  double stroke_width, GT_Color stroke_color);
/* Draws a caret (^) style glyph. */
void      gt_graphics_draw_caret(GT_Graphics*, double x, double y, double width,
                              double height, ArrowStatus arrow_status,
                              double arrow_width,  double stroke_width,
                              GT_Color stroke_color);
void      gt_graphics_draw_rectangle(GT_Graphics*, double x, double y,
                                  bool filled, GT_Color fill_color,
                                  bool outlined, GT_Color outgt_line_color,
                                  double outgt_line_width, double width);
void      gt_graphics_draw_arrowhead(GT_Graphics*, double x, double y, GT_Color,
                                     ArrowStatus);
/* Write out the Graphic to the given file with <filename>. */
int       gt_graphics_save_to_file(const GT_Graphics*, const char *filename,
                                   GT_Error*);
/* Write out the Graphic to the given <stream>. */
void      gt_graphics_save_to_stream(const GT_Graphics*, GT_Str *stream);

void      gt_graphics_delete(GT_Graphics*);

#endif
