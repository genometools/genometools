/*
  Copyright (c) 2007      Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
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

typedef struct GtGraphicsClass GtGraphicsClass;
typedef struct GtGraphics GtGraphics;

#include "annotationsketch/color.h"
#include "core/error.h"
#include "core/range.h"
#include "core/str.h"

void   gt_graphics_draw_text(GtGraphics*, double x, double y, const char*);
#define gt_graphics_draw_text_left(g,x,y,t) \
        gt_graphics_draw_text(g,x,y,t);
void   gt_graphics_draw_text_centered(GtGraphics*, double x, double y,
                                      const char*);
void   gt_graphics_draw_text_right(GtGraphics*, double x, double y,
                                   const char*);
void   gt_graphics_draw_colored_text(GtGraphics*, double x, double y,
                                     GtColor, const char*);
double gt_graphics_get_text_height(GtGraphics*);
double gt_graphics_get_text_width(GtGraphics*, const char *text);
void   gt_graphics_set_font(GtGraphics *g, const char *family,
                            FontSlant slant, FontWeight weight);
double gt_graphics_get_image_width(GtGraphics*);
double gt_graphics_get_image_height(GtGraphics*);
/* Set margins (space to the image boundaries that are clear of elements)
   in the graphics.
   <margin_x> denotes the Margin to the left and right, in pixels.
   <margin_y> denotes the Margin to the top and bottom, in pixels. */
void   gt_graphics_set_margins(GtGraphics*, double margin_x,
                                  double margin_y);
double gt_graphics_get_xmargins(GtGraphics*);
double gt_graphics_get_ymargins(GtGraphics*);
/* Draws a horizontal line beginning at the given coordinates downwards. */
void   gt_graphics_draw_horizontal_line(GtGraphics *gg, double x, double y,
                                        GtColor color, double width,
                                        double stroke_width);
/* Draws a vertical line beginning at the given coordinates downwards. */
void   gt_graphics_draw_vertical_line(GtGraphics *gg, double x, double y,
                                      GtColor color, double length,
                                      double stroke_width);
void   gt_graphics_draw_box(GtGraphics*, double x, double y, double width,
                            double height, GtColor fill_color,
                            ArrowStatus arrow_status, double arrow_width,
                            double stroke_width, GtColor stroke_color,
                            bool dashed);
void   gt_graphics_draw_dashes(GtGraphics*, double x, double y,
                                  double width, double height,
                                  ArrowStatus arrow_status, double arrow_width,
                                  double stroke_width, GtColor stroke_color);
/* Draws a caret (^) style glyph. */
void   gt_graphics_draw_caret(GtGraphics*, double x, double y, double width,
                              double height, ArrowStatus arrow_status,
                              double arrow_width,  double stroke_width,
                              GtColor stroke_color);
void   gt_graphics_draw_rectangle(GtGraphics*, double x, double y,
                                  bool filled, GtColor fill_color,
                                  bool outlined, GtColor outgt_line_color,
                                  double outgt_line_width, double width,
                                  double height);
void   gt_graphics_draw_arrowhead(GtGraphics*, double x, double y, GtColor,
                                  ArrowStatus);
void   gt_graphics_draw_curve_data(GtGraphics *g, double x, double y,
                                   GtColor color,
                                   double data[], unsigned long ndata,
                                   GtRange valrange, unsigned long height);
/* Write out the Graphic to the given file with <filename>. */
int    gt_graphics_save_to_file(const GtGraphics*, const char *filename,
                                GtError*);
/* Write out the Graphic to the given <stream>. */
void   gt_graphics_save_to_stream(const GtGraphics*, GtStr *stream);

void   gt_graphics_delete(GtGraphics*);

#endif
