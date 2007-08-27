/*
  Copyright (c) 2007 Sascha Steinbiss <ssteinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>,
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/
/**
 * \if INTERNAL \file graphics.h \endif
 * \author Sascha Steinbiss <ssteinbiss@stud.zbh.uni-hamburg.de>
 * \author Gordon Gremme <gremme@zbh.uni-hamburg.de>
 */

#ifndef GRAPHICS_H
#define GRAPHICS_H

#include "libgtview/color.h"
#include "libgtview/element.h"

/* specifies arrowhead directions on elements */
enum
{
  ARROW_LEFT,
  ARROW_RIGHT,
  ARROW_BOTH,
  ARROW_NONE
} ArrowStatus;

/* the graphics class */
typedef struct Graphics Graphics;

void      graphics_draw_horizontal_line(Graphics *g, double x, double y,
                                        double width);
void      graphics_draw_text(Graphics *g, double x, double y, const char *text);

/* new interface functions, ssteinbiss */

/*!
Creates a new Graphics object, which is an abstraction of a
drawing surface on which several relevant primitives can be drawn.
This constructor creates a Graphics that can be written out as a PNG file.
\param filename Filename of the resulting PNG.
\param width Surface width (in pixels).
\param height Surface height (in pixels).
\param env Pointer to Environment object.
\return Pointer to new Graphics object.
*/
Graphics* graphics_new_png(const char *filename, unsigned int width,
                           unsigned int height, Env*);

/*!
Sets margins (space to the image boundaries that are clear of elements)
in the graphics.
\param g Graphics object to modify.
\param margin_x Margin to the left and right, in pixels.
\param margin_y Margin to the top and bottom, in pixels.
\param width Surface width (in pixels).
\param height Surface height (in pixels).
\param env Pointer to Environment object.
\return Pointer to new Graphics object.
*/
void      graphics_set_margins(Graphics *g, double margin_x, double margin_y,
                               double width, double height);

/*!
Draws a vertical line beginning at the given coordinates downwards.
\param g Graphics object to modify.
\param x Horizontal coordinate to draw at.
\param y Vertical coordinate to draw at.
\param color Color of the line.
\param length Length of the line, in pixels.
*/
void      graphics_draw_vertical_line(Graphics *g, double x, double y,
                                      Color color, double length);

/*!
Draws a box style glyph at the given coordinates with the given style
parameters.
\param g Graphics object to modify.
\param x Horizontal coordinate to draw at.
\param y Vertical coordinate to draw at.
\param width Width of the box including arrowhead, in pixels.
\param height Height of the box, in pixels.
\param fill_color Color of the box.
\param arrow_status Arrow status, of ArrowStatus enum type.
\param arrow_width Width of the arrowhead, in pixels.
\param stroke_width Outline width, in pixels
                    (can be a decimal subpixel value).
\param stroke_color Outline color
*/
void      graphics_draw_box(Graphics *g, double x, double y, double width,
                            double height, Color fill_color, int arrow_status,
                            double arrow_width, double stroke_width,
                            Color stroke_color);

/*!
Draws a dashed line glyph at the given coordinates with the given style
parameters.
\param g Graphics object to modify.
\param x Horizontal coordinate to draw at.
\param y Vertical coordinate to draw at.
\param width Width of the box including arrowhead, in pixels.
\param height Height of the box with arrowhead, in pixels.
\param arrow_status Arrow status, of ArrowStatus enum type.
\param arrow_width Width of the arrowhead, in pixels.
\param stroke_width Outline width, in pixels
                    (can be a decimal subpixel value).
\param stroke_color Outline color
*/
void      graphics_draw_dashes(Graphics *g, double x, double y, double width,
                               double height, int arrow_status,
                               double arrow_width, double stroke_width,
                               Color stroke_color);

/*!
Draws a caret (^) style glyph at the given coordinates with the given style
parameters.
\param g Graphics object to modify.
\param x Horizontal coordinate to draw at.
\param y Vertical coordinate to draw at.
\param width Width of the box including arrowhead, in pixels.
\param height Height of the box, in pixels.
\param arrow_status Arrow status, of ArrowStatus enum type.
\param arrow_width Width of the arrowhead, in pixels.
\param stroke_width Outline width, in pixels
                    (can be a decimal subpixel value).
\param stroke_color Outline color
*/
void      graphics_draw_caret(Graphics *g, double x, double y, double width,
                              double height, int arrow_status,
                              double arrow_width,  double stroke_width,
                              Color stroke_color);

/*!
Draws a rectangle.
\param g Graphics object to modify.
\param x Horizontal coordinate to draw at.
\param y Vertical coordinate to draw at.
\param filled Defines whether the rectangle should be filled with a color.
\param fill_color Color to fill with.
\param stroked Defines whether the rectangle should be outlined.
\param stroke_color Outline color
\param stroke_width Outline width, in pixels
                    (can be a decimal subpixel value).
\param width Width of the rectangle, in pixels.
\param height Height of the rectangle, in pixels.
*/
void graphics_draw_rectangle(Graphics *g, double x, double y,
                             bool filled, Color fill_color, bool stroked,
                             Color stroke_color, double stroke_width,
                             double width, double height);

/*!
Draws the given string to the right of the given coordinates.
\param g Graphics object to modify.
\param x Horizontal coordinate to draw at.
\param y Vertical coordinate to draw at.
\param text String to draw.
*/
#define   graphics_draw_text_left(g,x,y,t) \
          graphics_draw_text(g,x,y,t);

/*!
Draws the given string centered over the given coordinates.
\param g Graphics object to modify.
\param x Horizontal coordinate to draw at.
\param y Vertical coordinate to draw at.
\param text String to draw.
*/
void      graphics_draw_text_centered(Graphics *g, double x, double y,
                                      const char *text);

/*!
Draws the given string to the left of the given coordinates.
\param g Graphics object to modify.
\param x Horizontal coordinate to draw at.
\param y Vertical coordinate to draw at.
\param text String to draw.
*/
void      graphics_draw_text_right(Graphics *g, double x, double y,
                                      const char *text);

/*!
Draws the given string to the right of the given coordinates.
Accepts a color for the text.
\param g Graphics object to modify.
\param x Horizontal coordinate to draw at.
\param y Vertical coordinate to draw at.
\param color Color for the text.
\param text String to draw.
*/
void      graphics_draw_colored_text(Graphics *g, double x, double y,
                                     Color color, const char *text);

/*!
Draws a single arrowhead.
\param g Graphics object to modify.
\param x Horizontal coordinate to draw at.
\param y Vertical coordinate to draw at.
\param color Color for the arrowhead.
\param arrow_status Arrow status, of ArrowStatus enum type.
*/
void      graphics_draw_arrowhead(Graphics *g, double x, double y,
                                  Color color, int arrow_status);

/*!
Returns text height.
\param g Graphics object to query.
\return Text height in pixels.
*/
double    graphics_get_text_height(Graphics *g);

/*!
Returns text width for a given string.
\param g Graphics object to query.
\param text String to get extents for.
\return Text height in pixels.
*/
double    graphics_get_text_width(Graphics *g, const char *text);

/*!
Writes out the Graphics to the chosen source.
\param g Graphics object to save.
*/
bool     graphics_save(const Graphics*);

/*!
Deletes a Graphics object.
\param g Graphics object to delete.
\param env Pointer to Environment object.
*/
void      graphics_delete(Graphics*, Env*);

#endif
