/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef GRAPHICS_H
#define GRAPHICS_H

#include <gtcore.h>
#include <libgtext/color.h>
#include <libgtext/element.h>

enum 
{
  ARROW_LEFT,
  ARROW_RIGHT,
  ARROW_BOTH,
  ARROW_NONE
} ArrowStatus;

/* the graphics class */
typedef struct Graphics Graphics;

Graphics* graphics_new(unsigned int width, unsigned int height, Env*);
void      graphics_draw_exon_box(Graphics *g, double x, double y, double width,
                                 double height, Strand strand);
void      graphics_draw_horizontal_line(Graphics *g, double x, double y,
                                        double width);
void      graphics_draw_text(Graphics *g, double x, double y, const char *text);
void      graphics_save_as_png(const Graphics*, const char *path);

/* new interface functions, ssteinbiss */
Graphics* graphics_new_png(const char *fname, unsigned int width,
                           unsigned int height, Env*);
void      graphics_set_margins(Graphics *g,
                               double margin_x,
															 double margin_y,
															 double width,
															 double height);
void      graphics_draw_vertical_line(Graphics *g, double x, double y,
                                      Color color, double length);
void      graphics_draw_box(Graphics *g, double x, double y, double width,
                            double height, Color color, int arrow_status,
                            double arrow_width, double stroke_width,
                            Color stroke_color);
void      graphics_draw_dashes(Graphics *g, double x, double y, double width,
                               double height, int arrow_status,
                               double arrow_width, double stroke_width,
                               Color stroke_color);
void      graphics_draw_caret(Graphics *g, double x, double y, double width,
                              double height, int arrow_status,
                              double arrow_width,  double stroke_width,
                              Color stroke_color);
#define   graphics_draw_text_left(g,x,y,t) \
          graphics_draw_text(g,x,y,t);
void      graphics_draw_text_centered(Graphics *g, double x, double y,
                                      const char *text);
void      graphics_draw_text_right(Graphics *g, double x, double y,
                                      const char *text);
void      graphics_draw_colored_text(Graphics *g,
                                     double x,
								      	 						 double y,
																     Color color,
																     const char *text);
void      graphics_draw_arrowhead(Graphics *g, double x, double y,
                                  Color color, int arrow_status);
double    graphics_get_text_height(Graphics *g);
double    graphics_get_text_width(Graphics *g, const char *text);
void      graphics_save(const Graphics*);
void      graphics_delete(Graphics*, Env*);

#endif
