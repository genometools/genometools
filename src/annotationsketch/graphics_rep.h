/*
  Copyright (c) 2008 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
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

typedef void   (*GtGraphicsDrawTextFunc)(GtGraphics*, double, double,
                                         const char*);
typedef void   (*GtGraphicsDrawColoredTextFunc)(GtGraphics*, double, double,
                                                GtColor,const char*);
typedef double (*GtGraphicsGetSingleExtentFunc)(GtGraphics*);
typedef double (*GtGraphicsGetTextWidthFunc)(GtGraphics*, const char*);
typedef void   (*GtGraphicsSetMarginsFunc)(GtGraphics*, double, double);
typedef int    (*GtGraphicsSetColorFunc)(GtGraphics*, GtColor);
typedef void   (*GtGraphicsSetFontFunc)(GtGraphics*, const char*, FontSlant,
                                        FontWeight, double);
typedef void   (*GtGraphicsDrawLineFunc)(GtGraphics*, double, double,
                                         GtColor, double, double);
typedef void   (*GtGraphicsDrawLineToFunc)(GtGraphics*, double, double,
                                           double, double, GtColor, double);
typedef void   (*GtGraphicsDrawBoxFunc)(GtGraphics*, double, double, double,
                                        double, GtColor, ArrowStatus, double,
                                        double, GtColor, bool);
typedef void    (*GtGraphicsDrawSimpleFunc)(GtGraphics*, double, double,
                                            double, double, ArrowStatus,
                                            double, double, GtColor);
typedef void    (*GtGraphicsDrawRectFunc)(GtGraphics*, double, double, bool,
                                          GtColor, bool, GtColor, double,
                                          double, double);
typedef void    (*GtGraphicsDrawArrowheadFunc)(GtGraphics*, double, double,
                                               GtColor, ArrowStatus);
typedef void    (*GtGraphicsDrawCurveDataFunc)(GtGraphics *g,
                                               double x, double y,
                                               GtColor color,
                                               double data[],
                                               unsigned long ndata,
                                               GtRange,
                                               unsigned long height);
typedef int     (*GtGraphicsSaveToFileFunc)(const GtGraphics*, const char*,
                                            GtError*);
typedef void    (*GtGraphicsSaveToStreamFunc)(const GtGraphics*, GtStr*);
typedef void    (*GtGraphicsFreeFunc)(GtGraphics*);

typedef struct GtGraphicsMembers GtGraphicsMembers;

struct GtGraphics {
  const GtGraphicsClass *c_class;
  GtGraphicsMembers *pvt;
};

const GtGraphicsClass* gt_graphics_class_new(size_t size,
                                         GtGraphicsDrawTextFunc draw_text,
                                         GtGraphicsDrawTextFunc draw_text_clip,
                                         GtGraphicsDrawTextFunc
                                                     draw_text_centered,
                                         GtGraphicsDrawTextFunc draw_text_right,
                                         GtGraphicsDrawColoredTextFunc
                                                     draw_colored_text,
                                         GtGraphicsGetSingleExtentFunc
                                                     get_text_height,
                                         GtGraphicsGetTextWidthFunc
                                                     get_text_width,
                                         GtGraphicsSetColorFunc
                                                     set_background_color,
                                         GtGraphicsSetFontFunc
                                                     set_font,
                                         GtGraphicsGetSingleExtentFunc
                                                     get_image_width,
                                         GtGraphicsGetSingleExtentFunc
                                                     get_image_height,
                                         GtGraphicsGetSingleExtentFunc
                                                     get_xmargins,
                                         GtGraphicsGetSingleExtentFunc
                                                     get_ymargins,
                                         GtGraphicsSetMarginsFunc
                                                     set_margins,
                                         GtGraphicsDrawLineToFunc
                                                     draw_line,
                                         GtGraphicsDrawLineFunc
                                                     draw_horizontal_line,
                                         GtGraphicsDrawLineFunc
                                                     draw_vertical_line,
                                         GtGraphicsDrawBoxFunc draw_box,
                                         GtGraphicsDrawSimpleFunc draw_dashes,
                                         GtGraphicsDrawSimpleFunc draw_caret,
                                         GtGraphicsDrawRectFunc draw_rectangle,
                                         GtGraphicsDrawArrowheadFunc
                                                     draw_arrowhead,
                                         GtGraphicsDrawCurveDataFunc
                                                     draw_curve,
                                         GtGraphicsSaveToFileFunc save_to_file,
                                         GtGraphicsSaveToStreamFunc
                                                     save_to_stream,
                                         GtGraphicsFreeFunc free);
GtGraphics* gt_graphics_create(const GtGraphicsClass*);
void*       gt_graphics_cast(const GtGraphicsClass*, GtGraphics*);

#endif
