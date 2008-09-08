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

#ifndef CANVAS_REP_H
#define CANVAS_REP_H

#include <stdio.h>
#include "core/bittab.h"
#include "core/range.h"
#include "annotationsketch/canvas.h"
#include "annotationsketch/graphics.h"
#include "annotationsketch/image_info.h"

struct GT_CanvasClass {
  size_t size;
  int           (*visit_diagram_pre)(GT_Canvas*, GT_Diagram*);
  int           (*visit_diagram_post)(GT_Canvas*, GT_Diagram*);
  void          (*free)(GT_Canvas*);
};

struct GT_Canvas {
  const GT_CanvasClass *c_class;
  GT_Range viewrange;
  double factor, y, margins;
  unsigned long width, height;
  GT_Style *sty;
  bool show_track_captions;
  Bittab *bt;
  Graphics *g;
  GT_ImageInfo *ii;
};

GT_Canvas* gt_canvas_create(const GT_CanvasClass*);
void*      gt_canvas_cast(const GT_CanvasClass*, GT_Canvas*);

#endif
