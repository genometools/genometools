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

#ifndef CANVAS_H
#define CANVAS_H

typedef struct Canvas Canvas;

#include "libgtview/block.h"
#include "libgtview/config.h"
#include "libgtview/diagram.h"
#include "libgtview/element.h"
#include "libgtview/line.h"
#include "libgtview/imageinfo.h"
#include "libgtview/track.h"

Canvas*      canvas_new(Config*, unsigned int width, ImageInfo*);
unsigned int canvas_get_height(Canvas*);
int          canvas_visit_diagram(Canvas*, Diagram*);
int          canvas_visit_track(Canvas*, Track*);
int          canvas_visit_line_pre(Canvas*, Line*);
int          canvas_visit_line_post(Canvas*, Line*);
int          canvas_visit_block(Canvas*, Block*);
int          canvas_visit_element(Canvas*, Element*);
int          canvas_to_png(Canvas*, const char*, Error*);
void         canvas_delete(Canvas*);

#endif
