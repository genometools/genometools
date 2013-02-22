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

#ifndef CANVAS_CAIRO_CONTEXT_H
#define CANVAS_CAIRO_CONTEXT_H

#include <cairo.h>
#include "annotationsketch/canvas.h"
#include "annotationsketch/canvas_cairo_context_api.h"

const GtCanvasClass* gt_canvas_cairo_context_class(void);

#define canvas_cairo_context_cast(C)\
        gt_canvas_cast(gt_canvas_cairo_context_class(), C)

#define canvas_cairo_context_try_cast(C)\
        gt_canvas_try_cast(gt_canvas_cairo_context_class(), C)

#endif
