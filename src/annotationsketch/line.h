/*
  Copyright (c) 2007      Christin Schaerfer <schaerfer@zbh.uni-hamburg.de>
  Copyright (c)      2008 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
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

#ifndef LINE_H
#define LINE_H

/* A GtLine contains GtBlock objects. */
typedef struct GtLine GtLine;

#include "annotationsketch/block.h"
#include "annotationsketch/canvas.h"
#include "annotationsketch/style.h"
#include "core/array.h"
#include "core/error.h"

GtLine*   gt_line_new(void);
void      gt_line_insert_block(GtLine*, GtBlock*); /* takes ownership */
bool      gt_line_has_captions(const GtLine*);
GtArray*  gt_line_get_blocks(GtLine*);
int       gt_line_sketch(GtLine*, GtCanvas*, GtError*);
int       gt_line_get_height(GtLine *line, double *height, const GtStyle *sty,
                             GtError *err);

int       gt_line_unit_test(GtError*);
void      gt_line_delete(GtLine*);

#endif
