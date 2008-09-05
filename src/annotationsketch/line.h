/*
  Copyright (c) 2007      Christin Schaerfer <cschaerfer@zbh.uni-hamburg.de>
  Copyright (c)      2008 Sascha Steinbiss <ssteinbiss@zbh.uni-hamburg.de>
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

/* A line contains block objects. */
typedef struct Line Line;

#include "annotationsketch/block.h"
#include "annotationsketch/canvas.h"
#include "annotationsketch/drawing_range.h"
#include "core/array.h"
#include "core/range.h"
#include "extended/genome_node.h"

Line*  line_new(void);
void   line_insert_block(Line*, Block*); /* takes ownership */
bool   line_has_captions(const Line*);
Array* line_get_blocks(Line*);
int    line_sketch(Line*, Canvas*);
int    line_unit_test(Error*);
void   line_delete(Line*);

#endif
