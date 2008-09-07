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

#ifndef LINE_BREAKER_H
#define LINE_BREAKER_H

/* The GT_LineBreaker interface. */
typedef struct GT_LineBreakerClass GT_LineBreakerClass;
typedef struct GT_LineBreaker GT_LineBreaker;

#include "annotationsketch/block.h"
#include "annotationsketch/line.h"

GT_LineBreaker* gt_line_breaker_ref(GT_LineBreaker*);
bool         gt_line_breaker_gt_line_is_occupied(GT_LineBreaker *lb,
                                           GT_Line *line, GT_Block *block);
void         gt_line_breaker_register_block(GT_LineBreaker *lb,
                                         GT_Line *line, GT_Block *block);
void         gt_line_breaker_delete(GT_LineBreaker*);

#endif
