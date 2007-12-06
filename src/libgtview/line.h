/*
  Copyright (c) 2007 Christin Schaerfer <cschaerfer@stud.zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

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

#include "libgtcore/array.h"
#include "libgtext/genome_node.h"
#include "libgtview/config.h"
#include "libgtview/block.h"

/* A line contains block objects. */
typedef struct Line Line;

Line*  line_new(void);
void   line_insert_block(Line*, Block*); /* takes ownership */
bool   line_is_occupied(const Line*, Range);
/* Returns Array containing Pointers to Block objects. */
Array* line_get_blocks(Line*);
int    line_unit_test(Error*);
void   line_delete(Line*);

#endif
