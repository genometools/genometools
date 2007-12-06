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

#ifndef TRACK_H
#define TRACK_H

#include "libgtcore/str.h"
#include "libgtcore/array.h"
#include "libgtext/genome_node.h"
#include "libgtview/config.h"
#include "libgtview/line.h"

/* A track has a title and a type und contains line objects. */
typedef struct Track Track;

Track* track_new(Str* title);
void   track_insert_block(Track*, Block*);
Str*   track_get_title(const Track*);
/* Returns Array containing Pointers to Line objects. */
Array* track_get_lines(const Track*);
int    track_get_number_of_lines(const Track*);
int    track_unit_test(Error*);
void   track_delete(Track*);

#endif
