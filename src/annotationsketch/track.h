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

#ifndef TRACK_H
#define TRACK_H

/* A track has a title and a type und contains line objects. */
typedef struct Track Track;

#include "annotationsketch/canvas.h"
#include "annotationsketch/line.h"
#include "annotationsketch/line_breaker.h"
#include "annotationsketch.h"
#include "core/array.h"
#include "extended/genome_feature_type.h"
#include "extended/genome_node.h"

Track*        track_new(Str *title, unsigned long max_num_lines,
                        bool split_lines, LineBreaker*);
void          track_insert_block(Track*, GT_Block*);
Str*          track_get_title(const Track*);
unsigned long track_get_number_of_lines(const Track*);
unsigned long track_get_number_of_lines_with_captions(const Track *track);
unsigned long track_get_number_of_discarded_blocks(Track *track);
int           track_sketch(Track*, Canvas*);
int           track_unit_test(Error*);
void          track_delete(Track*);

#endif
