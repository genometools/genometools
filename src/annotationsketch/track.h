/*
  Copyright (c) 2007      Christin Schaerfer <schaerfer@zbh.uni-hamburg.de>
  Copyright (c) 2008-2009 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2007-2009 Center for Bioinformatics, University of Hamburg

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

/* A <GtTrack> acts as a container for <GtLine> objects. */
typedef struct GtTrack GtTrack;

#include "annotationsketch/canvas.h"
#include "annotationsketch/line_breaker.h"
#include "core/array.h"

GtTrack*      gt_track_new(GtStr *title, unsigned long max_num_lines,
                           bool split_lines, GtLineBreaker *lb);
int           gt_track_insert_block(GtTrack*, GtBlock*, GtError*);
GtStr*        gt_track_get_title(const GtTrack*);
unsigned long gt_track_get_number_of_discarded_blocks(GtTrack *track);
int           gt_track_sketch(GtTrack*, GtCanvas*, GtError*);
int           gt_track_get_height(const GtTrack *track, double *height,
                                  const GtStyle *sty, GtError *err);
unsigned long gt_track_get_y_index(const GtTrack *track);
void          gt_track_set_y_index(GtTrack *track, unsigned long y_index);
void          gt_track_delete(GtTrack*);

int           gt_track_unit_test(GtError*);

#endif
