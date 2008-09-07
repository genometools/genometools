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

#include "core/ensure.h"
#include "core/hashtable.h"
#include "core/ma.h"
#include "core/undef.h"
#include "core/unused.h"
#include "annotationsketch/line_breaker_bases.h"
#include "annotationsketch/style.h"
#include "annotationsketch/track.h"

struct Track {
  Str *title;
  unsigned long max_num_lines, discarded_blocks;
  LineBreaker *lb;
  bool split;
  GT_Array *lines;
};

Track* track_new(Str *title, unsigned long max_num_lines, bool split,
                 LineBreaker *lb)
{
  Track *track;
  assert(title && lb);
  track = ma_calloc(1, sizeof (Track));
  assert(track);
  track->title = str_ref(title);
  track->lines = gt_array_new(sizeof (Line*));
  track->max_num_lines = max_num_lines;
  track->split = split;
  track->lb = lb;
  return track;
}

static Line* get_next_free_line(Track *track, GT_Block *block)
{
  unsigned long i;
  Line* line;

  assert(track);

  /* find unoccupied line -- may need optimisation */
  for (i = 0; i < gt_array_size(track->lines); i++) {
    line = *(Line**) gt_array_get(track->lines, i);
    if (!line_breaker_line_is_occupied(track->lb, line, block))
      return line;
  }
  /* if line limit is hit, do not create any more lines! */
  if (track->max_num_lines != UNDEF_ULONG
       && gt_array_size(track->lines) == track->max_num_lines)
  {
    track->discarded_blocks++;
    return NULL;
  }
  /* make sure there is only one line if 'split_lines' is set to false */
  if (!track->split)
  {
    if (gt_array_size(track->lines) < 1)
    {
      line = line_new();
      gt_array_add(track->lines, line);
    }
    else
      line = *(Line**) gt_array_get(track->lines, 0);
    assert(gt_array_size(track->lines) == 1);
  }
  else
  {
    line = line_new();
    gt_array_add(track->lines, line);
  }
  assert(line);
  return line;
}

unsigned long track_get_number_of_discarded_blocks(Track *track)
{
  assert(track);
  return track->discarded_blocks;
}

void track_insert_block(Track *track, GT_Block *block)
{
  Line *line;

  assert(track && block);
  line = get_next_free_line(track, block);
  block = gt_block_ref(block);
  if (line)
  {
    line_insert_block(line, block);
    line_breaker_register_block(track->lb, line, block);
  } else gt_block_delete(block);
}

Str* track_get_title(const Track *track)
{
  assert(track && track->title);
  return track->title;
}

unsigned long track_get_number_of_lines(const Track *track)
{
  assert(track);
  return gt_array_size(track->lines);
}

unsigned long track_get_number_of_lines_with_captions(const Track *track)
{
  unsigned long i = 0, nof_tracks = 0;
  assert(track);
  for (i = 0; i < gt_array_size(track->lines); i++) {
    if (line_has_captions(*(Line**) gt_array_get(track->lines, i)))
      nof_tracks++;
  }
  return nof_tracks;
}

int track_sketch(Track* track, GT_Canvas *canvas)
{
  int i = 0;
  assert(track && canvas);
  gt_canvas_visit_track_pre(canvas, track);
  for (i = 0; i < gt_array_size(track->lines); i++)
    line_sketch(*(Line**) gt_array_get(track->lines, i), canvas);
  gt_canvas_visit_track_post(canvas, track);
  return 0;
}

int track_unit_test(GT_Error *err)
{
  int had_err = 0;
  GT_Block *b1, *b2, *b3, *b4;
  Range r1, r2, r3, r4;
  Track *track;
  Str *title;
  error_check(err);
  LineBreaker *lb;

  title = str_new_cstr("test");

  r1.start=100UL;  r1.end=1000UL;
  r2.start=1001UL; r2.end=1500UL;
  r3.start=700UL;  r3.end=1200UL;
  r4.start=10UL;   r4.end=200UL;

  b1 = gt_block_new();
  gt_block_set_range(b1, r1);
  b2 = gt_block_new();
  gt_block_set_range(b2, r2);
  b3 = gt_block_new();
  gt_block_set_range(b3, r3);
  b4 = gt_block_new();
  gt_block_set_range(b4, r4);

  lb = line_breaker_bases_new();

  track = track_new(title, UNDEF_ULONG, true, lb);
  ensure(had_err, track);
  ensure(had_err, track_get_title(track) == title);

  ensure(had_err, track_get_number_of_lines(track) == 0);
  track_insert_block(track, b1);
  ensure(had_err, track_get_number_of_lines(track) == 1);
  track_insert_block(track, b2);
  ensure(had_err, track_get_number_of_lines(track) == 1);
  track_insert_block(track, b3);
  ensure(had_err, track_get_number_of_lines(track) == 2);
  track_insert_block(track, b4);
  ensure(had_err, track_get_number_of_lines(track) == 2);

  track_delete(track);
  str_delete(title);
  gt_block_delete(b1);
  gt_block_delete(b2);
  gt_block_delete(b3);
  gt_block_delete(b4);

  return had_err;
}

void track_delete(Track *track)
{
  unsigned long i;
  if (!track) return;
  if (track->lb)
    line_breaker_delete(track->lb);
  for (i = 0; i < gt_array_size(track->lines); i++)
    line_delete(*(Line**) gt_array_get(track->lines, i));
  gt_array_delete(track->lines);
  str_delete(track->title);
  ma_free(track);
}
