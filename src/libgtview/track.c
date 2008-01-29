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

#include "libgtcore/ensure.h"
#include "libgtcore/hashtable.h"
#include "libgtcore/ma.h"
#include "libgtview/track.h"

struct Track {
  Str *title;
  Array *lines;
};

Track* track_new(Str *title)
{
  Track *track;
  assert(title);
  track = ma_malloc(sizeof (Track));
  track->title = str_ref(title);
  track->lines = array_new(sizeof (Line*));
  assert(track);
  return track;
}

static Line* get_next_free_line(Track *track, Range r)
{
  unsigned long i;
  Line* line;

  assert(track);

  for (i = 0; i < array_size(track->lines); i++) {
    line = *(Line**) array_get(track->lines, i);
    if (!line_is_occupied(line, r))
      return line;
  }
  line = line_new();
  array_add(track->lines, line);

  assert(line);
  return line;
}

void track_insert_block(Track *track, Block *block)
{
  Range r;
  Line *line;

  assert(track && block);
  r = block_get_range(block);
  line = get_next_free_line(track, r);
  line_insert_block(line, block);
}

Str* track_get_title(const Track *track)
{
  assert(track && track->title);
  return track->title;
}

Array* track_get_lines(const Track *track)
{
  assert(track && track->lines);
  return track->lines;
}

int track_get_number_of_lines(const Track *track)
{
  int nof_tracks;
  assert(track);
  nof_tracks = (int) array_size(track->lines);
  return nof_tracks;
}

int track_unit_test(Error *err)
{
  int had_err = 0;
  Block *b1, *b2, *b3, *b4;
  Range r1, r2, r3, r4;
  Track *track;
  Str *title;
  error_check(err);

  title = str_new_cstr("test");

  r1.start=100UL;  r1.end=1000UL;
  r2.start=1001UL; r2.end=1500UL;
  r3.start=700UL;  r3.end=1200UL;
  r4.start=10UL;   r4.end=200UL;

  b1 = block_new();
  block_set_range(b1, r1);
  b2 = block_new();
  block_set_range(b2, r2);
  b3 = block_new();
  block_set_range(b3, r3);
  b4 = block_new();
  block_set_range(b4, r4);

  track = track_new(title);
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

  return had_err;
}

void track_delete(Track *track)
{
  unsigned long i;
  if (!track) return;
  for (i = 0; i < array_size(track->lines); i++)
    line_delete(*(Line**) array_get(track->lines, i));
  array_delete(track->lines);
  str_delete(track->title);
  ma_free(track);
}
