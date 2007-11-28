/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#include "libgtcore/ma.h"
#include "libgtext/transcript_bittabs.h"

struct TranscriptBittabs {
  Bittab *bittab_all,
         *bittab_single,
         *bittab_initial,
         *bittab_internal,
         *bittab_terminal;
};

TranscriptBittabs* transcript_bittabs_new(unsigned long size_all,
                                          unsigned long size_single,
                                          unsigned long size_initial,
                                          unsigned long size_internal,
                                          unsigned long size_terminal)
{
  TranscriptBittabs *tb = ma_calloc(1, sizeof (TranscriptBittabs));
  if (size_all) tb->bittab_all = bittab_new(size_all);
  if (size_single) tb->bittab_single = bittab_new(size_single);
  if (size_initial) tb->bittab_initial = bittab_new(size_initial);
  if (size_internal) tb->bittab_internal = bittab_new(size_internal);
  if (size_terminal) tb->bittab_terminal = bittab_new(size_all);
  return tb;
}

Bittab* transcript_bittabs_get_all(const TranscriptBittabs *tb)
{
  assert(tb);
  return tb->bittab_all;
}

Bittab* transcript_bittabs_get_single(const TranscriptBittabs *tb)
{
  assert(tb);
  return tb->bittab_single;
}

Bittab* transcript_bittabs_get_initial(const TranscriptBittabs *tb)
{
  assert(tb);
  return tb->bittab_initial;
}

Bittab* transcript_bittabs_get_internal(const TranscriptBittabs *tb)
{
  assert(tb);
  return tb->bittab_internal;
}

Bittab* transcript_bittabs_get_terminal(const TranscriptBittabs *tb)
{
  assert(tb);
  return tb->bittab_terminal;
}

void transcript_bittabs_delete(TranscriptBittabs *tb)
{
  if (!tb) return;
  bittab_delete(tb->bittab_all);
  bittab_delete(tb->bittab_single);
  bittab_delete(tb->bittab_initial);
  bittab_delete(tb->bittab_internal);
  bittab_delete(tb->bittab_terminal);
  ma_free(tb);
}
