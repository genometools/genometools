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

#include "core/ma.h"
#include "extended/transcript_bittabs.h"

struct TranscriptGtBittabs {
  GtBittab *gt_bittab_all,
         *gt_bittab_single,
         *gt_bittab_initial,
         *gt_bittab_internal,
         *gt_bittab_terminal;
};

TranscriptGtBittabs* transcript_bittabs_new(unsigned long size_all,
                                          unsigned long size_single,
                                          unsigned long size_initial,
                                          unsigned long size_internal,
                                          unsigned long size_terminal)
{
  TranscriptGtBittabs *tb = gt_calloc(1, sizeof (TranscriptGtBittabs));
  if (size_all) tb->gt_bittab_all = gt_bittab_new(size_all);
  if (size_single) tb->gt_bittab_single = gt_bittab_new(size_single);
  if (size_initial) tb->gt_bittab_initial = gt_bittab_new(size_initial);
  if (size_internal) tb->gt_bittab_internal = gt_bittab_new(size_internal);
  if (size_terminal) tb->gt_bittab_terminal = gt_bittab_new(size_all);
  return tb;
}

GtBittab* transcript_bittabs_get_all(const TranscriptGtBittabs *tb)
{
  assert(tb);
  return tb->gt_bittab_all;
}

GtBittab* transcript_bittabs_get_single(const TranscriptGtBittabs *tb)
{
  assert(tb);
  return tb->gt_bittab_single;
}

GtBittab* transcript_bittabs_get_initial(const TranscriptGtBittabs *tb)
{
  assert(tb);
  return tb->gt_bittab_initial;
}

GtBittab* transcript_bittabs_get_internal(const TranscriptGtBittabs *tb)
{
  assert(tb);
  return tb->gt_bittab_internal;
}

GtBittab* transcript_bittabs_get_terminal(const TranscriptGtBittabs *tb)
{
  assert(tb);
  return tb->gt_bittab_terminal;
}

void transcript_bittabs_delete(TranscriptGtBittabs *tb)
{
  if (!tb) return;
  gt_bittab_delete(tb->gt_bittab_all);
  gt_bittab_delete(tb->gt_bittab_single);
  gt_bittab_delete(tb->gt_bittab_initial);
  gt_bittab_delete(tb->gt_bittab_internal);
  gt_bittab_delete(tb->gt_bittab_terminal);
  gt_free(tb);
}
