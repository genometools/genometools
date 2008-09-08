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

struct TranscriptGT_Bittabs {
  GT_Bittab *gt_bittab_all,
         *gt_bittab_single,
         *gt_bittab_initial,
         *gt_bittab_internal,
         *gt_bittab_terminal;
};

TranscriptGT_Bittabs* transcript_bittabs_new(unsigned long size_all,
                                          unsigned long size_single,
                                          unsigned long size_initial,
                                          unsigned long size_internal,
                                          unsigned long size_terminal)
{
  TranscriptGT_Bittabs *tb = gt_calloc(1, sizeof (TranscriptGT_Bittabs));
  if (size_all) tb->gt_bittab_all = gt_bittab_new(size_all);
  if (size_single) tb->gt_bittab_single = gt_bittab_new(size_single);
  if (size_initial) tb->gt_bittab_initial = gt_bittab_new(size_initial);
  if (size_internal) tb->gt_bittab_internal = gt_bittab_new(size_internal);
  if (size_terminal) tb->gt_bittab_terminal = gt_bittab_new(size_all);
  return tb;
}

GT_Bittab* transcript_bittabs_get_all(const TranscriptGT_Bittabs *tb)
{
  assert(tb);
  return tb->gt_bittab_all;
}

GT_Bittab* transcript_bittabs_get_single(const TranscriptGT_Bittabs *tb)
{
  assert(tb);
  return tb->gt_bittab_single;
}

GT_Bittab* transcript_bittabs_get_initial(const TranscriptGT_Bittabs *tb)
{
  assert(tb);
  return tb->gt_bittab_initial;
}

GT_Bittab* transcript_bittabs_get_internal(const TranscriptGT_Bittabs *tb)
{
  assert(tb);
  return tb->gt_bittab_internal;
}

GT_Bittab* transcript_bittabs_get_terminal(const TranscriptGT_Bittabs *tb)
{
  assert(tb);
  return tb->gt_bittab_terminal;
}

void transcript_bittabs_delete(TranscriptGT_Bittabs *tb)
{
  if (!tb) return;
  gt_bittab_delete(tb->gt_bittab_all);
  gt_bittab_delete(tb->gt_bittab_single);
  gt_bittab_delete(tb->gt_bittab_initial);
  gt_bittab_delete(tb->gt_bittab_internal);
  gt_bittab_delete(tb->gt_bittab_terminal);
  gt_free(tb);
}
