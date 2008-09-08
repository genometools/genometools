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
#include "core/range.h"
#include "extended/transcript_exons.h"

struct TranscriptExons {
  GT_Array *exon_array_all,
        *exon_array_single,
        *exon_array_initial,
        *exon_array_internal,
        *exon_array_terminal;
};

TranscriptExons* transcript_exons_new(void)
{
  TranscriptExons *te = gt_malloc(sizeof (TranscriptExons));
  te->exon_array_all = gt_array_new(sizeof (GT_Range));
  te->exon_array_single = gt_array_new(sizeof (GT_Range));
  te->exon_array_initial = gt_array_new(sizeof (GT_Range));
  te->exon_array_internal = gt_array_new(sizeof (GT_Range));
  te->exon_array_terminal = gt_array_new(sizeof (GT_Range));
  return te;
}

GT_Array* transcript_exons_get_all(const TranscriptExons *te)
{
  assert(te);
  return te->exon_array_all;
}

GT_Array* transcript_exons_get_single(const TranscriptExons *te)
{
  assert(te);
  return te->exon_array_single;
}

GT_Array* transcript_exons_get_initial(const TranscriptExons *te)
{
  assert(te);
  return te->exon_array_initial;
}

GT_Array* transcript_exons_get_internal(const TranscriptExons *te)
{
  assert(te);
  return te->exon_array_internal;
}

GT_Array* transcript_exons_get_terminal(const TranscriptExons *te)
{
  assert(te);
  return te->exon_array_terminal;
}

void transcript_exons_sort(const TranscriptExons *te)
{
  assert(te);
  ranges_sort(te->exon_array_all);
  ranges_sort(te->exon_array_single);
  ranges_sort(te->exon_array_initial);
  ranges_sort(te->exon_array_internal);
  ranges_sort(te->exon_array_terminal);
}

TranscriptCounts* transcript_exons_uniq_in_place_count(TranscriptExons *te)
{
  TranscriptCounts *tc;
  GT_Array *counts;
  assert(te);
  tc = transcript_counts_new();
  counts = ranges_uniq_in_place_count(te->exon_array_all);
  transcript_counts_set_all(tc, counts);
  counts = ranges_uniq_in_place_count(te->exon_array_single);
  transcript_counts_set_single(tc, counts);
  counts = ranges_uniq_in_place_count(te->exon_array_initial);
  transcript_counts_set_initial(tc, counts);
  counts = ranges_uniq_in_place_count(te->exon_array_internal);
  transcript_counts_set_internal(tc, counts);
  counts = ranges_uniq_in_place_count(te->exon_array_terminal);
  transcript_counts_set_terminal(tc, counts);
  return tc;
}

bool transcript_exons_are_sorted(const TranscriptExons *te)
{
  assert(te);
  if (!ranges_are_sorted(te->exon_array_all)) return false;
  if (!ranges_are_sorted(te->exon_array_single)) return false;
  if (!ranges_are_sorted(te->exon_array_initial)) return false;
  if (!ranges_are_sorted(te->exon_array_internal)) return false;
  if (!ranges_are_sorted(te->exon_array_terminal)) return false;
  return true;
}

TranscriptBittabs* transcript_exons_create_bittabs(const TranscriptExons *te)
{
  assert(te);
  return transcript_bittabs_new(gt_array_size(te->exon_array_all),
                                gt_array_size(te->exon_array_single),
                                gt_array_size(te->exon_array_initial),
                                gt_array_size(te->exon_array_internal),
                                gt_array_size(te->exon_array_terminal));
}

void transcript_exons_delete(TranscriptExons *te)
{
  if (!te) return;
  gt_array_delete(te->exon_array_all);
  gt_array_delete(te->exon_array_single);
  gt_array_delete(te->exon_array_initial);
  gt_array_delete(te->exon_array_internal);
  gt_array_delete(te->exon_array_terminal);
  gt_free(te);
}
