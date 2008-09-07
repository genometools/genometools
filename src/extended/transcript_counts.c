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
#include "extended/transcript_counts.h"

struct TranscriptCounts {
  GT_Array *exon_array_all,
        *exon_array_single,
        *exon_array_initial,
        *exon_array_internal,
        *exon_array_terminal;
};

TranscriptCounts* transcript_counts_new(void)
{
  TranscriptCounts *tc = ma_calloc(1, sizeof (TranscriptCounts));
  return tc;
}

GT_Array* transcript_counts_get_all(const TranscriptCounts *tc)
{
  assert(tc);
  return tc->exon_array_all;
}

void transcript_counts_set_all(TranscriptCounts *tc, GT_Array *counts)
{
  assert(tc && counts);
  tc->exon_array_all = counts;
}

GT_Array* transcript_counts_get_single(const TranscriptCounts *tc)
{
  assert(tc);
  return tc->exon_array_single;
}

void transcript_counts_set_single(TranscriptCounts *tc, GT_Array *counts)
{
  assert(tc && counts);
  tc->exon_array_single = counts;
}

GT_Array* transcript_counts_get_initial(const TranscriptCounts *tc)
{
  assert(tc);
  return tc->exon_array_initial;
}

void transcript_counts_set_initial(TranscriptCounts *tc, GT_Array *counts)
{
  assert(tc && counts);
  tc->exon_array_initial = counts;
}

GT_Array* transcript_counts_get_internal(const TranscriptCounts *tc)
{
  assert(tc);
  return tc->exon_array_internal;
}

void transcript_counts_set_internal(TranscriptCounts *tc, GT_Array *counts)
{
  assert(tc && counts);
  tc->exon_array_internal = counts;
}

GT_Array* transcript_counts_get_terminal(const TranscriptCounts *tc)
{
  assert(tc);
  return tc->exon_array_terminal;
}

void transcript_counts_set_terminal(TranscriptCounts *tc, GT_Array *counts)
{
  assert(tc && counts);
  tc->exon_array_terminal = counts;
}

void transcript_counts_delete(TranscriptCounts *tc)
{
  if (!tc) return;
  gt_array_delete(tc->exon_array_all);
  gt_array_delete(tc->exon_array_single);
  gt_array_delete(tc->exon_array_initial);
  gt_array_delete(tc->exon_array_internal);
  gt_array_delete(tc->exon_array_terminal);
  ma_free(tc);
}
