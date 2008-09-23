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

struct GtTranscriptCounts {
  GtArray *exon_array_all,
        *exon_array_single,
        *exon_array_initial,
        *exon_array_internal,
        *exon_array_terminal;
};

GtTranscriptCounts* gt_transcript_counts_new(void)
{
  GtTranscriptCounts *tc = gt_calloc(1, sizeof (GtTranscriptCounts));
  return tc;
}

GtArray* gt_transcript_counts_get_all(const GtTranscriptCounts *tc)
{
  gt_assert(tc);
  return tc->exon_array_all;
}

void gt_transcript_counts_set_all(GtTranscriptCounts *tc, GtArray *counts)
{
  gt_assert(tc && counts);
  tc->exon_array_all = counts;
}

GtArray* gt_transcript_counts_get_single(const GtTranscriptCounts *tc)
{
  gt_assert(tc);
  return tc->exon_array_single;
}

void gt_transcript_counts_set_single(GtTranscriptCounts *tc, GtArray *counts)
{
  gt_assert(tc && counts);
  tc->exon_array_single = counts;
}

GtArray* gt_transcript_counts_get_initial(const GtTranscriptCounts *tc)
{
  gt_assert(tc);
  return tc->exon_array_initial;
}

void gt_transcript_counts_set_initial(GtTranscriptCounts *tc, GtArray *counts)
{
  gt_assert(tc && counts);
  tc->exon_array_initial = counts;
}

GtArray* gt_transcript_counts_get_internal(const GtTranscriptCounts *tc)
{
  gt_assert(tc);
  return tc->exon_array_internal;
}

void gt_transcript_counts_set_internal(GtTranscriptCounts *tc, GtArray *counts)
{
  gt_assert(tc && counts);
  tc->exon_array_internal = counts;
}

GtArray* gt_transcript_counts_get_terminal(const GtTranscriptCounts *tc)
{
  gt_assert(tc);
  return tc->exon_array_terminal;
}

void gt_transcript_counts_set_terminal(GtTranscriptCounts *tc, GtArray *counts)
{
  gt_assert(tc && counts);
  tc->exon_array_terminal = counts;
}

void gt_transcript_counts_delete(GtTranscriptCounts *tc)
{
  if (!tc) return;
  gt_array_delete(tc->exon_array_all);
  gt_array_delete(tc->exon_array_single);
  gt_array_delete(tc->exon_array_initial);
  gt_array_delete(tc->exon_array_internal);
  gt_array_delete(tc->exon_array_terminal);
  gt_free(tc);
}
