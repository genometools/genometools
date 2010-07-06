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

struct GtTranscriptExons {
  GtArray *exon_array_all,
          *exon_array_single,
          *exon_array_initial,
          *exon_array_internal,
          *exon_array_terminal;
};

GtTranscriptExons* gt_transcript_exons_new(void)
{
  GtTranscriptExons *te = gt_malloc(sizeof (GtTranscriptExons));
  te->exon_array_all = gt_array_new(sizeof (GtRange));
  te->exon_array_single = gt_array_new(sizeof (GtRange));
  te->exon_array_initial = gt_array_new(sizeof (GtRange));
  te->exon_array_internal = gt_array_new(sizeof (GtRange));
  te->exon_array_terminal = gt_array_new(sizeof (GtRange));
  return te;
}

GtArray* gt_transcript_exons_get_all(const GtTranscriptExons *te)
{
  gt_assert(te);
  return te->exon_array_all;
}

GtArray* gt_transcript_exons_get_single(const GtTranscriptExons *te)
{
  gt_assert(te);
  return te->exon_array_single;
}

GtArray* gt_transcript_exons_get_initial(const GtTranscriptExons *te)
{
  gt_assert(te);
  return te->exon_array_initial;
}

GtArray* gt_transcript_exons_get_internal(const GtTranscriptExons *te)
{
  gt_assert(te);
  return te->exon_array_internal;
}

GtArray* gt_transcript_exons_get_terminal(const GtTranscriptExons *te)
{
  gt_assert(te);
  return te->exon_array_terminal;
}

void gt_transcript_exons_sort(const GtTranscriptExons *te)
{
  gt_assert(te);
  gt_ranges_sort(te->exon_array_all);
  gt_ranges_sort(te->exon_array_single);
  gt_ranges_sort(te->exon_array_initial);
  gt_ranges_sort(te->exon_array_internal);
  gt_ranges_sort(te->exon_array_terminal);
}

GtTranscriptCounts* gt_transcript_exons_uniq_in_place_count(GtTranscriptExons
                                                                         *te)
{
  GtTranscriptCounts *tc;
  GtArray *counts;
  gt_assert(te);
  tc = gt_transcript_counts_new();
  counts = gt_ranges_uniq_in_place_count(te->exon_array_all);
  gt_transcript_counts_set_all(tc, counts);
  counts = gt_ranges_uniq_in_place_count(te->exon_array_single);
  gt_transcript_counts_set_single(tc, counts);
  counts = gt_ranges_uniq_in_place_count(te->exon_array_initial);
  gt_transcript_counts_set_initial(tc, counts);
  counts = gt_ranges_uniq_in_place_count(te->exon_array_internal);
  gt_transcript_counts_set_internal(tc, counts);
  counts = gt_ranges_uniq_in_place_count(te->exon_array_terminal);
  gt_transcript_counts_set_terminal(tc, counts);
  return tc;
}

bool gt_transcript_exons_are_sorted(const GtTranscriptExons *te)
{
  gt_assert(te);
  if (!gt_ranges_are_sorted(te->exon_array_all)) return false;
  if (!gt_ranges_are_sorted(te->exon_array_single)) return false;
  if (!gt_ranges_are_sorted(te->exon_array_initial)) return false;
  if (!gt_ranges_are_sorted(te->exon_array_internal)) return false;
  if (!gt_ranges_are_sorted(te->exon_array_terminal)) return false;
  return true;
}

GtTranscriptBittabs* gt_transcript_exons_create_bittabs(const
                                                        GtTranscriptExons *te)
{
  gt_assert(te);
  return gt_transcript_bittabs_new(gt_array_size(te->exon_array_all),
                                   gt_array_size(te->exon_array_single),
                                   gt_array_size(te->exon_array_initial),
                                   gt_array_size(te->exon_array_internal),
                                   gt_array_size(te->exon_array_terminal));
}

void gt_transcript_exons_delete(GtTranscriptExons *te)
{
  if (!te) return;
  gt_array_delete(te->exon_array_all);
  gt_array_delete(te->exon_array_single);
  gt_array_delete(te->exon_array_initial);
  gt_array_delete(te->exon_array_internal);
  gt_array_delete(te->exon_array_terminal);
  gt_free(te);
}
