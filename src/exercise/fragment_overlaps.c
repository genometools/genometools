/*
  Copyright (c) 2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2008 Center for Bioinformatics, University of Hamburg

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

#include "core/array.h"
#include "core/countingsort.h"
#include "core/ma.h"
#include "core/unused_api.h"
#include "exercise/fragment_overlaps.h"
#include "extended/string_matching.h"

struct FragmentOverlaps{
  GtArray *overlaps;
};

static void determine_overlaps(GtArray *overlaps, GT_Bioseq *fragments,
                               unsigned long i, unsigned long j,
                               unsigned long minlength)
{
  unsigned long cpl; /* common prefix length */
  Overlap overlap;
  assert(overlaps && fragments && i != j);
  cpl = string_matching_kmp(gt_bioseq_get_sequence(fragments, i),
                            gt_bioseq_get_sequence_length(fragments, i),
                            gt_bioseq_get_sequence(fragments, j),
                            gt_bioseq_get_sequence_length(fragments, j),
                            NULL, NULL);
  if (cpl >= minlength) {
    overlap.start = i;
    overlap.weight = cpl;
    overlap.end = j;
    gt_array_add(overlaps, overlap);
  }
  cpl = string_matching_kmp(gt_bioseq_get_sequence(fragments, j),
                            gt_bioseq_get_sequence_length(fragments, j),
                            gt_bioseq_get_sequence(fragments, i),
                            gt_bioseq_get_sequence_length(fragments, i),
                            NULL, NULL);
  if (cpl >= minlength) {
    overlap.start = j;
    overlap.weight = cpl;
    overlap.end = i;
    gt_array_add(overlaps, overlap);
  }
}

FragmentOverlaps* fragment_overlaps_new(GT_Bioseq *fragments,
                                        unsigned long minlength)
{
  FragmentOverlaps *fo;
  unsigned long i, j;
  assert(fragments);
  fo = gt_malloc(sizeof *fo);
  fo->overlaps = gt_array_new(sizeof (Overlap));
  for (i = 0; i < gt_bioseq_number_of_sequences(fragments); i++) {
    for (j = i + 1; j < gt_bioseq_number_of_sequences(fragments); j++)
      determine_overlaps(fo->overlaps, fragments, i, j, minlength);
  }
  return fo;
}

void fragment_overlaps_delete(FragmentOverlaps *fo)
{
  if (!fo) return;
  gt_array_delete(fo->overlaps);
  gt_free(fo);
}

static unsigned long get_weight(const void *elem, GT_UNUSED void *data)
{
  const Overlap *overlap = elem;
  assert(overlap);
  return overlap->weight;
}

void fragment_overlaps_sort(FragmentOverlaps *fo)
{
  Overlap *sorted_overlaps;
  unsigned long i, num_of_overlaps;
  assert(fo);
  num_of_overlaps = gt_array_size(fo->overlaps);
  sorted_overlaps = gt_malloc(sizeof (Overlap) * num_of_overlaps);
  /* sort overlaps by weight */
  gt_countingsort(sorted_overlaps, gt_array_get_space(fo->overlaps),
                  gt_array_elem_size(fo->overlaps), gt_array_size(fo->overlaps),
                  gt_countingsort_get_max(gt_array_get_space(fo->overlaps),
                                          gt_array_elem_size(fo->overlaps),
                                          gt_array_size(fo->overlaps), NULL,
                                          get_weight),
                  NULL, get_weight);
  gt_array_reset(fo->overlaps);
  for (i = 0; i < num_of_overlaps; i++)
    gt_array_add(fo->overlaps, sorted_overlaps[i]);
  gt_free(sorted_overlaps);
}

bool fragment_overlaps_are_sorted(const FragmentOverlaps *fo)
{
  Overlap *overlap_a, *overlap_b;
  unsigned long i;
  assert(fo);
  for (i = 1; i < gt_array_size(fo->overlaps); i++) {
    overlap_a = gt_array_get(fo->overlaps, i-1);
    overlap_b = gt_array_get(fo->overlaps, i);
    if (overlap_b->weight < overlap_a->weight) return false;
  }
  return true;
}

void fragment_overlaps_show(const FragmentOverlaps *fo)
{
  Overlap *overlap;
  unsigned long i;
  assert(fo);
  for (i = 0; i < gt_array_size(fo->overlaps); i++) {
    overlap = gt_array_get(fo->overlaps, i);
    printf("%lu %lu %lu\n", overlap->start+1, overlap->weight, overlap->end+1);
  }
}

const Overlap* fragment_overlaps_get(const FragmentOverlaps *fo,
                                     unsigned long fragnum)
{
  assert(fo && fragnum < gt_array_size(fo->overlaps));
  return gt_array_get(fo->overlaps, fragnum);
}

unsigned long fragment_overlaps_size(const FragmentOverlaps *fo)
{
  assert(fo);
  return gt_array_size(fo->overlaps);
}
