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

#include "libgtcore/array.h"
#include "libgtcore/countingsort.h"
#include "libgtcore/ma.h"
#include "libgtcore/unused.h"
#include "libgtexercise/fragment_overlaps.h"
#include "libgtext/string_matching.h"

struct FragmentOverlaps{
  Array *overlaps;
};

static void determine_overlaps(Array *overlaps, Bioseq *fragments,
                               unsigned long i, unsigned long j,
                               unsigned long minlength)
{
  unsigned long cpl; /* common prefix length */
  Overlap overlap;
  assert(overlaps && fragments && i != j);
  cpl = string_matching_kmp(bioseq_get_sequence(fragments, i),
                            bioseq_get_sequence_length(fragments, i),
                            bioseq_get_sequence(fragments, j),
                            bioseq_get_sequence_length(fragments, j),
                            NULL, NULL);
  if (cpl >= minlength) {
    overlap.start = i;
    overlap.weight = cpl;
    overlap.end = j;
    array_add(overlaps, overlap);
  }
  cpl = string_matching_kmp(bioseq_get_sequence(fragments, j),
                            bioseq_get_sequence_length(fragments, j),
                            bioseq_get_sequence(fragments, i),
                            bioseq_get_sequence_length(fragments, i),
                            NULL, NULL);
  if (cpl >= minlength) {
    overlap.start = j;
    overlap.weight = cpl;
    overlap.end = i;
    array_add(overlaps, overlap);
  }
}

FragmentOverlaps* fragment_overlaps_new(Bioseq *fragments,
                                        unsigned long minlength)
{
  FragmentOverlaps *fo;
  unsigned long i, j;
  assert(fragments);
  fo = ma_malloc(sizeof *fo);
  fo->overlaps = array_new(sizeof (Overlap));
  for (i = 0; i < bioseq_number_of_sequences(fragments); i++) {
    for (j = i + 1; j < bioseq_number_of_sequences(fragments); j++)
      determine_overlaps(fo->overlaps, fragments, i, j, minlength);
  }
  return fo;
}

void fragment_overlaps_delete(FragmentOverlaps *fo)
{
  if (!fo) return;
  array_delete(fo->overlaps);
  ma_free(fo);
}

static unsigned long get_weight(const void *elem, UNUSED void *data)
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
  num_of_overlaps = array_size(fo->overlaps);
  sorted_overlaps = ma_malloc(sizeof (Overlap) * num_of_overlaps);
  /* sort overlaps by weight */
  countingsort(sorted_overlaps, array_get_space(fo->overlaps),
               array_elem_size(fo->overlaps), array_size(fo->overlaps),
               countingsort_get_max(array_get_space(fo->overlaps),
                                    array_elem_size(fo->overlaps),
                                    array_size(fo->overlaps), NULL, get_weight),
               NULL, get_weight);
  array_reset(fo->overlaps);
  for (i = 0; i < num_of_overlaps; i++)
    array_add(fo->overlaps, sorted_overlaps[i]);
  ma_free(sorted_overlaps);
}

bool fragment_overlaps_are_sorted(const FragmentOverlaps *fo)
{
  Overlap *overlap_a, *overlap_b;
  unsigned long i;
  assert(fo);
  for (i = 1; i < array_size(fo->overlaps); i++) {
    overlap_a = array_get(fo->overlaps, i-1);
    overlap_b = array_get(fo->overlaps, i);
    if (overlap_b->weight < overlap_a->weight) return false;
  }
  return true;
}

void fragment_overlaps_show(const FragmentOverlaps *fo)
{
  Overlap *overlap;
  unsigned long i;
  assert(fo);
  for (i = 0; i < array_size(fo->overlaps); i++) {
    overlap = array_get(fo->overlaps, i);
    printf("%lu %lu %lu\n", overlap->start+1, overlap->weight, overlap->end+1);
  }
}

const Overlap* fragment_overlaps_get(const FragmentOverlaps *fo,
                                     unsigned long fragnum)
{
  assert(fo && fragnum < array_size(fo->overlaps));
  return array_get(fo->overlaps, fragnum);
}

unsigned long fragment_overlaps_size(const FragmentOverlaps *fo)
{
  assert(fo);
  return array_size(fo->overlaps);
}
