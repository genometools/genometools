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

#include "core/cstr_api.h"
#include "core/ma.h"
#include "extended/sspliced_alignment.h"

struct GtSSplicedAlignment {
  char *id;
  bool forward;
  GtArray *exons; /* the exon ranges */
};

GtSSplicedAlignment* gt_sspliced_alignment_new(const char *id, bool forward)
{
  GtSSplicedAlignment *sa;
  gt_assert(id);
  sa = gt_malloc(sizeof *sa);
  sa->id = gt_cstr_dup(id);
  sa->forward = forward;
  sa->exons = gt_array_new(sizeof (GtRange));
  return sa;
}

void gt_sspliced_alignment_delete(GtSSplicedAlignment *sa)
{
  if (!sa) return;
  gt_array_delete(sa->exons);
  gt_free(sa->id);
  gt_free(sa);
}

bool gt_sspliced_alignment_is_forward(const GtSSplicedAlignment *sa)
{
  gt_assert(sa);
  return sa->forward;
}

void gt_sspliced_alignment_add_exon(GtSSplicedAlignment *sa, GtRange exon)
{
  gt_assert(sa);
  gt_array_add(sa->exons, exon);
}

unsigned long gt_sspliced_alignment_num_of_exons(const GtSSplicedAlignment *sa)
{
  gt_assert(sa);
  return gt_array_size(sa->exons);
}

GtRange gt_sspliced_alignment_get_exon(const GtSSplicedAlignment *sa,
                                  unsigned long exon_number)
{
  gt_assert(sa);
  return *(GtRange*) gt_array_get(sa->exons, exon_number);
}

GtRange gt_sspliced_alignment_genomic_range(const GtSSplicedAlignment *sa)
{
  GtRange range;
  gt_assert(sa);
  gt_assert(gt_array_size(sa->exons));
  range.start = ((GtRange*) gt_array_get_first(sa->exons))->start;
  range.end   = ((GtRange*) gt_array_get_last(sa->exons))->end;
  return range;
}

static int range_compare_long_first(GtRange range_a, GtRange range_b)
{
  gt_assert(range_a.start <= range_a.end && range_b.start <= range_b.end);

  if ((range_a.start == range_b.start) && (range_a.end == range_b.end))
    return 0; /* range_a == range_b */

  if ((range_a.start < range_b.start) ||
      ((range_a.start == range_b.start) && (range_a.end > range_b.end)))
    return -1; /* range_a < range_b */

  return 1; /* range_a > range_b */
}

int gt_sspliced_alignment_compare_ptr(const GtSSplicedAlignment **sa_a,
                                   const GtSSplicedAlignment **sa_b)
{
  GtRange range_a, range_b;
  range_a = gt_sspliced_alignment_genomic_range(*sa_a);
  range_b = gt_sspliced_alignment_genomic_range(*sa_b);
  return range_compare_long_first(range_a, range_b);
}
