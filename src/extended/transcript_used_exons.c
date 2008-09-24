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
#include "extended/transcript_used_exons.h"

struct GtTranscriptUsedExons {
  GtDlist *used_exons_all,
        *used_exons_single,
        *used_exons_initial,
        *used_exons_internal,
        *used_exons_terminal;
};

GtTranscriptUsedExons* gt_transcript_used_exons_new(void)
{
  GtTranscriptUsedExons *tue = gt_malloc(sizeof (GtTranscriptUsedExons));
  tue->used_exons_all = gt_dlist_new((GtCompare) gt_range_compare);
  tue->used_exons_single = gt_dlist_new((GtCompare) gt_range_compare);
  tue->used_exons_initial = gt_dlist_new((GtCompare) gt_range_compare);
  tue->used_exons_internal = gt_dlist_new((GtCompare) gt_range_compare);
  tue->used_exons_terminal = gt_dlist_new((GtCompare) gt_range_compare);
  return tue;
}

GtDlist* gt_transcript_used_exons_get_all(GtTranscriptUsedExons *tue)
{
  gt_assert(tue);
  return tue->used_exons_all;
}

GtDlist* gt_transcript_used_exons_get_single(GtTranscriptUsedExons *tue)
{
  gt_assert(tue);
  return tue->used_exons_single;
}

GtDlist* gt_transcript_used_exons_get_initial(GtTranscriptUsedExons *tue)
{
  gt_assert(tue);
  return tue->used_exons_initial;
}

GtDlist* gt_transcript_used_exons_get_internal(GtTranscriptUsedExons *tue)
{
  gt_assert(tue);
  return tue->used_exons_internal;
}

GtDlist* gt_transcript_used_exons_get_terminal(GtTranscriptUsedExons *tue)
{
  gt_assert(tue);
  return tue->used_exons_terminal;
}

static void used_gt_dlist_delete(GtDlist *used_list)
{
  GtDlistelem *dlistelem;
  for (dlistelem = gt_dlist_first(used_list); dlistelem != NULL;
       dlistelem = gt_dlistelem_next(dlistelem)) {
    gt_free(gt_dlistelem_get_data(dlistelem));
  }
  gt_dlist_delete(used_list);
}

void gt_transcript_used_exons_delete(GtTranscriptUsedExons *tue)
{
  if (!tue) return;
  used_gt_dlist_delete(tue->used_exons_all);
  used_gt_dlist_delete(tue->used_exons_single);
  used_gt_dlist_delete(tue->used_exons_initial);
  used_gt_dlist_delete(tue->used_exons_internal);
  used_gt_dlist_delete(tue->used_exons_terminal);
  gt_free(tue);
}
