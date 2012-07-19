/*
  Copyright (c) 2012 Giorgio Gonnella <gonnella@zbh.uni-hamburg.de>
  Copyright (c) 2012 Center for Bioinformatics, University of Hamburg

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
#include "core/log_api.h"
#include "core/unused_api.h"
#include "core/undef_api.h"
#include "extended/aligned_segments_pile.h"

struct GtAlignedSegmentsPile
{
  GtSamfileIterator *sfi;
  GtSamfileEncseqMapping *sem;
  GtDlist *set;
  GtAlignedSegment *next_as;
  bool all_consumed;
  GtAlignedSegmentsPileProcessFunc process_complete;
  void *process_complete_data;
  GtAlignedSegmentsPileProcessFunc process_unmapped;
  void *process_unmapped_data;
  GtAlignedSegmentsPileProcessFunc process_skipped;
  void *process_skipped_data;
  unsigned long position;
  bool enable_edit_tracking;
  bool delete_processed_segments;
};

int gt_aligned_segments_pile_compare_as(const void *a, const void *b)
{
  unsigned long r_a, r_b;
  r_a = gt_aligned_segment_refregion_endpos((GtAlignedSegment*)a);
  r_b = gt_aligned_segment_refregion_endpos((GtAlignedSegment*)b);
  return (int)((r_a > r_b) - (r_a < r_b));
}

GtAlignedSegmentsPile *gt_aligned_segments_pile_new(GtSamfileIterator *sfi,
    GtSamfileEncseqMapping *sem)
{
  GtAlignedSegmentsPile *asp;
  gt_assert(sfi != NULL);
  gt_assert(sem != NULL);
  asp = gt_malloc(sizeof (GtAlignedSegmentsPile));
  asp->sfi = sfi;
  asp->sem = sem;
  asp->set = gt_dlist_new(gt_aligned_segments_pile_compare_as);
  asp->all_consumed = false;
  asp->next_as = NULL;
  asp->process_complete = NULL;
  asp->process_complete_data = NULL;
  asp->process_unmapped = NULL;
  asp->process_unmapped_data = NULL;
  asp->process_skipped = NULL;
  asp->process_skipped_data = NULL;
  asp->position = GT_UNDEF_ULONG;
  asp->enable_edit_tracking = false;
  asp->delete_processed_segments = true;
  return asp;
}

void gt_aligned_segments_pile_register_process_complete(
    GtAlignedSegmentsPile *asp,
    GtAlignedSegmentsPileProcessFunc process_complete, void *pc_data)
{
  gt_assert(asp != NULL);
  asp->process_complete = process_complete;
  asp->process_complete_data = pc_data;
}

void gt_aligned_segments_pile_register_process_skipped(
    GtAlignedSegmentsPile *asp,
    GtAlignedSegmentsPileProcessFunc process_skipped, void *ps_data)
{
  gt_assert(asp != NULL);
  asp->process_skipped = process_skipped;
  asp->process_skipped_data = ps_data;
}

void gt_aligned_segments_pile_register_process_unmapped(
    GtAlignedSegmentsPile *asp,
    GtAlignedSegmentsPileProcessFunc process_unmapped, void *pu_data)
{
  gt_assert(asp != NULL);
  asp->process_unmapped = process_unmapped;
  asp->process_unmapped_data = pu_data;
}

static int gt_aligned_segments_pile_fetch_sa(GtAlignedSegmentsPile *asp)
{
  int retvalue;
  GtSamAlignment *sa;
  gt_assert(asp != NULL);
  gt_assert(asp->next_as == NULL);
  while (true)
  {
    retvalue = gt_samfile_iterator_next(asp->sfi, &sa);
    if (retvalue == -1 || !gt_sam_alignment_is_unmapped(sa))
      break;
    if (!gt_sam_alignment_is_secondary(sa))
    {
      if (asp->process_unmapped != NULL)
      {
        GtAlignedSegment *as = gt_aligned_segment_new_from_sa(sa, asp->sem);
        if (asp->enable_edit_tracking)
          gt_aligned_segment_enable_edit_tracking(as);
        asp->process_unmapped(as, asp->process_unmapped_data);
        if (asp->delete_processed_segments)
          gt_aligned_segment_delete(as);
      }
    }
  }
  if (retvalue != -1)
  {
    asp->next_as = gt_aligned_segment_new_from_sa(sa, asp->sem);
    if (asp->enable_edit_tracking)
      gt_aligned_segment_enable_edit_tracking(asp->next_as);
  }
  return retvalue;
}

static void gt_aligned_segments_pile_delete_finishing_before(
    GtAlignedSegmentsPile *asp, unsigned long position)
{
  GtDlistelem *dlistelem, *dlistelem_to_remove = NULL;
  for (dlistelem = gt_dlist_first(asp->set); dlistelem != NULL;
       dlistelem = gt_dlistelem_next(dlistelem))
  {
    GtAlignedSegment *as;
    if (dlistelem_to_remove != NULL)
    {
      gt_dlist_remove(asp->set, dlistelem_to_remove);
      dlistelem_to_remove = NULL;
    }
    as = gt_dlistelem_get_data(dlistelem);
    if (gt_aligned_segment_refregion_endpos(as) < position)
    {
      if (asp->process_complete != NULL)
      {
        asp->process_complete(as, asp->process_complete_data);
      }
      if (asp->delete_processed_segments)
        gt_aligned_segment_delete(as);
      /* delay elem removal to avoid disrupting list traversal */
      dlistelem_to_remove = dlistelem;
    }
    else
    {
      break;
    }
  }
  if (dlistelem_to_remove != NULL)
  {
    gt_dlist_remove(asp->set, dlistelem_to_remove);
    dlistelem_to_remove = NULL;
  }
}

void gt_aligned_segments_pile_disable_segment_deletion(
    GtAlignedSegmentsPile *asp)
{
  gt_assert(asp != NULL);
  gt_assert(asp->process_skipped != NULL);
  gt_assert(asp->process_complete != NULL);
  gt_assert(asp->process_unmapped != NULL);
  asp->delete_processed_segments = false;
}

void gt_aligned_segments_pile_enable_edit_tracking(GtAlignedSegmentsPile *asp)
{
  gt_assert(asp != NULL);
  asp->enable_edit_tracking = true;
}

unsigned long gt_aligned_segments_pile_size(GtAlignedSegmentsPile *asp)
{
  gt_assert(asp != NULL);
  return gt_dlist_size(asp->set);
}

void gt_aligned_segments_pile_move_over_position(
    GtAlignedSegmentsPile *asp, unsigned long position)
{
  int retvalue;
  if (asp->position != GT_UNDEF_ULONG)
  {
    gt_assert(position > asp->position);
    gt_aligned_segments_pile_delete_finishing_before(asp, position);
  }
  while (true)
  {
    if (asp->next_as == NULL && !asp->all_consumed)
    {
      retvalue = gt_aligned_segments_pile_fetch_sa(asp);
      if (retvalue == -1)
      {
        asp->all_consumed = true;
        gt_assert(asp->next_as == NULL);
      }
    }
    if (asp->next_as != NULL)
    {
      if (gt_aligned_segment_refregion_endpos(asp->next_as) < position)
      {
        if (asp->process_skipped != NULL)
        {
          asp->process_skipped(asp->next_as, asp->process_unmapped_data);
        }
        if (asp->delete_processed_segments)
          gt_aligned_segment_delete(asp->next_as);
        asp->next_as = NULL;
      }
      else
      {
        if (gt_aligned_segment_refregion_startpos(asp->next_as) <= position)
        {
          gt_dlist_add(asp->set, asp->next_as);
          asp->next_as = NULL;
        }
        else
        {
          break;
        }
      }
    }
    else
    {
      break;
    }
  }
  asp->position = position;
}

const GtDlist *gt_aligned_segments_pile_get(const GtAlignedSegmentsPile *asp)
{
  gt_assert(asp != NULL);
  return asp->set;
}

void gt_aligned_segments_pile_flush(GtAlignedSegmentsPile *asp,
    bool skip_remaining)
{
  gt_assert(asp != NULL);
  gt_aligned_segments_pile_delete_finishing_before(asp, ULONG_MAX);
  if (asp->next_as != NULL)
  {
    if (skip_remaining && asp->process_skipped != NULL)
    {
      asp->process_skipped(asp->next_as, asp->process_skipped_data);
    }
    if (asp->delete_processed_segments)
      gt_aligned_segment_delete(asp->next_as);
    asp->next_as = NULL;
  }
  if (skip_remaining && asp->process_skipped != NULL)
  {
    int retvalue = 0;
    while (retvalue != -1)
    {
      retvalue = gt_aligned_segments_pile_fetch_sa(asp);
      if (asp->next_as != NULL)
      {
        asp->process_skipped(asp->next_as, asp->process_skipped_data);
      }
      if (asp->delete_processed_segments)
        gt_aligned_segment_delete(asp->next_as);
      asp->next_as = NULL;
    }
  }
  gt_assert(asp->next_as == NULL);
}

void gt_aligned_segments_pile_delete(GtAlignedSegmentsPile *asp)
{
  if (asp != NULL)
  {
    gt_aligned_segments_pile_flush(asp, true);
    gt_dlist_delete(asp->set);
    gt_free(asp);
  }
}
