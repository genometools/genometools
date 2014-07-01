/*
  Copyright (c) 2006-2011 Gordon Gremme <gordon@gremme.org>
  Copyright (c) 2014      Sascha Steinbiss <ss34@sanger.ac.uk>
  Copyright (c) 2006-2008 Center for Bioinformatics, University of Hamburg
  Copyright (c) 2014      Genome Research Ltd

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

#include "core/assert_api.h"
#include "core/class_alloc_lock.h"
#include "core/ma.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "extended/eof_node_api.h"
#include "extended/genome_node.h"
#include "extended/merge_stream.h"
#include "extended/node_stream_api.h"
#include "extended/priority_queue.h"
#include "extended/region_node.h"

typedef struct {
  GtGenomeNode *gn;
  GtUword input_index;
} GtMergeStreamItem;

struct GtMergeStream {
  const GtNodeStream parent_instance;
  GtArray *node_streams;
  GtGenomeNode *first_node, *second_node;
  GtMergeStreamItem *items;
  GtPriorityQueue *pq;
  bool filled;
};

#define gt_merge_stream_cast(GS)\
        gt_node_stream_cast(gt_merge_stream_class(), GS)

static int gt_merge_stream_item_compare(const void *a, const void *b)
{
  GtMergeStreamItem *item1, *item2;
  gt_assert(a && b);
  item1 = (GtMergeStreamItem*) a;
  item2 = (GtMergeStreamItem*) b;
  gt_assert(item1->gn && item2->gn);
  return gt_genome_node_compare(&item1->gn, &item2->gn);
}

static int merge_stream_next_in_order(GtNodeStream *ns, GtGenomeNode **gn,
                                      GtError *err)
{
  GtMergeStream *ms;
  GtGenomeNode *min_node = NULL;
  GtUword i;
  int had_err = 0;
  gt_error_check(err);
  ms = gt_merge_stream_cast(ns);

  /* initially fill the queue with first nodes from input streams */
  if (!ms->filled) {
    for (i = 0; !had_err && i < gt_array_size(ms->node_streams); i++) {
      GtGenomeNode *firstnode = NULL;
      ms->items[i].input_index = i;
      had_err = gt_node_stream_next(*(GtNodeStream**)
                                      gt_array_get(ms->node_streams, i),
                                    &firstnode, err);
      /* only add start of non-empty input streams */
      if (!had_err && firstnode) {
        if (!gt_eof_node_try_cast(firstnode)) {
          ms->items[i].gn = firstnode;
          ms->items[i].input_index = i;
          gt_priority_queue_add(ms->pq, ms->items+i);
        } else {
          gt_genome_node_delete(firstnode);
        }
      }
    }
    ms->filled = true;
  }

  /* stream the queue contents */
  if (!gt_priority_queue_is_empty(ms->pq)) {
    GtMergeStreamItem *min_item = gt_priority_queue_extract_min(ms->pq);
    GtGenomeNode *nextnode = NULL;
    gt_assert(min_item && min_item->gn);
    min_node = min_item->gn;
    /* get next element from the last stream queried */
    had_err = gt_node_stream_next(*(GtNodeStream**)
                                    gt_array_get(ms->node_streams,
                                                 min_item->input_index),
                                  &nextnode, err);
    /* add node to queue if still non-EOF nodes left in that stream */
    if (!had_err && nextnode) {
      if (!gt_eof_node_try_cast(nextnode)) {
        min_item->gn = nextnode;
        gt_priority_queue_add(ms->pq, min_item);
      } else {
        gt_genome_node_delete(nextnode);
      }
    }
  }

  *gn = min_node;
  return had_err;
}

static bool merge_stream_consolidate(GtGenomeNode **first_node,
                                     GtGenomeNode **second_node)
{
  gt_assert(*first_node && *second_node);
  if (gt_genome_nodes_are_equal_region_nodes(*first_node, *second_node)) {
    gt_region_node_consolidate((GtRegionNode*) *first_node,
                               (GtRegionNode*) *second_node);
    /* keep first node */
    gt_genome_node_delete(*second_node);
    *second_node = NULL;
    return true;
  }
  return false;
}

/* make sure that region nodes for equal regions (which follow each other in a
   sorted stream) are consolidated into a single one to avoid duplicate region
   definitions */
static int merge_stream_next(GtNodeStream *ns, GtGenomeNode **gn, GtError *err)
{
  GtMergeStream *ms;
  int had_err = 0;
  gt_error_check(err);
  ms = gt_merge_stream_cast(ns);

  gt_assert(!ms->second_node); /* the second buffer is always empty when this
                                  function is called */
  if (!ms->first_node) {
    /* both buffers are empty */
    had_err = merge_stream_next_in_order(ns, &ms->first_node, err);
    if (had_err)
      return had_err;
    if (!ms->first_node) {
      *gn = NULL;
      return 0;
    }
  }

  for (;;) {
    gt_assert(ms->first_node && !ms->second_node);
    had_err = merge_stream_next_in_order(ns, &ms->second_node, err);
    if (!had_err && ms->second_node) {
      if (!merge_stream_consolidate(&ms->first_node, &ms->second_node))
        break;
    }
    else
      break;
  }

  if (!had_err) {
    gt_assert(ms->first_node);
    *gn = ms->first_node;
    ms->first_node = ms->second_node;
    ms->second_node = NULL;
  }

  return had_err;
}

static void merge_stream_free(GtNodeStream *ns)
{
  GtMergeStream *ms = gt_merge_stream_cast(ns);
  GtUword i;
  for (i = 0; i < gt_array_size(ms->node_streams); i++)
    gt_node_stream_delete(*(GtNodeStream**) gt_array_get(ms->node_streams, i));
  gt_array_delete(ms->node_streams);
  gt_free(ms->items);
  gt_priority_queue_delete(ms->pq);
}

const GtNodeStreamClass* gt_merge_stream_class(void)
{
  static const GtNodeStreamClass *nsc = NULL;
  gt_class_alloc_lock_enter();
  if (!nsc) {
    nsc = gt_node_stream_class_new(sizeof (GtMergeStream),
                                   merge_stream_free,
                                   merge_stream_next);
  }
  gt_class_alloc_lock_leave();
  return nsc;
}

GtNodeStream* gt_merge_stream_new(const GtArray *node_streams)
{
  GtNodeStream *in_stream,
               *ns = gt_node_stream_create(gt_merge_stream_class(), true);
  GtMergeStream *ms = gt_merge_stream_cast(ns);
  GtUword i;
#ifndef NDEBUG
  gt_assert(gt_array_size(node_streams)); /* at least one input stream given */
  /* each input stream is sorted */
  for (i = 0; i < gt_array_size(node_streams); i++) {
    gt_assert(gt_node_stream_is_sorted(*(GtNodeStream**)
                                       gt_array_get(node_streams, i)));
  }
#endif
  ms->items = gt_calloc(gt_array_size(node_streams),
                        sizeof (GtMergeStreamItem));
  ms->node_streams = gt_array_new(sizeof (GtNodeStream*));
  for (i = 0; i < gt_array_size(node_streams); i++) {
    in_stream = gt_node_stream_ref(*(GtNodeStream**)
                                   gt_array_get(node_streams, i));
    gt_array_add(ms->node_streams, in_stream);
  }
  ms->pq = gt_priority_queue_new(gt_merge_stream_item_compare,
                                 gt_array_size(node_streams));
  ms->filled = false;
  ms->first_node = ms->second_node = NULL;
  return ns;
}
