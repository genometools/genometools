/*
  Copyright (c) 2006-2011 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2008 Center for Bioinformatics, University of Hamburg

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
#include "extended/eof_node.h"
#include "extended/genome_node.h"
#include "extended/merge_stream.h"
#include "extended/node_stream_api.h"
#include "extended/region_node.h"

struct GtMergeStream {
  const GtNodeStream parent_instance;
  GtArray *node_streams;
  GtGenomeNode **buffer;
};

#define gt_merge_stream_cast(GS)\
        gt_node_stream_cast(gt_merge_stream_class(), GS)

static int merge_stream_next(GtNodeStream *ns, GtGenomeNode **gn, GtError *err)
{
  GtMergeStream *ms;
  GtGenomeNode *min_node = NULL;
  unsigned long i, j, min_i = GT_UNDEF_ULONG;
  unsigned int gt_genome_node_consolidated;
  int had_err = 0;

  gt_error_check(err);

  ms = gt_merge_stream_cast(ns);

  /* fill buffers */
  for (i = 0; !had_err && i < gt_array_size(ms->node_streams); i++) {
    while (!ms->buffer[i]) {
      had_err = gt_node_stream_next(*(GtNodeStream**)
                                    gt_array_get(ms->node_streams, i),
                                    ms->buffer + i, err);
      if (had_err || !ms->buffer[i])
        break;
      /* remove EOF nodes */
      if (gt_eof_node_try_cast(ms->buffer[i])) {
        gt_genome_node_delete(ms->buffer[i]);
        ms->buffer[i] = NULL;
      }
    }
  }

  /* consolidate sequence regions (to avoid duplicates) */
  if (!had_err) {
    for (;;) {
      gt_genome_node_consolidated = 0;
      for (i = 0; i < gt_array_size(ms->node_streams); i++) {
        for (j = i+1; j < gt_array_size(ms->node_streams); j++) {
          gt_assert(i != j);
          if (gt_genome_nodes_are_equal_region_nodes(ms->buffer[i],
                                                     ms->buffer[j])) {
            gt_region_node_consolidate(gt_region_node_cast(ms->buffer[i]),
                                       gt_region_node_cast(ms->buffer[j]));
            gt_genome_node_delete(ms->buffer[j]);
            ms->buffer[j] = NULL;
          }
        }
      }
      if (!gt_genome_node_consolidated)
        break;
    }
  }

  /* find minimal node */
  if (!had_err) {
    for (i = 0; i < gt_array_size(ms->node_streams); i++) {
      if (ms->buffer[i]) {
        if (min_i != GT_UNDEF_ULONG) {
          if (gt_genome_node_compare(ms->buffer + i, ms->buffer + min_i) < 0)
            min_i = i;
        }
        else min_i = i;
      }
    }
    if (min_i != GT_UNDEF_ULONG) {
      min_node = ms->buffer[min_i];
      ms->buffer[min_i] = NULL;
    }
  }

  *gn = min_node;
  return had_err;
}

static void merge_stream_free(GtNodeStream *ns)
{
  GtMergeStream *ms = gt_merge_stream_cast(ns);
  unsigned long i;
  for (i = 0; i < gt_array_size(ms->node_streams); i++)
    gt_node_stream_delete(*(GtNodeStream**) gt_array_get(ms->node_streams, i));
  gt_array_delete(ms->node_streams);
  gt_free(ms->buffer);
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
  unsigned long i;
#ifndef NDEBUG
  gt_assert(gt_array_size(node_streams)); /* at least on input stream given */
  /* each input stream is sorted */
  for (i = 0; i < gt_array_size(node_streams); i++) {
    gt_assert(gt_node_stream_is_sorted(*(GtNodeStream**)
                                       gt_array_get(node_streams, i)));
  }
#endif
  ms->node_streams = gt_array_new(sizeof (GtNodeStream*));
  for (i = 0; i < gt_array_size(node_streams); i++) {
    in_stream = gt_node_stream_ref(*(GtNodeStream**)
                                   gt_array_get(node_streams, i));
    gt_array_add(ms->node_streams, in_stream);
  }
  ms->buffer = gt_calloc(gt_array_size(node_streams), sizeof (GtGenomeNode*));
  return ns;
}
