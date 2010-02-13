/*
  Copyright (c) 2006-2010 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#include "core/array.h"
#include "core/assert_api.h"
#include "extended/genome_node.h"
#include "extended/node_stream_api.h"
#include "extended/sort_stream.h"

struct GtSortStream
{
  const GtNodeStream parent_instance;
  GtNodeStream *in_stream;
  unsigned long idx;
  GtArray *trees;
  bool sorted;
};

#define gt_sort_stream_cast(GS)\
        gt_node_stream_cast(gt_sort_stream_class(), GS);

static int gt_sort_stream_next(GtNodeStream *gs, GtGenomeNode **gn,
                               GtError *err)
{
  GtSortStream *sort_stream;
  GtGenomeNode *node;
  int had_err = 0;
  gt_error_check(err);
  sort_stream = gt_sort_stream_cast(gs);

  if (!sort_stream->sorted) {
    while (!(had_err = gt_node_stream_next(sort_stream->in_stream, &node,
                                           err)) && node) {
      gt_array_add(sort_stream->trees, node);
    }
    if (!had_err) {
      gt_genome_nodes_sort_stable(sort_stream->trees);
      sort_stream->sorted = true;
    }
  }

  if (!had_err) {
    gt_assert(sort_stream->sorted);
    if (sort_stream->idx < gt_array_size(sort_stream->trees)) {
      *gn = *(GtGenomeNode**) gt_array_get(sort_stream->trees,
                                           sort_stream->idx);
      sort_stream->idx++;
      return 0;
    }
  }

  if (!had_err) {
    gt_array_reset(sort_stream->trees);
    *gn = NULL;
  }

  return had_err;
}

static void gt_sort_stream_free(GtNodeStream *gs)
{
  unsigned long i;
  GtSortStream *sort_stream = gt_sort_stream_cast(gs);
  for (i = sort_stream->idx; i < gt_array_size(sort_stream->trees); i++) {
    gt_genome_node_delete(*(GtGenomeNode**)
                              gt_array_get(sort_stream->trees, i));
  }
  gt_array_delete(sort_stream->trees);
  gt_node_stream_delete(sort_stream->in_stream);
}

const GtNodeStreamClass* gt_sort_stream_class(void)
{
  static const GtNodeStreamClass *nsc = NULL;
  if (!nsc) {
    nsc = gt_node_stream_class_new(sizeof (GtSortStream),
                                   gt_sort_stream_free,
                                   gt_sort_stream_next);
  }
  return nsc;
}

GtNodeStream* gt_sort_stream_new(GtNodeStream *in_stream)
{
  GtNodeStream *gs = gt_node_stream_create(gt_sort_stream_class(), true);
  GtSortStream *sort_stream = gt_sort_stream_cast(gs);
  gt_assert(in_stream);
  sort_stream->in_stream = gt_node_stream_ref(in_stream);
  sort_stream->sorted = false;
  sort_stream->idx = 0;
  sort_stream->trees = gt_array_new(sizeof (GtGenomeNode*));
  return gs;
}
