/*
  Copyright (c) 2011 Sascha Kastens <sascha.kastens@studium.uni-hamburg.de>
  Copyright (c) 2011 Center for Bioinformatics, University of Hamburg

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
#include "core/class_alloc_lock.h"
#include "core/error_api.h"
#include "core/unused_api.h"
#include "extended/array_in_stream.h"

struct GtArrayInStream {
  const GtNodeStream parent_instance;
  GtArray *nodes;
  unsigned long next_index,
                *progress;
};

#define gt_array_in_stream_cast(GS)\
        gt_node_stream_cast(gt_array_in_stream_class(), GS);

static int gt_array_in_stream_next(GtNodeStream *gs, GtGenomeNode **gn,
                                   GT_UNUSED GtError *err)
{
  GtArrayInStream *ais;
  int had_err = 0;

  gt_error_check(err);
  ais = gt_array_in_stream_cast(gs);
  if (ais->next_index >= gt_array_size(ais->nodes))
    *gn = NULL;
  else {
    *gn = *(GtGenomeNode**) gt_array_get(ais->nodes, ais->next_index);
    ais->next_index++;
    if (ais->progress)
      *ais->progress = *ais->progress + 1;
  }

  return had_err;
}

const GtNodeStreamClass* gt_array_in_stream_class(void)
{
  static const GtNodeStreamClass *gsc = NULL;
  gt_class_alloc_lock_enter();
  if (!gsc) {
    gsc = gt_node_stream_class_new(sizeof (GtArrayInStream),
                                   NULL,
                                   gt_array_in_stream_next);
  }
  gt_class_alloc_lock_leave();
  return gsc;
}

GtNodeStream* gt_array_in_stream_new(GtArray *nodes,
                                     unsigned long *progress,
                                     GT_UNUSED GtError *err)
{
  GtNodeStream *gs;
  GtArrayInStream *ais;
  gt_assert(nodes);
  gs = gt_node_stream_create(gt_array_in_stream_class(), false);
  ais = gt_array_in_stream_cast(gs);
  ais->nodes = nodes;
  ais->progress = progress;
  ais->next_index = 0;

  return gs;
}
