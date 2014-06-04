/*
  Copyright (c) 2011 Sascha Kastens <mail@skastens.de>
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
#include "extended/array_out_stream.h"
#include "extended/feature_node.h"
#include "extended/node_stream_api.h"

struct GtArrayOutStream {
  const GtNodeStream parent_instance;
  GtNodeStream *in_stream;
  GtArray *nodes;
  bool store_all;
};

#define gt_array_out_stream_cast(GS)\
        gt_node_stream_cast(gt_array_out_stream_class(), GS);

static int gt_array_out_stream_next(GtNodeStream *gs, GtGenomeNode **gn,
                                    GtError *err)
{
  GtArrayOutStream *aos;
  GT_UNUSED GtFeatureNode *node;
  GtGenomeNode  *gn_ref;
  int had_err = 0;

  gt_error_check(err);
  aos = gt_array_out_stream_cast(gs);
  had_err = gt_node_stream_next(aos->in_stream, gn, err);

  if (!had_err && *gn) {
    if (aos->store_all || (node = gt_feature_node_try_cast(*gn))) {
      gn_ref = gt_genome_node_ref(*gn);
      gt_array_add(aos->nodes, gn_ref);
    }
  }

  return had_err;
}

static void gt_array_out_stream_free(GtNodeStream *gs)
{
  GtArrayOutStream *aos = gt_array_out_stream_cast(gs);
  gt_node_stream_delete(aos->in_stream);
}

const GtNodeStreamClass* gt_array_out_stream_class(void)
{
  static const GtNodeStreamClass *gsc = NULL;
  gt_class_alloc_lock_enter();
  if (!gsc) {
    gsc = gt_node_stream_class_new(sizeof (GtArrayOutStream),
                                   gt_array_out_stream_free,
                                   gt_array_out_stream_next);
  }
  gt_class_alloc_lock_leave();
  return gsc;
}

static GtNodeStream* gt_array_out_stream_new_generic(GtNodeStream *in_stream,
                                                     GtArray *nodes,
                                                     bool store_all,
                                                     GT_UNUSED GtError *err)
{
  GtNodeStream *gs;
  GtArrayOutStream *aos;
  gt_assert(in_stream && nodes);
  gs = gt_node_stream_create(gt_array_out_stream_class(), false);
  aos = gt_array_out_stream_cast(gs);
  aos->in_stream = gt_node_stream_ref(in_stream);
  aos->nodes = nodes;
  aos->store_all = store_all;
  return gs;
}

GtNodeStream* gt_array_out_stream_new(GtNodeStream *in_stream,
                                      GtArray *nodes, GtError *err)
{
  return gt_array_out_stream_new_generic(in_stream, nodes, false, err);
}

GtNodeStream* gt_array_out_stream_all_new(GtNodeStream *in_stream,
                                          GtArray *nodes, GtError *err)
{
  return gt_array_out_stream_new_generic(in_stream, nodes, true, err);
}
