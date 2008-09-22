/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg

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

#include <assert.h>
#include "extended/mergefeat_stream_sorted.h"
#include "extended/mergefeat_stream_unsorted.h"
#include "extended/node_stream_rep.h"
#include "extended/sort_stream.h"

struct GtMergefeatStreamSorted {
  const GtNodeStream parent_instance;
  GtNodeStream *mergefeat_stream_unsorted,
               *sort_stream;
};

#define gt_mergefeat_stream_sorted_cast(GS)\
        gt_node_stream_cast(gt_mergefeat_stream_sorted_class(), GS)

static int mergefeat_stream_sorted_next(GtNodeStream *gs, GtGenomeNode **gn,
                                        GtError *err)
{
  GtMergefeatStreamSorted *mfs;
  gt_error_check(err);
  mfs = gt_mergefeat_stream_sorted_cast(gs);
  return gt_node_stream_next(mfs->sort_stream, gn, err);
}

static void mergefeat_stream_sorted_free(GtNodeStream *gs)
{
  GtMergefeatStreamSorted *mfs = gt_mergefeat_stream_sorted_cast(gs);
  gt_node_stream_delete(mfs->mergefeat_stream_unsorted);
  gt_node_stream_delete(mfs->sort_stream);
}

const GtNodeStreamClass* gt_mergefeat_stream_sorted_class(void)
{
  static const GtNodeStreamClass *nsc = NULL;
  if (!nsc) {
   nsc = gt_node_stream_class_new(sizeof (GtMergefeatStreamSorted),
                                  mergefeat_stream_sorted_free,
                                  mergefeat_stream_sorted_next);
  }
  return nsc;
}

GtNodeStream* gt_mergefeat_stream_sorted_new(GtNodeStream *in_stream)
{
  GtNodeStream *gs = gt_node_stream_create(gt_mergefeat_stream_sorted_class(),
                                          true);
  GtMergefeatStreamSorted *mfs = gt_mergefeat_stream_sorted_cast(gs);
  assert(in_stream && gt_node_stream_is_sorted(in_stream));
  mfs->mergefeat_stream_unsorted = gt_mergefeat_stream_unsorted_new(in_stream);
  mfs->sort_stream = sort_stream_new(mfs->mergefeat_stream_unsorted);
  return gs;
}
