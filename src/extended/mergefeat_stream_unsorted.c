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
#include "extended/mergefeat_stream_unsorted.h"
#include "extended/mergefeat_visitor.h"
#include "extended/node_stream_rep.h"

struct MergefeatStreamUnsorted {
  const GtNodeStream parent_instance;
  GtNodeStream *in_stream;
  GenomeVisitor *mergefeat_visitor;
};

#define mergefeat_stream_unsorted_cast(GS)\
        gt_node_stream_cast(mergefeat_stream_unsorted_class(), GS)

static int mergefeat_stream_unsorted_next_tree(GtNodeStream *gs,
                                               GtGenomeNode **gn,
                                               GtError *err)
{
  MergefeatStreamUnsorted *mfs;
  int had_err;
  gt_error_check(err);
  mfs = mergefeat_stream_unsorted_cast(gs);
  had_err = gt_node_stream_next(mfs->in_stream, gn, err);
  if (!had_err && *gn)
    had_err = gt_genome_node_accept(*gn, mfs->mergefeat_visitor, err);
  return had_err;
}

static void mergefeat_stream_unsorted_free(GtNodeStream *gs)
{
  MergefeatStreamUnsorted *mfs = mergefeat_stream_unsorted_cast(gs);
  genome_visitor_delete(mfs->mergefeat_visitor);
}

const GtNodeStreamClass* mergefeat_stream_unsorted_class(void)
{
  static const GtNodeStreamClass gsc = { sizeof (MergefeatStreamUnsorted),
                                         mergefeat_stream_unsorted_next_tree,
                                         mergefeat_stream_unsorted_free };
  return &gsc;
}

GtNodeStream* mergefeat_stream_unsorted_new(GtNodeStream *in_stream)
{
  GtNodeStream *gs = gt_node_stream_create(mergefeat_stream_unsorted_class(),
                                          false);
  MergefeatStreamUnsorted *mfs = mergefeat_stream_unsorted_cast(gs);
  mfs->in_stream = in_stream;
  mfs->mergefeat_visitor = mergefeat_visitor_new();
  return gs;
}
