/*
  Copyright (c) 2007-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007-2008 Center for Bioinformatics, University of Hamburg

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
#include "extended/node_stream_rep.h"
#include "extended/add_introns_stream.h"
#include "extended/add_introns_visitor.h"

struct GtAddIntronsStream{
  const GtNodeStream parent_instance;
  GtNodeStream *in_stream;
  GtNodeVisitor *add_introns_visitor;
};

#define gt_add_introns_stream_cast(GS)\
        gt_node_stream_cast(gt_add_introns_stream_class(), GS)

static int gt_add_introns_stream_next_tree(GtNodeStream *gs, GtGenomeNode **gn,
                                        GtError *err)
{
  GtAddIntronsStream *ais;
  int had_err;
  gt_error_check(err);
  ais = gt_add_introns_stream_cast(gs);
  had_err = gt_node_stream_next(ais->in_stream, gn, err);
  if (!had_err && *gn)
    had_err = gt_genome_node_accept(*gn, ais->add_introns_visitor, err);
  return had_err;
}

static void gt_add_introns_stream_free(GtNodeStream *gs)
{
  GtAddIntronsStream *ais = gt_add_introns_stream_cast(gs);
  gt_node_visitor_delete(ais->add_introns_visitor);
  gt_node_stream_delete(ais->in_stream);
}

const GtNodeStreamClass* gt_add_introns_stream_class(void)
{
  static const GtNodeStreamClass gsc = { sizeof (GtAddIntronsStream),
                                         gt_add_introns_stream_next_tree,
                                         gt_add_introns_stream_free };
  return &gsc;
}

GtNodeStream* gt_add_introns_stream_new(GtNodeStream *in_stream)
{
  GtNodeStream *gs = gt_node_stream_create(gt_add_introns_stream_class(), true);
  GtAddIntronsStream *ais = gt_add_introns_stream_cast(gs);
  assert(in_stream);
  ais->in_stream = gt_node_stream_ref(in_stream);
  ais->add_introns_visitor = gt_add_introns_visitor_new();
  return gs;
}
