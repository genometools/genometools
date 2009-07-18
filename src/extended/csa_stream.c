/*
  Copyright (c) 2006-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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
#include "extended/csa_stream.h"
#include "extended/csa_visitor.h"
#include "extended/consensus_sa.h"
#include "extended/node_stream_api.h"

struct GtCSAStream {
  const GtNodeStream parent_instance;
  GtNodeStream *in_stream;
  GtNodeVisitor *csa_visitor; /* the actual work is done in the visitor */
};

#define csa_stream_cast(GS)\
        gt_node_stream_cast(gt_csa_stream_class(), GS)

static int csa_stream_next(GtNodeStream *gs, GtGenomeNode **gn, GtError *err)
{
  GtCSAStream *cs;
  int had_err;
  gt_error_check(err);
  cs = csa_stream_cast(gs);

  /* we have still nodes in the buffer */
  if (gt_csa_visitor_node_buffer_size(cs->csa_visitor)) {
    *gn = gt_csa_visitor_get_node(cs->csa_visitor); /* return one of them */
    return 0;
  }

  /* no nodes in the buffer -> get new nodes */
  while (!(had_err = gt_node_stream_next(cs->in_stream, gn, err)) && *gn) {
    gt_assert(*gn && !had_err);
    had_err = gt_genome_node_accept(*gn, cs->csa_visitor, err);
    if (had_err)
      break;
    if (gt_csa_visitor_node_buffer_size(cs->csa_visitor)) {
      *gn = gt_csa_visitor_get_node(cs->csa_visitor);
      return 0;
    }
  }

  /* either we have an error or no new node */
  gt_assert(had_err || !*gn);

  /* if we have no error, process the last cluster */
  if (!had_err) {
    gt_csa_visitor_process_cluster(cs->csa_visitor, true);
    if (gt_csa_visitor_node_buffer_size(cs->csa_visitor)) {
      *gn = gt_csa_visitor_get_node(cs->csa_visitor);
      return 0;
    }
  }
  return had_err;
}

static void csa_stream_free(GtNodeStream *gs)
{
  GtCSAStream *cs = csa_stream_cast(gs);
  gt_node_visitor_delete(cs->csa_visitor);
  gt_node_stream_delete(cs->in_stream);
}

const GtNodeStreamClass* gt_csa_stream_class(void)
{
  static const GtNodeStreamClass *nsc= NULL;
  if (!nsc) {
    nsc = gt_node_stream_class_new(sizeof (GtCSAStream),
                                   csa_stream_free,
                                   csa_stream_next);
  }
  return nsc;
}

GtNodeStream* gt_csa_stream_new(GtNodeStream *in_stream,
                                unsigned long join_length)
{
  GtNodeStream *gs = gt_node_stream_create(gt_csa_stream_class(),
                                          gt_node_stream_is_sorted(in_stream));
  GtCSAStream *cs = csa_stream_cast(gs);
  cs->in_stream = gt_node_stream_ref(in_stream);
  cs->csa_visitor = gt_csa_visitor_new(join_length);
  return gs;
}
