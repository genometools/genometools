/*
  Copyright (c) 2010 Gordon Gremme <gremme@zbh.uni-hamburg.de>

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

#include "extended/seqids_to_md5s_stream.h"
#include "extended/seqids_to_md5s_visitor.h"
#include "extended/node_stream_api.h"

struct GtSeqidsToMD5Stream {
  const GtNodeStream parent_instance;
  GtNodeStream *in_stream;
  GtNodeVisitor *seqids_to_md5s_visitor;
};

#define seqids_to_md5s_stream_cast(GS)\
        gt_node_stream_cast(gt_seqids_to_md5s_stream_class(), GS)

static int seqids_to_md5s_stream_next(GtNodeStream *ns, GtGenomeNode **gn,
                                     GtError *err)
{
  GtSeqidsToMD5Stream *seqids_to_md5s_stream;
  int had_err;
  gt_error_check(err);
  seqids_to_md5s_stream = seqids_to_md5s_stream_cast(ns);
  had_err = gt_node_stream_next(seqids_to_md5s_stream->in_stream, gn, err);
  if (!had_err && *gn) {
    had_err = gt_genome_node_accept(*gn,
                                    seqids_to_md5s_stream
                                    ->seqids_to_md5s_visitor, err);
  }
  if (had_err) {
    /* we own the node -> delete it */
    gt_genome_node_delete(*gn);
    *gn = NULL;
  }
  return had_err;
}

static void seqids_to_md5s_stream_free(GtNodeStream *ns)
{
  GtSeqidsToMD5Stream *seqids_to_md5s_stream = seqids_to_md5s_stream_cast(ns);
  gt_node_visitor_delete(seqids_to_md5s_stream->seqids_to_md5s_visitor);
  gt_node_stream_delete(seqids_to_md5s_stream->in_stream);
}

const GtNodeStreamClass* gt_seqids_to_md5s_stream_class(void)
{
  static const GtNodeStreamClass *nsc = NULL;
  if (!nsc) {
    nsc = gt_node_stream_class_new(sizeof (GtSeqidsToMD5Stream),
                                   seqids_to_md5s_stream_free,
                                   seqids_to_md5s_stream_next);
  }
  return nsc;
}

GtNodeStream* gt_seqids_to_md5s_stream_new(GtNodeStream *in_stream,
                                          GtRegionMapping *rm)
{
  GtSeqidsToMD5Stream *seqids_to_md5s_stream;
  GtNodeStream *ns;
  ns = gt_node_stream_create(gt_seqids_to_md5s_stream_class(), true);
  seqids_to_md5s_stream = seqids_to_md5s_stream_cast(ns);
  seqids_to_md5s_stream->in_stream = gt_node_stream_ref(in_stream);
  seqids_to_md5s_stream->seqids_to_md5s_visitor =
    gt_seqids_to_md5s_visitor_new(rm);
  return ns;
}
