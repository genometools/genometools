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

#include "extended/gtf_out_stream.h"
#include "extended/gtf_visitor.h"
#include "extended/node_stream_api.h"

struct GtGTFOutStream {
  const GtNodeStream parent_instance;
  GtNodeStream *in_stream;
  GtNodeVisitor *gtf_visitor;
};

#define gtf_out_stream_cast(GS)\
        gt_node_stream_cast(gt_gtf_out_stream_class(), GS);

static int gtf_out_stream_next(GtNodeStream *gs, GtGenomeNode **gn,
                               GtError *err)
{
  GtGTFOutStream *gtf_out_stream;
  int had_err;
  gt_error_check(err);
  gtf_out_stream = gtf_out_stream_cast(gs);
  had_err = gt_node_stream_next(gtf_out_stream->in_stream, gn, err);
  if (!had_err && *gn)
    had_err = gt_genome_node_accept(*gn, gtf_out_stream->gtf_visitor, err);
  return had_err;
}

static void gtf_out_stream_free(GtNodeStream *gs)
{
  GtGTFOutStream *gtf_out_stream = gtf_out_stream_cast(gs);
  gt_node_visitor_delete(gtf_out_stream->gtf_visitor);
  gt_node_stream_delete(gtf_out_stream->in_stream);
}

const GtNodeStreamClass* gt_gtf_out_stream_class(void)
{
  static const GtNodeStreamClass *nsc = NULL;
  if (!nsc) {
   nsc = gt_node_stream_class_new(sizeof (GtGTFOutStream),
                                  gtf_out_stream_free,
                                  gtf_out_stream_next);
  }
  return nsc;
}

GtNodeStream* gt_gtf_out_stream_new(GtNodeStream *in_stream, GtGenFile *outfp)
{
  GtNodeStream *gs = gt_node_stream_create(gt_gtf_out_stream_class(),
                                           gt_node_stream_is_sorted(in_stream));
  GtGTFOutStream *gtf_out_stream = gtf_out_stream_cast(gs);
  gtf_out_stream->in_stream = gt_node_stream_ref(in_stream);
  gtf_out_stream->gtf_visitor = gt_gtf_visitor_new(outfp);
  return gs;
}
