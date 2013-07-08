/*
  Copyright (c) 2010-2011 Gordon Gremme <gordon@gremme.org>

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

#include "extended/genome_node.h"
#include "extended/node_stream_api.h"
#include "gth/bssm_train_stream.h"
#include "gth/bssm_train_visitor.h"

struct GthBSSMTrainStream {
  const GtNodeStream parent_instance;
  GtNodeStream *in_stream;
  GtNodeVisitor *bssm_train_visitor;
};

#define bssm_train_stream_cast(GS)\
        gt_node_stream_cast(gth_bssm_train_stream_class(), GS)

static int bssm_train_stream_next(GtNodeStream *ns, GtGenomeNode **gn,
                                  GtError *err)
{
  GthBSSMTrainStream *bssm_train_stream;
  int had_err;
  gt_error_check(err);
  bssm_train_stream = bssm_train_stream_cast(ns);
  had_err = gt_node_stream_next(bssm_train_stream->in_stream, gn, err);
  if (!had_err && *gn) {
    had_err = gt_genome_node_accept(*gn, bssm_train_stream->bssm_train_visitor,
                                    err);
  }
  if (had_err) {
    /* we own the node -> delete it */
    gt_genome_node_delete(*gn);
    *gn = NULL;
  }
  return had_err;
}

static void bssm_train_stream_free(GtNodeStream *ns)
{
  GthBSSMTrainStream *bssm_train_stream = bssm_train_stream_cast(ns);
  gt_node_visitor_delete(bssm_train_stream->bssm_train_visitor);
  gt_node_stream_delete(bssm_train_stream->in_stream);
}

const GtNodeStreamClass* gth_bssm_train_stream_class(void)
{
  static const GtNodeStreamClass *nsc = NULL;
  if (!nsc) {
    nsc = gt_node_stream_class_new(sizeof (GthBSSMTrainStream),
                                   bssm_train_stream_free,
                                   bssm_train_stream_next);
  }
  return nsc;
}

GtNodeStream* gth_bssm_train_stream_new(GtNodeStream *in_stream,
                                        GtRegionMapping *region_mapping,
                                        GthBSSMSeqProcessor *bsp,
                                        const char *filter_type,
                                        const char *extract_type,
                                        unsigned int good_exon_count,
                                        double cutoff)
{
  GtNodeStream *ns;
  GthBSSMTrainStream *bts;
  gt_assert(in_stream && region_mapping && filter_type && extract_type);
  ns = gt_node_stream_create(gth_bssm_train_stream_class(), true);
  bts = bssm_train_stream_cast(ns);
  bts->in_stream = gt_node_stream_ref(in_stream);
  bts->bssm_train_visitor =
    gth_bssm_train_visitor_new(region_mapping, bsp, filter_type, extract_type,
                               good_exon_count, cutoff);
  return ns;
}
