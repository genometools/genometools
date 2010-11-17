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
#include "extended/feature_node.h"
#include "extended/filter_stream.h"
#include "extended/filter_visitor.h"
#include "extended/node_stream_api.h"

struct GtFilterStream
{
  const GtNodeStream parent_instance;
  GtNodeStream *in_stream;
  GtNodeVisitor *filter_visitor; /* the actual work is done in the visitor */
};

#define gt_filter_stream_cast(GS)\
        gt_node_stream_cast(gt_filter_stream_class(), GS);

static int filter_stream_next(GtNodeStream *ns, GtGenomeNode **gn, GtError *err)
{
  GtFilterStream *fs;
  int had_err;
  gt_error_check(err);
  fs = gt_filter_stream_cast(ns);

  /* we still have nodes in the buffer */
  if (gt_filter_visitor_node_buffer_size(fs->filter_visitor)) {
    /* return one of them */
    *gn = gt_filter_visitor_get_node(fs->filter_visitor);
    return 0;
  }

  /* no nodes in the buffer -> get new nodes */
  while (!(had_err = gt_node_stream_next(fs->in_stream, gn, err)) && *gn) {
    gt_assert(*gn && !had_err);
    had_err = gt_genome_node_accept(*gn, fs->filter_visitor, err);
    if (had_err) {
      /* we own the node -> delete it */
      gt_genome_node_delete(*gn);
      *gn = NULL;
      break;
    }
    if (gt_filter_visitor_node_buffer_size(fs->filter_visitor)) {
      *gn = gt_filter_visitor_get_node(fs->filter_visitor);
      return 0;
    }
  }

  /* either we have an error or no new node */
  gt_assert(had_err || !*gn);
  return had_err;
}

static void filter_stream_free(GtNodeStream *ns)
{
  GtFilterStream *fs = gt_filter_stream_cast(ns);
  gt_node_visitor_delete(fs->filter_visitor);
  gt_node_stream_delete(fs->in_stream);
}

const GtNodeStreamClass* gt_filter_stream_class(void)
{
  static const GtNodeStreamClass *nsc = NULL;
  if (!nsc) {
    nsc = gt_node_stream_class_new(sizeof (GtFilterStream),
                                   filter_stream_free,
                                   filter_stream_next);
  }
  return nsc;
}

GtNodeStream* gt_filter_stream_new(GtNodeStream *in_stream,
                                GtStr *seqid, GtStr *typefilter,
                                GtRange contain_range, GtRange overlap_range,
                                GtStrand strand, GtStrand targetstrand,
                                bool has_CDS, unsigned long max_gene_length,
                                unsigned long max_gene_num,
                                double min_gene_score, double max_gene_score,
                                double min_average_splice_site_prob,
                                unsigned long feature_num)
{
  GtNodeStream *ns = gt_node_stream_create(gt_filter_stream_class(),
                                          gt_node_stream_is_sorted(in_stream));
  GtFilterStream *filter_stream = gt_filter_stream_cast(ns);
  gt_assert(in_stream);
  filter_stream->in_stream = gt_node_stream_ref(in_stream);
  filter_stream->filter_visitor = gt_filter_visitor_new(seqid, typefilter,
                                                     contain_range,
                                                     overlap_range, strand,
                                                     targetstrand, has_CDS,
                                                     max_gene_length,
                                                     max_gene_num,
                                                     min_gene_score,
                                                     max_gene_score,
                                                  min_average_splice_site_prob,
                                                     feature_num);
  return ns;
}
