/*
  Copyright (c) 2015 Sascha Steinbiss <sascha@steinbiss.name>

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
#include "core/cstr_table.h"
#include "core/queue.h"
#include "extended/genome_node.h"
#include "extended/gff3_numsorted_out_stream.h"
#include "extended/gff3_visitor_api.h"
#include "extended/node_stream_api.h"

struct GtGFF3NumsortedOutStream {
  const GtNodeStream parent_instance;
  GtNodeStream *in_stream;
  GtArray *buffer;
  GtQueue *outqueue;
  GtNodeVisitor *gff3_visitor;
};

#define gff3_numsorted_out_stream_cast(GS)\
        gt_node_stream_cast(gt_gff3_numsorted_out_stream_class(), GS);

static int gff3_numsorted_out_stream_next(GtNodeStream *ns, GtGenomeNode **gn,
                                           GtError *err)
{
  GtGFF3NumsortedOutStream *gff3_out_stream;
  int had_err = 0;
  GtUword i = 0;
  gt_error_check(err);
  gff3_out_stream = gff3_numsorted_out_stream_cast(ns);
  if (!gff3_out_stream->outqueue) {
    gff3_out_stream->outqueue = gt_queue_new();
    while (!(had_err =
                    gt_node_stream_next(gff3_out_stream->in_stream, gn, err))) {
      if (!*gn) break;
      gt_array_add(gff3_out_stream->buffer, *gn);
    }
    if (!had_err) {
      gt_genome_nodes_sort_stable_with_func(gff3_out_stream->buffer,
                             (GtCompare) gt_genome_node_compare_numeric_seqids);
      for (i = 0; !had_err && i < gt_array_size(gff3_out_stream->buffer); i++) {
        GtGenomeNode *mygn = *(GtGenomeNode**)
                                       gt_array_get(gff3_out_stream->buffer, i);
        gt_queue_add(gff3_out_stream->outqueue, mygn);
      }
    }
  }
  if (gff3_out_stream->outqueue && !had_err) {
    if (gt_queue_size(gff3_out_stream->outqueue) > 0) {
      GtGenomeNode *mygn = (GtGenomeNode*)
                                        gt_queue_get(gff3_out_stream->outqueue);
      gt_assert(mygn);
      had_err = gt_genome_node_accept(mygn, gff3_out_stream->gff3_visitor, err);
      if (!had_err)
        *gn = mygn;
    }
  }
  return had_err;
}

static void gff3_numsorted_out_stream_free(GtNodeStream *ns)
{
  GtGFF3NumsortedOutStream *gff3_out_stream =
                                             gff3_numsorted_out_stream_cast(ns);
  gt_node_stream_delete(gff3_out_stream->in_stream);
  gt_node_visitor_delete(gff3_out_stream->gff3_visitor);
  gt_array_delete(gff3_out_stream->buffer);
  gt_queue_delete(gff3_out_stream->outqueue);
}

const GtNodeStreamClass* gt_gff3_numsorted_out_stream_class(void)
{
  static const GtNodeStreamClass *nsc = NULL;
  gt_class_alloc_lock_enter();
  if (!nsc) {
    nsc = gt_node_stream_class_new(sizeof (GtGFF3NumsortedOutStream),
                                   gff3_numsorted_out_stream_free,
                                   gff3_numsorted_out_stream_next);
  }
  gt_class_alloc_lock_leave();
  return nsc;
}

GtNodeStream* gt_gff3_numsorted_out_stream_new(GtNodeStream *in_stream,
                                               GtFile *outfp)
{
  GtNodeStream *ns = gt_node_stream_create(gt_gff3_numsorted_out_stream_class(),
                                           false);
  GtGFF3NumsortedOutStream *gff3_out_stream =
                                             gff3_numsorted_out_stream_cast(ns);
  gff3_out_stream->in_stream = gt_node_stream_ref(in_stream);
  gff3_out_stream->buffer = gt_array_new(sizeof (GtGenomeNode*));
  gff3_out_stream->outqueue = NULL;
  gff3_out_stream->gff3_visitor = gt_gff3_visitor_new(outfp);
  return ns;
}

void gt_gff3_numsorted_out_stream_set_fasta_width(GtGFF3NumsortedOutStream
                                                               *gff3_out_stream,
                                                  GtUword fasta_width)
{
  gt_assert(gff3_out_stream);
  gt_gff3_visitor_set_fasta_width((GtGFF3Visitor*)
                                  gff3_out_stream->gff3_visitor, fasta_width);
}

void gt_gff3_numsorted_out_stream_retain_id_attributes(GtGFF3NumsortedOutStream
                                                               *gff3_out_stream)
{
  gt_assert(gff3_out_stream);
  gt_gff3_visitor_retain_id_attributes((GtGFF3Visitor*)
                                       gff3_out_stream->gff3_visitor);
}
