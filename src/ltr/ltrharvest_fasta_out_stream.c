/*
  Copyright (c) 2010 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2010 Center for Bioinformatics, University of Hamburg

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
#include "core/unused_api.h"
#include "extended/node_stream_api.h"
#include "extended/node_visitor.h"
#include "ltr/ltrharvest_fasta_out_stream.h"
#include "ltr/ltrharvest_fasta_out_visitor.h"

struct GtLTRharvestFastaOutStream
{
  const GtNodeStream parent_instance;
  GtNodeStream *in_stream;
  GtNodeVisitor *v;
  bool inner;
};

#define gt_ltrharvest_fasta_out_stream_cast(GS)\
        gt_node_stream_cast(gt_ltrharvest_fasta_out_stream_class(), GS)

static int ltrharvest_fasta_out_stream_next(GtNodeStream *gs, GtGenomeNode **gn,
                                            GtError *err)
{
  GtLTRharvestFastaOutStream *fos;
  int had_err;
  gt_error_check(err);
  fos = gt_ltrharvest_fasta_out_stream_cast(gs);
  had_err = gt_node_stream_next(fos->in_stream, gn, err);
  gt_assert(fos->v);
  if (!had_err && *gn != NULL) {
      had_err = gt_genome_node_accept(*gn, fos->v, err);
  }
  return had_err;
}

static void ltrharvest_fasta_out_stream_free(GtNodeStream *gs)
{
  GtLTRharvestFastaOutStream *fos = gt_ltrharvest_fasta_out_stream_cast(gs);
  gt_node_visitor_delete(fos->v);
  gt_node_stream_delete(fos->in_stream);
}

const GtNodeStreamClass* gt_ltrharvest_fasta_out_stream_class(void)
{
  static const GtNodeStreamClass *nsc = NULL;
  if (!nsc) {
    nsc = gt_node_stream_class_new(sizeof (GtLTRharvestFastaOutStream),
                                   ltrharvest_fasta_out_stream_free,
                                   ltrharvest_fasta_out_stream_next);
  }
  return nsc;
}

GtNodeStream* gt_ltrharvest_fasta_out_stream_new(GtNodeStream *in_stream,
                                                 bool inner,
                                                 const GtEncseq *encseq,
                                                 unsigned long width,
                                                 GtFile *outfp)
{
  GtNodeStream *gs;
  GtLTRharvestFastaOutStream *fos;
  gs = gt_node_stream_create(gt_ltrharvest_fasta_out_stream_class(), false);
  fos = gt_ltrharvest_fasta_out_stream_cast(gs);
  fos->in_stream = gt_node_stream_ref(in_stream);
  fos->v = gt_ltrharvest_fasta_out_visitor_new(encseq, inner, width, outfp);
  return gs;
}
