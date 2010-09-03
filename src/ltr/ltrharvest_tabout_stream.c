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

#include "extended/node_stream_api.h"
#include "ltr/ltrharvest_tabout_stream.h"
#include "ltr/ltrharvest_tabout_visitor.h"

struct GtLTRharvestTaboutStream {
  const GtNodeStream parent_instance;
  GtNodeStream *in_stream;
  GtNodeVisitor *v;
};

#define ltrharvest_tabout_stream_cast(GS)\
        gt_node_stream_cast(gt_ltrharvest_tabout_stream_class(), GS)

static int ltrharvest_tabout_stream_next(GtNodeStream *ns, GtGenomeNode **gn,
                                     GtError *err)
{
  GtLTRharvestTaboutStream *ltrharvest_tabout_stream;
  int had_err;
  gt_error_check(err);
  ltrharvest_tabout_stream = ltrharvest_tabout_stream_cast(ns);
  had_err = gt_node_stream_next(ltrharvest_tabout_stream->in_stream, gn, err);
  gt_assert(ltrharvest_tabout_stream->v != NULL);
  if (!had_err && *gn != NULL) {
    had_err = gt_genome_node_accept(*gn, ltrharvest_tabout_stream->v, err);
  }
  return had_err;
}

static void ltrharvest_tabout_stream_free(GtNodeStream *ns)
{
  GtLTRharvestTaboutStream *ltrharvest_tabout_stream =
                                              ltrharvest_tabout_stream_cast(ns);
  gt_node_stream_delete(ltrharvest_tabout_stream->in_stream);
}

const GtNodeStreamClass* gt_ltrharvest_tabout_stream_class(void)
{
  static const GtNodeStreamClass *nsc = NULL;
  if (!nsc) {
    nsc = gt_node_stream_class_new(sizeof (GtLTRharvestTaboutStream),
                                   ltrharvest_tabout_stream_free,
                                   ltrharvest_tabout_stream_next);
  }
  return nsc;
}

GtNodeStream* gt_ltrharvest_tabout_stream_new(GtNodeStream *in_stream,
                                              GtNodeVisitor *v)
{
  GtLTRharvestTaboutStream *ltrharvest_tabout_stream;
  GtNodeStream *ns;
  ns = gt_node_stream_create(gt_ltrharvest_tabout_stream_class(), false);
  ltrharvest_tabout_stream = ltrharvest_tabout_stream_cast(ns);
  ltrharvest_tabout_stream->in_stream = gt_node_stream_ref(in_stream);
  ltrharvest_tabout_stream->v = v;
  return ns;
}

void gt_ltrharvest_tabout_stream_printshortheader(void)
{
  printf("# predictions are reported in the following way\n");
  printf("# s(ret) e(ret) l(ret) s(lLTR) e(lLTR) l(lLTR)"
      " s(rLTR) e(rLTR) l(rLTR) sim(LTRs) seq-nr \n");
  printf("# where:\n");
  printf("# s = starting position\n");
  printf("# e = ending position\n");
  printf("# l = length\n");
  printf("# ret = LTR-retrotransposon\n");
  printf("# lLTR = left LTR\n");
  printf("# rLTR = right LTR\n");
  printf("# sim = similarity\n");
  printf("# seq-nr = sequence number\n");
}

void gt_ltrharvest_tabout_stream_printlongheader(bool withtsd, bool withmotif)
{
  printf("# predictions are reported in the following way\n");
  printf("# s(ret) e(ret) l(ret) ");
  printf("s(lLTR) e(lLTR) l(lLTR)");
  if (withtsd)
  {
    printf(" TSD l(TSD)");
  }
  if (withmotif)
  {
    printf(" m(lLTR)");
  }
  printf(" s(rLTR) e(rLTR) l(rLTR)");
  if (withtsd)
  {
    printf(" TSD l(TSD)");
  }
  if (withmotif)
  {
    printf(" m(rLTR)");
  }
  printf(" sim(LTRs)");
  printf(" seq-nr");
  printf("\n# where:\n");
  printf("# s = starting position\n");
  printf("# e = ending position\n");
  printf("# l = length\n");
  if (withmotif)
  {
    printf("# m = motif\n");
  }
  printf("# ret = LTR-retrotransposon\n");
  printf("# lLTR = left LTR\n");
  printf("# rLTR = right LTR\n");
  if (withtsd)
  {
    printf("# TSD = target site duplication\n");
  }
  printf("# sim = similarity\n");
  printf("# seq-nr = sequence number\n");
}
