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

#include "extended/visitor_stream_api.h"
#include "ltr/ltrharvest_tabout_stream.h"
#include "ltr/ltrharvest_tabout_visitor.h"

GtNodeStream* gt_ltrharvest_tabout_stream_new(GtNodeStream *in_stream,
                                              GtNodeVisitor *v)
{
  return gt_visitor_stream_new(in_stream, v);
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
