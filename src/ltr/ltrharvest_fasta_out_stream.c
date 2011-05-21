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
#include "extended/node_visitor_api.h"
#include "extended/visitor_stream_api.h"
#include "ltr/ltrharvest_fasta_out_stream.h"
#include "ltr/ltrharvest_fasta_out_visitor.h"

GtNodeStream* gt_ltrharvest_fasta_out_stream_new(GtNodeStream *in_stream,
                                                 bool inner,
                                                 const GtEncseq *encseq,
                                                 unsigned long width,
                                                 GtFile *outfp)
{
  GtNodeVisitor *nv;
  nv = gt_ltrharvest_fasta_out_visitor_new(encseq, inner, width, outfp);
  return gt_visitor_stream_new(in_stream, nv);
}
