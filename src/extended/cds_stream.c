/*
  Copyright (c) 2006-2011 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#include "extended/cds_stream_api.h"
#include "extended/cds_visitor.h"
#include "extended/visitor_stream_api.h"

GtNodeStream* gt_cds_stream_new(GtNodeStream *in_stream, GtRegionMapping *rm,
                                unsigned int minorflen, const char *source,
                                bool start_codon, bool final_stop_codon,
                                bool generic_start_codons)
{
  GtNodeVisitor *nv;
  GtStr *source_str;
  source_str = gt_str_new_cstr(source);
  nv = gt_cds_visitor_new(rm, minorflen, source_str, start_codon,
                          final_stop_codon, generic_start_codons);
  gt_str_delete(source_str);
  return gt_visitor_stream_new(in_stream, nv);
}
