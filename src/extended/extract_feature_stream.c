/*
  Copyright (c) 2006-2011 Gordon Gremme <gordon@gremme.org>
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
#include "core/trans_table_api.h"
#include "extended/extract_feature_stream_api.h"
#include "extended/extract_feature_visitor.h"
#include "extended/visitor_stream.h"

GtNodeStream* gt_extract_feature_stream_new(GtNodeStream *in_stream,
                                            GtRegionMapping *rm,
                                            const char *type, bool join,
                                            bool translate, bool seqid,
                                            bool target, GtUword width,
                                            GtFile *outfp)
{
  GtNodeVisitor *nv = gt_extract_feature_visitor_new(rm, type, join, translate,
                                                     seqid, target, width,
                                                     outfp);
  return gt_visitor_stream_new(in_stream, nv);
}

void gt_extract_feature_stream_retain_id_attributes(GtExtractFeatureStream *es)
{
  gt_assert(es);
  gt_extract_feature_visitor_retain_id_attributes((GtExtractFeatureVisitor*)
                         gt_visitor_stream_get_visitor((GtVisitorStream*) es));
}

void gt_extract_feature_stream_set_trans_table(GtExtractFeatureStream *es,
                                               GtTransTable *table)
{
  gt_assert(es);
  gt_extract_feature_visitor_set_trans_table((GtExtractFeatureVisitor*)
                         gt_visitor_stream_get_visitor((GtVisitorStream*) es),
                         table);
}

void gt_extract_feature_stream_show_coords(GtExtractFeatureStream *es)
{
  gt_assert(es);
  gt_extract_feature_visitor_show_coords((GtExtractFeatureVisitor*)
                         gt_visitor_stream_get_visitor((GtVisitorStream*) es));
}