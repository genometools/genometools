/*
  Copyright (c) 2009-2011 Gordon Gremme <gordon@gremme.org>

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
#include "extended/dup_feature_stream_api.h"
#include "extended/dup_feature_visitor.h"
#include "extended/visitor_stream_api.h"

GtNodeStream* gt_dup_feature_stream_new(GtNodeStream *in_stream,
                                        const char *dest_type,
                                        const char *source_type)
{
  GtNodeVisitor *nv;
  gt_assert(in_stream);
  nv = gt_dup_feature_visitor_new(dest_type, source_type);
  return gt_visitor_stream_new(in_stream, nv);
}
