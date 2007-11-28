/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg

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

#include "libgtview/feature_stream.h"
#include "libgtview/feature_visitor.h"
#include "libgtview/feature_index.h"
#include "libgtext/genome_stream_rep.h"

struct FeatureStream {
  const GenomeStream parent_instance;
  GenomeStream *in_stream;
  GenomeVisitor *feature_visitor;
};

#define feature_stream_cast(GS)\
        genome_stream_cast(feature_stream_class(), GS)

static int feature_stream_next_tree(GenomeStream *gs, GenomeNode **gn, Error *e)
{
  FeatureStream *feature_stream;
  int had_err;
  error_check(e);
  feature_stream = feature_stream_cast(gs);
  had_err = genome_stream_next_tree(feature_stream->in_stream, gn, e);
  if (!had_err && *gn)
    had_err = genome_node_accept(*gn, feature_stream->feature_visitor, e);
  return had_err;
}

static void feature_stream_free(GenomeStream *gs)
{
  FeatureStream *feature_stream = feature_stream_cast(gs);
  genome_stream_delete(feature_stream->in_stream);
  genome_visitor_delete(feature_stream->feature_visitor);
}

const GenomeStreamClass* feature_stream_class(void)
{
  static const GenomeStreamClass gsc = { sizeof (FeatureStream),
                                         feature_stream_next_tree,
                                         feature_stream_free };
  return &gsc;
}

GenomeStream* feature_stream_new(GenomeStream *in_stream, FeatureIndex *fi)
{
  GenomeStream *gs;
  FeatureStream *feature_stream;
  gs = genome_stream_create(feature_stream_class(), true);
  feature_stream = feature_stream_cast(gs);
  feature_stream->in_stream = genome_stream_ref(in_stream);
  feature_stream->feature_visitor = feature_visitor_new(fi);
  return gs;
}
