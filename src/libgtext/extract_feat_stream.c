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

#include <assert.h>
#include "libgtext/extract_feat_stream.h"
#include "libgtext/extract_feat_visitor.h"
#include "libgtext/genome_stream_rep.h"

struct ExtractFeatStream
{
  const GenomeStream parent_instance;
  GenomeStream *in_stream;
  GenomeVisitor *extract_feat_visitor;
};

#define extract_feat_stream_cast(GS)\
        genome_stream_cast(extract_feat_stream_class(), GS)

static int extract_feat_stream_next_tree(GenomeStream *gs, GenomeNode **gn,
                                         Error *err)
{
  ExtractFeatStream *efs;
  int had_err;
  error_check(err);
  efs = extract_feat_stream_cast(gs);
  had_err = genome_stream_next_tree(efs->in_stream, gn, err);
  if (!had_err) {
    assert(efs->extract_feat_visitor);
    if (*gn) {
      had_err = genome_node_accept(*gn, efs->extract_feat_visitor, err);
      if (had_err) {
        /* we own the node -> delete it */
        genome_node_rec_delete(*gn);
        *gn = NULL;
      }
    }
  }
  return had_err;
}

static void extract_feat_stream_free(GenomeStream *gs)
{
  ExtractFeatStream *extract_feat_stream = extract_feat_stream_cast(gs);
  genome_visitor_delete(extract_feat_stream->extract_feat_visitor);
  genome_stream_delete(extract_feat_stream->in_stream);
}

const GenomeStreamClass* extract_feat_stream_class(void)
{
  static const GenomeStreamClass gsc = { sizeof (ExtractFeatStream),
                                         extract_feat_stream_next_tree,
                                         extract_feat_stream_free };
  return &gsc;
}

GenomeStream* extract_feat_stream_new(GenomeStream *in_stream,
                                      RegionMapping *rm,
                                      GenomeFeatureType type,
                                      bool join, bool translate)
{
  GenomeStream *gs = genome_stream_create(extract_feat_stream_class(), true);
  ExtractFeatStream *efs = extract_feat_stream_cast(gs);
  efs->in_stream = genome_stream_ref(in_stream);
  efs->extract_feat_visitor = extract_feat_visitor_new(rm, type, join,
                                                       translate);
  return gs;
}
