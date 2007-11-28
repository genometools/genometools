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

#include <assert.h>
#include "libgtext/extractfeat_stream.h"
#include "libgtext/extractfeat_visitor.h"
#include "libgtext/genome_stream_rep.h"

struct ExtractFeatStream
{
  const GenomeStream parent_instance;
  GenomeStream *in_stream;
  GenomeVisitor *extractfeat_visitor;
};

#define extractfeat_stream_cast(GS)\
        genome_stream_cast(extractfeat_stream_class(), GS)

static int extractfeat_stream_next_tree(GenomeStream *gs, GenomeNode **gn,
                                        Error *e)
{
  ExtractFeatStream *extractfeat_stream;
  int had_err;
  error_check(e);
  extractfeat_stream = extractfeat_stream_cast(gs);
  had_err = genome_stream_next_tree(extractfeat_stream->in_stream, gn, e);
  if (!had_err) {
    assert(extractfeat_stream->extractfeat_visitor);
    if (*gn) {
      had_err = genome_node_accept(*gn, extractfeat_stream->extractfeat_visitor,
                                   e);
      if (had_err) {
        /* we own the node -> delete it */
        genome_node_rec_delete(*gn);
        *gn = NULL;
      }
    }
  }
  return had_err;
}

static void extractfeat_stream_free(GenomeStream *gs)
{
  ExtractFeatStream *extractfeat_stream = extractfeat_stream_cast(gs);
  genome_visitor_delete(extractfeat_stream->extractfeat_visitor);
}

const GenomeStreamClass* extractfeat_stream_class(void)
{
  static const GenomeStreamClass gsc = { sizeof (ExtractFeatStream),
                                         extractfeat_stream_next_tree,
                                         extractfeat_stream_free };
  return &gsc;
}

GenomeStream* extractfeat_stream_new(GenomeStream *in_stream, RegionMapping *rm,
                                     GenomeFeatureType type,
                                     bool join, bool translate)
{
  GenomeStream *gs = genome_stream_create(extractfeat_stream_class(), true);
  ExtractFeatStream *efs = extractfeat_stream_cast(gs);
  efs->in_stream = in_stream;
  efs->extractfeat_visitor = extractfeat_visitor_new(rm, type, join, translate);
  return gs;
}
