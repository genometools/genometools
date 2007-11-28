/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

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
#include "libgtext/splicesiteinfo_stream.h"
#include "libgtext/splicesiteinfo_visitor.h"
#include "libgtext/genome_stream_rep.h"

struct SpliceSiteInfoStream
{
  const GenomeStream parent_instance;
  GenomeStream *in_stream;
  GenomeVisitor *splicesiteinfo_visitor;
};

#define splicesiteinfo_stream_cast(GS)\
        genome_stream_cast(splicesiteinfo_stream_class(), GS)

static int splicesiteinfo_stream_next_tree(GenomeStream *gs, GenomeNode **gn,
                                           Error *e)
{
  SpliceSiteInfoStream *splicesiteinfo_stream;
  int had_err;
  error_check(e);
  splicesiteinfo_stream = splicesiteinfo_stream_cast(gs);
  had_err = genome_stream_next_tree(splicesiteinfo_stream->in_stream, gn, e);
  if (!had_err) {
    assert(splicesiteinfo_stream->splicesiteinfo_visitor);
    if (*gn) {
      had_err = genome_node_accept(*gn, splicesiteinfo_stream
                                        ->splicesiteinfo_visitor, e);
      if (had_err) {
        /* we own the node -> delete it */
        genome_node_rec_delete(*gn);
        *gn = NULL;
      }
    }
  }
  return had_err;
}

static void splicesiteinfo_stream_free(GenomeStream *gs)
{
  SpliceSiteInfoStream *splicesiteinfo_stream = splicesiteinfo_stream_cast(gs);
  genome_visitor_delete(splicesiteinfo_stream->splicesiteinfo_visitor);
}

const GenomeStreamClass* splicesiteinfo_stream_class(void)
{
  static const GenomeStreamClass gsc = { sizeof (SpliceSiteInfoStream),
                                         splicesiteinfo_stream_next_tree,
                                         splicesiteinfo_stream_free };
  return &gsc;
}

GenomeStream* splicesiteinfo_stream_new(GenomeStream *in_stream,
                                        RegionMapping *rm)
{
  GenomeStream *gs = genome_stream_create(splicesiteinfo_stream_class(), false);
  SpliceSiteInfoStream *ssis = splicesiteinfo_stream_cast(gs);
  ssis->in_stream = in_stream;
  ssis->splicesiteinfo_visitor = splicesiteinfo_visitor_new(rm);
  return gs;
}

bool splicesiteinfo_stream_show(GenomeStream *gs)
{
  SpliceSiteInfoStream *ssis;
  assert(gs);
  ssis = splicesiteinfo_stream_cast(gs);
  return splicesiteinfo_visitor_show(ssis->splicesiteinfo_visitor);
}
