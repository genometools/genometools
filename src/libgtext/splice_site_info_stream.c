/*
  Copyright (c) 2007-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007-2008 Center for Bioinformatics, University of Hamburg

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
#include "libgtext/splice_site_info_stream.h"
#include "libgtext/splicesiteinfo_visitor.h"
#include "libgtext/genome_stream_rep.h"

struct SpliceSiteInfoStream
{
  const GenomeStream parent_instance;
  GenomeStream *in_stream;
  GenomeVisitor *splice_site_info_visitor;
};

#define splice_site_info_stream_cast(GS)\
        genome_stream_cast(splice_site_info_stream_class(), GS)

static int splice_site_info_stream_next_tree(GenomeStream *gs, GenomeNode **gn,
                                             Error *err)
{
  SpliceSiteInfoStream *ssis;
  int had_err;
  error_check(err);
  ssis = splice_site_info_stream_cast(gs);
  had_err = genome_stream_next_tree(ssis->in_stream, gn, err);
  if (!had_err) {
    assert(ssis->splice_site_info_visitor);
    if (*gn) {
      had_err = genome_node_accept(*gn, ssis->splice_site_info_visitor, err);
      if (had_err) {
        /* we own the node -> delete it */
        genome_node_rec_delete(*gn);
        *gn = NULL;
      }
    }
  }
  return had_err;
}

static void splice_site_info_stream_free(GenomeStream *gs)
{
  SpliceSiteInfoStream *ssis = splice_site_info_stream_cast(gs);
  genome_visitor_delete(ssis->splice_site_info_visitor);
  genome_stream_delete(ssis->in_stream);
}

const GenomeStreamClass* splice_site_info_stream_class(void)
{
  static const GenomeStreamClass gsc = { sizeof (SpliceSiteInfoStream),
                                         splice_site_info_stream_next_tree,
                                         splice_site_info_stream_free };
  return &gsc;
}

GenomeStream* splice_site_info_stream_new(GenomeStream *in_stream,
                                          RegionMapping *rm)
{
  GenomeStream *gs = genome_stream_create(splice_site_info_stream_class(),
                                          false);
  SpliceSiteInfoStream *ssis = splice_site_info_stream_cast(gs);
  ssis->in_stream = genome_stream_ref(in_stream);
  ssis->splice_site_info_visitor = splicesiteinfo_visitor_new(rm);
  return gs;
}

bool splice_site_info_stream_show(GenomeStream *gs)
{
  SpliceSiteInfoStream *ssis;
  assert(gs);
  ssis = splice_site_info_stream_cast(gs);
  return splicesiteinfo_visitor_show(ssis->splice_site_info_visitor);
}
