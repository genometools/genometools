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

#include "libgtext/genome_stream_rep.h"
#include "libgtext/gff3_out_stream.h"
#include "libgtext/gff3_visitor.h"

struct GFF3OutStream {
  const GenomeStream parent_instance;
  GenomeStream *in_stream;
  GenomeVisitor *gff3_visitor;
};

#define gff3_out_stream_cast(GS)\
        genome_stream_cast(gff3_out_stream_class(), GS);

static int gff3_out_stream_next_tree(GenomeStream *gs, GenomeNode **gn,
                                     Error *e)
{
  GFF3OutStream *gff3_out_stream;
  int had_err;
  error_check(e);
  gff3_out_stream = gff3_out_stream_cast(gs);
  had_err = genome_stream_next_tree(gff3_out_stream->in_stream, gn, e);
  if (!had_err && *gn)
    had_err = genome_node_accept(*gn, gff3_out_stream->gff3_visitor, e);
  return had_err;
}

static void gff3_out_stream_free(GenomeStream *gs)
{
  GFF3OutStream *gff3_out_stream = gff3_out_stream_cast(gs);
  genome_stream_delete(gff3_out_stream->in_stream);
  genome_visitor_delete(gff3_out_stream->gff3_visitor);
}

const GenomeStreamClass* gff3_out_stream_class(void)
{
  static const GenomeStreamClass gsc = { sizeof (GFF3OutStream),
                                         gff3_out_stream_next_tree,
                                         gff3_out_stream_free };
  return &gsc;
}

GenomeStream* gff3_out_stream_new(GenomeStream *in_stream, GenFile *outfp)
{
  GenomeStream *gs = genome_stream_create(gff3_out_stream_class(),
                                          genome_stream_is_sorted(in_stream));
  GFF3OutStream *gff3_out_stream = gff3_out_stream_cast(gs);
  gff3_out_stream->in_stream = genome_stream_ref(in_stream);
  gff3_out_stream->gff3_visitor = gff3_visitor_new(outfp);
  return gs;
}
