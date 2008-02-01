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

#include "libgtext/genome_stream_rep.h"
#include "libgtext/gtf_out_stream.h"
#include "libgtext/gtf_visitor.h"

struct GTFOutStream {
  const GenomeStream parent_instance;
  GenomeStream *in_stream;
  GenomeVisitor *gtf_visitor;
};

#define gtf_out_stream_cast(GS)\
        genome_stream_cast(gtf_out_stream_class(), GS);

static int gtf_out_stream_next_tree(GenomeStream *gs, GenomeNode **gn, Error *e)
{
  GTFOutStream *gtf_out_stream;
  int had_err;
  error_check(e);
  gtf_out_stream = gtf_out_stream_cast(gs);
  had_err = genome_stream_next_tree(gtf_out_stream->in_stream, gn, e);
  if (!had_err && *gn)
    had_err = genome_node_accept(*gn, gtf_out_stream->gtf_visitor, e);
  return had_err;
}

static void gtf_out_stream_free(GenomeStream *gs)
{
  GTFOutStream *gtf_out_stream = gtf_out_stream_cast(gs);
  genome_visitor_delete(gtf_out_stream->gtf_visitor);
  genome_stream_delete(gtf_out_stream->in_stream);
}

const GenomeStreamClass* gtf_out_stream_class(void)
{
  static const GenomeStreamClass gsc = { sizeof (GTFOutStream),
                                         gtf_out_stream_next_tree,
                                         gtf_out_stream_free };
  return &gsc;
}

GenomeStream* gtf_out_stream_new(GenomeStream *in_stream, GenFile *outfp)
{
  GenomeStream *gs = genome_stream_create(gtf_out_stream_class(),
                                          genome_stream_is_sorted(in_stream));
  GTFOutStream *gtf_out_stream = gtf_out_stream_cast(gs);
  gtf_out_stream->in_stream = genome_stream_ref(in_stream);
  gtf_out_stream->gtf_visitor = gtf_visitor_new(outfp);
  return gs;
}
