/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "genome_stream_rep.h"
#include "gff3_out_stream.h"
#include "gff3_visitor.h"

struct Gff3_out_stream {
  const Genome_stream parent_instance;
  Genome_stream *in_stream;
  bool first_genome_feature;
  Genome_visitor *gff3_visitor;
};

#define gff3_out_stream_cast(GS)\
        genome_stream_cast(gff3_out_stream_class(), GS);

static Genome_node* gff3_out_stream_next_tree(Genome_stream *gs, Log *l)
{
  Gff3_out_stream *gff3_out_stream;
  Genome_node *gn;
  gff3_out_stream = gff3_out_stream_cast(gs);
  gn = genome_stream_next_tree(gff3_out_stream->in_stream, l);
  if (gn)
    genome_node_accept(gn, gff3_out_stream->gff3_visitor, l);
  return gn;
}

static void gff3_out_stream_free(Genome_stream *gs)
{
  Gff3_out_stream *gff3_out_stream = gff3_out_stream_cast(gs);
  genome_visitor_free(gff3_out_stream->gff3_visitor);
}

const Genome_stream_class* gff3_out_stream_class(void)
{
  static const Genome_stream_class gsc = { sizeof (Gff3_out_stream),
                                           gff3_out_stream_next_tree,
                                           gff3_out_stream_free };
  return &gsc;
}

Genome_stream* gff3_out_stream_new(Genome_stream *in_stream, FILE *outfp)
{
  Genome_stream *gs = genome_stream_create(gff3_out_stream_class(),
                                           genome_stream_is_sorted(in_stream));
  Gff3_out_stream *gff3_out_stream = gff3_out_stream_cast(gs);
  gff3_out_stream->in_stream = in_stream;
  gff3_out_stream->first_genome_feature = true;
  gff3_out_stream->gff3_visitor = gff3_visitor_new(outfp);
  return gs;
}
