/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "genome_stream_rep.h"
#include "gff3_out_stream.h"
#include "gff3_visitor.h"

struct GFF3OutStream {
  const GenomeStream parent_instance;
  GenomeStream *in_stream;
  bool first_genome_feature;
  GenomeVisitor *gff3_visitor;
};

#define gff3_out_stream_cast(GS)\
        genome_stream_cast(gff3_out_stream_class(), GS);

static int gff3_out_stream_next_tree(GenomeStream *gs, GenomeNode **gn, Log *l,
                                     Error *err)
{
  GFF3OutStream *gff3_out_stream;
  int has_err;
  error_check(err);
  gff3_out_stream = gff3_out_stream_cast(gs);
  has_err = genome_stream_next_tree(gff3_out_stream->in_stream, gn, l, err);
  if (!has_err && *gn)
    has_err = genome_node_accept(*gn, gff3_out_stream->gff3_visitor, l, err);
  return has_err;
}

static void gff3_out_stream_free(GenomeStream *gs)
{
  GFF3OutStream *gff3_out_stream = gff3_out_stream_cast(gs);
  genome_visitor_free(gff3_out_stream->gff3_visitor);
}

const GenomeStreamClass* gff3_out_stream_class(void)
{
  static const GenomeStreamClass gsc = { sizeof (GFF3OutStream),
                                         gff3_out_stream_next_tree,
                                         gff3_out_stream_free };
  return &gsc;
}

GenomeStream* gff3_out_stream_new(GenomeStream *in_stream, FILE *outfp)
{
  GenomeStream *gs = genome_stream_create(gff3_out_stream_class(),
                                          genome_stream_is_sorted(in_stream));
  GFF3OutStream *gff3_out_stream = gff3_out_stream_cast(gs);
  gff3_out_stream->in_stream = in_stream;
  gff3_out_stream->first_genome_feature = true;
  gff3_out_stream->gff3_visitor = gff3_visitor_new(outfp);
  return gs;
}
