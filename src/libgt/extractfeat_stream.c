/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "extractfeat_stream.h"
#include "extractfeat_visitor.h"
#include "genome_stream_rep.h"

struct Extractfeat_stream
{
  const Genome_stream parent_instance;
  Genome_stream *in_stream;
  Genome_visitor *extractfeat_visitor;
};

#define extractfeat_stream_cast(GS)\
        genome_stream_cast(extractfeat_stream_class(), GS)

static Genome_node* extractfeat_stream_next_tree(Genome_stream *gs, Log *l)
{
  Extractfeat_stream *extractfeat_stream = extractfeat_stream_cast(gs);
  Genome_node *gn = genome_stream_next_tree(extractfeat_stream->in_stream, l);
  if (gn)
    genome_node_accept(gn, extractfeat_stream->extractfeat_visitor, l);
  return gn;
}

static void extractfeat_stream_free(Genome_stream *gs)
{
  Extractfeat_stream *extractfeat_stream = extractfeat_stream_cast(gs);
  genome_visitor_free(extractfeat_stream->extractfeat_visitor);
}

const Genome_stream_class* extractfeat_stream_class(void)
{
  static const Genome_stream_class gsc = { sizeof(Extractfeat_stream),
                                           extractfeat_stream_next_tree,
                                           extractfeat_stream_free };
  return &gsc;
}

Genome_stream* extractfeat_stream_new(Genome_stream *in_stream,
                                      Str *sequence_file,
                                      Genome_feature_type type,
                                      unsigned int join, unsigned int translate)
{
  Genome_stream *gs = genome_stream_create(extractfeat_stream_class(), 1);
  Extractfeat_stream *extractfeat_stream = extractfeat_stream_cast(gs);

  extractfeat_stream->in_stream = in_stream;
  extractfeat_stream->extractfeat_visitor =
    extractfeat_visitor_new(sequence_file, type, join, translate);

  return gs;
}
