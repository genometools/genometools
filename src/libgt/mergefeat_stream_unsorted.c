/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include "genome_stream_rep.h"
#include "mergefeat_stream_unsorted.h"
#include "mergefeat_visitor.h"

struct Mergefeat_stream_unsorted {
  const Genome_stream parent_instance;
  Genome_stream *in_stream;
  Genome_visitor *mergefeat_visitor;
};

#define mergefeat_stream_unsorted_cast(GS)\
        genome_stream_cast(mergefeat_stream_unsorted_class(), GS)

Genome_node* mergefeat_stream_unsorted_next_tree(Genome_stream *gs, Log *l)
{
  Mergefeat_stream_unsorted *mfs = mergefeat_stream_unsorted_cast(gs);
  Genome_node *gn = genome_stream_next_tree(mfs->in_stream, l);
  if (gn)
    genome_node_accept(gn, mfs->mergefeat_visitor, l);
  return gn;
}

static void mergefeat_stream_unsorted_free(Genome_stream *gs)
{
  Mergefeat_stream_unsorted *mfs = mergefeat_stream_unsorted_cast(gs);
  genome_visitor_free(mfs->mergefeat_visitor);
}

const Genome_stream_class* mergefeat_stream_unsorted_class(void)
{
  static const Genome_stream_class gsc = { sizeof(Mergefeat_stream_unsorted),
                                           mergefeat_stream_unsorted_next_tree,
                                           mergefeat_stream_unsorted_free };
  return &gsc;
}

Genome_stream* mergefeat_stream_unsorted_new(Genome_stream *in_stream)
{
  Genome_stream *gs = genome_stream_create(mergefeat_stream_unsorted_class(),0);
  Mergefeat_stream_unsorted *mfs = mergefeat_stream_unsorted_cast(gs);
  mfs->in_stream = in_stream;
  mfs->mergefeat_visitor = mergefeat_visitor_new();
  return gs;
}
