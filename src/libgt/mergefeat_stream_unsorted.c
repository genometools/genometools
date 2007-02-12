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
  const GenomeStream parent_instance;
  GenomeStream *in_stream;
  GenomeVisitor *mergefeat_visitor;
};

#define mergefeat_stream_unsorted_cast(GS)\
        genome_stream_cast(mergefeat_stream_unsorted_class(), GS)

GenomeNode* mergefeat_stream_unsorted_next_tree(GenomeStream *gs, Log *l)
{
  Mergefeat_stream_unsorted *mfs = mergefeat_stream_unsorted_cast(gs);
  GenomeNode *gn = genome_stream_next_tree(mfs->in_stream, l);
  if (gn)
    genome_node_accept(gn, mfs->mergefeat_visitor, l);
  return gn;
}

static void mergefeat_stream_unsorted_free(GenomeStream *gs)
{
  Mergefeat_stream_unsorted *mfs = mergefeat_stream_unsorted_cast(gs);
  genome_visitor_free(mfs->mergefeat_visitor);
}

const GenomeStreamClass* mergefeat_stream_unsorted_class(void)
{
  static const GenomeStreamClass gsc = { sizeof(Mergefeat_stream_unsorted),
                                         mergefeat_stream_unsorted_next_tree,
                                         mergefeat_stream_unsorted_free };
  return &gsc;
}

GenomeStream* mergefeat_stream_unsorted_new(GenomeStream *in_stream)
{
  GenomeStream *gs = genome_stream_create(mergefeat_stream_unsorted_class(),
                                          false);
  Mergefeat_stream_unsorted *mfs = mergefeat_stream_unsorted_cast(gs);
  mfs->in_stream = in_stream;
  mfs->mergefeat_visitor = mergefeat_visitor_new();
  return gs;
}
