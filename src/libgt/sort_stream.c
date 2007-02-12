/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include "fptr.h"
#include "genome_stream_rep.h"
#include "sort_stream.h"

struct Sort_stream
{
  const Genome_stream parent_instance;
  Genome_stream *in_stream;
  unsigned long idx;
  Array *trees;
  bool sorted;
};

#define sort_stream_cast(GS)\
        genome_stream_cast(sort_stream_class(), GS);

static Genome_node* sort_stream_next_tree(Genome_stream *gs, Log *l)
{
  Sort_stream *sort_stream;
  Genome_node *gn;

  sort_stream = sort_stream_cast(gs);

  if (!sort_stream->sorted) {
    while ((gn = genome_stream_next_tree(sort_stream->in_stream, l)) != NULL)
      array_add(sort_stream->trees, gn);
    genome_nodes_sort_stable(sort_stream->trees);
    sort_stream->sorted = true;
  }

  assert(sort_stream->sorted);
  if (sort_stream->idx < array_size(sort_stream->trees)) {
    gn = *(Genome_node**) array_get(sort_stream->trees, sort_stream->idx);
    sort_stream->idx++;
    return gn;
  }
  array_set_size(sort_stream->trees, 0);
  return NULL;
}

static void sort_stream_free(Genome_stream *gs)
{
  Sort_stream *sort_stream = sort_stream_cast(gs);
  assert(!array_size(sort_stream->trees));
  array_free(sort_stream->trees);
}

const Genome_stream_class* sort_stream_class(void)
{
  static const Genome_stream_class gsc = { sizeof (Sort_stream),
                                           sort_stream_next_tree,
                                           sort_stream_free };
  return &gsc;
}

Genome_stream* sort_stream_new(Genome_stream *in_stream)
{
  Genome_stream *gs = genome_stream_create(sort_stream_class(), 1);
  Sort_stream *sort_stream = sort_stream_cast(gs);
  assert(in_stream);
  sort_stream->in_stream = in_stream;
  sort_stream->sorted = false;
  sort_stream->idx = 0;
  sort_stream->trees = array_new(sizeof(Genome_node*));
  return gs;
}
