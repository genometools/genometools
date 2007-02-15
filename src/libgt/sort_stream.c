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
  const GenomeStream parent_instance;
  GenomeStream *in_stream;
  unsigned long idx;
  Array *trees;
  bool sorted;
};

#define sort_stream_cast(GS)\
        genome_stream_cast(sort_stream_class(), GS);

static int sort_stream_next_tree(GenomeStream *gs, GenomeNode **gn, Log *l,
                                 Error *err)
{
  Sort_stream *sort_stream;
  int has_err = 0;
  error_check(err);
  sort_stream = sort_stream_cast(gs);

  if (!sort_stream->sorted) {
    while (!(has_err = genome_stream_next_tree(sort_stream->in_stream, gn, l,
             err)) && *gn) {
      array_add(sort_stream->trees, *gn);
    }
    if (!has_err) {
      genome_nodes_sort_stable(sort_stream->trees);
      sort_stream->sorted = true;
    }
  }

  if (!has_err) {
    assert(sort_stream->sorted);
    if (sort_stream->idx < array_size(sort_stream->trees)) {
      *gn = *(GenomeNode**) array_get(sort_stream->trees, sort_stream->idx);
      sort_stream->idx++;
      return 0;
    }
  }

  if (!has_err) {
    array_set_size(sort_stream->trees, 0);
    *gn = NULL;
  }

  return has_err;
}

static void sort_stream_free(GenomeStream *gs)
{
  Sort_stream *sort_stream = sort_stream_cast(gs);
  assert(!array_size(sort_stream->trees));
  array_free(sort_stream->trees);
}

const GenomeStreamClass* sort_stream_class(void)
{
  static const GenomeStreamClass gsc = { sizeof (Sort_stream),
                                         sort_stream_next_tree,
                                         sort_stream_free };
  return &gsc;
}

GenomeStream* sort_stream_new(GenomeStream *in_stream)
{
  GenomeStream *gs = genome_stream_create(sort_stream_class(), true);
  Sort_stream *sort_stream = sort_stream_cast(gs);
  assert(in_stream);
  sort_stream->in_stream = in_stream;
  sort_stream->sorted = false;
  sort_stream->idx = 0;
  sort_stream->trees = array_new(sizeof(GenomeNode*));
  return gs;
}
