/*
  Copyright (c) 2006-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2008 Center for Bioinformatics, University of Hamburg

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
#include "libgtext/genome_stream_rep.h"
#include "libgtext/sort_stream.h"

struct SortStream
{
  const GenomeStream parent_instance;
  GenomeStream *in_stream;
  unsigned long idx;
  Array *trees;
  bool sorted;
};

#define sort_stream_cast(GS)\
        genome_stream_cast(sort_stream_class(), GS);

static int sort_stream_next_tree(GenomeStream *gs, GenomeNode **gn, Error *e)
{
  SortStream *sort_stream;
  int had_err = 0;
  error_check(e);
  sort_stream = sort_stream_cast(gs);

  if (!sort_stream->sorted) {
    while (!(had_err = genome_stream_next_tree(sort_stream->in_stream, gn, e))
           && *gn) {
      array_add(sort_stream->trees, *gn);
    }
    if (!had_err) {
      genome_nodes_sort_stable(sort_stream->trees);
      sort_stream->sorted = true;
    }
  }

  if (!had_err) {
    assert(sort_stream->sorted);
    if (sort_stream->idx < array_size(sort_stream->trees)) {
      *gn = *(GenomeNode**) array_get(sort_stream->trees, sort_stream->idx);
      sort_stream->idx++;
      return 0;
    }
  }

  if (!had_err) {
    array_reset(sort_stream->trees);
    *gn = NULL;
  }

  return had_err;
}

static void sort_stream_free(GenomeStream *gs)
{
  unsigned long i;
  SortStream *sort_stream = sort_stream_cast(gs);
  for (i = 0; i < array_size(sort_stream->trees); i++)
    genome_node_rec_delete(*(GenomeNode**) array_get(sort_stream->trees, i));
  array_delete(sort_stream->trees);
  genome_stream_delete(sort_stream->in_stream);
}

const GenomeStreamClass* sort_stream_class(void)
{
  static const GenomeStreamClass gsc = { sizeof (SortStream),
                                         sort_stream_next_tree,
                                         sort_stream_free };
  return &gsc;
}

GenomeStream* sort_stream_new(GenomeStream *in_stream)
{
  GenomeStream *gs = genome_stream_create(sort_stream_class(), true);
  SortStream *sort_stream = sort_stream_cast(gs);
  assert(in_stream);
  sort_stream->in_stream = genome_stream_ref(in_stream);
  sort_stream->sorted = false;
  sort_stream->idx = 0;
  sort_stream->trees = array_new(sizeof (GenomeNode*));
  return gs;
}
