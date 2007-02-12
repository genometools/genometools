/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include "genome_stream_rep.h"
#include "mergefeat_stream_sorted.h"
#include "mergefeat_stream_unsorted.h"
#include "sort_stream.h"

struct Mergefeat_stream_sorted {
  const GenomeStream parent_instance;
  GenomeStream *mergefeat_stream_unsorted,
                *sort_stream;
};

#define mergefeat_stream_sorted_cast(GS)\
        genome_stream_cast(mergefeat_stream_sorted_class(), GS)

GenomeNode* mergefeat_stream_sorted_next_tree(GenomeStream *gs, Log *l)
{
  Mergefeat_stream_sorted *mfs = mergefeat_stream_sorted_cast(gs);
  return genome_stream_next_tree(mfs->sort_stream, l);
}

static void mergefeat_stream_sorted_free(GenomeStream *gs)
{
  Mergefeat_stream_sorted *mfs = mergefeat_stream_sorted_cast(gs);
  genome_stream_free(mfs->mergefeat_stream_unsorted);
  genome_stream_free(mfs->sort_stream);
}

const GenomeStreamClass* mergefeat_stream_sorted_class(void)
{
  static const GenomeStreamClass gsc = { sizeof(Mergefeat_stream_sorted),
                                         mergefeat_stream_sorted_next_tree,
                                         mergefeat_stream_sorted_free };
  return &gsc;
}

GenomeStream* mergefeat_stream_sorted_new(GenomeStream *in_stream)
{
  GenomeStream *gs = genome_stream_create(mergefeat_stream_sorted_class(),
                                          true);
  Mergefeat_stream_sorted *mfs = mergefeat_stream_sorted_cast(gs);
  assert(in_stream && genome_stream_is_sorted(in_stream));
  mfs->mergefeat_stream_unsorted = mergefeat_stream_unsorted_new(in_stream);
  mfs->sort_stream = sort_stream_new(mfs->mergefeat_stream_unsorted);
  return gs;
}
