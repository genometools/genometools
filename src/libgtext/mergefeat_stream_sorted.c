/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg

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
#include "libgtext/mergefeat_stream_sorted.h"
#include "libgtext/mergefeat_stream_unsorted.h"
#include "libgtext/sort_stream.h"

struct MergefeatStreamSorted {
  const GenomeStream parent_instance;
  GenomeStream *mergefeat_stream_unsorted,
               *sort_stream;
};

#define mergefeat_stream_sorted_cast(GS)\
        genome_stream_cast(mergefeat_stream_sorted_class(), GS)

static int mergefeat_stream_sorted_next_tree(GenomeStream *gs, GenomeNode **gn,
                                             Error *e)
{
  MergefeatStreamSorted *mfs;
  error_check(e);
  mfs = mergefeat_stream_sorted_cast(gs);
  return genome_stream_next_tree(mfs->sort_stream, gn, e);
}

static void mergefeat_stream_sorted_free(GenomeStream *gs)
{
  MergefeatStreamSorted *mfs = mergefeat_stream_sorted_cast(gs);
  genome_stream_delete(mfs->mergefeat_stream_unsorted);
  genome_stream_delete(mfs->sort_stream);
}

const GenomeStreamClass* mergefeat_stream_sorted_class(void)
{
  static const GenomeStreamClass gsc = { sizeof (MergefeatStreamSorted),
                                         mergefeat_stream_sorted_next_tree,
                                         mergefeat_stream_sorted_free };
  return &gsc;
}

GenomeStream* mergefeat_stream_sorted_new(GenomeStream *in_stream)
{
  GenomeStream *gs = genome_stream_create(mergefeat_stream_sorted_class(),
                                          true);
  MergefeatStreamSorted *mfs = mergefeat_stream_sorted_cast(gs);
  assert(in_stream && genome_stream_is_sorted(in_stream));
  mfs->mergefeat_stream_unsorted = mergefeat_stream_unsorted_new(in_stream);
  mfs->sort_stream = sort_stream_new(mfs->mergefeat_stream_unsorted);
  return gs;
}
