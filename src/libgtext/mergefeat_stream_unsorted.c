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
#include "libgtext/mergefeat_stream_unsorted.h"
#include "libgtext/mergefeat_visitor.h"

struct MergefeatStreamUnsorted {
  const GenomeStream parent_instance;
  GenomeStream *in_stream;
  GenomeVisitor *mergefeat_visitor;
};

#define mergefeat_stream_unsorted_cast(GS)\
        genome_stream_cast(mergefeat_stream_unsorted_class(), GS)

static int mergefeat_stream_unsorted_next_tree(GenomeStream *gs,
                                               GenomeNode **gn, Error *e)
{
  MergefeatStreamUnsorted *mfs;
  int had_err;
  error_check(e);
  mfs = mergefeat_stream_unsorted_cast(gs);
  had_err = genome_stream_next_tree(mfs->in_stream, gn, e);
  if (!had_err && *gn)
    had_err = genome_node_accept(*gn, mfs->mergefeat_visitor, e);
  return had_err;
}

static void mergefeat_stream_unsorted_free(GenomeStream *gs)
{
  MergefeatStreamUnsorted *mfs = mergefeat_stream_unsorted_cast(gs);
  genome_visitor_delete(mfs->mergefeat_visitor);
}

const GenomeStreamClass* mergefeat_stream_unsorted_class(void)
{
  static const GenomeStreamClass gsc = { sizeof (MergefeatStreamUnsorted),
                                         mergefeat_stream_unsorted_next_tree,
                                         mergefeat_stream_unsorted_free };
  return &gsc;
}

GenomeStream* mergefeat_stream_unsorted_new(GenomeStream *in_stream)
{
  GenomeStream *gs = genome_stream_create(mergefeat_stream_unsorted_class(),
                                          false);
  MergefeatStreamUnsorted *mfs = mergefeat_stream_unsorted_cast(gs);
  mfs->in_stream = in_stream;
  mfs->mergefeat_visitor = mergefeat_visitor_new();
  return gs;
}
