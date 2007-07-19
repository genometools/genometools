/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
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
                                               GenomeNode **gn, Env *env)
{
  MergefeatStreamUnsorted *mfs;
  int had_err;
  env_error_check(env);
  mfs = mergefeat_stream_unsorted_cast(gs);
  had_err = genome_stream_next_tree(mfs->in_stream, gn, env);
  if (!had_err && *gn)
    had_err = genome_node_accept(*gn, mfs->mergefeat_visitor, env);
  return had_err;
}

static void mergefeat_stream_unsorted_free(GenomeStream *gs, Env *env)
{
  MergefeatStreamUnsorted *mfs = mergefeat_stream_unsorted_cast(gs);
  genome_visitor_delete(mfs->mergefeat_visitor, env);
}

const GenomeStreamClass* mergefeat_stream_unsorted_class(void)
{
  static const GenomeStreamClass gsc = { sizeof (MergefeatStreamUnsorted),
                                         mergefeat_stream_unsorted_next_tree,
                                         mergefeat_stream_unsorted_free };
  return &gsc;
}

GenomeStream* mergefeat_stream_unsorted_new(GenomeStream *in_stream, Env *env)
{
  GenomeStream *gs = genome_stream_create(mergefeat_stream_unsorted_class(),
                                          false, env);
  MergefeatStreamUnsorted *mfs = mergefeat_stream_unsorted_cast(gs);
  mfs->in_stream = in_stream;
  mfs->mergefeat_visitor = mergefeat_visitor_new(env);
  return gs;
}
