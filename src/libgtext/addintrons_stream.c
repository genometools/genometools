/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include "libgtext/genome_stream_rep.h"
#include "libgtext/addintrons_stream.h"
#include "libgtext/addintrons_visitor.h"

struct AddIntronsStream{
  const GenomeStream parent_instance;
  GenomeStream *in_stream;
  GenomeVisitor *addintrons_visitor;
};

#define addintrons_stream_cast(GS)\
        genome_stream_cast(addintrons_stream_class(), GS)

static int addintrons_stream_next_tree(GenomeStream *gs, GenomeNode **gn,
                                             Env *env)
{
  AddIntronsStream *ais;
  int had_err;
  env_error_check(env);
  ais = addintrons_stream_cast(gs);
  had_err = genome_stream_next_tree(ais->in_stream, gn, env);
  if (!had_err && *gn)
    had_err = genome_node_accept(*gn, ais->addintrons_visitor, env);
  return had_err;
}

static void addintrons_stream_free(GenomeStream *gs, Env *env)
{
  AddIntronsStream *ais = addintrons_stream_cast(gs);
  genome_visitor_delete(ais->addintrons_visitor, env);
}

const GenomeStreamClass* addintrons_stream_class(void)
{
  static const GenomeStreamClass gsc = { sizeof (AddIntronsStream),
                                         addintrons_stream_next_tree,
                                         addintrons_stream_free };
  return &gsc;
}

GenomeStream* addintrons_stream_new(GenomeStream *in_stream, Env *env)
{
  GenomeStream *gs = genome_stream_create(addintrons_stream_class(), true, env);
  AddIntronsStream *ais = addintrons_stream_cast(gs);
  assert(in_stream);
  ais->in_stream = in_stream;
  ais->addintrons_visitor = addintrons_visitor_new(env);
  return gs;
}
