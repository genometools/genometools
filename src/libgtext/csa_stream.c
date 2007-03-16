/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include <libgtext/csa_stream.h>
#include <libgtext/csa_visitor.h>
#include <libgtext/consensus_sa.h>
#include <libgtext/genome_stream_rep.h>

struct CSAStream {
  const GenomeStream parent_instance;
  GenomeStream *in_stream;
  GenomeVisitor *csa_visitor; /* the actual work is done in the visitor */
};

#define csa_stream_cast(GS)\
        genome_stream_cast(csa_stream_class(), GS)

int csa_stream_next_tree(GenomeStream *gs, GenomeNode **gn, Env *env)
{
  CSAStream *cs;
  int has_err;
  env_error_check(env);
  cs = csa_stream_cast(gs);

  /* we have still nodes in the buffer */
  if (csa_visitor_node_buffer_size(cs->csa_visitor)) {
    *gn = csa_visitor_get_node(cs->csa_visitor); /* return one of them */
    return 0;
  }

  /* no nodes in the buffer -> get new nodes */
  while (!(has_err = genome_stream_next_tree(cs->in_stream, gn, env)) && *gn) {
    assert(*gn && !has_err);
    has_err = genome_node_accept(*gn, cs->csa_visitor, env);
    if (has_err)
      break;
    if (csa_visitor_node_buffer_size(cs->csa_visitor)) {
      *gn = csa_visitor_get_node(cs->csa_visitor);
      return 0;
    }
  }

  /* either we have an error or no new node */
  assert(has_err || !*gn);

  /* if we have no error, process the last cluster */
  if (!has_err) {
    csa_visitor_process_cluster(cs->csa_visitor, true, env);
    if (csa_visitor_node_buffer_size(cs->csa_visitor)) {
      *gn = csa_visitor_get_node(cs->csa_visitor);
      return 0;
    }
  }
  return has_err;
}

static void csa_stream_free(GenomeStream *gs, Env *env)
{
  CSAStream *cs = csa_stream_cast(gs);
  genome_visitor_delete(cs->csa_visitor, env);
}

const GenomeStreamClass* csa_stream_class(void)
{
  static const GenomeStreamClass gsc = { sizeof (CSAStream),
                                         csa_stream_next_tree,
                                         csa_stream_free };
  return &gsc;
}

GenomeStream* csa_stream_new(GenomeStream *in_stream, unsigned long join_length,
                             Env *env)
{
  GenomeStream *gs = genome_stream_create(csa_stream_class(),
                                          genome_stream_is_sorted(in_stream),
                                          env);
  CSAStream *cs = csa_stream_cast(gs);
  cs->in_stream = in_stream;
  cs->csa_visitor = csa_visitor_new(join_length, env);
  return gs;
}
