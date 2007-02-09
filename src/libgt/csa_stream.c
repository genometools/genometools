/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include "csa_stream.h"
#include "csa_visitor.h"
#include "consensus_sa.h"
#include "genome_stream_rep.h"
#include "queue.h"

struct Csa_stream {
  const Genome_stream parent_instance;
  Genome_stream *in_stream;
  Genome_visitor *csa_visitor; /* the actual work is done in the visitor */
};

#define csa_stream_cast(GS)\
        genome_stream_cast(csa_stream_class(), GS)

Genome_node* csa_stream_next_tree(Genome_stream *gs, Log *l)
{
  Csa_stream *cs = csa_stream_cast(gs);
  Genome_node *gn;

  if (csa_visitor_node_buffer_size(cs->csa_visitor))
    return csa_visitor_get_node(cs->csa_visitor);

  while ((gn = genome_stream_next_tree(cs->in_stream, l))) {
    assert(gn);
    genome_node_accept(gn, cs->csa_visitor, l);
    if (csa_visitor_node_buffer_size(cs->csa_visitor))
      return csa_visitor_get_node(cs->csa_visitor);
  }
  assert(!gn);
  csa_visitor_process_cluster(cs->csa_visitor, 1, l);
  if (csa_visitor_node_buffer_size(cs->csa_visitor))
    return csa_visitor_get_node(cs->csa_visitor);
  return NULL;
}

static void csa_stream_free(Genome_stream *gs)
{
  Csa_stream *cs = csa_stream_cast(gs);
  genome_visitor_free(cs->csa_visitor);
}

const Genome_stream_class* csa_stream_class(void)
{
  static const Genome_stream_class gsc = { sizeof(Csa_stream),
                                           csa_stream_next_tree,
                                           csa_stream_free };
  return &gsc;
}

Genome_stream* csa_stream_new(Genome_stream *in_stream,
                              unsigned long join_length)
{
  Genome_stream *gs = genome_stream_create(csa_stream_class(),
                                           genome_stream_is_sorted(in_stream));
  Csa_stream *cs = csa_stream_cast(gs);
  cs->in_stream = in_stream;
  cs->csa_visitor = csa_visitor_new(join_length);
  return gs;
}
