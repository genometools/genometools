/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include "filter_stream.h"
#include "filter_visitor.h"
#include "genome_feature.h"
#include "genome_stream_rep.h"
#include "undef.h"

struct FilterStream
{
  const GenomeStream parent_instance;
  GenomeStream *in_stream;
  GenomeVisitor *filter_visitor; /* the actual work is done in the visitor */
};

#define filter_stream_cast(GS)\
        genome_stream_cast(filter_stream_class(), GS);

static int filter_stream_next_tree(GenomeStream *gs, GenomeNode **gn, Env *env)
{
  FilterStream *fs;
  int has_err;
  env_error_check(env);
  fs = filter_stream_cast(gs);

  /* we still have nodes in the buffer */
  if (filter_visitor_node_buffer_size(fs->filter_visitor)) {
    *gn = filter_visitor_get_node(fs->filter_visitor); /* return one of them */
    return 0;
  }

  /* no nodes in the buffer -> get new nodes */
  while (!(has_err = genome_stream_next_tree(fs->in_stream, gn, env)) && *gn) {
    assert(*gn && !has_err);
    has_err = genome_node_accept(*gn, fs->filter_visitor, env);
    if (has_err)
      break;
    if (filter_visitor_node_buffer_size(fs->filter_visitor)) {
      *gn = filter_visitor_get_node(fs->filter_visitor);
      return 0;
    }
  }

  /* either we have an error or no new node */
  assert(has_err || !*gn);
  return has_err;
}

static void filter_stream_free(GenomeStream *gs, Env *env)
{
  FilterStream *fs = filter_stream_cast(gs);
  genome_visitor_delete(fs->filter_visitor, env);
}

const GenomeStreamClass* filter_stream_class(void)
{
  static const GenomeStreamClass gsc = { sizeof (FilterStream),
                                         filter_stream_next_tree,
                                         filter_stream_free };
  return &gsc;
}

GenomeStream* filter_stream_new(GenomeStream *in_stream,
                                Str *seqid, Str *typefilter,
                                unsigned long max_gene_length,
                                unsigned long max_gene_num,
                                double min_gene_score, Env *env)
{
  GenomeStream *gs = genome_stream_create(filter_stream_class(),
                                          genome_stream_is_sorted(in_stream),
                                          env);
  FilterStream *filter_stream = filter_stream_cast(gs);
  assert(in_stream);
  filter_stream->in_stream = in_stream;
  filter_stream->filter_visitor = filter_visitor_new(seqid, typefilter,
                                                     max_gene_length,
                                                     max_gene_num,
                                                     min_gene_score, env);
  return gs;
}
