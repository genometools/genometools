/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <gtcore.h>
#include <libgtext/feature_stream.h>
#include <libgtext/feature_visitor.h>
#include <libgtext/feature_index.h>
#include <libgtext/genome_stream_rep.h>

struct FeatureStream
{
  const GenomeStream parent_instance;
  GenomeStream *in_stream;
  GenomeVisitor *feature_visitor;
};

#define feature_stream_cast(GS)\
        genome_stream_cast(feature_stream_class(), GS)

static int feature_stream_next_tree(GenomeStream *gs, GenomeNode **gn, Env *env)
{
  FeatureStream *feature_stream;
  int has_err;
  env_error_check(env);
  feature_stream = feature_stream_cast(gs);
  has_err = genome_stream_next_tree(feature_stream->in_stream, gn, env);
  if (!has_err && *gn)
    has_err = genome_node_accept(*gn, feature_stream->feature_visitor, env);
  return has_err;
}

static void feature_stream_free(GenomeStream *gs, Env *env)
{
  FeatureStream *feature_stream = feature_stream_cast(gs);
  genome_visitor_delete(feature_stream->feature_visitor, env);
}

const GenomeStreamClass* feature_stream_class(void)
{
  static const GenomeStreamClass gsc = { sizeof (FeatureStream),
                                         feature_stream_next_tree,
                                         feature_stream_free };
  return &gsc;
}

GenomeStream* feature_stream_new(GenomeStream *in_stream,
                             FeatureIndex *fi, Env *env)
{
  GenomeStream *gs;
  FeatureStream *feature_stream;
  int has_err = 0;
  env_error_check(env);
  gs = genome_stream_create(feature_stream_class(), true, env);
  feature_stream = feature_stream_cast(gs);
  feature_stream->in_stream = in_stream;
	
  feature_stream->feature_visitor = feature_visitor_new(fi, env);
  
	if (!feature_stream->feature_visitor)
    has_err = -1;
  if (has_err) {
    feature_stream_free(gs, env);
    return NULL;
  }
  return gs;
}
