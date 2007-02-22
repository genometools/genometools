/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include "extractfeat_stream.h"
#include "extractfeat_visitor.h"
#include "genome_stream_rep.h"

struct ExtractFeatStream
{
  const GenomeStream parent_instance;
  GenomeStream *in_stream;
  GenomeVisitor *extractfeat_visitor;
  GenomeFeatureType type;
  bool join,
       translate;
};

#define extractfeat_stream_cast(GS)\
        genome_stream_cast(extractfeat_stream_class(), GS)

static int extractfeat_stream_next_tree(GenomeStream *gs, GenomeNode **gn,
                                        Env *env)
{
  ExtractFeatStream *extractfeat_stream;
  int has_err;
  env_error_check(env);
  extractfeat_stream = extractfeat_stream_cast(gs);
  has_err = genome_stream_next_tree(extractfeat_stream->in_stream, gn, env);
  if (!has_err) {
  assert(extractfeat_stream->extractfeat_visitor);
  if (*gn)
    has_err = genome_node_accept(*gn, extractfeat_stream->extractfeat_visitor,
                                 env);
  }
  return has_err;
}

static void extractfeat_stream_free(GenomeStream *gs, Env *env)
{
  ExtractFeatStream *extractfeat_stream = extractfeat_stream_cast(gs);
  genome_visitor_delete(extractfeat_stream->extractfeat_visitor, env);
}

const GenomeStreamClass* extractfeat_stream_class(void)
{
  static const GenomeStreamClass gsc = { sizeof (ExtractFeatStream),
                                         extractfeat_stream_next_tree,
                                         extractfeat_stream_free };
  return &gsc;
}

GenomeStream* extractfeat_stream_new(GenomeStream *in_stream,
                                     GenomeFeatureType type,
                                     bool join, bool translate)
{
  GenomeStream *gs = genome_stream_create(extractfeat_stream_class(), true);
  ExtractFeatStream *efs = extractfeat_stream_cast(gs);
  efs->in_stream = in_stream;
  efs->type = type;
  efs->join = join;
  efs->translate = translate;
  return gs;
}

int extractfeat_stream_use_sequence_file(GenomeStream *gs, Str *seqfile,
                                         Env *env)
{
  ExtractFeatStream *efs;
  env_error_check(env);
  efs = extractfeat_stream_cast(gs);
  efs->extractfeat_visitor = extractfeat_visitor_new_seqfile(seqfile, efs->type,
                                                             efs->join,
                                                             efs->translate,
                                                             env);
  if (!efs->extractfeat_visitor)
    return -1;
  return 0;
}

void extractfeat_stream_use_region_mapping(GenomeStream *gs, RegionMapping *rm,
                                           Env *env)
{
  ExtractFeatStream *efs = extractfeat_stream_cast(gs);
  efs->extractfeat_visitor = extractfeat_visitor_new_regionmapping(rm,
                                                                   efs->type,
                                                                   efs->join,
                                                                   efs
                                                                   ->translate,
                                                                   env);
}
