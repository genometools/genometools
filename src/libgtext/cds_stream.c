/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <libgtext/cds_stream.h>
#include <libgtext/cds_visitor.h>
#include <libgtext/genome_stream_rep.h>

struct CDSStream
{
  const GenomeStream parent_instance;
  GenomeStream *in_stream;
  GenomeVisitor *cds_visitor;
};

#define cds_stream_cast(GS)\
        genome_stream_cast(cds_stream_class(), GS)

static int cds_stream_next_tree(GenomeStream *gs, GenomeNode **gn, Env *env)
{
  CDSStream *cds_stream;
  int had_err;
  env_error_check(env);
  cds_stream = cds_stream_cast(gs);
  had_err = genome_stream_next_tree(cds_stream->in_stream, gn, env);
  if (!had_err && *gn)
    had_err = genome_node_accept(*gn, cds_stream->cds_visitor, env);
  return had_err;
}

static void cds_stream_free(GenomeStream *gs, Env *env)
{
  CDSStream *cds_stream = cds_stream_cast(gs);
  genome_visitor_delete(cds_stream->cds_visitor, env);
}

const GenomeStreamClass* cds_stream_class(void)
{
  static const GenomeStreamClass gsc = { sizeof (CDSStream),
                                         cds_stream_next_tree,
                                         cds_stream_free };
  return &gsc;
}

GenomeStream* cds_stream_new(GenomeStream *in_stream, RegionMapping *rm,
                             const char *source, Env *env)
{
  GenomeStream *gs;
  CDSStream *cds_stream;
  Str *source_str;
  int had_err = 0;
  env_error_check(env);
  gs = genome_stream_create(cds_stream_class(), true, env);
  cds_stream = cds_stream_cast(gs);
  source_str = str_new_cstr(source, env);
  cds_stream->in_stream = in_stream;
  cds_stream->cds_visitor = cds_visitor_new(rm, source_str, env);
  if (!cds_stream->cds_visitor)
    had_err = -1;
  str_delete(source_str, env);
  if (had_err) {
    cds_stream_free(gs, env);
    return NULL;
  }
  return gs;
}
