/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <libgtext/genome_stream_rep.h>
#include <libgtext/gtf_out_stream.h>
#include <libgtext/gtf_visitor.h>

struct GTFOutStream {
  const GenomeStream parent_instance;
  GenomeStream *in_stream;
  GenomeVisitor *gtf_visitor;
};

#define gtf_out_stream_cast(GS)\
        genome_stream_cast(gtf_out_stream_class(), GS);

static int gtf_out_stream_next_tree(GenomeStream *gs, GenomeNode **gn,
                                     Env *env)
{
  GTFOutStream *gtf_out_stream;
  int had_err;
  env_error_check(env);
  gtf_out_stream = gtf_out_stream_cast(gs);
  had_err = genome_stream_next_tree(gtf_out_stream->in_stream, gn, env);
  if (!had_err && *gn)
    had_err = genome_node_accept(*gn, gtf_out_stream->gtf_visitor, env);
  return had_err;
}

static void gtf_out_stream_free(GenomeStream *gs, Env *env)
{
  GTFOutStream *gtf_out_stream = gtf_out_stream_cast(gs);
  genome_visitor_delete(gtf_out_stream->gtf_visitor, env);
}

const GenomeStreamClass* gtf_out_stream_class(void)
{
  static const GenomeStreamClass gsc = { sizeof (GTFOutStream),
                                         gtf_out_stream_next_tree,
                                         gtf_out_stream_free };
  return &gsc;
}

GenomeStream* gtf_out_stream_new(GenomeStream *in_stream, GenFile *outfp,
                                  Env *env)
{
  GenomeStream *gs = genome_stream_create(gtf_out_stream_class(),
                                          genome_stream_is_sorted(in_stream),
                                          env);
  GTFOutStream *gtf_out_stream = gtf_out_stream_cast(gs);
  gtf_out_stream->in_stream = in_stream;
  gtf_out_stream->gtf_visitor = gtf_visitor_new(outfp, env);
  return gs;
}
