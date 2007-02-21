/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "cds_stream.h"
#include "cds_visitor.h"
#include "genome_stream_rep.h"
#include "str.h"

struct CDSStream
{
  const GenomeStream parent_instance;
  GenomeStream *in_stream;
  GenomeVisitor *cds_visitor;
};

#define cds_stream_cast(GS)\
        genome_stream_cast(cds_stream_class(), GS)

static int cds_stream_next_tree(GenomeStream *gs, GenomeNode **gn, Log *l,
                                Error *err)
{
  CDSStream *cds_stream;
  int has_err;
  error_check(err);
  cds_stream = cds_stream_cast(gs);
  has_err = genome_stream_next_tree(cds_stream->in_stream, gn, l, err);
  if (!has_err && *gn)
    has_err = genome_node_accept(*gn, cds_stream->cds_visitor, l, err);
  return has_err;
}

static void cds_stream_free(GenomeStream *gs)
{
  CDSStream *cds_stream = cds_stream_cast(gs);
  genome_visitor_delete(cds_stream->cds_visitor);
}

const GenomeStreamClass* cds_stream_class(void)
{
  static const GenomeStreamClass gsc = { sizeof (CDSStream),
                                         cds_stream_next_tree,
                                         cds_stream_free };
  return &gsc;
}

GenomeStream* cds_stream_new(GenomeStream *in_stream, const char *sequence_file,
                             const char *source, Error *err)
{
  GenomeStream *gs;
  CDSStream *cds_stream;
  Str *sequence_file_str, *source_str;
  int has_err = 0;
  error_check(err);
  gs = genome_stream_create(cds_stream_class(), true);
  cds_stream = cds_stream_cast(gs);
  sequence_file_str = str_new_cstr(sequence_file),
  source_str = str_new_cstr(source);
  cds_stream->in_stream = in_stream;
  cds_stream->cds_visitor = cds_visitor_new(sequence_file_str, source_str, err);
  if (!cds_stream->cds_visitor)
    has_err = -1;
  str_delete(sequence_file_str);
  str_delete(source_str);
  if (has_err) {
    cds_stream_free(gs);
    return NULL;
  }
  return gs;
}
