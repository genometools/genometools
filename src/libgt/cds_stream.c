/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "cds_stream.h"
#include "cds_visitor.h"
#include "genome_stream_rep.h"
#include "str.h"

struct CDS_stream
{
  const GenomeStream parent_instance;
  GenomeStream *in_stream;
  GenomeVisitor *cds_visitor;
};

#define cds_stream_cast(GS)\
        genome_stream_cast(cds_stream_class(), GS)

static GenomeNode* cds_stream_next_tree(GenomeStream *gs, Log *l)
{
  CDS_stream *cds_stream = cds_stream_cast(gs);
  GenomeNode *gn = genome_stream_next_tree(cds_stream->in_stream, l);
  if (gn)
    genome_node_accept(gn, cds_stream->cds_visitor, l);
  return gn;
}

static void cds_stream_free(GenomeStream *gs)
{
  CDS_stream *cds_stream = cds_stream_cast(gs);
  genome_visitor_free(cds_stream->cds_visitor);
}

const GenomeStreamClass* cds_stream_class(void)
{
  static const GenomeStreamClass gsc = { sizeof(CDS_stream),
                                           cds_stream_next_tree,
                                           cds_stream_free };
  return &gsc;
}

GenomeStream* cds_stream_new(GenomeStream *in_stream,
                              const char *sequence_file, const char *source)
{
  GenomeStream *gs = genome_stream_create(cds_stream_class(), 1);
  CDS_stream *cds_stream = cds_stream_cast(gs);
  Str *sequence_file_str = str_new_cstr(sequence_file),
      *source_str = str_new_cstr(source);

  cds_stream->in_stream = in_stream;
  cds_stream->cds_visitor = cds_visitor_new(sequence_file_str, source_str);

  str_free(sequence_file_str);
  str_free(source_str);

  return gs;
}
