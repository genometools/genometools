/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include "genome_stream_rep.h"
#include "gtf_in_stream.h"
#include "gtf_parser.h"
#include "queue.h"
#include "xansi.h"

struct GTFInStream
{
  const GenomeStream parent_instance;
  Queue *genome_node_buffer;
};

#define gtf_in_stream_cast(GS)\
        genome_stream_cast(gtf_in_stream_class(), GS)

static int gtf_in_stream_next_tree(GenomeStream *gs, GenomeNode **gn,
                                   /*@unused@*/ Log *l, Error *err)
{
  GTFInStream *is;
  error_check(err);
  is = gtf_in_stream_cast(gs);
  if (queue_size(is->genome_node_buffer)) {
    /* we still have a node in the buffer -> serve it from there */
    *gn = *(GenomeNode**) queue_get(is->genome_node_buffer);
    return 0;
  }
  /* the buffer is empty */
  assert(!queue_size(is->genome_node_buffer));
  *gn = NULL;
  return 0;
}

static void gtf_in_stream_free(GenomeStream *gs)
{
  GTFInStream *gtf_in_stream = gtf_in_stream_cast(gs);
  queue_free(gtf_in_stream->genome_node_buffer);
}

const GenomeStreamClass* gtf_in_stream_class(void)
{
  static const GenomeStreamClass gsc = { sizeof (GTFInStream),
                                         gtf_in_stream_next_tree,
                                         gtf_in_stream_free };
  return &gsc;
}

GenomeStream* gtf_in_stream_new(const char *filename, bool be_tolerant,
                                Error *err)
{
  GenomeStream *gs;
  GTFInStream *gtf_in_stream;
  GTF_parser *gtf_parser;
  int has_err;
  FILE *fpin;

  error_check(err);

  gs = genome_stream_create(gtf_in_stream_class(), false);
  gtf_in_stream = gtf_in_stream_cast(gs);
  gtf_parser = gtf_parser_new();

  gtf_in_stream->genome_node_buffer = queue_new(sizeof (GenomeNode*));

  /* open input file */
  if (filename)
    fpin = xfopen(filename, "r");
  else
    fpin = stdin;

  /* parse input file */
  has_err = gtf_parser_parse(gtf_parser, gtf_in_stream->genome_node_buffer,
                             filename ? filename : "stdin", fpin, be_tolerant,
                             err);

  /* close input file, if necessary */
  if (filename)
    xfclose(fpin);

  /* free */
  gtf_parser_free(gtf_parser);

  if (has_err) {
    genome_stream_free(gs);
    return NULL;
  }
  return gs;
}
