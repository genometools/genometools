/*
  Copyright (c) 2006-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2008 Center for Bioinformatics, University of Hamburg

  Permission to use, copy, modify, and distribute this software for any
  purpose with or without fee is hereby granted, provided that the above
  copyright notice and this permission notice appear in all copies.

  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/

#include <assert.h>
#include "core/fa.h"
#include "core/unused_api.h"
#include "extended/genome_stream_rep.h"
#include "extended/gtf_in_stream.h"
#include "extended/gtf_parser.h"
#include "extended/type_factory_builtin.h"

struct GTFInStream
{
  const GenomeStream parent_instance;
  GT_Queue *gt_genome_node_buffer;
  GT_TypeFactory *feature_type_factory;
};

#define gtf_in_stream_cast(GS)\
        genome_stream_cast(gtf_in_stream_class(), GS)

static int gtf_in_stream_next_tree(GenomeStream *gs, GT_GenomeNode **gn,
                                   GT_UNUSED GT_Error *err)
{
  GTFInStream *is;
  gt_error_check(err);
  is = gtf_in_stream_cast(gs);
  if (gt_queue_size(is->gt_genome_node_buffer)) {
    /* we still have a node in the buffer -> serve it from there */
    *gn = gt_queue_get(is->gt_genome_node_buffer);
    return 0;
  }
  /* the buffer is empty */
  assert(!gt_queue_size(is->gt_genome_node_buffer));
  *gn = NULL;
  return 0;
}

static void gtf_in_stream_free(GenomeStream *gs)
{
  GTFInStream *gtf_in_stream = gtf_in_stream_cast(gs);
  gt_type_factory_delete(gtf_in_stream->feature_type_factory);
  gt_queue_delete(gtf_in_stream->gt_genome_node_buffer);
}

const GenomeStreamClass* gtf_in_stream_class(void)
{
  static const GenomeStreamClass gsc = { sizeof (GTFInStream),
                                         gtf_in_stream_next_tree,
                                         gtf_in_stream_free };
  return &gsc;
}

GenomeStream* gtf_in_stream_new(const char *filename, bool be_tolerant,
                                GT_Error *err)
{
  GenomeStream *gs;
  GTFInStream *gtf_in_stream;
  GTF_parser *gtf_parser;
  GT_Str *filenamestr;
  int had_err;
  FILE *fpin;

  gt_error_check(err);

  gs = genome_stream_create(gtf_in_stream_class(), false);
  gtf_in_stream = gtf_in_stream_cast(gs);
  gtf_in_stream->gt_genome_node_buffer = gt_queue_new();
  gtf_in_stream->feature_type_factory = gt_type_factory_builtin_new();

  gtf_parser = gtf_parser_new(gtf_in_stream->feature_type_factory);

  /* open input file */
  if (filename)
    fpin = gt_xfopen(filename, "r");
  else
    fpin = stdin;

  /* parse input file */
  filenamestr = gt_str_new_cstr(filename ? filename : "stdin");
  had_err = gtf_parser_parse(gtf_parser, gtf_in_stream->gt_genome_node_buffer,
                             filenamestr, fpin, be_tolerant, err);
  gt_str_delete(filenamestr);

  /* close input file, if necessary */
  if (filename)
    gt_xfclose(fpin);

  /* free */
  gtf_parser_delete(gtf_parser);

  if (had_err) {
    genome_stream_delete(gs);
    return NULL;
  }
  return gs;
}
