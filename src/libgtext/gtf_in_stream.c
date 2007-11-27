/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg

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
#include "libgtcore/fa.h"
#include "libgtext/genome_stream_rep.h"
#include "libgtext/gtf_in_stream.h"
#include "libgtext/gtf_parser.h"

struct GTFInStream
{
  const GenomeStream parent_instance;
  Queue *genome_node_buffer;
};

#define gtf_in_stream_cast(GS)\
        genome_stream_cast(gtf_in_stream_class(), GS)

static int gtf_in_stream_next_tree(GenomeStream *gs, GenomeNode **gn, Env *env)
{
  GTFInStream *is;
  env_error_check(env);
  is = gtf_in_stream_cast(gs);
  if (queue_size(is->genome_node_buffer)) {
    /* we still have a node in the buffer -> serve it from there */
    *gn = queue_get(is->genome_node_buffer, env);
    return 0;
  }
  /* the buffer is empty */
  assert(!queue_size(is->genome_node_buffer));
  *gn = NULL;
  return 0;
}

static void gtf_in_stream_free(GenomeStream *gs, Env *env)
{
  GTFInStream *gtf_in_stream = gtf_in_stream_cast(gs);
  queue_delete(gtf_in_stream->genome_node_buffer, env);
}

const GenomeStreamClass* gtf_in_stream_class(void)
{
  static const GenomeStreamClass gsc = { sizeof (GTFInStream),
                                         gtf_in_stream_next_tree,
                                         gtf_in_stream_free };
  return &gsc;
}

GenomeStream* gtf_in_stream_new(const char *filename, bool be_tolerant,
                                Env *env)
{
  GenomeStream *gs;
  GTFInStream *gtf_in_stream;
  GTF_parser *gtf_parser;
  Str *filenamestr;
  int had_err;
  FILE *fpin;

  env_error_check(env);

  gs = genome_stream_create(gtf_in_stream_class(), false, env);
  gtf_in_stream = gtf_in_stream_cast(gs);
  gtf_parser = gtf_parser_new(env);

  gtf_in_stream->genome_node_buffer = queue_new(env);

  /* open input file */
  if (filename)
    fpin = fa_xfopen(filename, "r");
  else
    fpin = stdin;

  /* parse input file */
  filenamestr = str_new_cstr(filename ? filename : "stdin", env);
  had_err = gtf_parser_parse(gtf_parser, gtf_in_stream->genome_node_buffer,
                             filenamestr, fpin, be_tolerant, env);
  str_delete(filenamestr, env);

  /* close input file, if necessary */
  if (filename)
    fa_xfclose(fpin);

  /* free */
  gtf_parser_delete(gtf_parser, env);

  if (had_err) {
    genome_stream_delete(gs, env);
    return NULL;
  }
  return gs;
}
