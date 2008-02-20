/*
  Copyright (c) 2007-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007-2008 Center for Bioinformatics, University of Hamburg

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
#include "libgtext/stat_stream.h"
#include "libgtext/stat_visitor.h"
#include "libgtext/genome_stream_rep.h"

struct StatStream
{
  const GenomeStream parent_instance;
  GenomeStream *in_stream;
  GenomeVisitor *stat_visitor;
  unsigned long number_of_trees;
};

#define stat_stream_cast(GS)\
        genome_stream_cast(stat_stream_class(), GS)

static int stat_stream_next_tree(GenomeStream *gs, GenomeNode **gn, Error *e)
{
  StatStream *stat_stream;
  int had_err;
  error_check(e);
  stat_stream = stat_stream_cast(gs);
  had_err = genome_stream_next_tree(stat_stream->in_stream, gn, e);
  if (!had_err) {
    assert(stat_stream->stat_visitor);
    if (*gn) {
      stat_stream->number_of_trees++;
      had_err = genome_node_accept(*gn, stat_stream->stat_visitor, e);
      assert(!had_err); /* the status visitor is sane */
    }
  }
  return had_err;
}

static void stat_stream_free(GenomeStream *gs)
{
  StatStream *stat_stream = stat_stream_cast(gs);
  genome_visitor_delete(stat_stream->stat_visitor);
  genome_stream_delete(stat_stream->in_stream);
}

const GenomeStreamClass* stat_stream_class(void)
{
  static const GenomeStreamClass gsc = { sizeof (StatStream),
                                         stat_stream_next_tree,
                                         stat_stream_free };
  return &gsc;
}

GenomeStream* stat_stream_new(GenomeStream *in_stream,
                              bool gene_length_distri,
                              bool gene_score_distri,
                              bool exon_length_distri,
                              bool exon_number_distri,
                              bool intron_length_distri)
{
  GenomeStream *gs = genome_stream_create(stat_stream_class(), false);
  StatStream *ss = stat_stream_cast(gs);
  ss->in_stream = genome_stream_ref(in_stream);
  ss->stat_visitor = stat_visitor_new(gene_length_distri, gene_score_distri,
                                      exon_length_distri, exon_number_distri,
                                      intron_length_distri);
  return gs;
}

void stat_stream_show_stats(GenomeStream *gs)
{
  StatStream *ss = stat_stream_cast(gs);
  printf("parsed feature trees: %lu\n", ss->number_of_trees);
  stat_visitor_show_stats(ss->stat_visitor);
}
