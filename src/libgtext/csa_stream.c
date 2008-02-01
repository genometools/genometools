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
#include "libgtext/csa_stream.h"
#include "libgtext/csa_visitor.h"
#include "libgtext/consensus_sa.h"
#include "libgtext/genome_stream_rep.h"

struct CSAStream {
  const GenomeStream parent_instance;
  GenomeStream *in_stream;
  GenomeVisitor *csa_visitor; /* the actual work is done in the visitor */
};

#define csa_stream_cast(GS)\
        genome_stream_cast(csa_stream_class(), GS)

int csa_stream_next_tree(GenomeStream *gs, GenomeNode **gn, Error *e)
{
  CSAStream *cs;
  int had_err;
  error_check(e);
  cs = csa_stream_cast(gs);

  /* we have still nodes in the buffer */
  if (csa_visitor_node_buffer_size(cs->csa_visitor)) {
    *gn = csa_visitor_get_node(cs->csa_visitor); /* return one of them */
    return 0;
  }

  /* no nodes in the buffer -> get new nodes */
  while (!(had_err = genome_stream_next_tree(cs->in_stream, gn, e)) && *gn) {
    assert(*gn && !had_err);
    had_err = genome_node_accept(*gn, cs->csa_visitor, e);
    if (had_err)
      break;
    if (csa_visitor_node_buffer_size(cs->csa_visitor)) {
      *gn = csa_visitor_get_node(cs->csa_visitor);
      return 0;
    }
  }

  /* either we have an error or no new node */
  assert(had_err || !*gn);

  /* if we have no error, process the last cluster */
  if (!had_err) {
    csa_visitor_process_cluster(cs->csa_visitor, true);
    if (csa_visitor_node_buffer_size(cs->csa_visitor)) {
      *gn = csa_visitor_get_node(cs->csa_visitor);
      return 0;
    }
  }
  return had_err;
}

static void csa_stream_free(GenomeStream *gs)
{
  CSAStream *cs = csa_stream_cast(gs);
  genome_visitor_delete(cs->csa_visitor);
  genome_stream_delete(cs->in_stream);
}

const GenomeStreamClass* csa_stream_class(void)
{
  static const GenomeStreamClass gsc = { sizeof (CSAStream),
                                         csa_stream_next_tree,
                                         csa_stream_free };
  return &gsc;
}

GenomeStream* csa_stream_new(GenomeStream *in_stream, unsigned long join_length)
{
  GenomeStream *gs = genome_stream_create(csa_stream_class(),
                                          genome_stream_is_sorted(in_stream));
  CSAStream *cs = csa_stream_cast(gs);
  cs->in_stream = genome_stream_ref(in_stream);
  cs->csa_visitor = csa_visitor_new(join_length);
  return gs;
}
