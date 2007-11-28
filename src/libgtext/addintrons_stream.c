/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

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
#include "libgtext/genome_stream_rep.h"
#include "libgtext/addintrons_stream.h"
#include "libgtext/addintrons_visitor.h"

struct AddIntronsStream{
  const GenomeStream parent_instance;
  GenomeStream *in_stream;
  GenomeVisitor *addintrons_visitor;
};

#define addintrons_stream_cast(GS)\
        genome_stream_cast(addintrons_stream_class(), GS)

static int addintrons_stream_next_tree(GenomeStream *gs, GenomeNode **gn,
                                       Error *e)
{
  AddIntronsStream *ais;
  int had_err;
  error_check(e);
  ais = addintrons_stream_cast(gs);
  had_err = genome_stream_next_tree(ais->in_stream, gn, e);
  if (!had_err && *gn)
    had_err = genome_node_accept(*gn, ais->addintrons_visitor, e);
  return had_err;
}

static void addintrons_stream_free(GenomeStream *gs)
{
  AddIntronsStream *ais = addintrons_stream_cast(gs);
  genome_visitor_delete(ais->addintrons_visitor);
}

const GenomeStreamClass* addintrons_stream_class(void)
{
  static const GenomeStreamClass gsc = { sizeof (AddIntronsStream),
                                         addintrons_stream_next_tree,
                                         addintrons_stream_free };
  return &gsc;
}

GenomeStream* addintrons_stream_new(GenomeStream *in_stream)
{
  GenomeStream *gs = genome_stream_create(addintrons_stream_class(), true);
  AddIntronsStream *ais = addintrons_stream_cast(gs);
  assert(in_stream);
  ais->in_stream = in_stream;
  ais->addintrons_visitor = addintrons_visitor_new();
  return gs;
}
