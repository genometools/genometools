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
#include "libgtext/genome_stream_rep.h"
#include "libgtext/add_introns_stream.h"
#include "libgtext/add_introns_visitor.h"

struct AddIntronsStream{
  const GenomeStream parent_instance;
  GenomeStream *in_stream;
  GenomeVisitor *add_introns_visitor;
};

#define add_introns_stream_cast(GS)\
        genome_stream_cast(add_introns_stream_class(), GS)

static int add_introns_stream_next_tree(GenomeStream *gs, GenomeNode **gn,
                                        Error *err)
{
  AddIntronsStream *ais;
  int had_err;
  error_check(err);
  ais = add_introns_stream_cast(gs);
  had_err = genome_stream_next_tree(ais->in_stream, gn, err);
  if (!had_err && *gn)
    had_err = genome_node_accept(*gn, ais->add_introns_visitor, err);
  return had_err;
}

static void add_introns_stream_free(GenomeStream *gs)
{
  AddIntronsStream *ais = add_introns_stream_cast(gs);
  genome_visitor_delete(ais->add_introns_visitor);
  genome_stream_delete(ais->in_stream);
}

const GenomeStreamClass* add_introns_stream_class(void)
{
  static const GenomeStreamClass gsc = { sizeof (AddIntronsStream),
                                         add_introns_stream_next_tree,
                                         add_introns_stream_free };
  return &gsc;
}

GenomeStream* add_introns_stream_new(GenomeStream *in_stream)
{
  GenomeStream *gs = genome_stream_create(add_introns_stream_class(), true);
  AddIntronsStream *ais = add_introns_stream_cast(gs);
  assert(in_stream);
  ais->in_stream = genome_stream_ref(in_stream);
  ais->add_introns_visitor = add_introns_visitor_new();
  return gs;
}
