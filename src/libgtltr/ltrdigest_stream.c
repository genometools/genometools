/*
  Copyright (c) 2008 Sascha Steinbiss <ssteinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2008 Center for Bioinformatics, University of Hamburg

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

#include "libgtext/genome_stream_rep.h"
#include "libgtltr/ltrdigest_stream.h"

struct LTRdigestStream {
  const GenomeStream parent_instance;
  GenomeStream *in_stream;
};

#define ltrdigest_stream_cast(GS)\
        genome_stream_cast(ltrdigest_stream_class(), GS)

static int ltrdigest_stream_next_tree(GenomeStream *gs, GenomeNode **gn, Error *e)
{
  LTRdigestStream *ls;
  int had_err;
  error_check(e);
  ls = ltrdigest_stream_cast(gs);
  had_err = genome_stream_next_tree(ls->in_stream, gn, e);


  return had_err;
}

static void ltrdigest_stream_free(GenomeStream *gs)
{
  LTRdigestStream *ls = ltrdigest_stream_cast(gs);
  genome_stream_delete(ls->in_stream);
}

const GenomeStreamClass* ltrdigest_stream_class(void)
{
  static const GenomeStreamClass gsc = { sizeof (LTRdigestStream),
                                         ltrdigest_stream_next_tree,
                                         ltrdigest_stream_free };
  return &gsc;
}

GenomeStream* ltrdigest_stream_new(GenomeStream *in_stream)
{
  GenomeStream *gs;
  LTRdigestStream *ls;
  gs = genome_stream_create(ltrdigest_stream_class(), true);
  ls = ltrdigest_stream_cast(gs);
  ls->in_stream = genome_stream_ref(in_stream);
  return gs;
}
