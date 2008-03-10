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

#include <string.h>
#include "libgtcore/ma.h"
#include "libgtcore/range.h"
#include "libgtext/genome_stream_rep.h"
#include "libgtext/genome_feature.h"
#include "libgtext/genome_feature_type.h"
#include "libgtext/genome_node_iterator.h"
#include "libgtltr/ltrdigest_stream.h"
#include "libgtltr/ltrdigest_visitor.h"
#include "libgtltr/pbs.h"
#include "libgtltr/ppt.h"
#include "libgtltr/pdom.h"

struct LTRdigestStream {
  const GenomeStream parent_instance;
  GenomeStream *in_stream;
  LTRdigestVisitor *lv;
  LTRElement element;
};

void range_print(Range r)
{
  fprintf(stderr, "(%lu-%lu)\n", r.start, r.end);
}

#define ltrdigest_stream_cast(GS)\
        genome_stream_cast(ltrdigest_stream_class(), GS)

static int ltrdigest_stream_next_tree(GenomeStream *gs, GenomeNode **gn,
                                      Error *e)
{
  LTRdigestStream *ls;
  int had_err;
  GenomeNodeIterator* gni;

  error_check(e);
  ls = ltrdigest_stream_cast(gs);

  /* initialize this element */
  memset(&ls->element, 0, sizeof (LTRElement));

  /* get annotations from parser */
  had_err = genome_stream_next_tree(ls->in_stream, gn, e);
  if (!had_err && *gn)
  {
    GenomeNode *mygn;
    /* fill LTRElement structure from GFF3 subgraph */
    gni = genome_node_iterator_new(*gn);
    for(mygn = *gn; mygn; mygn = genome_node_iterator_next(gni))
      genome_node_accept(mygn, (GenomeVisitor*) ls->lv, e);
    genome_node_iterator_delete(gni);
  }

  if(ls->element.mainnode)
    fprintf(stderr,"%lu %lu %lu %lu %lu %lu\n", ls->element.leftLTR_5,
                                                ls->element.rightLTR_3,
                                                ls->element.leftLTR_5,
                                                ls->element.leftLTR_3,
                                                ls->element.rightLTR_5,
                                                ls->element.rightLTR_3);
  return had_err;
}

static void ltrdigest_stream_free(GenomeStream *gs)
{
  LTRdigestStream *ls = ltrdigest_stream_cast(gs);
  genome_visitor_delete((GenomeVisitor*) ls->lv);
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
    ls->lv = (LTRdigestVisitor*) ltrdigest_visitor_new(&ls->element);
  return gs;
}
