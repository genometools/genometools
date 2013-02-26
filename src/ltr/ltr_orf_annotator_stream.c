/*
  Copyright (c) 2012 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2012 Center for Bioinformatics, University of Hamburg

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

#include "core/class_alloc_lock.h"
#include "core/hashmap_api.h"
#include "extended/feature_type.h"
#include "extended/orf_finder_stream.h"
#include "ltr/ltr_orf_annotator_stream_api.h"

struct GtLTRORFAnnotatorStream {
  const GtNodeStream parent_instance;
  GtNodeStream *orf_stream;
  unsigned long *progress_loc;
  GtHashmap *types;
};

const GtNodeStreamClass* gt_ltr_orf_annotator_stream_class(void);

#define ltr_orf_annotator_stream_cast(NS)\
        gt_node_stream_cast(gt_ltr_orf_annotator_stream_class(), NS)

static int ltr_orf_annotator_stream_next(GtNodeStream *ns, GtGenomeNode **gn,
                                         GtError *err)
{
  GtLTRORFAnnotatorStream *bs;
  gt_error_check(err);
  bs = ltr_orf_annotator_stream_cast(ns);
  if (bs->progress_loc != NULL)
    (*bs->progress_loc)++;
  return gt_node_stream_next(bs->orf_stream, gn, err);
}

static void ltr_orf_annotator_stream_free(GtNodeStream *ns)
{
  GtLTRORFAnnotatorStream *bs = ltr_orf_annotator_stream_cast(ns);
  gt_hashmap_delete(bs->types);
  gt_node_stream_delete(bs->orf_stream);
}

const GtNodeStreamClass* gt_ltr_orf_annotator_stream_class(void)
{
  static const GtNodeStreamClass *nsc = NULL;
  gt_class_alloc_lock_enter();
  if (!nsc) {
    nsc = gt_node_stream_class_new(sizeof (GtLTRORFAnnotatorStream),
                                   ltr_orf_annotator_stream_free,
                                   ltr_orf_annotator_stream_next);
  }
  gt_class_alloc_lock_leave();
  return nsc;
}

GtNodeStream* gt_ltr_orf_annotator_stream_new(GtNodeStream *in_stream,
                                              GtEncseq *encseq,
                                              unsigned int min,
                                              unsigned int max,
                                              bool all,
                                              GtError *err)
{
  GtLTRORFAnnotatorStream *bs;
  GtNodeStream *ns;
  gt_assert(in_stream);
  ns = gt_node_stream_create(gt_ltr_orf_annotator_stream_class(), false);
  bs = ltr_orf_annotator_stream_cast(ns);
  bs->types = gt_hashmap_new(GT_HASH_STRING, NULL, NULL);
  bs->progress_loc = NULL;
  gt_hashmap_add(bs->types, gt_ft_LTR_retrotransposon, (void*) 1);
  bs->orf_stream = gt_orf_finder_stream_new(in_stream, encseq, bs->types, min,
                                            max, all, err);
  return ns;
}

void gt_ltr_orf_annotator_stream_set_progress_location(
                                                     GtLTRORFAnnotatorStream *s,
                                                     unsigned long *loc)
{
  gt_assert(s);
  s->progress_loc = loc;
}
