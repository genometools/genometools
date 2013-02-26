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

#include "core/assert_api.h"
#include "core/class_alloc_lock.h"
#include "core/undef_api.h"
#include "extended/feature_node.h"
#include "extended/feature_node_iterator_api.h"
#include "extended/node_stream_api.h"
#include "extended/uniq_stream.h"

struct GtUniqStream{
  const GtNodeStream parent_instance;
  GtNodeStream *in_stream;
  GtGenomeNode *first_node,
               *second_node;
};

#define uniq_stream_cast(GS)\
        gt_node_stream_cast(gt_uniq_stream_class(), GS)

static bool nodes_are_equal_feature_trees(GtGenomeNode *first_node,
                                          GtGenomeNode *second_node)
{
  bool equal = false;
  GtFeatureNodeIterator *fni_a, *fni_b;
  GtFeatureNode *fn_a, *fn_b;
  fn_a = gt_feature_node_try_cast(first_node);
  fn_b = gt_feature_node_try_cast(second_node);
  if (fn_a && fn_b) {
    fni_a = gt_feature_node_iterator_new((GtFeatureNode*) first_node);
    fni_b = gt_feature_node_iterator_new((GtFeatureNode*) second_node);
    for (fn_a = gt_feature_node_iterator_next(fni_a),
         fn_b = gt_feature_node_iterator_next(fni_b);
         fn_a && fn_b;
         fn_a = gt_feature_node_iterator_next(fni_a),
         fn_b = gt_feature_node_iterator_next(fni_b)) {
      if (!fn_b || !gt_feature_node_is_similar(fn_a, fn_b))
        break;
    }
    fn_b = (GtFeatureNode*) gt_feature_node_iterator_next(fni_b);
    if (!fn_a && !fn_b)
      equal = true;
    gt_feature_node_iterator_delete(fni_a);
    gt_feature_node_iterator_delete(fni_b);
    return equal;
  }
  return false;
}

static bool uniq(GtGenomeNode **first_node, GtGenomeNode **second_node)
{
  bool first_score_is_defined, second_score_is_defined;
  float first_score = 0.0, second_score = 0.0;
  gt_assert(*first_node && *second_node);
  if (nodes_are_equal_feature_trees(*first_node, *second_node)) {
    if ((first_score_is_defined =
           gt_feature_node_score_is_defined((GtFeatureNode*) *first_node))) {
      first_score = gt_feature_node_get_score((GtFeatureNode*) *first_node);
    }
    if ((second_score_is_defined =
           gt_feature_node_score_is_defined((GtFeatureNode*) *second_node))) {
      second_score = gt_feature_node_get_score((GtFeatureNode*) *second_node);
    }
    if ((!first_score_is_defined && !second_score_is_defined) ||
        (first_score_is_defined && !second_score_is_defined) ||
        (first_score_is_defined && second_score_is_defined &&
         first_score >= second_score)) {
      /* keep first node */
      gt_genome_node_delete(*second_node);
    }
    else {
      /* keep second node */
      gt_genome_node_delete(*first_node);
      *first_node = *second_node;
    }
    *second_node = NULL;
    return true;
  }
  return false;
}

static int uniq_stream_next(GtNodeStream *ns, GtGenomeNode **gn, GtError *err)
{
  GtUniqStream *us;
  int had_err;
  gt_error_check(err);
  us = uniq_stream_cast(ns);

  gt_assert(!us->second_node); /* the second buffer is always empty when this
                               function is called */
  if (!us->first_node) {
    /* both buffers are empty */
    had_err = gt_node_stream_next(us->in_stream, &us->first_node, err);
    if (had_err)
      return had_err;
    if (!us->first_node) {
      *gn = NULL;
      return 0;
    }
  }

  /* uniq loop */
  for (;;) {
    gt_assert(us->first_node && !us->second_node);
    had_err = gt_node_stream_next(us->in_stream, &us->second_node, err);
    if (!had_err && us->second_node) {
      if (!uniq(&us->first_node, &us->second_node))
        break; /* no uniq possible */
    }
    else
      break;
  }

  /* serve node */
  if (!had_err) {
    gt_assert(us->first_node);
    *gn = us->first_node;
    us->first_node = us->second_node;
    us->second_node = NULL;
  }

  return had_err;
}

static void uniq_stream_free(GtNodeStream *ns)
{
  GtUniqStream *us = uniq_stream_cast(ns);
  gt_genome_node_delete(us->first_node);
  gt_genome_node_delete(us->second_node);
  gt_node_stream_delete(us->in_stream);
}

const GtNodeStreamClass* gt_uniq_stream_class(void)
{
  static const GtNodeStreamClass *nsc = NULL;
  gt_class_alloc_lock_enter();
  if (!nsc) {
    nsc = gt_node_stream_class_new(sizeof (GtUniqStream),
                                   uniq_stream_free,
                                   uniq_stream_next);
  }
  gt_class_alloc_lock_leave();
  return nsc;
}

GtNodeStream* gt_uniq_stream_new(GtNodeStream *in_stream)
{
  GtNodeStream *ns;
  GtUniqStream *us;
  gt_assert(in_stream && gt_node_stream_is_sorted(in_stream));
  ns = gt_node_stream_create(gt_uniq_stream_class(), true);
  us = uniq_stream_cast(ns);
  us->in_stream = gt_node_stream_ref(in_stream);
  return ns;
}
