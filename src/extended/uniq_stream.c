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
#include "core/undef.h"
#include "extended/genome_node_iterator.h"
#include "extended/node_stream_rep.h"
#include "extended/uniq_stream.h"

struct UniqStream{
  const GtNodeStream parent_instance;
  GtNodeStream *in_stream;
  GtGenomeNode *first_node,
             *second_node;
};

#define uniq_stream_cast(GS)\
        gt_node_stream_cast(uniq_stream_class(), GS)

static bool nodes_are_equal_feature_trees(GtGenomeNode *first_node,
                                          GtGenomeNode *second_node)
{
  bool equal = false;
  GtGenomeNodeIterator *gni_a, *gni_b;
  GtFeatureNode *gf_a, *gf_b;
  gf_a = gt_genome_node_cast(gt_feature_node_class(), first_node);
  gf_b = gt_genome_node_cast(gt_feature_node_class(), second_node);
  if (gf_a && gf_b) {
    gni_a = gt_genome_node_iterator_new(first_node);
    gni_b = gt_genome_node_iterator_new(second_node);
    for (gf_a = (GtFeatureNode*) gt_genome_node_iterator_next(gni_a),
         gf_b = (GtFeatureNode*) gt_genome_node_iterator_next(gni_b);
         gf_a && gf_b;
         gf_a = (GtFeatureNode*) gt_genome_node_iterator_next(gni_a),
         gf_b = (GtFeatureNode*) gt_genome_node_iterator_next(gni_b)) {
      if (!gf_b || !gt_genome_features_are_similar(gf_a, gf_b))
        break;
    }
    gf_b = (GtFeatureNode*) gt_genome_node_iterator_next(gni_b);
    if (!gf_a && !gf_b)
      equal = true;
    gt_genome_node_iterator_delete(gni_a);
    gt_genome_node_iterator_delete(gni_b);
    return equal;
  }
  return false;
}

static bool uniq(GtGenomeNode **first_node, GtGenomeNode **second_node)
{
  bool first_score_is_defined, second_score_is_defined;
  float first_score = 0.0, second_score;
  assert(*first_node && *second_node);
  if (nodes_are_equal_feature_trees(*first_node, *second_node)) {
    if ((first_score_is_defined =
           gt_feature_node_score_is_defined((GtFeatureNode*)
                                              *first_node))) {
      first_score = gt_feature_node_get_score((GtFeatureNode*)
                                                *first_node);
    }
    if ((second_score_is_defined =
           gt_feature_node_score_is_defined((GtFeatureNode*)
                                              *second_node))) {
      second_score = gt_feature_node_get_score((GtFeatureNode*)
                                                 *second_node);
    }
    if ((!first_score_is_defined && !second_score_is_defined) ||
        (first_score_is_defined && !second_score_is_defined) ||
        (first_score_is_defined && second_score_is_defined &&
         first_score >= second_score)) {
      /* keep first node */
      gt_genome_node_rec_delete(*second_node);
    }
    else {
      /* keep second node */
      gt_genome_node_rec_delete(*first_node);
      *first_node = *second_node;
    }
    *second_node = NULL;
    return true;
  }
  return false;
}

static int uniq_stream_next_tree(GtNodeStream *gs, GtGenomeNode **gn,
                                 GtError *err)
{
  UniqStream *us;
  int had_err;
  gt_error_check(err);
  us = uniq_stream_cast(gs);

  assert(!us->second_node); /* the second buffer is always empty when this
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
    assert(us->first_node && !us->second_node);
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
    assert(us->first_node);
    *gn = us->first_node;
    us->first_node = us->second_node;
    us->second_node = NULL;
  }

  return had_err;
}

static void uniq_stream_free(GtNodeStream *gs)
{
  UniqStream *us = uniq_stream_cast(gs);
  gt_genome_node_rec_delete(us->first_node);
  gt_genome_node_rec_delete(us->second_node);
  gt_node_stream_delete(us->in_stream);
}

const GtNodeStreamClass* uniq_stream_class(void)
{
  static const GtNodeStreamClass gsc = { sizeof (UniqStream),
                                         uniq_stream_next_tree,
                                         uniq_stream_free };
  return &gsc;
}

GtNodeStream* uniq_stream_new(GtNodeStream *in_stream)
{
  GtNodeStream *gs;
  UniqStream *us;
  assert(in_stream && gt_node_stream_is_sorted(in_stream));
  gs = gt_node_stream_create(uniq_stream_class(), true);
  us = uniq_stream_cast(gs);
  us->in_stream = gt_node_stream_ref(in_stream);
  return gs;
}
