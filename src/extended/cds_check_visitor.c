/*
  Copyright (c) 2010, 2012 Gordon Gremme <gremme@zbh.uni-hamburg.de>

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
#include "core/unused_api.h"
#include "core/warning_api.h"
#include "extended/cds_check_visitor.h"
#include "extended/feature_node_iterator_api.h"
#include "extended/feature_node.h"
#include "extended/feature_type.h"
#include "extended/gff3_defines.h"
#include "extended/node_visitor_api.h"

struct GtCDSCheckVisitor {
  const GtNodeVisitor parent_instance;
  GtHashmap *cds_features,
            *cds_features_to_split;
  bool tidy,
       splitting_is_necessary;
};

#define cds_check_visitor_cast(NS)\
        gt_node_visitor_cast(gt_cds_check_visitor_class(), NS)

static int check_cds_phases(GtArray *cds_features, GtCDSCheckVisitor *v,
                            bool is_multi, bool second_pass, GtError *err)
{
  GtPhase current_phase, correct_phase = GT_PHASE_ZERO;
  GtFeatureNode *fn;
  GtStrand strand;
  unsigned long i, current_length;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(cds_features);
  gt_assert(gt_array_size(cds_features));
  fn = *(GtFeatureNode**) gt_array_get_first(cds_features);
  strand = gt_feature_node_get_strand(fn);
  if (strand == GT_STRAND_REVERSE)
    gt_array_reverse(cds_features);
  for (i = 0; !had_err && i < gt_array_size(cds_features); i++) {
    fn = *(GtFeatureNode**) gt_array_get(cds_features, i);
    /* the first phase can be anything (except being undefined), because the
       GFF3 spec says:

       NOTE 4 - CDS features MUST have have a defined phase field. Otherwise it
       is not possible to infer the correct polypeptides corresponding to
       partially annotated genes. */
    if ((!i && gt_feature_node_get_phase(fn) == GT_PHASE_UNDEFINED) ||
        (i && gt_feature_node_get_phase(fn) != correct_phase)) {
      if (gt_hashmap_get(v->cds_features, fn)) {
        if (v->tidy && !is_multi && !gt_feature_node_has_children(fn)) {
          /* we can split the feature */
          gt_warning("%s feature on line %u in file \"%s\" has multiple "
                     "parents which require different phases; split feature",
                     gt_ft_CDS,
                     gt_genome_node_get_line_number((GtGenomeNode*) fn),
                     gt_genome_node_get_filename((GtGenomeNode*) fn));
          gt_hashmap_add(v->cds_features_to_split, fn, fn);
          v->splitting_is_necessary = true; /* split later */
        }
        else {
          gt_error_set(err, "%s feature on line %u in file \"%s\" has multiple "
                       "parents which require different phases",
                       gt_ft_CDS,
                       gt_genome_node_get_line_number((GtGenomeNode*) fn),
                       gt_genome_node_get_filename((GtGenomeNode*) fn));
          had_err = -1;
        }
      }
      else {
        if (v->tidy) {
          if (!second_pass) {
            gt_warning("%s feature on line %u in file \"%s\" has the wrong "
                       "phase %c -> correcting it to %c", gt_ft_CDS,
                       gt_genome_node_get_line_number((GtGenomeNode*) fn),
                       gt_genome_node_get_filename((GtGenomeNode*) fn),
                       GT_PHASE_CHARS[gt_feature_node_get_phase(fn)],
                       GT_PHASE_CHARS[correct_phase]);
          }
          gt_feature_node_set_phase(fn, correct_phase);
        }
        else {
          gt_error_set(err, "%s feature on line %u in file \"%s\" has the "
                       "wrong phase %c (should be %c)", gt_ft_CDS,
                       gt_genome_node_get_line_number((GtGenomeNode*) fn),
                       gt_genome_node_get_filename((GtGenomeNode*) fn),
                       GT_PHASE_CHARS[gt_feature_node_get_phase(fn)],
                       GT_PHASE_CHARS[correct_phase]);
          had_err = -1;
        }
      }
    }
    if (!had_err) {
      current_phase = gt_feature_node_get_phase(fn);
      current_length = gt_genome_node_get_length((GtGenomeNode*) fn);
      correct_phase = (3 - (current_length - current_phase) % 3) % 3;
      gt_hashmap_add(v->cds_features, fn, fn); /* record CDS feature */
    }
  }
  return had_err;
}

static int check_cds_phases_hm(GT_UNUSED void *key, void *value, void *data,
                               GtError *err)
{
  GtArray *cds_features = value;
  gt_error_check(err);
  gt_assert(cds_features && data);
  return check_cds_phases(cds_features, data, true, false, err);
}

static int check_cds_phases_if_necessary(GtFeatureNode *fn,
                                         GtCDSCheckVisitor *v,
                                         bool second_pass, GtError *err)
{
  GtFeatureNodeIterator *fni;
  GtFeatureNode *node;
  GtArray *cds_features = NULL;
  GtHashmap *multi_features = NULL;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(fn);
  fni = gt_feature_node_iterator_new_direct(fn);
  while ((node = gt_feature_node_iterator_next(fni))) {
    if (gt_feature_node_has_type(node, gt_ft_CDS)) {
      if (gt_feature_node_is_multi(node)) {
        GtArray *features;
        if (!multi_features)
          multi_features = gt_hashmap_new(GT_HASH_DIRECT, NULL,
                                          (GtFree) gt_array_delete);
        if ((features =
                gt_hashmap_get(multi_features,
                             gt_feature_node_get_multi_representative(node)))) {
          gt_array_add(features, node);
        }
        else {
          GtFeatureNode *representative;
          features = gt_array_new(sizeof (GtFeatureNode*));
          representative = gt_feature_node_get_multi_representative(node);
          gt_array_add(features, representative);
          gt_hashmap_add(multi_features, representative, features);
        }
      }
      else {
        if (!cds_features)
          cds_features = gt_array_new(sizeof (GtFeatureNode*));
        gt_array_add(cds_features, node);
      }
    }
  }
  if (cds_features)
    had_err = check_cds_phases(cds_features, v, false, second_pass, err);
  if (!had_err && multi_features)
    had_err = gt_hashmap_foreach(multi_features, check_cds_phases_hm, v, err);
  gt_array_delete(cds_features);
  gt_hashmap_delete(multi_features);
  gt_feature_node_iterator_delete(fni);
  return had_err;
}

static int collect_cds_feature(void *key, GT_UNUSED void *value, void *data,
                               GT_UNUSED GtError *err)
{
  GtArray *cds_features = data;
  gt_error_check(err);
  gt_array_add(cds_features, key);
  return 0;
}

static GtArray* find_cds_parents(GtFeatureNode *cds_feature, GtFeatureNode *fn)
{
  GtFeatureNodeIterator *fni, *di;
  GtFeatureNode *parent, *child;
  GtArray *parents;
  gt_assert(cds_feature && fn);
  parents = gt_array_new(sizeof (GtFeatureNode*));
  fni = gt_feature_node_iterator_new(fn);
  while ((parent = gt_feature_node_iterator_next(fni))) {
    di = gt_feature_node_iterator_new_direct(parent);
    while ((child = gt_feature_node_iterator_next(di))) {
      if (child == cds_feature)
        gt_array_add(parents, parent);
    }
    gt_feature_node_iterator_delete(di);
  }
  gt_feature_node_iterator_delete(fni);
  return parents;
}

static void split_cds_feature(GtFeatureNode *cds_feature, GtFeatureNode *fn)
{
  GtArray *parents;
  unsigned long i;
  gt_assert(cds_feature && fn);

  /* find parents */
  parents = find_cds_parents(cds_feature, fn);

  /* remove CDS feature */
  gt_feature_node_remove_leaf(fn, cds_feature);

  /* add CDS feature to all parents */
  for (i = 0; i < gt_array_size(parents); i++) {
    GtFeatureNode *parent = *(GtFeatureNode**) gt_array_get(parents, i);
    const char *id = gt_feature_node_get_attribute(parent, GT_GFF_ID);
    if (!i) {
      gt_feature_node_set_attribute(cds_feature, GT_GFF_PARENT, id);
      gt_feature_node_add_child(parent, cds_feature);
    }
    else {
      GtFeatureNode *new_cds = gt_feature_node_clone(cds_feature);
      gt_feature_node_set_attribute(new_cds, GT_GFF_PARENT, id);
      gt_feature_node_add_child(parent, new_cds);
      gt_genome_node_delete((GtGenomeNode*) cds_feature);
    }
  }

  gt_array_delete(parents);
}

static void split_cds_features(GtHashmap *features_to_split, GtFeatureNode *fn)
{
  GtArray *cds_features;
  unsigned long i;
  GT_UNUSED int had_err;
  gt_assert(features_to_split && fn);
  cds_features = gt_array_new(sizeof (GtFeatureNode*));
  had_err = gt_hashmap_foreach(features_to_split, collect_cds_feature,
                               cds_features, NULL);
  gt_assert(!had_err); /* collect_cds_feature() is sane */
  gt_assert(gt_array_size(cds_features));
  for (i = 0; i < gt_array_size(cds_features); i++) {
    GtFeatureNode *cds_feature;
    cds_feature = *(GtFeatureNode**) gt_array_get(cds_features, i);
    split_cds_feature(cds_feature, fn);
  }
  gt_array_delete(cds_features);
}

static int cds_check_visitor_feature_node(GtNodeVisitor *nv, GtFeatureNode *fn,
                                          GtError *err)
{
  GtCDSCheckVisitor *v = cds_check_visitor_cast(nv);
  GtFeatureNodeIterator *fni;
  GtFeatureNode *node;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(v && fn);
  fni = gt_feature_node_iterator_new(fn);
  while (!had_err && (node = gt_feature_node_iterator_next(fni)))
    had_err = check_cds_phases_if_necessary(node, v, false, err);
  gt_feature_node_iterator_delete(fni);
  gt_hashmap_reset(v->cds_features);
  while (v->splitting_is_necessary) {
    split_cds_features(v->cds_features_to_split, fn);
    gt_hashmap_reset(v->cds_features_to_split);
    v->splitting_is_necessary = false;
    /* perform second pass to correct phases */
    fni = gt_feature_node_iterator_new(fn);
    while (!had_err && (node = gt_feature_node_iterator_next(fni)))
      had_err = check_cds_phases_if_necessary(node, v, false, err);
    gt_feature_node_iterator_delete(fni);
    gt_hashmap_reset(v->cds_features);
  }
  return had_err;
}

static void cds_check_visitor_free(GtNodeVisitor *nv)
{
  GtCDSCheckVisitor *v = cds_check_visitor_cast(nv);
  gt_hashmap_delete(v->cds_features_to_split);
  gt_hashmap_delete(v->cds_features);
}

const GtNodeVisitorClass* gt_cds_check_visitor_class()
{
  static const GtNodeVisitorClass *nvc = NULL;
  gt_class_alloc_lock_enter();
  if (!nvc) {
    nvc = gt_node_visitor_class_new(sizeof (GtCDSCheckVisitor),
                                    cds_check_visitor_free,
                                    NULL,
                                    cds_check_visitor_feature_node,
                                    NULL,
                                    NULL,
                                    NULL);
  }
  gt_class_alloc_lock_leave();
  return nvc;
}

GtNodeVisitor* gt_cds_check_visitor_new(void)
{
  GtNodeVisitor *nv = gt_node_visitor_create(gt_cds_check_visitor_class());
  GtCDSCheckVisitor *v = cds_check_visitor_cast(nv);
  v->cds_features = gt_hashmap_new(GT_HASH_DIRECT, NULL, NULL);
  v->cds_features_to_split = gt_hashmap_new(GT_HASH_DIRECT, NULL, NULL);
  return nv;
}

void gt_cds_check_visitor_enable_tidy_mode(GtCDSCheckVisitor *v)
{
  gt_assert(v);
  v->tidy = true;
}
