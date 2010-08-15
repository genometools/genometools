/*
  Copyright (c) 2010 Gordon Gremme <gremme@zbh.uni-hamburg.de>

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

#include "core/hashmap_api.h"
#include "core/unused_api.h"
#include "core/warning_api.h"
#include "extended/cds_check_visitor.h"
#include "extended/feature_node_iterator_api.h"
#include "extended/node_visitor_rep.h"

struct GtCDSCheckVisitor {
  const GtNodeVisitor parent_instance;
  bool tidy;
};

#define cds_check_visitor_cast(NS)\
        gt_node_visitor_cast(gt_cds_check_visitor_class(), NS)

static int check_cds_phases(GtArray *cds_features, bool tidy, GtError *err)
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
    if (gt_feature_node_get_phase(fn) != correct_phase) {
      if (!tidy) {
        gt_error_set(err, "%s feature on line %u in file \"%s\" has the wrong "
                     "phase %c (should be %c)", gt_ft_CDS,
                     gt_genome_node_get_line_number((GtGenomeNode*) fn),
                     gt_genome_node_get_filename((GtGenomeNode*) fn),
                     GT_PHASE_CHARS[gt_feature_node_get_phase(fn)],
                     GT_PHASE_CHARS[correct_phase]);
        had_err = -1;
      }
      else {
        gt_warning("%s feature on line %u in file \"%s\" has the wrong phase "
                   "%c -> correcting it to %c", gt_ft_CDS,
                   gt_genome_node_get_line_number((GtGenomeNode*) fn),
                   gt_genome_node_get_filename((GtGenomeNode*) fn),
                   GT_PHASE_CHARS[gt_feature_node_get_phase(fn)],
                   GT_PHASE_CHARS[correct_phase]);
        gt_feature_node_set_phase(fn, correct_phase);
      }
    }
    if (!had_err) {
      current_phase = gt_feature_node_get_phase(fn);
      current_length = gt_genome_node_get_length((GtGenomeNode*) fn);
      correct_phase = (3 - (current_length - current_phase) % 3) % 3;
    }
  }
  return had_err;
}

static int check_cds_phases_hm(GT_UNUSED void *key, void *value, void *data,
                               GtError *err)
{
  GtArray *cds_features = value;
  bool *tidy = data;
  gt_error_check(err);
  gt_assert(cds_features && tidy);
  return check_cds_phases(cds_features, *tidy, err);
}

static int check_cds_phases_if_necessary(GtFeatureNode *fn, bool tidy,
                                         GtError *err)
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
    had_err = check_cds_phases(cds_features, tidy, err);
  if (!had_err && multi_features) {
    had_err = gt_hashmap_foreach(multi_features, check_cds_phases_hm, &tidy,
                                 err);
  }
  gt_array_delete(cds_features);
  gt_hashmap_delete(multi_features);
  gt_feature_node_iterator_delete(fni);
  return had_err;
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
    had_err = check_cds_phases_if_necessary(node, v->tidy, err);
  gt_feature_node_iterator_delete(fni);
  return had_err;
}

const GtNodeVisitorClass* gt_cds_check_visitor_class()
{
  static const GtNodeVisitorClass *nvc = NULL;
  if (!nvc) {
    nvc = gt_node_visitor_class_new(sizeof (GtCDSCheckVisitor),
                                    NULL,
                                    NULL,
                                    cds_check_visitor_feature_node,
                                    NULL,
                                    NULL);
  }
  return nvc;
}

GtNodeVisitor* gt_cds_check_visitor_new(void)
{
  return gt_node_visitor_create(gt_cds_check_visitor_class());
}

void gt_cds_check_visitor_enable_tidy_mode(GtCDSCheckVisitor *v)
{
  gt_assert(v);
  v->tidy = true;
}
