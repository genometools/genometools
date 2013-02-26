/*
  Copyright (c) 2011 Sascha Kastens <sascha.kastens@studium.uni-hamburg.de>
  Copyright (c) 2011 Center for Bioinformatics, University of Hamburg

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
#include "core/array.h"
#include "core/bittab.h"
#include "core/class_alloc_lock.h"
#include "core/ma.h"
#include "core/cstr_api.h"
#include "core/hashmap.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "extended/feature_node.h"
#include "extended/feature_node_iterator_api.h"
#include "extended/feature_type.h"
#include "extended/genome_node.h"
#include "ltr/ltr_classify_stream.h"

struct GtLTRClassifyStream {
  const GtNodeStream parent_instance;
  GtNodeStream *in_stream;
  GtArray *nodes;
  GtHashmap *features;
  const char *famprefix;
  char **current_state;
  bool first_next;
  unsigned long next_index,
                *progress;
};

#define gt_ltr_classify_stream_cast(CS)\
        gt_node_stream_cast(gt_ltr_classify_stream_class(), CS);

static bool ltr_candidates_compatible(GtGenomeNode *candidate1,
                                      GtGenomeNode *candidate2,
                                      GtHashmap *features,
                                      GT_UNUSED GtError *err)
{
  GtGenomeNode *candidate;
  GtFeatureNode *curnode1,
                *curnode2;
  GtFeatureNodeIterator *fni1,
                        *fni2;
  GtHashmap *fnmap;
  const char *clid1,
             *clid2,
             *fnt;
  bool compatible = false,
       seen_a = false,
       first_ltr = true;
  unsigned long clnum1,
                clnum2;

  gt_error_check(err);
  gt_assert(candidate1 && candidate2);

  fni1 = gt_feature_node_iterator_new((GtFeatureNode*) candidate1);
  curnode1 = gt_feature_node_iterator_next(fni1);
  fni2 = gt_feature_node_iterator_new((GtFeatureNode*) candidate2);
  curnode2 = gt_feature_node_iterator_next(fni2);

  if (gt_feature_node_number_of_children(curnode1) <=
      gt_feature_node_number_of_children(curnode2)) {
    gt_feature_node_iterator_delete(fni2);
    candidate = candidate2;
  } else {
    gt_feature_node_iterator_delete(fni1);
    fni1 = fni2;
    candidate = candidate1;
  }
  fnmap = (GtHashmap*) gt_genome_node_get_user_data(candidate, "fnmap");

  while ((curnode1 = gt_feature_node_iterator_next(fni1)) != NULL) {
    fnt = gt_feature_node_get_type(curnode1);
    if ((strcmp(fnt, gt_ft_LTR_retrotransposon) == 0))
      continue;
    else if (strcmp(fnt, gt_ft_long_terminal_repeat) == 0) {
      if (first_ltr) {
        fnt = "lLTR";
        first_ltr = false;
      } else
        fnt = "rLTR";
    } else if (strcmp(fnt, gt_ft_protein_match) == 0) {
      fnt = gt_feature_node_get_attribute(curnode1, "name");
    }
    if (!fnt)
      continue;
    if (features != NULL)
      if (gt_hashmap_get(features, fnt) == NULL)
        continue;
    curnode2 = (GtFeatureNode*) gt_hashmap_get(fnmap, fnt);
    if (curnode2 != NULL) {
      clid1 = gt_feature_node_get_attribute(curnode1, "clid");
      if (clid1 != NULL)
        (void) sscanf(clid1, "%lu", &clnum1);
      else
        clnum1 = GT_UNDEF_ULONG;
      clid2 = gt_feature_node_get_attribute(curnode2, "clid");
      if (clid2 != NULL)
        (void) sscanf(clid2, "%lu", &clnum2);
      else
        clnum2 = GT_UNDEF_ULONG;
      if (clnum1 == clnum2) {
        if ((clnum1 != GT_UNDEF_ULONG) && (clnum2 != GT_UNDEF_ULONG))
          seen_a = true;
        compatible = true;
      }
      else if ((clnum1 == GT_UNDEF_ULONG && clnum2 != GT_UNDEF_ULONG) ||
               (clnum1 != GT_UNDEF_ULONG && clnum2 == GT_UNDEF_ULONG))
        compatible = true;
      else if (clnum1 != clnum2) {
        compatible = false;
        break;
      }
    }
  }
  gt_feature_node_iterator_delete(fni1);

  return (compatible && seen_a);
}

static bool ltr_group_compatible(GtArray *candidates, GtGenomeNode *candidate1,
                                 GtBittab *group, GtHashmap *features,
                                 GtError *err)
{
  GtArray *group_member;
  GtGenomeNode *candidate2;
  bool compatible = true;
  unsigned long i, index;

  gt_assert(candidates && candidate1 && group);
  gt_error_check(err);

  group_member = gt_array_new(sizeof(unsigned long));
  gt_bittab_get_all_bitnums(group, group_member);

  for (i = 0; i < gt_array_size(group_member); i++) {
    index = *(unsigned long*) gt_array_get(group_member, i);
    candidate2 = *(GtGenomeNode**) gt_array_get(candidates, index);
    if (!(compatible = ltr_candidates_compatible(candidate1, candidate2,
                                                 features, err)))
      break;
  }
  gt_array_delete(group_member);

  return compatible;
}

static int check_ambiguous_candidates(GtArray *candidates, GtArray *groups,
                                      GtHashmap *features,
                                      unsigned long *progress,
                                      GtError *err)
{
  GtGenomeNode *candidate;
  GtBittab *group;
  GtArray *compat_groups;
  int had_err = 0;
  unsigned long i,
                j;

  gt_error_check(err);
  gt_assert(candidates && groups);

  for (i = 0; i < gt_array_size(candidates); i++) {
    candidate = *(GtGenomeNode**) gt_array_get(candidates, i);
    if (gt_feature_node_try_cast(candidate) == NULL)
      continue;
    compat_groups = gt_array_new(sizeof(GtBittab*));
    for (j = 0; j < gt_array_size(groups); j++) {
      group = *(GtBittab**) gt_array_get(groups, j);
      if (ltr_group_compatible(candidates, candidate, group, features, err))
        gt_array_add(compat_groups, group);
    }
    if (gt_array_size(compat_groups) > 1UL) {
      for (j = 0; j < gt_array_size(groups); j++) {
        group = *(GtBittab**) gt_array_get(groups, j);
        gt_bittab_unset_bit(group, i);
        if (gt_bittab_count_set_bits(group) == 0)
          gt_array_rem(groups, j);
      }
    }
    gt_array_delete(compat_groups);
    if (progress != NULL)
      *progress = *progress + 1;
  }
  return had_err;
}

static int annotate_nodes(GtArray *candidates, GtArray *groups,
                          const char *famprefix, GtError *err)
{
  GtArray *group_member;
  GtBittab *group;
  GtFeatureNode *curnode;
  GtFeatureNodeIterator *fni = NULL;
  GtGenomeNode *gn;
  int had_err = 0;
  unsigned long i, j, index, famno = 0;

  gt_assert(candidates && groups);
  gt_error_check(err);

  for (i = 0; i < gt_array_size(groups); i++) {
    group = *(GtBittab**) gt_array_get(groups, i);
    group_member = gt_array_new(sizeof(unsigned long));
    gt_bittab_get_all_bitnums(group, group_member);
    if (gt_array_size(group_member) < 2UL) {
      gt_array_delete(group_member);
      continue;
    }
    for (j = 0; j < gt_array_size(group_member); j++) {
      index = *(unsigned long*) gt_array_get(group_member, j);
      gn = *(GtGenomeNode**) gt_array_get(candidates, index);
      fni = gt_feature_node_iterator_new((GtFeatureNode*) gn);
      curnode = gt_feature_node_iterator_next(fni);
      if (strcmp(gt_feature_node_get_type(curnode), gt_ft_repeat_region) == 0) {
        char fam[BUFSIZ];
        if (famprefix != NULL)
          (void) snprintf(fam, BUFSIZ, "%s%lu", famprefix, famno);
        else
          (void) snprintf(fam, BUFSIZ, "ltrfam_%lu", famno);
        gt_feature_node_set_attribute(curnode, "ltrfam", fam);
      } else {
        gt_feature_node_iterator_delete(fni);
        gt_error_set(err, "repeat_region is not root node");
        had_err = -1;
        break;
      }
      gt_feature_node_iterator_delete(fni);
    }
    gt_array_delete(group_member);
    if (had_err)
      break;
    famno++;
  }
  return had_err;
}

static void delete_hash(void *hash)
{
  gt_hashmap_delete((GtHashmap*) hash);
}

static void add_fnmap_to_candidates(GtArray *candidates)
{
  GtGenomeNode *gn;
  GtFeatureNode *curnode;
  GtFeatureNodeIterator *fni;
  GtHashmap *fnmap;
  bool first_ltr;
  const char *fnt;
  unsigned long i;

  for (i = 0; i < gt_array_size(candidates); i++) {
    gn = *(GtGenomeNode**) gt_array_get(candidates, i);
    if (gt_feature_node_try_cast(gn) == NULL)
      continue;
    first_ltr = true;
    fnmap = gt_hashmap_new(GT_HASH_STRING, gt_free_func, NULL);
    fni = gt_feature_node_iterator_new((GtFeatureNode*) gn);
    while ((curnode = gt_feature_node_iterator_next(fni)) != NULL) {
      fnt = gt_feature_node_get_type(curnode);
      if ((strcmp(fnt, gt_ft_repeat_region) == 0) ||
          (strcmp(fnt, gt_ft_LTR_retrotransposon) == 0))
        continue;
      else if (strcmp(fnt, gt_ft_long_terminal_repeat) == 0) {
        if (first_ltr) {
          fnt = "lLTR";
          first_ltr = false;
        } else
          fnt = "rLTR";
      } else if (strcmp(fnt, gt_ft_protein_match) == 0)
        fnt = gt_feature_node_get_attribute(curnode, "name");
      if (!fnt)
        continue;
      gt_hashmap_add(fnmap, (void*) gt_cstr_dup(fnt), (void*) curnode);
    }
    gt_genome_node_add_user_data(gn, "fnmap", (void*) fnmap, delete_hash);
    gt_feature_node_iterator_delete(fni);
    fni = NULL;
  }
}

static void remove_fnmap_from_candidates(GtArray *candidates)
{
  GtGenomeNode *gn;
  unsigned long i;

  for (i = 0; i < gt_array_size(candidates); i++) {
    gn = *(GtGenomeNode**) gt_array_get(candidates, i);
    if (gt_feature_node_try_cast(gn) == NULL)
      continue;
    gt_genome_node_release_user_data(gn, "fnmap");
  }
}

static int classify_ltrs(GtArray *candidates, GtHashmap *features,
                         const char *famprefix, char **current_state,
                         unsigned long *progress, GtError *err)
{
  GtBittab *group,
           *new_group;
  GtArray *groups;
  GtGenomeNode *candidate;
  bool sorted;
  int had_err = 0;
  unsigned long i, j, num_of_cands;

  gt_error_check(err);
  gt_assert(candidates);

  add_fnmap_to_candidates(candidates);
  num_of_cands = gt_array_size(candidates);
  groups = gt_array_new(sizeof(GtBittab*));

  if (current_state != NULL) {
    gt_free(*current_state);
    *current_state = gt_cstr_dup("Assigning candidates to families");
  }

  for (i = 0; i < num_of_cands; i++) {
    sorted = false;
    candidate = *(GtGenomeNode**) gt_array_get(candidates, i);
    if (gt_feature_node_try_cast(candidate) == NULL)
      continue;
    for (j = 0; j < gt_array_size(groups); j++) {
      group = *(GtBittab**) gt_array_get(groups, j);
      if (!sorted && ltr_group_compatible(candidates, candidate, group,
                                          features, err)) {
        gt_bittab_set_bit(group, i);
        sorted = true;
        break;
      }
    }
    if (!sorted) {
      new_group = gt_bittab_new(num_of_cands);
      gt_bittab_set_bit(new_group, i);
      gt_array_add(groups, new_group);
    }
    if (progress != NULL)
      *progress = *progress + 1;
  }
  if (current_state != NULL) {
    gt_free(*current_state);
    *current_state = gt_cstr_dup("Checking for ambiguity");
  }
  had_err = check_ambiguous_candidates(candidates, groups, features,
                                       progress, err);

  if (!had_err)
    had_err = annotate_nodes(candidates, groups, famprefix, err);

  for (i = 0; i < gt_array_size(groups); i++) {
    gt_bittab_delete(*(GtBittab**) gt_array_get(groups, i));
  }
  gt_array_delete(groups);
  remove_fnmap_from_candidates(candidates);

  return had_err;
}

static int gt_ltr_classify_stream_next(GtNodeStream *ns,
                                       GtGenomeNode **gn,
                                       GtError *err)
{
  GtLTRClassifyStream *lcs;
  GtGenomeNode *ref_gn;
  int had_err = 0;

  gt_error_check(err);
  lcs = gt_ltr_classify_stream_cast(ns);
  if (lcs->first_next) {
    while (!(had_err = gt_node_stream_next(lcs->in_stream, gn, err)) && *gn) {
      gt_assert(*gn && !had_err);
      ref_gn = gt_genome_node_ref(*gn);
      gt_array_add(lcs->nodes, ref_gn);
    }
    if (!had_err)
      had_err = classify_ltrs(lcs->nodes, lcs->features, lcs->famprefix,
                              lcs->current_state, lcs->progress, err);
    if (!had_err) {
      *gn = *(GtGenomeNode**) gt_array_get(lcs->nodes, lcs->next_index);
      lcs->next_index++;
      lcs->first_next = false;
      return 0;
    }
  } else {
    if (lcs->next_index >= gt_array_size(lcs->nodes))
      *gn = NULL;
    else {
      *gn = *(GtGenomeNode**) gt_array_get(lcs->nodes, lcs->next_index);
      lcs->next_index++;
    }
    return 0;
  }
  return had_err;
}

static void gt_ltr_classify_stream_free(GtNodeStream *ns)
{
  unsigned long i;
  GtLTRClassifyStream *lcs = gt_ltr_classify_stream_cast(ns);
  for (i = 0; i < gt_array_size(lcs->nodes); i++)
    gt_genome_node_delete(*(GtGenomeNode**) gt_array_get(lcs->nodes, i));
  gt_array_delete(lcs->nodes);
  gt_node_stream_delete(lcs->in_stream);
}

const GtNodeStreamClass* gt_ltr_classify_stream_class(void)
{
  static const GtNodeStreamClass *nsc;
  gt_class_alloc_lock_enter();
  if (!nsc) {
    nsc = gt_node_stream_class_new(sizeof(GtLTRClassifyStream),
                                   gt_ltr_classify_stream_free,
                                   gt_ltr_classify_stream_next);
  }
  gt_class_alloc_lock_leave();
  return nsc;
}

GtNodeStream* gt_ltr_classify_stream_new(GtNodeStream *in_stream,
                                         GtHashmap *features,
                                         const char *famprefix,
                                         char **current_state,
                                         unsigned long *progress,
                                         GT_UNUSED GtError *err)
{
  GtNodeStream *ns;
  GtLTRClassifyStream *lcs;
  ns = gt_node_stream_create(gt_ltr_classify_stream_class(), false);
  lcs = gt_ltr_classify_stream_cast(ns);
  lcs->in_stream = gt_node_stream_ref(in_stream);
  lcs->nodes = gt_array_new(sizeof(GtGenomeNode*));
  lcs->features = features;
  lcs->first_next = true;
  lcs->next_index = 0;
  lcs->famprefix = famprefix;
  lcs->progress = progress;
  lcs->current_state = current_state;
  return ns;
}
