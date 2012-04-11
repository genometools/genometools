/*
  Copyright (c) 2008-2009 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2008      Center for Bioinformatics, University of Hamburg

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
#include "core/cstr_api.h"
#include "core/hashmap.h"
#include "core/ma.h"
#include "extended/feature_info.h"
#include "extended/genome_node.h"
#include "extended/gff3_defines.h"

struct GtFeatureInfo {
  GtHashmap *id_to_genome_node,
            *id_to_pseudo_parent;
};

GtFeatureInfo* gt_feature_info_new(void)
{
  GtFeatureInfo *fi = gt_malloc(sizeof *fi);
  fi->id_to_genome_node = gt_hashmap_new(GT_HASH_STRING, gt_free_func,
                                         (GtFree) gt_genome_node_delete);
  fi->id_to_pseudo_parent = gt_hashmap_new(GT_HASH_STRING, gt_free_func,
                                           (GtFree) gt_genome_node_delete);
  return fi;
}

void gt_feature_info_delete(GtFeatureInfo *fi)
{
  if (!fi) return;
  gt_hashmap_delete(fi->id_to_genome_node);
  gt_hashmap_delete(fi->id_to_pseudo_parent);
  gt_free(fi);
}

void gt_feature_info_reset(GtFeatureInfo *fi)
{
  gt_assert(fi);
  gt_hashmap_reset(fi->id_to_genome_node);
  gt_hashmap_reset(fi->id_to_pseudo_parent);
}

GtFeatureNode* gt_feature_info_get(const GtFeatureInfo *fi, const char *id)
{
  gt_assert(fi && id);
  return gt_hashmap_get(fi->id_to_genome_node, id);
}

void gt_feature_info_add(GtFeatureInfo *fi, const char *id, GtFeatureNode *fn)
{
  gt_assert(fi && id && fn);
  gt_assert(!gt_feature_node_is_pseudo((GtFeatureNode*) fn));
  gt_hashmap_add(fi->id_to_genome_node, gt_cstr_dup(id),
                 gt_genome_node_ref((GtGenomeNode*) fn));
}

GtFeatureNode* gt_feature_info_get_pseudo_parent(const GtFeatureInfo *fi,
                                                 const char *id)
{
  gt_assert(fi && id);
  return gt_hashmap_get(fi->id_to_pseudo_parent, id);
}

void gt_feature_info_add_pseudo_parent(GtFeatureInfo *fi, const char *id,
                                       GtFeatureNode *pseudo_parent)
{
  gt_assert(fi && id && pseudo_parent);
  gt_assert(gt_feature_node_is_pseudo((GtFeatureNode*) pseudo_parent));
  gt_hashmap_add(fi->id_to_pseudo_parent, gt_cstr_dup(id),
                 gt_genome_node_ref((GtGenomeNode*) pseudo_parent));
}

void gt_feature_info_replace_pseudo_parent(GtFeatureInfo *fi,
                                           GtFeatureNode *child,
                                           GtFeatureNode *new_pseudo_parent)
{
  const char *id;
  gt_assert(fi && child && new_pseudo_parent);
  gt_assert(gt_feature_node_is_pseudo((GtFeatureNode*) new_pseudo_parent));
  id = gt_feature_node_get_attribute(child, GT_GFF_ID);
  gt_assert(id);
  gt_hashmap_remove(fi->id_to_pseudo_parent, id);
  gt_feature_info_add_pseudo_parent(fi, id, new_pseudo_parent);
}

static GtFeatureNode* find_root(const GtFeatureInfo *fi, const char *id)
{
  const char *delim, *parents;
  GtFeatureNode *this_feature, *parent_pseudo_feature;
  gt_assert(fi && id);
  /* get feature */
  delim = strchr(id, ';');
  if (delim) {
    char *first_parent = gt_cstr_dup_nt(id, delim - id);
    this_feature = gt_hashmap_get(fi->id_to_genome_node, first_parent);
    parent_pseudo_feature = gt_hashmap_get(fi->id_to_pseudo_parent,
                                           first_parent);
    gt_free(first_parent);
  }
  else {
    this_feature = gt_hashmap_get(fi->id_to_genome_node, id);
    parent_pseudo_feature = gt_hashmap_get(fi->id_to_pseudo_parent, id);
  }
  gt_assert(this_feature);
  /* recursion */
  parents = gt_feature_node_get_attribute(this_feature, GT_GFF_PARENT);
  if (parents)
    return find_root(fi, parents);
  else if (parent_pseudo_feature)
    return parent_pseudo_feature;
  return this_feature;
}

GtFeatureNode* gt_feature_info_find_root(const GtFeatureInfo *fi,
                                         const char *id)
{
  gt_assert(fi && id);
  gt_assert(gt_feature_info_get(fi, id));
  return find_root(fi, id);
}
