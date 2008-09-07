/*
  Copyright (c) 2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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
#include "core/cstr.h"
#include "core/hashmap.h"
#include "core/ma.h"
#include "extended/feature_info.h"
#include "extended/genome_node.h"
#include "extended/gff3_parser.h"

struct FeatureInfo {
  Hashmap *id_to_genome_node,
          *id_to_pseudo_parent;
};

FeatureInfo* feature_info_new(void)
{
  FeatureInfo *fi = ma_malloc(sizeof *fi);
  fi->id_to_genome_node = hashmap_new(HASH_STRING, ma_free_func,
                                      (FreeFunc) genome_node_delete);
  fi->id_to_pseudo_parent = hashmap_new(HASH_STRING, ma_free_func,
                                        (FreeFunc) genome_node_delete);
  return fi;
}

void feature_info_delete(FeatureInfo *fi)
{
  if (!fi) return;
  hashmap_delete(fi->id_to_genome_node);
  hashmap_delete(fi->id_to_pseudo_parent);
  ma_free(fi);
}

void feature_info_reset(FeatureInfo *fi)
{
  assert(fi);
  hashmap_reset(fi->id_to_genome_node);
  hashmap_reset(fi->id_to_pseudo_parent);
}

GT_GenomeNode* feature_info_get(const FeatureInfo *fi, const char *id)
{
  assert(fi && id);
  return hashmap_get(fi->id_to_genome_node, id);
}

void feature_info_add(FeatureInfo *fi, const char *id, GT_GenomeNode *gf)
{
  assert(fi && id && gf);
  assert(!genome_feature_is_pseudo((GenomeFeature*) gf));
  hashmap_add(fi->id_to_genome_node, cstr_dup(id), genome_node_ref(gf));
}

GT_GenomeNode* feature_info_get_pseudo_parent(const FeatureInfo *fi,
                                           const char *id)
{
  assert(fi && id);
  return hashmap_get(fi->id_to_pseudo_parent, id);
}

void feature_info_add_pseudo_parent(FeatureInfo *fi, const char *id,
                                    GT_GenomeNode *pseudo_parent)
{
  assert(fi && id && pseudo_parent);
  assert(genome_feature_is_pseudo((GenomeFeature*) pseudo_parent));
  hashmap_add(fi->id_to_pseudo_parent, cstr_dup(id),
              genome_node_ref(pseudo_parent));
}

void feature_info_replace_pseudo_parent(FeatureInfo *fi, GT_GenomeNode *child,
                                        GT_GenomeNode *new_pseudo_parent)
{
  const char *id;
  assert(fi && child && new_pseudo_parent);
  assert(genome_feature_is_pseudo((GenomeFeature*) new_pseudo_parent));
  id = genome_feature_get_attribute(child, ID_STRING);
  assert(id);
  hashmap_remove(fi->id_to_pseudo_parent, id);
  feature_info_add_pseudo_parent(fi, id, new_pseudo_parent);
}

static GT_GenomeNode* find_root(const FeatureInfo *fi, const char *id)
{
  const char *delim, *parents;
  GT_GenomeNode *this_feature, *parent_pseudo_feature;
  assert(fi && id);
  /* get feature */
  delim = strchr(id, ';');
  if (delim) {
    char *first_parent = cstr_dup_nt(id, delim - id);
    this_feature = hashmap_get(fi->id_to_genome_node, first_parent);
    parent_pseudo_feature = hashmap_get(fi->id_to_pseudo_parent, first_parent);
    ma_free(first_parent);
  }
  else {
    this_feature = hashmap_get(fi->id_to_genome_node, id);
    parent_pseudo_feature = hashmap_get(fi->id_to_pseudo_parent, id);
  }
  assert(this_feature);
  /* recursion */
  parents = genome_feature_get_attribute(this_feature, PARENT_STRING);
  if (parents)
    return find_root(fi, parents);
  else if (parent_pseudo_feature)
    return parent_pseudo_feature;
  return this_feature;
}

GT_GenomeNode* feature_info_find_root(const FeatureInfo *fi, const char *id)
{
  assert(fi && id);
  assert(feature_info_get(fi, id));
  return find_root(fi, id);
}
