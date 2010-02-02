/*
  Copyright (c) 2007-2009 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2007      Malte Mader <mader@zbh.uni-hamburg.de>
  Copyright (c) 2007      Christin Schaerfer <schaerfer@zbh.uni-hamburg.de>
  Copyright (c) 2007-2009 Center for Bioinformatics, University of Hamburg

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

#include "annotationsketch/canvas.h"
#include "annotationsketch/canvas_cairo_file.h"
#include "annotationsketch/diagram.h"
#include "annotationsketch/feature_index_memory_api.h"
#include "annotationsketch/line_breaker_captions.h"
#include "annotationsketch/style.h"
#include "annotationsketch/track.h"
#include "core/basename_api.h"
#include "core/cstr.h"
#include "core/ensure.h"
#include "core/hashmap.h"
#include "core/ma.h"
#include "core/msort.h"
#include "core/log.h"
#include "core/str.h"
#include "core/thread.h"
#include "core/undef.h"
#include "core/unused_api.h"
#include "extended/feature_node.h"
#include "extended/genome_node.h"

/* used to index non-multiline-feature blocks
   undefined pointer -- never dereference! */
#define GT_UNDEF_REPR               (void*)~0
/* used to separate a filename from the type in a track name */
#define FILENAME_TYPE_SEPARATOR  '|'

struct GtDiagram {
  /* GtBlock lists indexed by track keys */
  GtHashmap *blocks;
  /* Reverse lookup structure (per node) */
  GtHashmap *nodeinfo;
  /* Cache tables for configuration data */
  GtHashmap *collapsingtypes, *caption_display_status, *groupedtypes;
  GtStyle *style;
  GtArray *features,
          *custom_tracks;
  GtRange range;
  void *ptr;
  GtTrackSelectorFunc select_func;
  GtRWLock *lock;
};

typedef enum {
  GT_DO_NOT_GROUP_BY_PARENT,
  GT_GROUP_BY_PARENT,
  GT_UNDEFINED_GROUPING
} GtShouldGroupByParent;

/* holds a GtBlock with associated type */
typedef struct {
  const char *gft;
  GtFeatureNode *rep;
  GtBlock *block;
} GtBlockTuple;

/* a node in the reverse lookup structure used for collapsing */
typedef struct {
  GtFeatureNode *parent;
  GtHashmap *type_index;
  GtStrArray *types;
} NodeInfoElement;

typedef struct {
  GtHashmap *rep_index;
  GtArray *blocktuples;
  bool must_merge;
} PerTypeInfo;

typedef struct {
  GtFeatureNode *parent;
  GtDiagram *diagram;
} NodeTraverseInfo;

static inline GtBlockTuple* blocktuple_new(const char *gft,
                                           GtFeatureNode *rep,
                                           GtBlock *block)
{
  GtBlockTuple *bt;
  gt_assert(block);
  bt = gt_calloc(1, sizeof (GtBlockTuple));
  bt->gft = gft;
  bt->rep = rep;
  bt->block = block;
  return bt;
}

static inline NodeInfoElement* nodeinfo_get(GtDiagram *d,
                                            GtFeatureNode *node)
{
  NodeInfoElement *ni;
  gt_assert(d && node);
  if (!(ni = gt_hashmap_get(d->nodeinfo, node))) {
    ni = gt_calloc(1, sizeof (NodeInfoElement));
    ni->type_index  = gt_hashmap_new(HASH_STRING, NULL,
                                     gt_free_func);
    ni->types       = gt_str_array_new();
    gt_hashmap_add(d->nodeinfo, node, ni);
  }
  return ni;
}

static inline GtBlock* nodeinfo_find_block(NodeInfoElement* ni,
                                           const char *gft,
                                           GtFeatureNode *gf)
{
  PerTypeInfo *type_struc = NULL;
  GtBlockTuple *bt = NULL;
  gt_assert(ni);
  if (!(type_struc = gt_hashmap_get(ni->type_index, gft)))
    return NULL;
  if (!(bt = gt_hashmap_get(type_struc->rep_index, gf)))
    return NULL;
  gt_assert(bt);
  return bt->block;
}

static inline void nodeinfo_add_block(NodeInfoElement *ni,
                                      const char *gft,
                                      GtFeatureNode *rep,
                                      GtBlock *block)
{
  GtBlockTuple *bt;
  PerTypeInfo *type_struc = NULL;
  gt_assert(ni);
  bt = blocktuple_new(gft, rep, block);
  if (!(ni->type_index))
  {
    ni->type_index  = gt_hashmap_new(HASH_STRING, NULL,
                                     gt_free_func);
  }
  if (!(type_struc = gt_hashmap_get(ni->type_index, gft)))
  {
    type_struc = gt_calloc(1, sizeof (PerTypeInfo));
    type_struc->rep_index = gt_hashmap_new(HASH_DIRECT, NULL, NULL);
    type_struc->blocktuples = gt_array_new(sizeof (GtBlockTuple*));
    gt_hashmap_add(ni->type_index, (char*) gft, type_struc);
    gt_str_array_add_cstr(ni->types, gft);
  }
  gt_hashmap_add(type_struc->rep_index, rep, bt);
  if (rep != GT_UNDEF_REPR)
    type_struc->must_merge = true;
  gt_array_add(type_struc->blocktuples, bt);
}

static inline const char* get_node_name_or_id(GtFeatureNode *gn)
{
  const char *ret;
  if (!gn) return NULL;
  if (!(ret = gt_feature_node_get_attribute(gn, "Name"))) {
    if (!(ret = gt_feature_node_get_attribute(gn, "ID")))
      ret = NULL;
  }
  return ret;
}

static inline bool get_caption_display_status(GtDiagram *d, const char *gft)
{
  bool *status;
  gt_assert(d && gft);
  status = (bool*) gt_hashmap_get(d->caption_display_status, gft);
  if (!status)
  {
    unsigned long threshold = GT_UNDEF_ULONG;
    double tmp;
    status = gt_malloc(sizeof (bool*));
    if (!gt_style_get_bool(d->style, "format", "show_block_captions", status,
                           NULL))
      *status = true;
    if (*status)
    {
      if (gt_style_get_num(d->style, gft, "max_capt_show_width", &tmp, NULL))
        threshold = tmp;
      if (threshold == GT_UNDEF_ULONG)
        *status = true;
      else
        *status = (gt_range_length(&d->range) <= threshold);
    }
    gt_hashmap_add(d->caption_display_status, (void*) gft, status);
  }
  return *status;
}

static inline void assign_block_caption(GtDiagram *d,
                                        GtFeatureNode *node,
                                        GtFeatureNode *parent,
                                        GtBlock *block)
{
  const char *nnid_p = NULL, *nnid_n = NULL;
  GtStr *caption = NULL;
  nnid_p = get_node_name_or_id(parent);
  nnid_n = get_node_name_or_id(node);
  if ((nnid_p || nnid_n) && get_caption_display_status(d,
                                          gt_feature_node_get_type(node)))
  {
    caption = gt_str_new();
    if (parent)
    {
      if (nnid_p && gt_genome_node_has_children((GtGenomeNode*) parent))
        gt_str_append_cstr(caption, nnid_p);
      else
        gt_str_append_cstr(caption, "-");
      gt_str_append_cstr(caption, "/");
    }
    if (nnid_n)
      gt_str_append_cstr(caption, nnid_n);
  }
  gt_block_set_caption(block, caption);
}

static void add_to_current(GtDiagram *d, GtFeatureNode *node,
                           GtFeatureNode *parent)
{
  GtBlock *block;
  NodeInfoElement *ni;
  GtStr *caption = NULL;
  const char *nnid_p = NULL, *nnid_n = NULL;
  gt_assert(d && node);
  /* Get nodeinfo element and set itself as parent */
  ni = nodeinfo_get(d, node);
  gt_log_log("adding %s to self", gt_feature_node_get_type(node));
  ni->parent = node;
  /* create new GtBlock tuple and add to node info */
  block = gt_block_new_from_node(node);
  caption = gt_str_new();
  if (!gt_style_get_str(d->style,
                        gt_feature_node_get_type(node),
                        "block_caption",
                        caption,
                        node))
  {
    nnid_p = get_node_name_or_id(parent);
    nnid_n = get_node_name_or_id(node);
    if ((nnid_p || nnid_n) && get_caption_display_status(d,
                                            gt_feature_node_get_type(node)))
    {
      if (parent) {
        if (nnid_p && gt_genome_node_has_children((GtGenomeNode*) parent))
          gt_str_append_cstr(caption, nnid_p);
        else
          gt_str_append_cstr(caption, "-");
        gt_str_append_cstr(caption, "/");
      }
      if (nnid_n)
        gt_str_append_cstr(caption, nnid_n);
    } else {
      gt_str_delete(caption);
      caption = NULL;
    }
  }
  gt_block_set_caption(block, caption);
  gt_block_insert_element(block, node);
  nodeinfo_add_block(ni,
                     gt_feature_node_get_type(node),
                     GT_UNDEF_REPR, block);
}

static void add_to_parent(GtDiagram *d, GtFeatureNode *node,
                          GtFeatureNode *parent)
{
  GtBlock *block = NULL;
  NodeInfoElement *par_ni, *ni;
  gt_assert(d && node);
  if (!parent) return;
  par_ni = nodeinfo_get(d, parent);
  ni = nodeinfo_get(d, node);
  gt_log_log("adding %s to parent %p", gt_feature_node_get_type(node), parent);
  ni->parent = parent;
  block = nodeinfo_find_block(par_ni,
                              gt_feature_node_get_type(node),
                              parent);
  if (!block) {
    block = gt_block_new_from_node(parent);
    gt_block_set_type(block, gt_feature_node_get_type(node));
    assign_block_caption(d, node, parent, block);
    nodeinfo_add_block(par_ni,
                     gt_feature_node_get_type((GtFeatureNode*) node),
                     parent,
                     block);
  }
  gt_assert(block);
  gt_block_insert_element(block, node);
}

static void add_to_rep(GtDiagram *d, GtFeatureNode *node,
                       GtFeatureNode* parent)
{
  GtBlock *block = NULL;
  GtFeatureNode *rep = GT_UNDEF_REPR;
  NodeInfoElement *ni;
  gt_assert(d && node && gt_feature_node_is_multi(node));
  rep = gt_feature_node_get_multi_representative(node);
  gt_log_log("adding %s to representative %p", gt_feature_node_get_type(node),
                                               rep);
  ni = nodeinfo_get(d, rep);

  block = nodeinfo_find_block(ni,
                              gt_feature_node_get_type(node),
                              rep);
  if (!block) {
    block = gt_block_new_from_node(parent);
    gt_block_set_type(block, gt_feature_node_get_type(node));
    /* if parent is a pseudonode, then we have a multiline feature without
       a parent. we must not access the parent in this case! */
    if (gt_feature_node_is_pseudo(parent))
      assign_block_caption(d, node, NULL, block);
    else
      assign_block_caption(d, node, parent, block);
    nodeinfo_add_block(ni, gt_feature_node_get_type(node),
                       rep, block);
  }
  gt_assert(block);
  gt_block_insert_element(block, node);
}

static void add_recursive(GtDiagram *d, GtFeatureNode *node,
                          GtFeatureNode* parent,
                          GtFeatureNode *original_node)
{
  NodeInfoElement *ni;
  GtFeatureNode *rep = GT_UNDEF_REPR;
  gt_assert(d && node && original_node);
  if (!parent) return;
  ni = nodeinfo_get(d, node);
  if (gt_feature_node_is_multi(original_node))
  {
    rep = gt_feature_node_get_multi_representative(original_node);
  }
    /* end of recursion, insert into target block */
  if (parent == node) {
    GtBlock *block ;
    block = nodeinfo_find_block(ni,
                                gt_feature_node_get_type(node),
                                rep);
    if (!block) {
      block = gt_block_new_from_node(node);
      nodeinfo_add_block(ni,
                         gt_feature_node_get_type(node),
                         rep,
                         block);
    }
    gt_block_insert_element(block, original_node);
    gt_log_log("add %s to target %s", gt_feature_node_get_type(original_node),
                                      gt_block_get_type(block));
  }
  else {
    /* not at target type block yet, set up reverse entry and follow */
    NodeInfoElement *parent_ni;
    /* set up reverse entry */
    ni->parent = parent;
    parent_ni = gt_hashmap_get(d->nodeinfo, parent);
    if (parent_ni) {
      gt_log_log("recursion: %s -> %s", gt_feature_node_get_type(node),
                                        gt_feature_node_get_type(parent));
      add_recursive(d, parent, parent_ni->parent, original_node);
    }
  }
}

static void process_node(GtDiagram *d, GtFeatureNode *node,
                         GtFeatureNode *parent)
{
  GtRange elem_range;
  bool *collapse;
  GtShouldGroupByParent *group;
  const char *feature_type = NULL,
             *parent_gft = NULL;
  double tmp;
  unsigned long max_show_width = GT_UNDEF_ULONG,
                par_max_show_width = GT_UNDEF_ULONG;

  gt_assert(d && node);

  gt_log_log(">> getting '%s'", gt_feature_node_get_type(node));

  /* skip pseudonodes */
  if (gt_feature_node_is_pseudo(node))
    return;

  feature_type = gt_feature_node_get_type(node);
  gt_assert(feature_type);

  /* discard elements that do not overlap with visible range */
  elem_range = gt_genome_node_get_range((GtGenomeNode*) node);
  if (!gt_range_overlap(&d->range, &elem_range))
    return;

  /* get maximal view widths in nucleotides to show this type */
  if (gt_style_get_num(d->style, feature_type, "max_show_width", &tmp, NULL))
    max_show_width = tmp;
  else
    max_show_width = GT_UNDEF_ULONG;

  /* for non-root nodes, get maximal view with to show parent */
  if (parent)
  {
    if (!gt_feature_node_is_pseudo(parent))
    {
      parent_gft = gt_feature_node_get_type(parent);
      if (gt_style_get_num(d->style, parent_gft, "max_show_width", &tmp, NULL))
        par_max_show_width = tmp;
      else
        par_max_show_width = GT_UNDEF_ULONG;
    } else par_max_show_width = GT_UNDEF_ULONG;
  }

  /* check if this type is to be displayed at all */
  if (max_show_width != GT_UNDEF_ULONG &&
      gt_range_length(&d->range) > max_show_width)
  {
    return;
  }

  /* disregard parent node if it is configured not to be shown */
  if (parent
        && par_max_show_width != GT_UNDEF_ULONG
        && gt_range_length(&d->range) > par_max_show_width)
  {
    parent = NULL;
  }

  /* check if this is a collapsing type, cache result */
  if ((collapse = (bool*) gt_hashmap_get(d->collapsingtypes,
                                         feature_type)) == NULL)
  {
    collapse = gt_malloc(sizeof (bool));
    if (!gt_style_get_bool(d->style, feature_type, "collapse_to_parent",
                           collapse, NULL))
      *collapse = false;
    gt_hashmap_add(d->collapsingtypes, (void*) feature_type, collapse);
  }

  /* check if type should be grouped by parent, cache result */
  if ((group = (GtShouldGroupByParent*) gt_hashmap_get(d->groupedtypes,
                                                       feature_type)) == NULL)
  {
    bool tmp;
    group = gt_malloc(sizeof (GtShouldGroupByParent));
    if (!gt_style_get_bool(d->style, feature_type, "group_by_parent",
                           &tmp, NULL)) {
      *group = GT_UNDEFINED_GROUPING;
    } else {
      if (tmp)
        *group = GT_GROUP_BY_PARENT;
      else
        *group = GT_DO_NOT_GROUP_BY_PARENT;
    }
    gt_hashmap_add(d->groupedtypes, (void*) feature_type, group);
  }

  /* decide where to place this feature: */
  if (*collapse)
  {
    /* user has specified collapsing to parent for this type */
    if (parent && !gt_feature_node_is_pseudo(parent)) {
      /* collapsing child nodes are added to upwards blocks,
         but never collapse into pseudo nodes */
      add_recursive(d, node, parent, node);
    } else {
      /* if no parent or only pseudo-parent, do not collapse */
      add_to_current(d, node, parent);
    }
  }
  else  /* (!*collapse) */
  {
    if (parent) {
      bool do_not_overlap = false;
      do_not_overlap =
        gt_genome_node_direct_children_do_not_overlap_st((GtGenomeNode*) parent,
                                                         (GtGenomeNode*) node);
      if (*group == GT_GROUP_BY_PARENT
          || (do_not_overlap && *group == GT_UNDEFINED_GROUPING))
      {
        if (gt_feature_node_is_pseudo(parent)
              && gt_feature_node_is_multi(node))
        {
          add_to_rep(d, node, parent);
        } else if
            (gt_genome_node_number_of_children((GtGenomeNode*) parent) > 1)
        {
          add_to_parent(d, node, parent);
        } else {
          add_to_current(d, node, parent);
        }
      } else {
        add_to_current(d, node, parent);
      }
    } else {
      /* root nodes always get their own block */
      add_to_current(d, node, parent);
    }
  }

  /* we can now assume that this node (or its representative)
     has been processed into the reverse lookup structure */
  if (gt_feature_node_is_multi(node))
  {
    GtFeatureNode *rep;
    rep = gt_feature_node_get_multi_representative((GtFeatureNode*) node);
    gt_assert(gt_hashmap_get(d->nodeinfo, rep));
  }
  else
    gt_assert(gt_hashmap_get(d->nodeinfo, node));
}

static int visit_child(GtGenomeNode* gn, void *nti,
                       GtError *err)
{
  NodeTraverseInfo* gt_genome_node_info;
   int had_err;
  gt_genome_node_info = (NodeTraverseInfo*) nti;
  gt_error_check(err);

  if (gt_genome_node_has_children(gn))
  {
    GtFeatureNode *oldparent = gt_genome_node_info->parent;
    process_node(gt_genome_node_info->diagram, (GtFeatureNode*) gn,
                 gt_genome_node_info->parent);
    gt_genome_node_info->parent = (GtFeatureNode*) gn;
    had_err = gt_genome_node_traverse_direct_children(gn,
                                                      gt_genome_node_info,
                                                      visit_child,
                                                      err);
    gt_assert(!had_err); /* visit_child() is sane */
    gt_genome_node_info->parent = oldparent;
  }
  else
    process_node(gt_genome_node_info->diagram, (GtFeatureNode*)gn,
                 gt_genome_node_info->parent);
  return 0;
}

static void default_track_selector(GtBlock *block, GtStr *result,
                                   GT_UNUSED void *data)
{
  GtGenomeNode *top;
  char *basename;
  gt_assert(block && result);
  gt_str_reset(result);
  top = (GtGenomeNode*) gt_block_get_top_level_feature(block);
  /* we take the basename of the filename to have nicer output in the
     generated graphic. this might lead to ``collapsed'' tracks, if two files
     with different paths have the same basename. */
  basename = gt_basename(gt_genome_node_get_filename(top));
  gt_str_append_cstr(result, basename);
  gt_free(basename);
  gt_str_append_char(result, FILENAME_TYPE_SEPARATOR);
  gt_str_append_cstr(result, gt_block_get_type(block));
}

/* Create lists of all GtBlocks in the diagram. */
static int collect_blocks(GT_UNUSED void *key, void *value, void *data,
                          GT_UNUSED GtError *err)
{
  NodeInfoElement *ni = (NodeInfoElement*) value;
  GtDiagram *diagram = (GtDiagram*) data;
  GtBlock *block = NULL;
  GtStr *trackid_str;
  unsigned long i = 0;
  trackid_str = gt_str_new();
  for (i = 0; i < gt_str_array_size(ni->types); i++) {
    const char *type;
    unsigned long j;
    GtArray *list;
    PerTypeInfo *type_struc = NULL;
    GtBlock* mainblock = NULL;
    type = gt_str_array_get(ni->types, i);
    type_struc = gt_hashmap_get(ni->type_index, type);
    gt_assert(type_struc);
    for (j=0; j<gt_array_size(type_struc->blocktuples); j++)
    {
      GtBlockTuple *bt;
      bt  = *(GtBlockTuple**) gt_array_get(type_struc->blocktuples, j);
      if (bt->rep == GT_UNDEF_REPR && type_struc->must_merge)
      {
        block = mainblock = gt_block_ref(bt->block);
        gt_block_delete(mainblock);
        gt_free(bt);
        continue;
      }
      else
      {
        if (mainblock)
        {
          block = gt_block_clone(mainblock);
          gt_block_merge(block, bt->block);
          gt_block_delete(bt->block);
        } else block = bt->block;
      }
      gt_assert(block);
      gt_str_reset(trackid_str);
      /* execute hook for track selector function */
      diagram->select_func(block, trackid_str, diagram->ptr);

      if (!(list = (GtArray*) gt_hashmap_get(diagram->blocks,
                                             gt_str_get(trackid_str))))
      {
        list = gt_array_new(sizeof (GtBlock*));
        gt_hashmap_add(diagram->blocks, gt_cstr_dup(gt_str_get(trackid_str)),
                       list);
      };
      gt_assert(list);
      gt_array_add(list, block);
      gt_free(bt);
    }
    gt_array_delete(type_struc->blocktuples);
    gt_hashmap_delete(type_struc->rep_index);
    gt_block_delete(mainblock);
  }
  gt_hashmap_delete(ni->type_index);
  gt_str_array_delete(ni->types);
  gt_free(ni);
  gt_str_delete(trackid_str);
  return 0;
}

/* Traverse a genome node graph with depth first search. */
static void traverse_genome_nodes(GtFeatureNode *gn,
                                  void *nti)
{
  NodeTraverseInfo* gt_genome_node_info;
  int had_err;
  gt_assert(nti);
  gt_genome_node_info = (NodeTraverseInfo*) nti;
  gt_genome_node_info->parent = gn;
  /* handle root nodes */
  process_node(gt_genome_node_info->diagram, (GtFeatureNode*)gn, NULL);
  if (gt_genome_node_has_children((GtGenomeNode*) gn)) {
    had_err = gt_genome_node_traverse_direct_children((GtGenomeNode*)gn,
                                                      gt_genome_node_info,
                                                      visit_child, NULL);
    gt_assert(!had_err); /* visit_child() is sane */
  }
}

static void blocklist_delete(void *value)
{
  unsigned long i;
  GtArray *a;
  if (!value) return;
  a = (GtArray*) value;
  for (i = 0; i < gt_array_size(a); i++)
    gt_block_delete(*(GtBlock**) gt_array_get(a, i));
  gt_array_delete(a);
}

int gt_diagram_build(GtDiagram *diagram)
{
  unsigned long i = 0;
  int had_err = 0;
  NodeTraverseInfo nti;
  gt_assert(diagram);

  nti.diagram = diagram;
  /* clear caches */
  gt_hashmap_reset(diagram->collapsingtypes);
  gt_hashmap_reset(diagram->groupedtypes);
  gt_hashmap_reset(diagram->caption_display_status);

  if (!diagram->blocks)
  {
    gt_hashmap_reset(diagram->nodeinfo);
    /* do node traversal for each root feature */
    for (i = 0; i < gt_array_size(diagram->features); i++)
    {
      GtFeatureNode *current_root = *(GtFeatureNode**)
                                           gt_array_get(diagram->features,i);
      traverse_genome_nodes(current_root, &nti);
    }
    diagram->blocks = gt_hashmap_new(HASH_STRING, gt_free_func,
                                    (GtFree) blocklist_delete);
    /* collect blocks from nodeinfo structures */
    had_err = gt_hashmap_foreach_ordered(diagram->nodeinfo,
                                         collect_blocks,
                                         diagram,
                                         (GtCompare) gt_genome_node_cmp,
                                         NULL);
    gt_assert(!had_err); /* collect_blocks() is sane */
  }

  return had_err;
}

static GtDiagram* gt_diagram_new_generic(GtArray *features,
                                         const GtRange *range,
                                         GtStyle *style,
                                         bool ref_features)
{
  GtDiagram *diagram;
  diagram = gt_calloc(1, sizeof (GtDiagram));
  diagram->nodeinfo = gt_hashmap_new(HASH_DIRECT, NULL, NULL);
  diagram->style = style;
  diagram->lock = gt_rwlock_new();
  diagram->range = *range;
  if (ref_features)
    diagram->features = gt_array_ref(features);
  else
    diagram->features = features;
  diagram->select_func = default_track_selector;
  diagram->custom_tracks = gt_array_new(sizeof (GtCustomTrack*));
  /* init caches */
  diagram->collapsingtypes = gt_hashmap_new(HASH_STRING, NULL, gt_free_func);
  diagram->groupedtypes = gt_hashmap_new(HASH_STRING, NULL, gt_free_func);
  diagram->caption_display_status = gt_hashmap_new(HASH_DIRECT, NULL,
                                                   gt_free_func);
  return diagram;
}

GtDiagram* gt_diagram_new(GtFeatureIndex *feature_index, const char *seqid,
                          const GtRange *range, GtStyle *style,
                          GtError *err)
{
  GtDiagram *diagram;
  int had_err = 0;
  GtArray *features = NULL;
  gt_assert(seqid && range && style);
  if (range->start == range->end)
  {
    gt_error_set(err, "range start must not be equal to range end");
    return NULL;
  }
  features = gt_array_new(sizeof (GtGenomeNode*));
  had_err = gt_feature_index_get_features_for_range(feature_index, features,
                                                    seqid, range, err);
  if (had_err)
  {
    gt_array_delete(features);
    return NULL;
  }
  diagram = gt_diagram_new_generic(features, range, style, false);
  return diagram;
}

GtDiagram* gt_diagram_new_from_array(GtArray *features, const GtRange *range,
                                     GtStyle *style)
{
  gt_assert(features && range && style);
  return gt_diagram_new_generic(features, range, style, true);
}

GtRange gt_diagram_get_range(const GtDiagram *diagram)
{
  GtRange rng;
  gt_assert(diagram);
  gt_rwlock_rdlock(diagram->lock);
  rng = diagram->range;
  gt_rwlock_unlock(diagram->lock);
  return rng;
}

void gt_diagram_set_track_selector_func(GtDiagram *diagram,
                                        GtTrackSelectorFunc bsfunc,
                                        void *ptr)
{
  gt_assert(diagram);
  gt_rwlock_wrlock(diagram->lock);
  /* register selector function and attached pointer */
  diagram->select_func = bsfunc;
  diagram->ptr = ptr;
  /* this could change track assignment -> discard current blocks and requeue */
  gt_hashmap_delete(diagram->blocks);
  diagram->blocks = NULL;
  gt_rwlock_unlock(diagram->lock);
}

void gt_diagram_reset_track_selector_func(GtDiagram *diagram)
{
  gt_assert(diagram);
  gt_rwlock_wrlock(diagram->lock);
  diagram->select_func = default_track_selector;
  gt_hashmap_delete(diagram->blocks);
  diagram->blocks = NULL;
  gt_rwlock_unlock(diagram->lock);
}

GtHashmap* gt_diagram_get_blocks(const GtDiagram *diagram)
{
  GtHashmap *ret;
  gt_assert(diagram);
  gt_rwlock_wrlock(diagram->lock);
  (void) gt_diagram_build((GtDiagram*) diagram);
  ret = diagram->blocks;
  gt_rwlock_unlock(diagram->lock);
  return ret;
}

GtArray* gt_diagram_get_custom_tracks(const GtDiagram *diagram)
{
  GtArray *ret;
  gt_assert(diagram);
  gt_rwlock_rdlock(diagram->lock);
  ret = diagram->custom_tracks;
  gt_rwlock_unlock(diagram->lock);
  return ret;
}

void gt_diagram_add_custom_track(GtDiagram *diagram, GtCustomTrack* ctrack)
{
  gt_assert(diagram && ctrack);
  gt_rwlock_wrlock(diagram->lock);
  gt_array_add(diagram->custom_tracks, ctrack);
  gt_rwlock_unlock(diagram->lock);
}

typedef struct {
  GtFeatureIndex *fi;
  GtError *err;
  GtDiagram *d;
  GtStyle *sty;
  int errstatus;
} GtDiagramTestShared;

void* gt_diagram_unit_test_sketch_func(void *data)
{
  int had_err = 0;
  GtLayout *l = NULL;
  GtCanvas *c = NULL;
  GtDiagramTestShared *sh = (GtDiagramTestShared*) data;

  l = gt_layout_new(sh->d, 1000, sh->sty, sh->err);
  if (!l)
    had_err = -1;
  if (!had_err) {
    c = gt_canvas_cairo_file_new(sh->sty, GT_GRAPHICS_PNG, 1000,
                                 gt_layout_get_height(l), NULL);
    if (!c)
      had_err = -1;
  }
  if (!had_err) {
    had_err = gt_layout_sketch(l, c, sh->err);
  }
  if (had_err)
    sh->errstatus = 1;
  gt_layout_delete(l);
  gt_canvas_delete(c);
  return NULL;
}

int gt_diagram_unit_test(GtError *err)
{
  int had_err = 0;
  GtGenomeNode *gn;
  GtDiagramTestShared sh;
  GtRange testrng = {100, 10000};
  gt_error_check(err);

  gn = gt_feature_node_new_standard_gene();
  sh.fi = gt_feature_index_memory_new();
  sh.sty = gt_style_new(err);
  sh.err = err;
  sh.errstatus = 0;
  gt_feature_index_add_feature_node(sh.fi, gt_feature_node_cast(gn));
  gt_genome_node_delete(gn);
  sh.d = gt_diagram_new(sh.fi, "ctg123", &testrng, sh.sty, err);

  gt_multithread(gt_diagram_unit_test_sketch_func, &sh, err);
  ensure(had_err, sh.errstatus == 0);

  gt_style_delete(sh.sty);
  gt_diagram_delete(sh.d);
  gt_feature_index_delete(sh.fi);

  return had_err;
}

void gt_diagram_delete(GtDiagram *diagram)
{
  if (!diagram) return;
  gt_rwlock_wrlock(diagram->lock);
  gt_array_delete(diagram->features);
  if (diagram->blocks)
    gt_hashmap_delete(diagram->blocks);
  gt_hashmap_delete(diagram->nodeinfo);
  gt_hashmap_delete(diagram->collapsingtypes);
  gt_hashmap_delete(diagram->groupedtypes);
  gt_hashmap_delete(diagram->caption_display_status);
  gt_array_delete(diagram->custom_tracks);
  gt_rwlock_unlock(diagram->lock);
  gt_rwlock_delete(diagram->lock);
  gt_free(diagram);
}
