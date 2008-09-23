/*
  Copyright (c) 2007-2008 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2007      Malte Mader <mmader@zbh.uni-hamburg.de>
  Copyright (c) 2007      Christin Schaerfer <cschaerfer@zbh.uni-hamburg.de>
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

#include "annotationsketch/canvas.h"
#include "annotationsketch/canvas_cairo_file.h"
#include "annotationsketch/diagram.h"
#include "annotationsketch/feature_index_memory_api.h"
#include "annotationsketch/line_breaker_captions.h"
#include "annotationsketch/style.h"
#include "annotationsketch/track.h"
#include "core/basename.h"
#include "core/cstr.h"
#include "core/ensure.h"
#include "core/hashmap.h"
#include "core/ma.h"
#include "core/msort.h"
#include "core/str.h"
#include "core/undef.h"
#include "core/unused_api.h"
#include "core/warning_api.h"
#include "extended/feature_node.h"
#include "extended/genome_node.h"

/* used to separate a filename from the type in a track name */
#define FILENAME_TYPE_SEPARATOR  '|'

/* used to index non-multiline-feature blocks */
#define UNDEF_REPR               (void*)"UNDEF"

struct GtDiagram {
  /* GT_Tracks indexed by track keys */
  GtHashmap *tracks;
  /* GtBlock lists indexed by track keys */
  GtHashmap *blocks;
  /* Reverse lookup structure (per node) */
  GtHashmap *nodeinfo;
  /* Cache tables for configuration data */
  GtHashmap *collapsingtypes, *caption_display_status;
  int nof_tracks;
  GtStyle *style;
  GtRange range;
};

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

typedef struct {
  GtCanvas *canvas;
  GtDiagram *dia;
} GtTrackTraverseInfo;

static GtBlockTuple* blocktuple_new(const char *gft,
                                     GtFeatureNode *rep,
                                     GtBlock *block)
{
  GtBlockTuple *bt;
  assert(block);
  bt = gt_calloc(1, sizeof (GtBlockTuple));
  bt->gft = gft;
  bt->rep = rep;
  bt->block = block;
  return bt;
}

static NodeInfoElement* nodeinfo_get(GtDiagram *d,
                                     GtFeatureNode *node)
{
  NodeInfoElement *ni;
  assert(d && node);
  if (!(ni = gt_hashmap_get(d->nodeinfo, node))) {
    ni = gt_calloc(1, sizeof (NodeInfoElement));
    ni->type_index  = gt_hashmap_new(HASH_STRING, NULL,
                                  gt_free_func);
    ni->types       = gt_strarray_new();
    gt_hashmap_add(d->nodeinfo, node, ni);
  }
  return ni;
}

static GtBlock* nodeinfo_find_block(NodeInfoElement* ni,
                                    const char *gft,
                                    GtFeatureNode *gf)
{
  PerTypeInfo *type_struc = NULL;
  GtBlockTuple *bt = NULL;
  assert(ni);
  if (!(type_struc = gt_hashmap_get(ni->type_index, gft)))
    return NULL;
  if (!(bt = gt_hashmap_get(type_struc->rep_index, gf)))
    return NULL;
  assert(bt);
  return bt->block;
}

static void nodeinfo_add_block(NodeInfoElement *ni,
                               const char *gft,
                               GtFeatureNode *rep,
                               GtBlock *block)
{
  GtBlockTuple *bt;
  PerTypeInfo *type_struc = NULL;
  assert(ni && !nodeinfo_find_block(ni, gft, rep));
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
    gt_strarray_add_cstr(ni->types, gft);
  }
  gt_hashmap_add(type_struc->rep_index, rep, bt);
  if (rep != UNDEF_REPR)
    type_struc->must_merge = true;
  gt_array_add(type_struc->blocktuples, bt);
}

static const char* get_node_name_or_id(GtFeatureNode *gn)
{
  const char *ret;
  if (!gn) return NULL;
  if (!(ret = gt_feature_node_get_attribute(gn, "Name"))) {
    if (!(ret = gt_feature_node_get_attribute(gn, "ID")))
      ret = NULL;
  }
  return ret;
}

static bool get_caption_display_status(GtDiagram *d, const char *gft)
{
  assert(d && gft);
  bool *status;
  status = (bool*) gt_hashmap_get(d->caption_display_status, gft);
  if (!status)
  {
    unsigned long threshold;
    double tmp;
    status = gt_malloc(sizeof (bool*));
    if (!gt_style_get_bool(d->style, "format", "show_block_captions", status,
                           NULL))
      *status = true;
    if (*status)
    {
      if (gt_style_get_num(d->style, gft, "max_capt_show_width", &tmp, NULL))
        threshold = tmp;
      else
        threshold = UNDEF_ULONG;
      if (threshold == UNDEF_ULONG)
        *status = true;
      else
        *status = (gt_range_length(d->range) <= threshold);
    }
    gt_hashmap_add(d->caption_display_status, (void*) gft, status);
  }
  return *status;
}

static void assign_block_caption(GtDiagram *d,
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
    caption = gt_str_new_cstr("");
    if (parent) {
      if (gt_genome_node_has_children((GtGenomeNode*) parent))
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
  GtStr *caption;
  const char *nnid_p = NULL, *nnid_n = NULL;
  assert(d && node);
  /* Get nodeinfo element and set itself as parent */
  ni = nodeinfo_get(d, node);
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
        if (gt_genome_node_has_children((GtGenomeNode*) parent))
          gt_str_append_cstr(caption, nnid_p);
        else
          gt_str_append_cstr(caption, "-");
        gt_str_append_cstr(caption, "/");
      }
      if (nnid_n)
        gt_str_append_cstr(caption, nnid_n);
    }
  }
  gt_block_set_caption(block, caption);
  gt_block_insert_element(block, node);
  nodeinfo_add_block(ni,
                     gt_feature_node_get_type(node),
                     UNDEF_REPR, block);
}

static void add_to_parent(GtDiagram *d, GtFeatureNode *node,
                          GtFeatureNode *parent)
{
  GtBlock *block = NULL;
  NodeInfoElement *par_ni, *ni;
  assert(d && node);
  if (!parent) return;
  par_ni = nodeinfo_get(d, parent);
  ni = nodeinfo_get(d, node);
  ni->parent = parent;
  block = nodeinfo_find_block(par_ni,
                              gt_feature_node_get_type(node),
                              parent);
  if (!block) {
    block = gt_block_new_from_node(parent);
    assign_block_caption(d, node, parent, block);
    nodeinfo_add_block(par_ni,
                     gt_feature_node_get_type((GtFeatureNode*) node),
                     parent,
                     block);
  }
  assert(block);
  gt_block_insert_element(block, node);
}

static void add_to_rep(GtDiagram *d, GtFeatureNode *node,
                       GtFeatureNode* parent)
{
  GtBlock *block = NULL;
  GtFeatureNode *rep = UNDEF_REPR;
  NodeInfoElement *ni;
  assert(d && node && gt_feature_node_is_multi(node));
  rep = gt_feature_node_get_multi_representative(node);
  ni = nodeinfo_get(d, rep);

  block = nodeinfo_find_block(ni,
                              gt_feature_node_get_type(node),
                              rep);
  if (!block) {
    block = gt_block_new_from_node(parent);
    assign_block_caption(d, node, parent, block);
    nodeinfo_add_block(ni, gt_feature_node_get_type(node),
                     rep, block);
  }
  assert(block);
  gt_block_insert_element(block, node);
}

static void add_recursive(GtDiagram *d, GtFeatureNode *node,
                          GtFeatureNode* parent,
                          GtFeatureNode *original_node)
{
  NodeInfoElement *ni;
  GtFeatureNode *rep = UNDEF_REPR;
  assert(d && node && original_node);
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
  }
  else {
    /* not at target type block yet, set up reverse entry and follow */
    NodeInfoElement *parent_ni;
    /* set up reverse entry */
    ni->parent = parent;
    parent_ni = gt_hashmap_get(d->nodeinfo, parent);
    if (parent_ni)
      add_recursive(d, parent, parent_ni->parent, original_node);
  }
}

static void process_node(GtDiagram *d, GtFeatureNode *node,
                         GtFeatureNode *parent)
{
  GtRange elem_range;
  bool *collapse;
  bool do_not_overlap=false;
  const char *feature_type = NULL, *parent_gft = NULL;
  double tmp;
  unsigned long max_show_width = UNDEF_ULONG,
                par_max_show_width = UNDEF_ULONG;

  assert(d && node);

  feature_type = gt_feature_node_get_type(node);
  if (parent)
    parent_gft = gt_feature_node_get_type(parent);

  /* discard elements that do not overlap with visible range */
  elem_range = gt_genome_node_get_range((GtGenomeNode*) node);
  if (!gt_range_overlap(d->range, elem_range))
    return;

  /* get maximal view widths in nucleotides to show this type */
  if (gt_style_get_num(d->style, feature_type, "max_show_width", &tmp, NULL))
    max_show_width = tmp;
  else
    max_show_width = UNDEF_ULONG;
  if (parent)
  {
    if (gt_style_get_num(d->style, parent_gft, "max_show_width", &tmp, NULL))
    par_max_show_width = tmp;
  else
    par_max_show_width = UNDEF_ULONG;

  }
  /* check if this type is to be displayed */
  if (max_show_width != UNDEF_ULONG &&
      gt_range_length(d->range) > max_show_width) {
    return;
  }
  if (parent && par_max_show_width != UNDEF_ULONG
        && gt_range_length(d->range) > par_max_show_width)
    parent = NULL;

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

   /* check if direct children overlap */
  if (parent)
    do_not_overlap =
      gt_genome_node_direct_children_do_not_overlap_st((GtGenomeNode*) parent,
                                                       (GtGenomeNode*) node);

  /* decide how to continue: */
  if (*collapse && parent)
  {
    add_recursive(d, node, parent, node);
  }
  else if (!(*collapse)
             && gt_feature_node_is_multi(node))
  {
    add_to_rep(d, node, parent);
  }
  else if (!(*collapse)
             && do_not_overlap
             && gt_genome_node_number_of_children((GtGenomeNode*) parent) > 1)
  {
    add_to_parent(d, node, parent);
  }
  else
  {
    add_to_current(d, node, parent);
  }

  /* we can now assume that this node has been processed into the reverse
     lookup structure */
  if (gt_feature_node_is_multi(node))
  {
    GtFeatureNode *rep;
    rep = gt_feature_node_get_multi_representative((GtFeatureNode*) node);
    assert(gt_hashmap_get(d->nodeinfo, rep));
  }
  else
    assert(gt_hashmap_get(d->nodeinfo, node));
}

static int gt_diagram_add_tracklines(GT_UNUSED void *key, void *value,
                                     void *data, GT_UNUSED GtError *err)
{
  GtTracklineInfo *add = (GtTracklineInfo*) data;
  add->total_lines += gt_track_get_number_of_lines((GtTrack*) value);
  add->total_captionlines += gt_track_get_number_of_lines_with_captions(
                                                             (GtTrack*) value);
  return 0;
}

static int visit_child(GtGenomeNode* gn, void *gt_genome_node_children,
                       GtError *err)
{
  NodeTraverseInfo* gt_genome_node_info;
  gt_genome_node_info = (NodeTraverseInfo*) gt_genome_node_children;
  gt_error_check(err);
  int had_err;
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
    assert(!had_err); /* visit_child() is sane */
    gt_genome_node_info->parent = oldparent;
  }
  else
    process_node(gt_genome_node_info->diagram, (GtFeatureNode*)gn,
                 gt_genome_node_info->parent);
  return 0;
}

static GtStr* gt_track_key_new(const char *filename, const char *type)
{
  GtStr *gt_track_key;
  gt_track_key = gt_str_new_cstr(filename);
  gt_str_append_char(gt_track_key, FILENAME_TYPE_SEPARATOR);
  gt_str_append_cstr(gt_track_key, type);
  return gt_track_key;
}

/* Create lists of all GtBlocks in the diagram. */
static int collect_blocks(GT_UNUSED void *key, void *value, void *data,
                          GT_UNUSED GtError *err)
{
  NodeInfoElement *ni = (NodeInfoElement*) value;
  GtDiagram *diagram = (GtDiagram*) data;
  GtBlock *block = NULL;
  unsigned long i = 0;
  for (i = 0; i < gt_strarray_size(ni->types); i++) {
    const char *type;
    unsigned long j;
    GtArray *list;
    PerTypeInfo *type_struc = NULL;
    GtBlock* mainblock = NULL;
    type = gt_strarray_get(ni->types, i);
    type_struc = gt_hashmap_get(ni->type_index, type);
    assert(type_struc);
    for (j=0; j<gt_array_size(type_struc->blocktuples); j++)
    {
      GtBlockTuple *bt;
      bt  = *(GtBlockTuple**) gt_array_get(type_struc->blocktuples, j);
      if (bt->rep == UNDEF_REPR && type_struc->must_merge)
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
        } else
            block = bt->block;
      }
      assert(block);
      list = (GtArray*) gt_hashmap_get(diagram->blocks, bt->gft);
      if (!list)
      {
        list = gt_array_new(sizeof (GtBlock*));
        gt_hashmap_add(diagram->blocks, (void*) bt->gft, list);
      }
      assert(list);
      gt_array_add(list, block);
      gt_free(bt);
    }
    gt_array_delete(type_struc->blocktuples);
    gt_hashmap_delete(type_struc->rep_index);
    gt_block_delete(mainblock);
  }
  gt_hashmap_delete(ni->type_index);
  gt_strarray_delete(ni->types);
  gt_free(ni);
  return 0;
}

/* Traverse a genome node graph with depth first search. */
static void traverse_genome_nodes(GtFeatureNode *gn,
                                  void *gt_genome_node_children)
{
  NodeTraverseInfo* gt_genome_node_info;
  int had_err;
  assert(gt_genome_node_children);
  gt_genome_node_info = (NodeTraverseInfo*) gt_genome_node_children;
  gt_genome_node_info->parent = gn;
  /* handle root nodes */
  process_node(gt_genome_node_info->diagram, (GtFeatureNode*)gn, NULL);
  if (gt_genome_node_has_children((GtGenomeNode*) gn)) {
    had_err = gt_genome_node_traverse_direct_children((GtGenomeNode*)gn,
                                                      gt_genome_node_info,
                                                      visit_child, NULL);
    assert(!had_err); /* visit_child() is sane */
  }
}

static void gt_diagram_build(GtDiagram *diagram, GtArray *features)
{
  unsigned long i = 0;
  int had_err;
  NodeTraverseInfo gt_genome_node_children;
  gt_genome_node_children.diagram = diagram;

  /* initialise caches */
  diagram->collapsingtypes = gt_hashmap_new(HASH_STRING, NULL, gt_free_func);
  diagram->caption_display_status = gt_hashmap_new(HASH_DIRECT,
                                                  NULL, gt_free_func);

  /* do node traversal for each root feature */
  for (i = 0; i < gt_array_size(features); i++) {
    GtFeatureNode *current_root = *(GtFeatureNode**)
                                                 gt_array_get(features,i);
    traverse_genome_nodes(current_root, &gt_genome_node_children);
  }
  /* collect blocks from nodeinfo structures and create the tracks */
  had_err = gt_hashmap_foreach_ordered(diagram->nodeinfo, collect_blocks,
                                       diagram, (GtCompare) gt_genome_node_cmp,
                                       NULL);
  assert(!had_err); /* collect_blocks() is sane */

  /* clear caches */
  gt_hashmap_delete(diagram->collapsingtypes);
  gt_hashmap_delete(diagram->caption_display_status);
}

static int blocklist_delete(void *value)
{
  unsigned long i;
  GtArray *a = (GtArray*) value;
  for (i = 0; i < gt_array_size(a); i++)
    gt_block_delete(*(GtBlock**) gt_array_get(a, i));
  gt_array_delete(a);
  return 0;
}

static GtDiagram* gt_diagram_new_generic(GtArray *features,
                                          const GtRange *range,
                                          GtStyle *style)
{
  GtDiagram *diagram;
  diagram = gt_malloc(sizeof (GtDiagram));
  diagram->tracks = gt_hashmap_new(HASH_STRING, gt_free_func,
                                (GtFree) gt_track_delete);
  diagram->blocks = gt_hashmap_new(HASH_DIRECT, NULL,
                                (GtFree) blocklist_delete);
  diagram->nodeinfo = gt_hashmap_new(HASH_DIRECT, NULL, NULL);
  diagram->nof_tracks = 0;
  diagram->style = style;
  diagram->range = *range;
  gt_diagram_build(diagram, features);
  return diagram;
}

GtDiagram* gt_diagram_new(GtFeatureIndex *fi, const char *seqid,
                           const GtRange *range, GtStyle *style)
{
  GtDiagram *diagram;
  int had_err = 0;
  GtArray *features = gt_array_new(sizeof (GtGenomeNode*));
  assert(features && seqid && range && style);
  had_err = gt_feature_index_get_features_for_range(fi, features, seqid, *range,
                                                 NULL);
  assert(!had_err); /* <fi> must contain <seqid> */
  diagram = gt_diagram_new_generic(features, range, style);
  gt_array_delete(features);
  return diagram;
}

GtDiagram* gt_diagram_new_from_array(GtArray *features, const GtRange *range,
                                GtStyle *style)
{
  assert(features && range && style);
  return gt_diagram_new_generic(features, range, style);
}

GtRange gt_diagram_get_range(GtDiagram* diagram)
{
  assert(diagram);
  return diagram->range;
}

GtHashmap* gt_diagram_get_tracks(const GtDiagram *diagram)
{
  assert(diagram);
  return diagram->tracks;
}

void gt_diagram_get_lineinfo(const GtDiagram *diagram, GtTracklineInfo *tli)
{
  int had_err;
  assert(diagram);
  had_err = gt_hashmap_foreach(diagram->tracks, gt_diagram_add_tracklines,
                              tli, NULL);
  assert(!had_err); /* gt_diagram_add_tracklines() is sane */
}

int gt_diagram_get_number_of_tracks(const GtDiagram *diagram)
{
  assert(diagram);
  return diagram->nof_tracks;
}

static int blocklist_block_compare(const void *item1, const void *item2)
{
  assert(item1 && item2);
  return gt_block_compare(*(GtBlock**) item1, *(GtBlock**) item2);
}

static int layout_tracks(void *key, void *value, void *data,
                         GT_UNUSED GtError *err)
{
  unsigned long i, max;
  GtTrack *track;
  GtTrackTraverseInfo *tti = (GtTrackTraverseInfo*) data;
  GtArray *list = (GtArray*) value;
  char *filename;
  GtStr *gt_track_key;
  const char *type = key;
  GtBlock *block;
  bool split;
  double tmp;
  assert(type && list);

  /* to get a deterministic layout, we sort the GtBlocks for each type */
  gt_array_sort_stable(list, blocklist_block_compare);
  /* we take the basename of the filename to have nicer output in the
     generated graphic. this might lead to ``collapsed'' tracks, if two files
     with different paths have the same basename. */
  block = *(GtBlock**) gt_array_get(list, 0);
  filename = gt_basename(gt_genome_node_get_filename(
                                     (GtGenomeNode*)
                                       gt_block_get_top_level_feature(block)));
  gt_track_key = gt_track_key_new(filename, type);
  gt_free(filename);

  if (!gt_style_get_bool(tti->dia->style, "format", "split_lines", &split,
                         NULL))
  {
    split = true;
  }
  if (split)
    if (!gt_style_get_bool(tti->dia->style, type, "split_lines", &split, NULL))
      split = true;
  if (gt_style_get_num(tti->dia->style, type, "max_num_lines", &tmp, NULL))
    max = tmp;
  else
    max = 50;

  /* For now, use the captions line breaker */
  track = gt_track_new(gt_track_key, max, split,
                    gt_line_breaker_captions_new(tti->canvas));
  tti->dia->nof_tracks++;
  for (i = 0; i < gt_array_size(list); i++) {
    block = *(GtBlock**) gt_array_get(list, i);
    gt_track_insert_block(track, block);
  }
  gt_hashmap_add(tti->dia->tracks, gt_cstr_dup(gt_str_get(gt_track_key)),
                 track);
  gt_str_delete(gt_track_key);
  return 0;
}

static int render_tracks(GT_UNUSED void *key, void *value, void *data,
                     GT_UNUSED GtError *err)
{
  GtTrackTraverseInfo *tti = (GtTrackTraverseInfo*) data;
  GtTrack *track = (GtTrack*) value;
  int had_err = 0;
  assert(tti && track);
  had_err = gt_track_sketch(track, tti->canvas);
  return had_err;
}

int gt_diagram_sketch(GtDiagram *dia, GtCanvas *canvas)
{
  int had_err = 0;
  GtTrackTraverseInfo tti;
  tti.dia = dia;
  tti.canvas = canvas;
  gt_canvas_visit_diagram_pre(canvas, dia);
  gt_hashmap_reset(dia->tracks);
  dia->nof_tracks = 0;
  (void) gt_hashmap_foreach(dia->blocks, layout_tracks, &tti, NULL);
  gt_canvas_visit_diagram_post(canvas, dia);
  had_err = gt_hashmap_foreach_in_key_order(dia->tracks, render_tracks,
                                         &tti, NULL);

  return had_err;
}

int gt_diagram_unit_test(GtError *err)
{
  GtGenomeNode *gn1, *gn2, *ex1, *ex2, *ex3, *cds1;
  GtFeatureIndex *fi;
  GtRange dr1, rs;
  GtStr *seqid1, *seqid2, *gt_track_key;
  GtRegionNode *rn1, *rn2;
  int had_err=0;
  GtStyle *sty = NULL;
  GtDiagram *dia = NULL, *dia2 = NULL, *dia3 = NULL;
  GtArray *features;
  GtCanvas *canvas = NULL;
  GtWarningHandler warning_handler;
  void *warning_data;
  gt_error_check(err);

  /* backup warning handler and disable warnings */
  warning_handler = gt_warning_get_handler();
  warning_data = gt_warning_get_data();
  gt_warning_disable();

  /* generating some ranges */
  rs.start=100; rs.end=1200;

  /* generating sequence IDs */
  seqid1 = gt_str_new_cstr("test1");
  seqid2 = gt_str_new_cstr("test2");

  rn1 = (GtRegionNode*) gt_region_node_new(seqid1, rs.start, rs.end);
  rn2 = (GtRegionNode*) gt_region_node_new(seqid2, rs.start, rs.end);

  gn1 = gt_feature_node_new(seqid1, gft_gene, 100, 1000, GT_STRAND_UNKNOWN);
  gn2 = gt_feature_node_new(seqid2, gft_gene, 600, 1200, GT_STRAND_UNKNOWN);
  ex1 = gt_feature_node_new(seqid1, gft_exon, 100, 300, GT_STRAND_UNKNOWN);
  ex2 = gt_feature_node_new(seqid1, gft_exon, 500, 1000, GT_STRAND_UNKNOWN);
  ex3 = gt_feature_node_new(seqid2, gft_exon, 600, 1200 , GT_STRAND_UNKNOWN);
  cds1 = gt_feature_node_new(seqid2, gft_CDS, 600, 1200, GT_STRAND_UNKNOWN);

  /* determine the structure of our feature tree */
  gt_feature_node_add_child((GtFeatureNode*) gn1, (GtFeatureNode*) ex1);
  gt_feature_node_add_child((GtFeatureNode*) gn1, (GtFeatureNode*) ex2);
  gt_feature_node_add_child((GtFeatureNode*) gn2, (GtFeatureNode*) ex3);
  gt_feature_node_add_child((GtFeatureNode*) gn2, (GtFeatureNode*) cds1);

  /* create a new feature index on which we can perform some tests */
  fi = gt_feature_index_memory_new();

  /* add features to every sequence region */
  gt_feature_index_add_feature_node(fi, (GtFeatureNode*) gn1);
  gt_feature_index_add_feature_node(fi, (GtFeatureNode*) gn2);

  /* set the GtRange for the diagram */
  dr1.start = 400UL;
  dr1.end   = 900UL;

  /* create a style object */
  if (!had_err) {
    if (!(sty = gt_style_new(err)))
      had_err = -1;
  }

  /* create a diagram object and test it */
  if (!had_err)
    dia = gt_diagram_new(fi, "test1", &dr1, sty);

  ensure(had_err, dia->style);
  ensure(had_err, dia->range.start == 400UL);
  ensure(had_err, dia->range.end == 900UL);

  if (!had_err)
  {
    canvas = gt_canvas_cairo_file_new(sty, GT_GRAPHICS_PNG, 600, NULL);
    gt_diagram_sketch(dia, canvas);
  }

  if (!had_err &&
      !gt_style_get_bool(dia->style, "gene", "collapse_to_parent", false, NULL))
  {
    gt_track_key = gt_track_key_new("generated", gft_gene);
    ensure(had_err, gt_hashmap_get(dia->tracks, gt_str_get(gt_track_key)));
    gt_str_delete(gt_track_key);
  }

  if (!had_err &&
      !gt_style_get_bool(dia->style, "exon", "collapse_to_parent", false, NULL))
  {
    gt_track_key = gt_track_key_new("generated", gft_exon);
    ensure(had_err, gt_hashmap_get(dia->tracks, gt_str_get(gt_track_key)));
    gt_str_delete(gt_track_key);
  }
  ensure(had_err, gt_range_compare(gt_diagram_get_range(dia),dr1) == 0);

  /* create a diagram object and test it */
  if (!had_err) {
    dia2 = gt_diagram_new(fi, "test2", &dr1, sty);
    ensure(had_err, dia->range.start == 400UL);
    ensure(had_err, dia->range.end == 900UL);
  }

  if (!had_err &&
      !gt_style_get_bool(dia2->style, "gene", "collapse_to_parent", false,
                         NULL))
  {
    gt_diagram_sketch(dia2, canvas);
    gt_track_key = gt_track_key_new("generated", gft_gene);
    ensure(had_err, gt_hashmap_get(dia2->tracks, gt_str_get(gt_track_key)));
    gt_str_delete(gt_track_key);
  }

  if (!had_err &&
      !gt_style_get_bool(dia2->style, "exon", "collapse_to_parent", false,
                         NULL))
  {
    gt_track_key = gt_track_key_new("generated", gft_exon);
    ensure(had_err, gt_hashmap_get(dia2->tracks, gt_str_get(gt_track_key)));
    gt_str_delete(gt_track_key);
  }

  if (!had_err &&
      !gt_style_get_bool(dia2->style, "CDS", "collapse_to_parent", false, NULL))
  {
    gt_track_key = gt_track_key_new("generated", gft_CDS);
    ensure(had_err, gt_hashmap_get(dia2->tracks, gt_str_get(gt_track_key)));
    gt_str_delete(gt_track_key);
  }
  ensure(had_err, gt_range_compare(gt_diagram_get_range(dia),dr1) == 0);

  features = gt_array_new(sizeof (GtGenomeNode*));
  gt_array_add(features, gn1);
  gt_array_add(features, gn2);
  dia3 = gt_diagram_new_from_array(features, &rs, sty);

  ensure(had_err, dia3->style);

  if (!had_err &&
      !gt_style_get_bool(dia3->style, "gene", "collapse_to_parent", false,
                         NULL))
  {
    gt_diagram_sketch(dia3, canvas);
    gt_track_key = gt_track_key_new("generated", gft_gene);
    ensure(had_err, gt_hashmap_get(dia3->tracks, gt_str_get(gt_track_key)));
    gt_str_delete(gt_track_key);
  }

  if (!had_err &&
      !gt_style_get_bool(dia3->style, "exon", "collapse_to_parent", false,
                         NULL))
  {
    gt_track_key = gt_track_key_new("generated", gft_exon);
    ensure(had_err, gt_hashmap_get(dia3->tracks, gt_str_get(gt_track_key)));
    gt_str_delete(gt_track_key);
  }
  ensure(had_err, gt_range_compare(gt_diagram_get_range(dia3),rs) == 0);

  /* delete all generated objects */
  gt_style_delete(sty);
  gt_array_delete(features);
  gt_diagram_delete(dia);
  gt_diagram_delete(dia2);
  gt_diagram_delete(dia3);
  gt_canvas_delete(canvas);
  gt_feature_index_delete(fi);
  gt_genome_node_rec_delete(gn1);
  gt_genome_node_rec_delete(gn2);
  gt_genome_node_rec_delete((GtGenomeNode*) rn1);
  gt_genome_node_rec_delete((GtGenomeNode*) rn2);
  gt_str_delete(seqid1);
  gt_str_delete(seqid2);

  /* restore warning handler */
  gt_warning_set_handler(warning_handler, warning_data);

  return had_err;
}

void gt_diagram_delete(GtDiagram *diagram)
{
  if (!diagram) return;
  gt_hashmap_delete(diagram->tracks);
  gt_hashmap_delete(diagram->blocks);
  gt_hashmap_delete(diagram->nodeinfo);
  gt_free(diagram);
}
