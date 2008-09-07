/*
  Copyright (c) 2007 Sascha Steinbiss <ssteinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2007 Malte Mader <mmader@zbh.uni-hamburg.de>
  Copyright (c) 2007 Christin Schaerfer <cschaerfer@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

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

#include "core/cstr.h"
#include "core/ensure.h"
#include "core/getbasename.h"
#include "core/hashmap.h"
#include "core/ma.h"
#include "core/str.h"
#include "core/undef.h"
#include "core/unused_api.h"
#include "extended/feature_type_factory_builtin.h"
#include "extended/genome_node.h"
#include "extended/genome_feature.h"
#include "extended/genome_feature_type.h"
#include "annotationsketch/canvas.h"
#include "annotationsketch/diagram.h"
#include "annotationsketch/feature_index.h"
#include "annotationsketch/line_breaker_captions.h"
#include "annotationsketch/line_breaker_bases.h"
#include "annotationsketch/track.h"

/* used to separate a filename from the type in a track name */
#define FILENAME_TYPE_SEPARATOR  '|'

struct GT_Diagram {
  /* GT_Tracks indexed by track keys */
  Hashmap *tracks;
  /* GT_Block lists indexed by track keys */
  Hashmap *blocks;
  /* Reverse lookup structure (per node) */
  Hashmap *nodeinfo;
  /* Cache tables for configuration data */
  Hashmap *collapsingtypes, *caption_display_status;
  int nof_tracks;
  GT_Style *style;
  GT_Range range;
};

/* holds a GT_Block with associated type */
typedef struct {
  GT_GenomeFeatureType *gft;
  GT_Block *block;
} GT_BlockTuple;

/* a node in the reverse lookup structure used for collapsing */
typedef struct {
  GT_GenomeNode *parent;
  GT_Array *blocktuples;
} NodeInfoGT_Element;

typedef struct {
  GT_GenomeNode *parent;
  GT_Diagram *diagram;
} NodeTraverseInfo;

typedef struct {
  GT_Canvas *canvas;
  GT_Diagram *dia;
} GT_TrackTraverseInfo;

static GT_BlockTuple* blocktuple_new(GT_GenomeFeatureType *gft, GT_Block *block)
{
  GT_BlockTuple *bt;
  assert(block);
  bt = ma_malloc(sizeof (GT_BlockTuple));
  bt->gft = gft;
  bt->block = block;
  return bt;
}

static NodeInfoGT_Element* get_or_create_node_info(GT_Diagram *d, GT_GenomeNode *node)
{
  NodeInfoGT_Element *ni;
  assert(d && node);
  ni = hashmap_get(d->nodeinfo, node);
  if (ni == NULL) {
    NodeInfoGT_Element *new_ni = ma_malloc(sizeof (NodeInfoGT_Element));
    new_ni->blocktuples = gt_array_new(sizeof (GT_BlockTuple*));
    hashmap_add(d->nodeinfo, node, new_ni);
    ni = new_ni;
  }
  return ni;
}

static GT_Block* find_block_for_type(NodeInfoGT_Element* ni,
                                     GT_GenomeFeatureType *gft)
{
  GT_Block *block = NULL;
  unsigned long i;
  assert(ni);
  for (i = 0; i < gt_array_size(ni->blocktuples); i++) {
    GT_BlockTuple *bt = *(GT_BlockTuple**) gt_array_get(ni->blocktuples, i);
    if (bt->gft == gft) {
      block = bt->block;
      break;
    }
  }
  return block;
}

static const char* get_node_name_or_id(GT_GenomeNode *gn)
{
  const char *ret;
  if (!gn) return NULL;
  if (!(ret = gt_genome_feature_get_attribute(gn, "Name"))) {
    if (!(ret = gt_genome_feature_get_attribute(gn, "ID")))
      ret = NULL;
  }
  return ret;
}

static bool get_caption_display_status(GT_Diagram *d, GT_GenomeFeatureType *gft)
{
  assert(d && gft);
  bool *status;

  status = (bool*) hashmap_get(d->caption_display_status, gft);
  if (!status)
  {
    unsigned long threshold;
    double tmp;
    status = ma_malloc(sizeof (bool*));
    if (!gt_style_get_bool(d->style, "format", "show_block_captions",
                       status, NULL))
      *status = true;
    if (*status)
    {
      if (gt_style_get_num(d->style, gt_genome_feature_type_get_cstr(gft),
                         "max_capt_show_width", &tmp, NULL))
        threshold = tmp;
      else
        threshold = UNDEF_ULONG;
      if (threshold == UNDEF_ULONG)
        *status = true;
      else
        *status = (gt_range_length(d->range) <= threshold);
    }
    hashmap_add(d->caption_display_status, gft, status);
  }
  return *status;
}

static void add_to_current(GT_Diagram *d, GT_GenomeNode *node, GT_GenomeNode *parent)
{
  NodeInfoGT_Element *ni;
  GT_Block *block;
  GT_BlockTuple *bt;
  GT_Str *caption = NULL;
  const char *nnid_p = NULL, *nnid_n = NULL;
  assert(d && node);

  /* Lookup node info and set itself as parent */
  ni = get_or_create_node_info(d, node);
  ni->parent = node;
  /* create new GT_Block tuple and add to node info */
  block = gt_block_new_from_node(node);
  /* assign block caption */

  caption = gt_str_new();
  if (!gt_style_get_str(d->style,
                     gt_genome_feature_type_get_cstr(
                         gt_genome_feature_get_type((GT_GenomeFeature*) node)),
                     "block_caption",
                     caption,
                     node))
  {
    nnid_p = get_node_name_or_id(parent);
    nnid_n = get_node_name_or_id(node);
    if ((nnid_p || nnid_n) && get_caption_display_status(d,
                  gt_genome_feature_get_type((GT_GenomeFeature*) node)))
    {
      if (parent) {
        if (gt_genome_node_has_children(parent))
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
  /* insert node into block */
  gt_block_insert_element(block, node);
  bt = blocktuple_new(gt_genome_feature_get_type((GT_GenomeFeature*) node), block);
  gt_array_add(ni->blocktuples, bt);
}

static void add_to_parent(GT_Diagram *d, GT_GenomeNode *node, GT_GenomeNode* parent)
{
  GT_Block *block = NULL;
  NodeInfoGT_Element *par_ni, *ni;
  const char *nnid_p = NULL, *nnid_n = NULL;

  assert(d && node);

  if (!parent) return;

  par_ni = get_or_create_node_info(d, parent);
  ni = get_or_create_node_info(d, node);
  ni->parent = parent;

  /* try to find the right block to insert */
  block = find_block_for_type(par_ni,
                              gt_genome_feature_get_type((GT_GenomeFeature*) node));
  /* no fitting block was found, create a new one */
  if (block == NULL) {
    GT_BlockTuple *bt;
    GT_Str *caption = NULL;
    block = gt_block_new_from_node(parent);
    /* assign block caption */
    nnid_p = get_node_name_or_id(parent);
    nnid_n = get_node_name_or_id(node);
    if ((nnid_p || nnid_n) && get_caption_display_status(d,
                  gt_genome_feature_get_type((GT_GenomeFeature*) node)))
    {
      caption = gt_str_new_cstr("");
      if (parent) {
        if (gt_genome_node_has_children(parent))
          gt_str_append_cstr(caption, nnid_p);
        else
          gt_str_append_cstr(caption, "-");
        gt_str_append_cstr(caption, "/");
      }
      if (nnid_n)
        gt_str_append_cstr(caption, nnid_n);
    }
    gt_block_set_caption(block, caption);
    /* add block to nodeinfo */
    bt = blocktuple_new(gt_genome_feature_get_type((GT_GenomeFeature*) node), block);
    gt_array_add(par_ni->blocktuples, bt);
  }
  /* now we have a block to insert into */
  gt_block_insert_element(block, node);
}

static void add_recursive(GT_Diagram *d, GT_GenomeNode *node,
                          GT_GenomeNode* parent, GT_GenomeNode *original_node)
{
  NodeInfoGT_Element *ni;

  assert(d && node && original_node);
  if (!parent) return;

  ni = get_or_create_node_info(d, node);

  /* end of recursion, insert into target block */
  if (parent == node) {
    GT_Block *block ;
    GT_BlockTuple *bt;
    /* try to find the right block to insert */
    block = find_block_for_type(ni,
                                gt_genome_feature_get_type((GT_GenomeFeature*) node));
    if (block == NULL) {
      block = gt_block_new_from_node(node);
      bt = blocktuple_new(gt_genome_feature_get_type((GT_GenomeFeature*) node),
                          block);
      gt_array_add(ni->blocktuples, bt);
    }
    gt_block_insert_element(block, original_node);
  }
  else {
    /* not at target type block yet, set up reverse entry and follow */
    NodeInfoGT_Element *parent_ni;
    /* set up reverse entry */
    ni->parent = parent;
    /* recursively call with parent node and its parent */
    parent_ni = hashmap_get(d->nodeinfo, parent);
    if (parent_ni)
      add_recursive(d, parent, parent_ni->parent, original_node);
  }
}

static void process_node(GT_Diagram *d, GT_GenomeNode *node, GT_GenomeNode *parent)
{
  GT_Range elem_range;
  bool *collapse, do_not_overlap=false;
  const char *feature_type = NULL, *parent_gft = NULL;
  double tmp;
  unsigned long max_show_width = UNDEF_ULONG,
                par_max_show_width = UNDEF_ULONG;

  assert(d && node);

  feature_type = gt_genome_feature_type_get_cstr(
                        gt_genome_feature_get_type((GT_GenomeFeature*) node));
  if (parent)
    parent_gft = gt_genome_feature_type_get_cstr(
                        gt_genome_feature_get_type((GT_GenomeFeature*) parent));

  /* discard elements that do not overlap with visible range */
  elem_range = gt_genome_node_get_range(node);
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
  if (max_show_width != UNDEF_ULONG && gt_range_length(d->range) > max_show_width)
    return;
  if (parent && par_max_show_width != UNDEF_ULONG
        && gt_range_length(d->range) > par_max_show_width)
    parent = NULL;

  /* check if this is a collapsing type, cache result */
  if ((collapse = (bool*) hashmap_get(d->collapsingtypes,
                                        feature_type)) == NULL)
  {
    collapse = ma_malloc(sizeof (bool));
    if (!gt_style_get_bool(d->style, feature_type, "collapse_to_parent",
                        collapse, NULL))
      *collapse = false;
    hashmap_add(d->collapsingtypes, (char*) feature_type, collapse);
  }

  /* check if direct children overlap */
  if (parent)
    do_not_overlap =
      gt_genome_node_direct_children_do_not_overlap_st(parent, node);

  /* decide how to continue: */
  if (*collapse && parent) {
    /* collapsing features recursively search their target blocks */
    add_recursive(d, node, parent, node);
  }
  else if (do_not_overlap
            && gt_genome_node_number_of_children(parent) > 1)
  {
    /* group non-overlapping child nodes of a non-collapsing type by parent */
    add_to_parent(d, node, parent);
  }
  else {
    /* nodes that belong into their own track and block */
    add_to_current(d, node, parent);
  }

  /* we can now assume that this node has been processed into the reverse
     lookup structure */
  assert(hashmap_get(d->nodeinfo, node));
}

static int gt_diagram_add_tracklines(GT_UNUSED void *key, void *value, void *data,
                                  GT_UNUSED GT_Error *err)
{
  GT_TracklineInfo *add = (GT_TracklineInfo*) data;
  add->total_lines += gt_track_get_number_of_lines((GT_Track*) value);
  add->total_captionlines += gt_track_get_number_of_lines_with_captions(
                                                               (GT_Track*) value);
  return 0;
}

static int visit_child(GT_GenomeNode* gn, void* gt_genome_node_children, GT_Error* e)
{
  NodeTraverseInfo* gt_genome_node_info;
  gt_genome_node_info = (NodeTraverseInfo*) gt_genome_node_children;
  int had_err;
  if (gt_genome_node_has_children(gn))
  {
    GT_GenomeNode *oldparent = gt_genome_node_info->parent;
    process_node(gt_genome_node_info->diagram, gn, gt_genome_node_info->parent);
    gt_genome_node_info->parent = gn;
    had_err = gt_genome_node_traverse_direct_children(gn, gt_genome_node_info,
                                                   visit_child, e);
    assert(!had_err); /* visit_child() is sane */
    gt_genome_node_info->parent = oldparent;
  }
  else
    process_node(gt_genome_node_info->diagram, gn, gt_genome_node_info->parent);
  return 0;
}

static GT_Str* gt_track_key_new(const char *filename, GT_GenomeFeatureType *type)
{
  GT_Str *gt_track_key;
  gt_track_key = gt_str_new_cstr(filename);
  gt_str_append_char(gt_track_key, FILENAME_TYPE_SEPARATOR);
  gt_str_append_cstr(gt_track_key, gt_genome_feature_type_get_cstr(type));
  return gt_track_key;
}

/* Create lists of all GT_Blocks in the diagram. */
static int collect_blocks(GT_UNUSED void *key, void *value, void *data,
                          GT_UNUSED GT_Error *err)
{
  NodeInfoGT_Element *ni = (NodeInfoGT_Element*) value;
  GT_Diagram *diagram = (GT_Diagram*) data;
  unsigned long i = 0;

  for (i = 0; i < gt_array_size(ni->blocktuples); i++) {
    GT_Array *list;
    GT_BlockTuple *bt = *(GT_BlockTuple**) gt_array_get(ni->blocktuples, i);
    list = (GT_Array*) hashmap_get(diagram->blocks, bt->gft);
    if (!list)
    {
      list = gt_array_new(sizeof (GT_Block*));
      hashmap_add(diagram->blocks, bt->gft, list);
    }
    assert(list);
    gt_array_add(list, bt->block);
    ma_free(bt);
  }
  gt_array_delete(ni->blocktuples);
  ma_free(ni);
  return 0;
}

/* Traverse a genome node graph with depth first search. */
static void traverse_genome_nodes(GT_GenomeNode *gn, void *gt_genome_node_children)
{
  NodeTraverseInfo* gt_genome_node_info;
  int had_err;
  assert(gt_genome_node_children);
  gt_genome_node_info = (NodeTraverseInfo*) gt_genome_node_children;
  gt_genome_node_info->parent = gn;
  /* handle root nodes */
  process_node(gt_genome_node_info->diagram, gn, NULL);
  if (gt_genome_node_has_children(gn)) {
    had_err = gt_genome_node_traverse_direct_children(gn, gt_genome_node_info,
                                                   visit_child, NULL);
    assert(!had_err); /* visit_child() is sane */
  }
}

static void gt_diagram_build(GT_Diagram *diagram, GT_Array *features)
{
  unsigned long i = 0;
  int had_err;
  NodeTraverseInfo gt_genome_node_children;
  gt_genome_node_children.diagram = diagram;

  /* initialise caches */
  diagram->collapsingtypes = hashmap_new(HASH_STRING, NULL, ma_free_func);
  diagram->caption_display_status = hashmap_new(HASH_DIRECT,
                                                  NULL, ma_free_func);

  /* do node traversal for each root feature */
  for (i = 0; i < gt_array_size(features); i++) {
    GT_GenomeNode *current_root = *(GT_GenomeNode**) gt_array_get(features,i);
    traverse_genome_nodes(current_root, &gt_genome_node_children);
  }
  /* collect blocks from nodeinfo structures and create the tracks */
  had_err = hashmap_foreach_ordered(diagram->nodeinfo, collect_blocks,
                                      diagram, (GT_Compare) gt_genome_node_cmp, NULL);
  assert(!had_err); /* collect_blocks() is sane */

  /* clear caches */
  hashmap_delete(diagram->collapsingtypes);
  hashmap_delete(diagram->caption_display_status);
}

static int blocklist_delete(void *value)
{
  unsigned long i;
  GT_Array *a = (GT_Array*) value;
  for (i = 0; i < gt_array_size(a); i++)
    gt_block_delete(*(GT_Block**) gt_array_get(a, i));
  gt_array_delete(a);
  return 0;
}

static GT_Diagram* gt_diagram_new_generic(GT_Array *features, const GT_Range *range,
                                    GT_Style *style)
{
  GT_Diagram *diagram;
  diagram = ma_malloc(sizeof (GT_Diagram));
  diagram->tracks = hashmap_new(HASH_STRING, ma_free_func,
                                (GT_FreeFunc) gt_track_delete);
  diagram->blocks = hashmap_new(HASH_DIRECT, NULL,
                                  (GT_FreeFunc) blocklist_delete);
  diagram->nodeinfo = hashmap_new(HASH_DIRECT, NULL, NULL);
  diagram->nof_tracks = 0;
  diagram->style = style;
  diagram->range = *range;
  gt_diagram_build(diagram, features);
  return diagram;
}

GT_Diagram* gt_diagram_new(GT_FeatureIndex *fi, const char *seqid,
                           const GT_Range *range, GT_Style *style)
{
  GT_Diagram *diagram;
  int had_err = 0;
  GT_Array *features = gt_array_new(sizeof (GT_GenomeNode*));
  assert(features && seqid && range && style);
  had_err = gt_feature_index_get_features_for_range(fi, features, seqid, *range,
                                                 NULL);
  assert(!had_err); /* <fi> must contain <seqid> */
  diagram = gt_diagram_new_generic(features, range, style);
  gt_array_delete(features);
  return diagram;
}

GT_Diagram* gt_diagram_new_from_array(GT_Array *features, const GT_Range *range,
                                GT_Style *style)
{
  assert(features && range && style);
  return gt_diagram_new_generic(features, range, style);
}

GT_Range gt_diagram_get_range(GT_Diagram* diagram)
{
  assert(diagram);
  return diagram->range;
}

Hashmap* gt_diagram_get_tracks(const GT_Diagram *diagram)
{
  assert(diagram);
  return diagram->tracks;
}

void gt_diagram_get_lineinfo(const GT_Diagram *diagram, GT_TracklineInfo *tli)
{
  int had_err;
  assert(diagram);
  had_err = hashmap_foreach(diagram->tracks, gt_diagram_add_tracklines,
                              tli, NULL);
  assert(!had_err); /* gt_diagram_add_tracklines() is sane */
}

int gt_diagram_get_number_of_tracks(const GT_Diagram *diagram)
{
  assert(diagram);
  return diagram->nof_tracks;
}

static int blocklist_block_compare(const void *item1, const void *item2)
{
  assert(item1 && item2);
  return gt_block_compare(*(GT_Block**) item1, *(GT_Block**) item2);
}

static int layout_tracks(void *key, void *value, void *data,
                         GT_UNUSED GT_Error *err)
{
  unsigned long i, max;
  GT_Track *track;
  GT_TrackTraverseInfo *tti = (GT_TrackTraverseInfo*) data;
  GT_GenomeFeatureType *gft = (GT_GenomeFeatureType*) key;
  GT_Array *list = (GT_Array*) value;
  char *filename;
  GT_Str *gt_track_key;
  const char *type;
  GT_Block *block;
  bool split;
  double tmp;
  assert(gft && list);

  /* to get a deterministic layout, we sort the GT_Blocks for each type */
  gt_array_sort(list, blocklist_block_compare);
  /* we take the basename of the filename to have nicer output in the
     generated graphic. this might lead to ``collapsed'' tracks, if two files
     with different paths have the same basename. */
  block = *(GT_Block**) gt_array_get(list, 0);
  filename = getbasename(gt_genome_node_get_filename(
                                        gt_block_get_top_level_feature(block)));
  gt_track_key = gt_track_key_new(filename, gft);
  ma_free(filename);
  type = gt_genome_feature_type_get_cstr(gft);

  if (!gt_style_get_bool(tti->dia->style, "format", "split_lines", &split,
                         NULL)) {
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
    block = *(GT_Block**) gt_array_get(list, i);
    gt_track_insert_block(track, block);
  }
  hashmap_add(tti->dia->tracks, cstr_dup(gt_str_get(gt_track_key)), track);
  gt_str_delete(gt_track_key);
  return 0;
}

static int render_tracks(GT_UNUSED void *key, void *value, void *data,
                     GT_UNUSED GT_Error *err)
{
  GT_TrackTraverseInfo *tti = (GT_TrackTraverseInfo*) data;
  GT_UNUSED GT_Track *track = (GT_Track*) value;
  int had_err = 0;
  assert(tti && track);
  had_err = gt_track_sketch((GT_Track*) value, tti->canvas);
  return had_err;
}

int gt_diagram_sketch(GT_Diagram *dia, GT_Canvas *canvas)
{
  int had_err = 0;
  GT_TrackTraverseInfo tti;
  tti.dia = dia;
  tti.canvas = canvas;
  gt_canvas_visit_gt_diagram_pre(canvas, dia);
  hashmap_reset(dia->tracks);
  dia->nof_tracks = 0;
  (void) hashmap_foreach(dia->blocks, layout_tracks, &tti, NULL);
  gt_canvas_visit_gt_diagram_post(canvas, dia);
  had_err = hashmap_foreach_in_key_order(dia->tracks, render_tracks,
                                         &tti, NULL);

  return had_err;
}

int gt_diagram_unit_test(GT_Error *err)
{
  GT_FeatureTypeFactory *feature_type_factory;
  GT_GenomeFeatureType *gene_type, *exon_type, *CDS_type;
  GT_GenomeNode *gn1, *gn2, *ex1, *ex2, *ex3, *cds1;
  GT_FeatureIndex *fi;
  GT_Range r1, r2, r3, r4, r5, dr1, rs;
  GT_Str *seqid1, *seqid2, *gt_track_key;
  GT_SequenceRegion *sr1, *sr2;
  int had_err=0;
  GT_Style *sty = NULL;
  GT_Diagram *dia = NULL, *dia2 = NULL, *dia3 = NULL;
  GT_Array *features;
  GT_Canvas *canvas = NULL;
  gt_error_check(err);

  feature_type_factory = gt_feature_type_factory_builtin_new();
  gene_type = gt_feature_type_factory_create_gft(feature_type_factory, gft_gene);
  exon_type = gt_feature_type_factory_create_gft(feature_type_factory, gft_exon);
  CDS_type = gt_feature_type_factory_create_gft(feature_type_factory, gft_CDS);

  /* generating some ranges */
  r1.start=100UL; r1.end=1000UL;
  r2.start=100UL; r2.end=300UL;
  r3.start=500UL; r3.end=1000UL;
  r4.start=600UL; r4.end=1200UL;
  r5.start=600UL; r5.end=1000UL;
  rs.start=100UL; rs.end=1200UL;

  /* generating sequence IDs */
  seqid1 = gt_str_new_cstr("test1");
  seqid2 = gt_str_new_cstr("test2");

  sr1 = (GT_SequenceRegion*) gt_sequence_regionnew(seqid1, rs);
  sr2 = (GT_SequenceRegion*) gt_sequence_regionnew(seqid2, rs);

  gn1 = gt_genome_feature_new(seqid1, gene_type, r1, GT_STRAND_UNKNOWN);

  gn2 = gt_genome_feature_new(seqid2, gene_type, r4, GT_STRAND_UNKNOWN);

  ex1 = gt_genome_feature_new(seqid1, exon_type, r2, GT_STRAND_UNKNOWN);

  ex2 = gt_genome_feature_new(seqid1, exon_type, r3, GT_STRAND_UNKNOWN);

  ex3 = gt_genome_feature_new(seqid2, exon_type, r4, GT_STRAND_UNKNOWN);

  cds1 = gt_genome_feature_new(seqid2, CDS_type, r5, GT_STRAND_UNKNOWN);

  /* determine the structure of our feature tree */
  gt_genome_node_is_part_of_genome_node(gn1, ex1);
  gt_genome_node_is_part_of_genome_node(gn1, ex2);
  gt_genome_node_is_part_of_genome_node(gn2, ex3);
  gt_genome_node_is_part_of_genome_node(gn2, cds1);

  /* create a new feature index on which we can perform some tests */
  fi = gt_feature_index_new();

  /* add features to every sequence region */
  gt_feature_index_add_genome_feature(fi, (GT_GenomeFeature*) gn1);
  gt_feature_index_add_genome_feature(fi, (GT_GenomeFeature*) gn2);

  /* set the GT_Range for the diagram */
  dr1.start = 400UL;
  dr1.end   = 900UL;

  /* create a style object */
  if (!had_err) {
    if (!(sty = gt_style_new(false, err)))
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
    canvas = gt_canvas_new(sty, GRAPHICS_PNG, 600, NULL);
    gt_diagram_sketch(dia, canvas);
  }

  if (!had_err &&
      !gt_style_get_bool(dia->style, "gene", "collapse_to_parent", false, NULL))
  {
    gt_track_key = gt_track_key_new("generated", gene_type);
    ensure(had_err, hashmap_get(dia->tracks, gt_str_get(gt_track_key)));
    gt_str_delete(gt_track_key);
  }

  if (!had_err &&
      !gt_style_get_bool(dia->style, "exon", "collapse_to_parent", false, NULL))
  {
    gt_track_key = gt_track_key_new("generated", exon_type);
    ensure(had_err, hashmap_get(dia->tracks, gt_str_get(gt_track_key)));
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
    gt_track_key = gt_track_key_new("generated", gene_type);
    ensure(had_err, hashmap_get(dia2->tracks, gt_str_get(gt_track_key)));
    gt_str_delete(gt_track_key);
  }

  if (!had_err &&
      !gt_style_get_bool(dia2->style, "exon", "collapse_to_parent", false,
                         NULL))
  {
    gt_track_key = gt_track_key_new("generated", exon_type);
    ensure(had_err, hashmap_get(dia2->tracks, gt_str_get(gt_track_key)));
    gt_str_delete(gt_track_key);
  }

  if (!had_err &&
      !gt_style_get_bool(dia2->style, "CDS", "collapse_to_parent", false, NULL))
  {
    gt_track_key = gt_track_key_new("generated", CDS_type);
    ensure(had_err, hashmap_get(dia2->tracks, gt_str_get(gt_track_key)));
    gt_str_delete(gt_track_key);
  }
  ensure(had_err, gt_range_compare(gt_diagram_get_range(dia),dr1) == 0);

  features = gt_array_new(sizeof (GT_GenomeNode*));
  gt_array_add(features, gn1);
  gt_array_add(features, gn2);
  dia3 = gt_diagram_new_from_array(features, &rs, sty);

  ensure(had_err, dia3->style);

  if (!had_err &&
      !gt_style_get_bool(dia3->style, "gene", "collapse_to_parent", false,
                         NULL))
  {
    gt_diagram_sketch(dia3, canvas);
    gt_track_key = gt_track_key_new("generated", gene_type);
    ensure(had_err, hashmap_get(dia3->tracks, gt_str_get(gt_track_key)));
    gt_str_delete(gt_track_key);
  }

  if (!had_err &&
      !gt_style_get_bool(dia3->style, "exon", "collapse_to_parent", false,
                         NULL))
  {
    gt_track_key = gt_track_key_new("generated", exon_type);
    ensure(had_err, hashmap_get(dia3->tracks, gt_str_get(gt_track_key)));
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
  gt_genome_node_rec_delete((GT_GenomeNode*) sr1);
  gt_genome_node_rec_delete((GT_GenomeNode*) sr2);
  gt_str_delete(seqid1);
  gt_str_delete(seqid2);
  gt_feature_type_factory_delete(feature_type_factory);

  return had_err;
}

void gt_diagram_delete(GT_Diagram *diagram)
{
  if (!diagram) return;
  hashmap_delete(diagram->tracks);
  hashmap_delete(diagram->blocks);
  hashmap_delete(diagram->nodeinfo);
  ma_free(diagram);
}
