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

#include "libgtcore/cstr.h"
#include "libgtcore/ensure.h"
#include "libgtcore/getbasename.h"
#include "libgtcore/log.h"
#include "libgtcore/ma.h"
#include "libgtcore/str.h"
#include "libgtcore/undef.h"
#include "libgtcore/unused.h"
#include "libgtcore/warning.h"
#include "libgtext/feature_type_factory_builtin.h"
#include "libgtext/genome_node.h"
#include "libgtext/genome_feature.h"
#include "libgtext/genome_feature_type.h"
#include "libgtview/canvas.h"
#include "libgtview/diagram.h"
#include "libgtview/feature_index.h"
#include "libgtview/line_breaker_captions.h"
#include "libgtview/line_breaker_bases.h"
#include "libgtview/track.h"

/* used to separate a filename from the type in a track name */
#define FILENAME_TYPE_SEPARATOR  '|'

struct Diagram {
  /* Tracks indexed by track keys */
  Hashtable *tracks;
  /* Block lists indexed by track keys */
  Hashtable *blocks;
  /* Reverse lookup structure (per node) */
  Hashtable *nodeinfo;
  /* Cache tables for configuration data */
  Hashtable *collapsingtypes, *caption_display_status;
  int nof_tracks;
  Config *config;
  Range range;
};

/* holds a Block with associated type */
typedef struct {
  GenomeFeatureType *gft;
  Block *block;
} BlockTuple;

/* a node in the reverse lookup structure used for collapsing */
typedef struct {
  GenomeNode *parent;
  Array *blocktuples;
} NodeInfoElement;

typedef struct {
  GenomeNode *parent;
  Diagram *diagram;
} NodeTraverseInfo;

typedef struct {
  Canvas *canvas;
  Diagram *dia;
} TrackTraverseInfo;

static BlockTuple* blocktuple_new(GenomeFeatureType *gft, Block *block)
{
  BlockTuple *bt;
  assert(block);
  bt = ma_malloc(sizeof (BlockTuple));
  bt->gft = gft;
  bt->block = block;
  return bt;
}

static NodeInfoElement* get_or_create_node_info(Diagram *d, GenomeNode *node)
{
  NodeInfoElement *ni;
  assert(d && node);
  ni = hashtable_get(d->nodeinfo, node);
  if (ni == NULL) {
    NodeInfoElement *new_ni = ma_malloc(sizeof (NodeInfoElement));
    new_ni->blocktuples = array_new(sizeof (BlockTuple*));
    hashtable_add(d->nodeinfo, node, new_ni);
    ni = new_ni;
  }
  return ni;
}

static Block* find_block_for_type(NodeInfoElement* ni, GenomeFeatureType *gft)
{
  Block *block = NULL;
  unsigned long i;
  assert(ni);
  for (i = 0; i < array_size(ni->blocktuples); i++) {
    BlockTuple *bt = *(BlockTuple**) array_get(ni->blocktuples, i);
    if (bt->gft == gft) {
      log_log("    found root block with type %s, inserting...",
              genome_feature_type_get_cstr(bt->gft));
      block = bt->block;
      break;
    }
  }
  return block;
}

static const char* get_node_name_or_id(GenomeNode *gn)
{
  const char *ret;
  if (!gn) return NULL;
  if (!(ret = genome_feature_get_attribute(gn, "Name"))) {
    if (!(ret = genome_feature_get_attribute(gn, "ID")))
      ret = NULL;
  }
  return ret;
}

static bool get_caption_display_status(Diagram *d, GenomeFeatureType *gft)
{
  assert(d && gft);
  bool *status;

  status = (bool*) hashtable_get(d->caption_display_status, gft);
  if (!status)
  {
    unsigned long threshold;
    double tmp;
    status = ma_malloc(sizeof (bool*));
    if (!config_get_bool(d->config, "format", "show_block_captions",
                       status))
      *status = true;
    if (*status)
    {
      if (config_get_num(d->config, genome_feature_type_get_cstr(gft),
                         "max_capt_show_width", &tmp))
        threshold = tmp;
      else
        threshold = UNDEF_ULONG;
      if (threshold == UNDEF_ULONG)
        *status = true;
      else
        *status = (range_length(d->range) <= threshold);
    }
    hashtable_add(d->caption_display_status, gft, status);
  }
  return *status;
}

static void add_to_current(Diagram *d, GenomeNode *node, GenomeNode *parent)
{
  NodeInfoElement *ni;
  Block *block;
  BlockTuple *bt;
  Str *caption = NULL;
  const char *nnid_p = NULL, *nnid_n = NULL;

  log_log("  calling add_to_current");

  assert(d && node);

  /* Lookup node info and set itself as parent */
  ni = get_or_create_node_info(d, node);
  ni->parent = node;
  /* create new Block tuple and add to node info */
  block = block_new_from_node(node);
  /* assign block caption */
  nnid_p = get_node_name_or_id(parent);
  nnid_n = get_node_name_or_id(node);
  if ((nnid_p || nnid_n) && get_caption_display_status(d,
                genome_feature_get_type((GenomeFeature*) node)))
  {
    caption = str_new_cstr("");
    if (parent) {
      if (genome_node_has_children(parent))
        str_append_cstr(caption, nnid_p);
      else
        str_append_cstr(caption, "-");
      str_append_cstr(caption, "/");
    }
    if (nnid_n)
      str_append_cstr(caption, nnid_n);
  }
  block_set_caption(block, caption);
  /* insert node into block */
  block_insert_element(block, node);
  bt = blocktuple_new(genome_feature_get_type((GenomeFeature*) node), block);
  array_add(ni->blocktuples, bt);
}

static void add_to_parent(Diagram *d, GenomeNode *node, GenomeNode* parent)
{
  Block *block = NULL;
  NodeInfoElement *par_ni, *ni;
  const char *nnid_p = NULL, *nnid_n = NULL;

  assert(d && node);

  if (!parent) return;

  log_log("  calling add_to_parent: %s -> %s",
          genome_feature_type_get_cstr(
            genome_feature_get_type((GenomeFeature*) node)),
          genome_feature_type_get_cstr(
            genome_feature_get_type((GenomeFeature*) parent)));

  par_ni = get_or_create_node_info(d, parent);
  ni = get_or_create_node_info(d, node);
  ni->parent = parent;

  /* try to find the right block to insert */
  block = find_block_for_type(par_ni,
                              genome_feature_get_type((GenomeFeature*) node));
  /* no fitting block was found, create a new one */
  if (block == NULL) {
    BlockTuple *bt;
    Str *caption = NULL;
    block = block_new_from_node(parent);
    /* assign block caption */
    nnid_p = get_node_name_or_id(parent);
    nnid_n = get_node_name_or_id(node);
    if ((nnid_p || nnid_n) && get_caption_display_status(d,
                  genome_feature_get_type((GenomeFeature*) node)))
    {
      caption = str_new_cstr("");
      if (parent) {
        if (genome_node_has_children(parent))
          str_append_cstr(caption, nnid_p);
        else
          str_append_cstr(caption, "-");
        str_append_cstr(caption, "/");
      }
      if (nnid_n)
        str_append_cstr(caption, nnid_n);
    }
    block_set_caption(block, caption);
    /* add block to nodeinfo */
    bt = blocktuple_new(genome_feature_get_type((GenomeFeature*) node), block);
    array_add(par_ni->blocktuples, bt);

    log_log("    added %s to (%s) block for node %s",
                genome_feature_type_get_cstr(
                  genome_feature_get_type((GenomeFeature*) node)),
                genome_feature_type_get_cstr(
                  genome_feature_get_type((GenomeFeature*) parent)),
                genome_feature_get_attribute(parent, "ID" ));
  }
  /* now we have a block to insert into */
  block_insert_element(block, node);
}

static void add_recursive(Diagram *d, GenomeNode *node,
                          GenomeNode* parent, GenomeNode *original_node)
{
  NodeInfoElement *ni;

  assert(d && node &&
          original_node);

  if (!parent) return;

  log_log("  calling add_recursive: %s (%s) -> %s (%s)",
              genome_feature_type_get_cstr(
                genome_feature_get_type((GenomeFeature*) node)),
              genome_feature_get_attribute(node,"ID"),
              genome_feature_type_get_cstr(
                genome_feature_get_type((GenomeFeature*) parent)),
              genome_feature_get_attribute(parent,"ID"));

  ni = get_or_create_node_info(d, node);

  /* end of recursion, insert into target block */
  if (parent == node) {
    Block *block ;
    BlockTuple *bt;
    /* try to find the right block to insert */
    block = find_block_for_type(ni,
                                genome_feature_get_type((GenomeFeature*) node));
    if (block == NULL) {
      block = block_new_from_node(node);
      bt = blocktuple_new(genome_feature_get_type((GenomeFeature*) node),
                          block);
      array_add(ni->blocktuples, bt);
    }
    block_insert_element(block, original_node);
  }
  else {
    /* not at target type block yet, set up reverse entry and follow */
    NodeInfoElement *parent_ni;
    /* set up reverse entry */
    ni->parent = parent;
    /* recursively call with parent node and its parent */
    parent_ni = hashtable_get(d->nodeinfo, parent);
    if (parent_ni)
      add_recursive(d, parent, parent_ni->parent, original_node);
  }
}

static void process_node(Diagram *d, GenomeNode *node, GenomeNode *parent)
{
  Range elem_range;
  bool *collapse, do_not_overlap=false;
  const char *feature_type = NULL, *parent_gft = NULL;
  double tmp;
  unsigned long max_show_width = UNDEF_ULONG,
                par_max_show_width = UNDEF_ULONG;

  assert(d && node);

  feature_type = genome_feature_type_get_cstr(
                        genome_feature_get_type((GenomeFeature*) node));
  if (parent)
    parent_gft = genome_feature_type_get_cstr(
                        genome_feature_get_type((GenomeFeature*) parent));

  /* discard elements that do not overlap with visible range */
  elem_range = genome_node_get_range(node);
  if (!range_overlap(d->range, elem_range))
    return;

  /* get maximal view widths in nucleotides to show this type */
  if (config_get_num(d->config, feature_type, "max_show_width", &tmp))
    max_show_width = tmp;
  else
    max_show_width = UNDEF_ULONG;
  if (parent)
  {
    if (config_get_num(d->config, parent_gft, "max_show_width", &tmp))
    par_max_show_width = tmp;
  else
    par_max_show_width = UNDEF_ULONG;

  }
  /* check if this type is to be displayed */
  if (max_show_width != UNDEF_ULONG && range_length(d->range) > max_show_width)
    return;
  if (parent && par_max_show_width != UNDEF_ULONG
        && range_length(d->range) > par_max_show_width)
    parent = NULL;

  /* check if this is a collapsing type, cache result */
  if ((collapse = (bool*) hashtable_get(d->collapsingtypes,
                                        feature_type)) == NULL)
  {
    collapse = ma_malloc(sizeof (bool));
    if (!config_get_bool(d->config, feature_type, "collapse_to_parent",
                        collapse))
      *collapse = false;
    hashtable_add(d->collapsingtypes, (char*) feature_type, collapse);
  }

  /* check if direct children overlap */
  if (parent)
    do_not_overlap =
      genome_node_direct_children_do_not_overlap_st(parent, node);

  log_log("processing node: %s (%s)", feature_type,
          genome_feature_get_attribute(node, "ID" ));

  /* decide how to continue: */
  if (*collapse && parent) {
    /* collapsing features recursively search their target blocks */
    if (!do_not_overlap)
      warning("collapsing %s features overlap "
              "and will be missing in the %s block!",
              feature_type,
              genome_feature_get_attribute(parent, "ID" ));
    add_recursive(d, node, parent, node);
  }
  else if (do_not_overlap
            && genome_node_number_of_children(parent) > 1)
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
  assert(hashtable_get(d->nodeinfo, node));
}

static int diagram_add_tracklines(UNUSED void *key, void *value, void *data,
                                  UNUSED Error *err)
{
  TracklineInfo *add = (TracklineInfo*) data;
  add->total_lines += track_get_number_of_lines((Track*) value);
  add->total_captionlines += track_get_number_of_lines_with_captions(
                                                               (Track*) value);
  return 0;
}

static int visit_child(GenomeNode* gn, void* genome_node_children, Error* e)
{
  NodeTraverseInfo* genome_node_info;
  genome_node_info = (NodeTraverseInfo*) genome_node_children;
  int had_err;
  if (genome_node_has_children(gn))
  {
    GenomeNode *oldparent = genome_node_info->parent;
    process_node(genome_node_info->diagram, gn, genome_node_info->parent);
    genome_node_info->parent = gn;
    had_err = genome_node_traverse_direct_children(gn, genome_node_info,
                                                   visit_child, e);
    assert(!had_err); /* visit_child() is sane */
    genome_node_info->parent = oldparent;
  }
  else
    process_node(genome_node_info->diagram, gn, genome_node_info->parent);
  return 0;
}

static Str* track_key_new(const char *filename, GenomeFeatureType *type)
{
  Str *track_key;
  track_key = str_new_cstr(filename);
  str_append_char(track_key, FILENAME_TYPE_SEPARATOR);
  str_append_cstr(track_key, genome_feature_type_get_cstr(type));
  return track_key;
}

/* Create lists of all Blocks in the diagram. */
static int collect_blocks(UNUSED void *key, void *value, void *data,
                          UNUSED Error *err)
{
  NodeInfoElement *ni = (NodeInfoElement*) value;
  Diagram *diagram = (Diagram*) data;
  unsigned long i = 0;

  for (i = 0; i < array_size(ni->blocktuples); i++) {
    Array *list;
    BlockTuple *bt = *(BlockTuple**) array_get(ni->blocktuples, i);
    list = (Array*) hashtable_get(diagram->blocks, bt->gft);
    if (!list)
    {
      list = array_new(sizeof (Block*));
      hashtable_add(diagram->blocks, bt->gft, list);
    }
    assert(list);
    array_add(list, bt->block);
    ma_free(bt);
  }
  array_delete(ni->blocktuples);
  ma_free(ni);
  return 0;
}

/* Traverse a genome node graph with depth first search. */
static void traverse_genome_nodes(GenomeNode *gn, void *genome_node_children)
{
  NodeTraverseInfo* genome_node_info;
  int had_err;
  assert(genome_node_children);
  genome_node_info = (NodeTraverseInfo*) genome_node_children;
  genome_node_info->parent = gn;
  /* handle root nodes */
  process_node(genome_node_info->diagram, gn, NULL);
  if (genome_node_has_children(gn)) {
    had_err = genome_node_traverse_direct_children(gn, genome_node_info,
                                                   visit_child, NULL);
    assert(!had_err); /* visit_child() is sane */
  }
}

static void diagram_build(Diagram *diagram, Array *features)
{
  unsigned long i = 0;
  int had_err;
  NodeTraverseInfo genome_node_children;
  genome_node_children.diagram = diagram;

  /* initialise caches */
  diagram->collapsingtypes = hashtable_new(HASH_STRING, NULL, ma_free_func);
  diagram->caption_display_status = hashtable_new(HASH_DIRECT,
                                                  NULL, ma_free_func);

  /* do node traversal for each root feature */
  for (i = 0; i < array_size(features); i++) {
    GenomeNode *current_root = *(GenomeNode**) array_get(features,i);
    traverse_genome_nodes(current_root, &genome_node_children);
  }
  /* collect blocks from nodeinfo structures and create the tracks */
  had_err = hashtable_foreach_ordered(diagram->nodeinfo, collect_blocks,
                                      diagram, (Compare) genome_node_cmp, NULL);
  assert(!had_err); /* collect_blocks() is sane */

  /* clear caches */
  hashtable_delete(diagram->collapsingtypes);
  hashtable_delete(diagram->caption_display_status);
}

static int blocklist_delete(void *value)
{
  unsigned long i;
  Array *a = (Array*) value;
  for (i=0;i<array_size(a);i++)
    block_delete(*(Block**) array_get(a, i));
  array_delete(a);
  return 0;
}

Diagram* diagram_new(FeatureIndex *fi, const char *seqid, const Range *range,
                     Config *config)
{
  Diagram *diagram;
  Array *features = array_new(sizeof (GenomeNode*));
  int had_err;
  diagram = ma_malloc(sizeof (Diagram));
  diagram->tracks = hashtable_new(HASH_STRING, ma_free_func,
                                  (FreeFunc) track_delete);
  diagram->blocks = hashtable_new(HASH_DIRECT, NULL,
                                  (FreeFunc) blocklist_delete);
  diagram->nodeinfo = hashtable_new(HASH_DIRECT, NULL, NULL);
  diagram->nof_tracks = 0;
  diagram->config = config;
  diagram->range = *range;
  had_err = feature_index_get_features_for_range(fi, features, seqid, *range,
                                                 NULL);
  assert(!had_err); /* <fi> must contain <seqid> */
  diagram_build(diagram, features);
  array_delete(features);
  return diagram;
}

Range diagram_get_range(Diagram* diagram)
{
  assert(diagram);
  return diagram->range;
}

void diagram_set_config(Diagram *diagram, Config *config)
{
  assert(diagram && config);
  diagram->config = config;
}

Hashtable* diagram_get_tracks(const Diagram *diagram)
{
  assert(diagram);
  return diagram->tracks;
}

void diagram_get_lineinfo(const Diagram *diagram, TracklineInfo *tli)
{
  int had_err;
  assert(diagram);
  had_err = hashtable_foreach(diagram->tracks, diagram_add_tracklines,
                              tli, NULL);
  assert(!had_err); /* diagram_add_tracklines() is sane */
}

int diagram_get_number_of_tracks(const Diagram *diagram)
{
  assert(diagram);
  return diagram->nof_tracks;
}

static int blocklist_block_compare(const void *item1, const void *item2)
{
  assert(item1 && item2);
  return block_compare(*(Block**) item1, *(Block**) item2);
}

static int layout_tracks(void *key, void *value, void *data,
                         UNUSED Error *err)
{
  unsigned long i, max;
  Track *track;
  TrackTraverseInfo *tti = (TrackTraverseInfo*) data;
  GenomeFeatureType *gft = (GenomeFeatureType*) key;
  Array *list = (Array*) value;
  char *filename;
  Str *track_key;
  const char *type;
  Block *block;
  bool split;
  double tmp;
  assert(gft && list);

  /* to get a deterministic layout, we sort the Blocks for each type */
  array_sort(list, blocklist_block_compare);
  /* we take the basename of the filename to have nicer output in the
     generated graphic. this might lead to ``collapsed'' tracks, if two files
     with different paths have the same basename. */
  block = *(Block**) array_get(list, 0);
  filename = getbasename(genome_node_get_filename(
                                 block_get_top_level_feature(block)));
  track_key = track_key_new(filename, gft);
  ma_free(filename);
  type = genome_feature_type_get_cstr(gft);
  log_log("layouting track %s", type);

  if (!config_get_bool(tti->dia->config, "format", "split_lines", &split))
    split = true;
  if (split)
    if (!config_get_bool(tti->dia->config, type, "split_lines", &split))
      split = true;
  if (config_get_num(tti->dia->config, type, "max_num_lines", &tmp))
    max = tmp;
  else
    max = 50;

  /* For now, use the captions line breaker */
  track = track_new(track_key, max, split,
                    line_breaker_captions_new(tti->canvas));
  tti->dia->nof_tracks++;
  log_log("created track: %s, diagram has now %d tracks",
          str_get(track_key), diagram_get_number_of_tracks(tti->dia));
  for (i=0;i<array_size(list);i++)
  {
    block = *(Block**) array_get(list, i);
    track_insert_block(track, block);
  }
  hashtable_add(tti->dia->tracks, cstr_dup(str_get(track_key)), track);
  str_delete(track_key);
  return 0;
}

static int render_tracks(UNUSED void *key, void *value, void *data,
                     UNUSED Error *err)
{
  TrackTraverseInfo *tti = (TrackTraverseInfo*) data;
  Track *track = (Track*) value;
  int had_err = 0;
  assert(tti && track);
  log_log("rendering track");
  had_err = track_render((Track*) value, tti->canvas);
  return had_err;
}

int diagram_render(Diagram *dia, Canvas *canvas)
{
  int had_err = 0;
  TrackTraverseInfo tti;
  tti.dia = dia;
  tti.canvas = canvas;
  canvas_visit_diagram_pre(canvas, dia);
  hashtable_reset(dia->tracks);
  dia->nof_tracks = 0;
  (void) hashtable_foreach(dia->blocks, layout_tracks,
                           &tti, NULL);
  canvas_visit_diagram_post(canvas, dia);
  had_err = hashtable_foreach_ao(dia->tracks, render_tracks,
                                 &tti, NULL);

  log_log("finished rendering!\n");
  return had_err;
}

int diagram_unit_test(Error *err)
{
  FeatureTypeFactory *feature_type_factory;
  GenomeFeatureType *gene_type, *exon_type, *CDS_type;
  GenomeNode *gn1, *gn2, *ex1, *ex2, *ex3, *cds1;
  FeatureIndex *fi;
  Range r1, r2, r3, r4, r5, dr1, rs;
  Str *seqid1, *seqid2, *track_key;
  SequenceRegion *sr1, *sr2;
  int had_err=0;
  Config *cfg = NULL;
  Diagram *dia = NULL, *dia2 = NULL;
  Canvas *canvas = NULL;
  error_check(err);

  feature_type_factory = feature_type_factory_builtin_new();
  gene_type = feature_type_factory_create_gft(feature_type_factory, gft_gene);
  exon_type = feature_type_factory_create_gft(feature_type_factory, gft_exon);
  CDS_type = feature_type_factory_create_gft(feature_type_factory, gft_CDS);

  /* generating some ranges */
  r1.start=100UL; r1.end=1000UL;
  r2.start=100UL; r2.end=300UL;
  r3.start=500UL; r3.end=1000UL;
  r4.start=600UL; r4.end=1200UL;
  r5.start=600UL; r5.end=1000UL;
  rs.start=100UL; rs.end=1200UL;

  /* generating sequence IDs */
  seqid1 = str_new_cstr("test1");
  seqid2 = str_new_cstr("test2");

  sr1 = (SequenceRegion*) sequence_region_new(seqid1, rs);
  sr2 = (SequenceRegion*) sequence_region_new(seqid2, rs);

  gn1 = genome_feature_new(gene_type, r1, STRAND_UNKNOWN);
  genome_node_set_seqid((GenomeNode*) gn1, seqid1);

  gn2 = genome_feature_new(gene_type, r4, STRAND_UNKNOWN);
  genome_node_set_seqid((GenomeNode*) gn2, seqid2);

  ex1 = genome_feature_new(exon_type, r2, STRAND_UNKNOWN);
  genome_node_set_seqid((GenomeNode*) ex1, seqid1);

  ex2 = genome_feature_new(exon_type, r3, STRAND_UNKNOWN);
  genome_node_set_seqid((GenomeNode*) ex2, seqid1);

  ex3 = genome_feature_new(exon_type, r4, STRAND_UNKNOWN);
  genome_node_set_seqid((GenomeNode*) ex3, seqid2);

  cds1 = genome_feature_new(CDS_type, r5, STRAND_UNKNOWN);
  genome_node_set_seqid((GenomeNode*) cds1, seqid2);

  /* determine the structure of our feature tree */
  genome_node_is_part_of_genome_node(gn1, ex1);
  genome_node_is_part_of_genome_node(gn1, ex2);
  genome_node_is_part_of_genome_node(gn2, ex3);
  genome_node_is_part_of_genome_node(gn2, cds1);

  /* create a new feature index on which we can perfom some tests */
  fi = feature_index_new();

  /* add a sequence region the feature index and test if it has really been
     added */
  feature_index_add_sequence_region(fi, sr1);
  feature_index_add_sequence_region(fi, sr2);

  /* add features to every sequence region */
  feature_index_add_genome_feature(fi, (GenomeFeature*) gn1);
  feature_index_add_genome_feature(fi, (GenomeFeature*) gn2);

  /* set the Range for the diagram */
  dr1.start = 400UL;
  dr1.end   = 900UL;

  /* create a config object */
  if (!had_err) {
    if (!(cfg = config_new(false, err)))
      had_err = -1;
  }

  /* create a diagram object and test it */
  if (!had_err)
    dia = diagram_new(fi, "test1", &dr1, cfg);

  ensure(had_err, dia->config);
  ensure(had_err, dia->range.start == 400UL);
  ensure(had_err, dia->range.end == 900UL);

  if (!had_err)
  {
    canvas = canvas_new(cfg, GRAPHICS_PDF, 600, NULL);
    diagram_render(dia, canvas);
  }

  if (!had_err &&
      !config_get_bool(dia->config, "gene", "collapse_to_parent", false)) {
    track_key = track_key_new("generated", gene_type);
    ensure(had_err, hashtable_get(dia->tracks, str_get(track_key)));
    str_delete(track_key);
  }

  if (!had_err &&
      !config_get_bool(dia->config, "exon", "collapse_to_parent", false)) {
    track_key = track_key_new("generated", exon_type);
    ensure(had_err, hashtable_get(dia->tracks, str_get(track_key)));
    str_delete(track_key);
  }
  ensure(had_err, range_compare(diagram_get_range(dia),dr1) == 0);

  /* create a diagram object and test it */
  if (!had_err) {
    dia2 = diagram_new(fi, "test2", &dr1, cfg);
    ensure(had_err, dia->range.start == 400UL);
    ensure(had_err, dia->range.end == 900UL);
  }

  if (!had_err &&
      !config_get_bool(dia2->config, "gene", "collapse_to_parent", false)) {
    diagram_render(dia2, canvas);
    track_key = track_key_new("generated", gene_type);
    ensure(had_err, hashtable_get(dia2->tracks, str_get(track_key)));
    str_delete(track_key);
  }

  if (!had_err &&
      !config_get_bool(dia2->config, "exon", "collapse_to_parent", false)) {
    track_key = track_key_new("generated", exon_type);
    ensure(had_err, hashtable_get(dia2->tracks, str_get(track_key)));
    str_delete(track_key);
  }

  if (!had_err &&
      !config_get_bool(dia2->config, "CDS", "collapse_to_parent", false)) {
    track_key = track_key_new("generated", CDS_type);
    ensure(had_err, hashtable_get(dia2->tracks, str_get(track_key)));
    str_delete(track_key);
  }
  ensure(had_err, range_compare(diagram_get_range(dia),dr1) == 0);

  /* delete all generated objects */
  config_delete(cfg);
  diagram_delete(dia);
  diagram_delete(dia2);
  canvas_delete(canvas);
  feature_index_delete(fi);
  genome_node_rec_delete(gn1);
  genome_node_rec_delete(gn2);
  genome_node_rec_delete((GenomeNode*) sr1);
  genome_node_rec_delete((GenomeNode*) sr2);
  str_delete(seqid1);
  str_delete(seqid2);
  feature_type_factory_delete(feature_type_factory);

  return had_err;
}

void diagram_delete(Diagram *diagram)
{
  if (!diagram) return;
  hashtable_delete(diagram->tracks);
  hashtable_delete(diagram->blocks);
  hashtable_delete(diagram->nodeinfo);
  ma_free(diagram);
}
