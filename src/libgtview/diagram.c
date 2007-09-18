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
#include "libgtcore/warning.h"
#include "libgtcore/str.h"
#include "libgtcore/undef.h"
#include "libgtext/genome_node.h"
#include "libgtext/genome_feature.h"
#include "libgtext/genome_feature_type.h"
#include "libgtview/diagram.h"
#include "libgtview/feature_index.h"
#include "libgtview/track.h"

/* used to separate a filename from the type in a track name */
#define FILENAME_TYPE_SEPARATOR  '|'

struct Diagram {
  Hashtable *tracks;
  Hashtable *nodeinfo;
  int nof_tracks;
  Config *config;
  Range range;
};

/* holds a Block with associated type */
typedef struct {
  GenomeFeatureType gft;
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

static BlockTuple* blocktuple_new(GenomeFeatureType gft, Block *block, Env* env)
{
  BlockTuple *bt;
  assert(block);
  bt = env_ma_malloc(env, sizeof (BlockTuple));
  bt->gft = gft;
  bt->block = block;
  return bt;
}

static NodeInfoElement* get_or_create_node_info(Diagram *d,
                                                GenomeNode *node,
                                                Env* env)
{
  NodeInfoElement *ni;
  assert(d && node);
  ni = hashtable_get(d->nodeinfo, node);
  if (ni == NULL) {
    NodeInfoElement *new_ni = env_ma_malloc(env, sizeof (NodeInfoElement));
    new_ni->blocktuples = array_new(sizeof (BlockTuple*), env);
    hashtable_add(d->nodeinfo, node, new_ni, env);
    ni = new_ni;
  }
  return ni;
}

static Block* find_block_for_type(NodeInfoElement* ni,
                                  GenomeFeatureType gft,
                                  Env *env)
{
  Block *block = NULL;
  unsigned long i;
  assert(ni);
  for (i = 0; i < array_size(ni->blocktuples); i++) {
    BlockTuple *bt = *(BlockTuple**) array_get(ni->blocktuples, i);
    if (bt->gft == gft)
    {
      env_log_log(env, "    found root block with type %s, inserting...\n",
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
  assert(gn);
  if (!(ret = genome_feature_get_attribute(gn, "Name"))) {
    if (!(ret = genome_feature_get_attribute(gn, "ID")))
      ret = "-";
  }
  return ret;
}

static void add_to_current(Diagram *d, GenomeNode *node,
                           GenomeNode *parent, Env *env)
{
  NodeInfoElement *ni;
  Block *block;
  BlockTuple *bt;
  Str *caption;

  env_log_log(env, "  calling add_to_current\n");

  assert(d && node);

  /* Lookup node info and set itself as parent */
  ni = get_or_create_node_info(d, node, env);
  ni->parent = node;
  /* create new Block tuple and add to node info */
  block = block_new_from_node(node,env);
  /* assign block caption */
  caption = str_new_cstr("",env);
  if (parent) {
    if (genome_node_has_children(parent))
      str_append_cstr(caption, get_node_name_or_id(parent), env);
    else
      str_append_cstr(caption, "-", env);
    str_append_cstr(caption, "/", env);
  }
  str_append_cstr(caption, get_node_name_or_id(node), env);
  block_set_caption(block, caption);
  /* insert node into block */
  block_insert_element(block, node, d->config, env);
  bt = blocktuple_new(genome_feature_get_type((GenomeFeature*) node),
                      block,
                      env);
  array_add(ni->blocktuples, bt, env);
}

static void add_to_parent(Diagram *d, GenomeNode *node,
                          GenomeNode* parent, Env *env)
{
  Block *block = NULL;
  NodeInfoElement *par_ni, *ni;

  assert(d && parent && node);

  env_log_log(env,"  calling add_to_parent: %s -> %s\n",
            genome_feature_type_get_cstr(
              genome_feature_get_type((GenomeFeature*) node)),
            genome_feature_type_get_cstr(
              genome_feature_get_type((GenomeFeature*) parent)));

  par_ni = get_or_create_node_info(d, parent, env);
  ni = get_or_create_node_info(d, node, env);
  ni->parent = parent;

  /* try to find the right block to insert */
  block = find_block_for_type(par_ni,
                              genome_feature_get_type((GenomeFeature*) node),
                              env);
  /* no fitting block was found, create a new one */
  if (block == NULL) {
    BlockTuple *bt;
    Str *caption;
    block = block_new_from_node(parent, env);
    /* assign caption */
    caption = str_new_cstr("",env);
    if (genome_node_has_children(parent))
      str_append_cstr(caption, get_node_name_or_id(parent), env);
    else
      str_append_cstr(caption, "-", env);
    str_append_cstr(caption, "/", env);
    str_append_cstr(caption, get_node_name_or_id(node), env);
    block_set_caption(block, caption);
    /* add block to nodeinfo */
    bt = blocktuple_new(genome_feature_get_type((GenomeFeature*) node),
                        block,
                        env);
    array_add(par_ni->blocktuples, bt, env);

    env_log_log(env, "    added %s to (%s) block for node %s\n",
                genome_feature_type_get_cstr(
                  genome_feature_get_type((GenomeFeature*) node)),
                genome_feature_type_get_cstr(
                  genome_feature_get_type((GenomeFeature*) parent)),
                genome_feature_get_attribute(parent, "ID" ));
  }
  /* now we have a block to insert into */
  block_insert_element(block, node, d->config, env);
}

static void add_recursive(Diagram *d, GenomeNode *node,
                          GenomeNode* parent, GenomeNode *original_node,
                          Env *env)
{
  NodeInfoElement *ni;

  assert(d && node &&
           parent && original_node);

  env_log_log(env ,"  calling add_recursive: %s (%s) -> %s (%s)\n",
              genome_feature_type_get_cstr(
                genome_feature_get_type((GenomeFeature*) node)),
              genome_feature_get_attribute(node,"ID"),
              genome_feature_type_get_cstr(
                genome_feature_get_type((GenomeFeature*) parent)),
              genome_feature_get_attribute(parent,"ID"));

  ni = get_or_create_node_info(d, node, env);

  /* end of recursion, insert into target block */
  if (parent == node) {
    Block *block ;
    BlockTuple *bt;
    /* try to find the right block to insert */
    block = find_block_for_type(ni,
                                genome_feature_get_type((GenomeFeature*) node),
                                env);
    if (block == NULL) {
      block = block_new_from_node(node,env);
      bt = blocktuple_new(genome_feature_get_type((GenomeFeature*) node),
                          block,
                          env);
      array_add(ni->blocktuples, bt, env);
    }
    block_insert_element(block, original_node, d->config, env);
  }
  else {
    /* not at target type block yet, set up reverse entry and follow */
    NodeInfoElement *parent_ni;
    /* set up reverse entry */
    ni->parent = parent;
    /* recursively call with parent node and its parent */
    parent_ni = hashtable_get(d->nodeinfo, parent);
    add_recursive(d, parent, parent_ni->parent, original_node, env);
  }
}

static void process_node(Diagram *d, GenomeNode *node, GenomeNode *parent,
                         Env *env)
{
  Range elem_range;
  bool collapse, do_not_overlap=false;
  const char *feature_type;

  assert(d && node);

  /* discard elements that do not overlap with visible range */
  elem_range = genome_node_get_range(node);
  if (!range_overlap(d->range, elem_range))
    return;

  /* check if this is a collapsing type */
  feature_type = genome_feature_type_get_cstr(
                   genome_feature_get_type((GenomeFeature*) node));
  collapse = config_cstr_in_list(d->config,"collapse","to_parent",feature_type,
                                 env);

  /* check if direct children overlap */
  if (parent)
    do_not_overlap =
      genome_node_direct_children_do_not_overlap_st(parent, node, env);

  env_log_log(env ,"processing node: %s (%s)\n",
              feature_type,
              genome_feature_get_attribute(node, "ID" ));

  /* decide how to continue: */
  if (collapse) {
    /* collapsing features recursively search their target blocks */
    if (!do_not_overlap)
      warning("collapsing %s features overlap "
              "and will be missing in the %s block!",
              feature_type,
              genome_feature_get_attribute(parent, "ID" ));
    add_recursive(d, node, parent, node, env);
  }
  else if (do_not_overlap && genome_node_number_of_children(parent) > 1)
  {
    /* group non-overlapping child nodes of a non-collapsing type by parent */
    add_to_parent(d, node, parent, env);
  }
  else {
    /* nodes that belong into their own track and block */
    add_to_current(d, node, parent, env);
  }

  /* we can now assume that this node has been processed into the reverse
     lookup structure */
  assert(hashtable_get(d->nodeinfo, node));
}

static int diagram_add_tracklines(void *key, void *value, void *data,
                                  Env* env)
{
  int *add = (int*) data;
  *add += track_get_number_of_lines((Track*) value);
  return 0;
}

static int diagram_track_delete(void *key, void *value, void *data, Env* env)
{
  track_delete((Track*) value, env);
  return 0;
}

static int visit_child(GenomeNode* gn, void* genome_node_children, Env* env)
{
  NodeTraverseInfo* genome_node_info;
  genome_node_info = (NodeTraverseInfo*) genome_node_children;
  if (genome_node_has_children(gn))
  {
    GenomeNode *oldparent = genome_node_info->parent;
    process_node(genome_node_info->diagram,
                 gn,
                 genome_node_info->parent,
                 env);
    genome_node_info->parent = gn;
    (void) genome_node_traverse_direct_children(gn, genome_node_info,
                                                visit_child, env);
    genome_node_info->parent = oldparent;
  }
  else
    process_node(genome_node_info->diagram, gn, genome_node_info->parent, env);
  return 0;
}

static Str* track_key_new(const char *filename, GenomeFeatureType type,
                          Env *env)
{
  Str *track_key;
  env_error_check(env);
  track_key = str_new_cstr(filename, env);
  str_append_char(track_key, FILENAME_TYPE_SEPARATOR, env);
  str_append_cstr(track_key, genome_feature_type_get_cstr(type), env);
  return track_key;
}

/* Create Tracks for all Blocks in the diagram. */
static int collect_blocks(void *key, void *value, void *data, Env *env)
{
  NodeInfoElement *ni = (NodeInfoElement*) value;
  Diagram *diagram = (Diagram*) data;
  unsigned long i = 0;

  if (array_size(ni->blocktuples) > 0)
    env_log_log(env, "collecting blocks from %lu tuples...\n",
                array_size(ni->blocktuples));

  for (i = 0; i < array_size(ni->blocktuples); i++) {
    Track *track;
    const char *filename;
    Str *track_key;
    BlockTuple *bt = *(BlockTuple**) array_get(ni->blocktuples, i);
    filename = genome_node_get_filename(ni->parent);
    track_key = track_key_new(filename, bt->gft, env);
    track = hashtable_get(diagram->tracks, str_get(track_key));

    if (track == NULL) {
      track = track_new(track_key, env);
      hashtable_add(diagram->tracks, cstr_dup(str_get(track_key), env), track,
                    env);
      diagram->nof_tracks++;
      env_log_log(env , "created track: %s, diagram has now %d tracks\n",
                  str_get(track_key), diagram_get_number_of_tracks(diagram));
    }
    track_insert_block(track, bt->block, env);
    env_log_log(env, "inserted block %s into track %s\n",
                str_get(block_get_caption(bt->block)), str_get(track_key));
    str_delete(track_key, env);
    env_ma_free(bt, env);
  }
  array_delete(ni->blocktuples, env);
  env_ma_free(ni, env);
  return 0;
}

/* Traverse a genome node tree with depth first search. */
static void traverse_genome_nodes(GenomeNode *gn, void *genome_node_children,
                                  Env *env)
{
  NodeTraverseInfo* genome_node_info;
  assert(genome_node_children);
  genome_node_info = (NodeTraverseInfo*) genome_node_children;
  genome_node_info->parent = gn;
  /* handle root nodes */
  process_node(genome_node_info->diagram, gn, NULL, env);
  if (genome_node_has_children(gn)) {
    (void) genome_node_traverse_direct_children(gn, genome_node_info,
                                                visit_child, env);
  }
}

static void diagram_build(Diagram *diagram, Array *features, Env *env)
{
  unsigned long i = 0;
  NodeTraverseInfo genome_node_children;
  genome_node_children.diagram = diagram;
  /* do node traversal for each root feature */
  for (i = 0; i < array_size(features); i++) {
    GenomeNode *current_root = *(GenomeNode**) array_get(features,i);
    traverse_genome_nodes(current_root, &genome_node_children, env);
  }
  /* collect blocks from nodeinfo structures and create the tracks */
  (void) hashtable_foreach(diagram->nodeinfo, collect_blocks, diagram, env);
}

Diagram* diagram_new(FeatureIndex *fi, Range range, const char *seqid,
                     Config *config, Env *env)
{
  Diagram *diagram;
  Array *features = array_new(sizeof (GenomeNode*), env);
  int had_err;
  env_error_check(env);
  diagram = env_ma_malloc(env, sizeof (Diagram));
  diagram->tracks = hashtable_new(HASH_STRING, env_ma_free_func, NULL, env);
  diagram->nodeinfo = hashtable_new(HASH_DIRECT, NULL, NULL, env);
  diagram->nof_tracks = 0;
  diagram->config = config;
  diagram->range = range;
  had_err = feature_index_get_features_for_range(fi,
                                                 features,
                                                 seqid,
                                                 range,
                                                 env);
  assert(!had_err); /* <fi> must contain <seqid> */
  diagram_build(diagram, features, env);
  array_delete(features, env);
  return diagram;
}

Range diagram_get_range(Diagram* diagram)
{
  assert(diagram);
  return diagram->range;
}

void diagram_set_config(Diagram *diagram, Config *config, Env *env)
{
  assert(diagram && config);
  diagram->config = config;
}

Hashtable* diagram_get_tracks(const Diagram *diagram)
{
  assert(diagram);
  return diagram->tracks;
}

int diagram_get_total_lines(const Diagram *diagram, Env *env)
{
  int total_lines = 0;
  assert(diagram);
  (void) hashtable_foreach(diagram->tracks, diagram_add_tracklines,
                           &total_lines, env);
  return total_lines;
}

int diagram_get_number_of_tracks(const Diagram *diagram)
{
  assert(diagram);
  return diagram->nof_tracks;
}

int diagram_unit_test(Env *env)
{
  GenomeNode *gn1, *gn2, *ex1, *ex2, *ex3, *cds1;
  FeatureIndex *fi;
  Range r1, r2, r3, r4, r5, dr1, rs;
  Str *seqid1, *seqid2, *track_key;
  SequenceRegion *sr1, *sr2;
  int had_err=0;
  Config *cfg = NULL;
  Diagram *dia = NULL, *dia2 = NULL;

  /* generating some ranges */
  r1.start=100UL; r1.end=1000UL;
  r2.start=100UL; r2.end=300UL;
  r3.start=500UL; r3.end=1000UL;
  r4.start=600UL; r4.end=1200UL;
  r5.start=600UL; r5.end=1000UL;
  rs.start=100UL; rs.end=1200UL;

  /* generating sequence IDs */
  seqid1 = str_new_cstr("test1", env);
  seqid2 = str_new_cstr("test2", env);

  sr1 = (SequenceRegion*) sequence_region_new(seqid1, rs, NULL, 0, env);
  sr2 = (SequenceRegion*) sequence_region_new(seqid2, rs, NULL, 0, env);

  gn1 = genome_feature_new(gft_gene, r1, STRAND_UNKNOWN, NULL,
                           UNDEF_ULONG, env);
  genome_node_set_seqid((GenomeNode*) gn1, seqid1, env);

  gn2 = genome_feature_new(gft_gene, r4, STRAND_UNKNOWN, NULL,
                           UNDEF_ULONG, env);
  genome_node_set_seqid((GenomeNode*) gn2, seqid2, env);

  ex1 = genome_feature_new(gft_exon, r2, STRAND_UNKNOWN, NULL,
                           UNDEF_ULONG, env);
  genome_node_set_seqid((GenomeNode*) ex1, seqid1, env);

  ex2 = genome_feature_new(gft_exon, r3, STRAND_UNKNOWN, NULL,
                           UNDEF_ULONG, env);
  genome_node_set_seqid((GenomeNode*) ex2, seqid1, env);

  ex3 = genome_feature_new(gft_exon, r4, STRAND_UNKNOWN, NULL,
                           UNDEF_ULONG, env);
  genome_node_set_seqid((GenomeNode*) ex3, seqid2, env);

  cds1 = genome_feature_new(gft_CDS, r5, STRAND_UNKNOWN, NULL,
                            UNDEF_ULONG, env);
  genome_node_set_seqid((GenomeNode*) cds1, seqid2, env);

  /* determine the structure of our feature tree */
  genome_node_is_part_of_genome_node(gn1, ex1, env);
  genome_node_is_part_of_genome_node(gn1, ex2, env);
  genome_node_is_part_of_genome_node(gn2, ex3, env);
  genome_node_is_part_of_genome_node(gn2, cds1, env);

  /* create a new feature index on which we can perfom some tests */
  fi = feature_index_new(env);

  /* add a sequence region the feature index and test if it has really been
     added */
  feature_index_add_sequence_region(fi, sr1, env);
  feature_index_add_sequence_region(fi, sr2, env);

  /* add features to every sequence region */
  feature_index_add_genome_feature(fi, (GenomeFeature*) gn1, env);
  feature_index_add_genome_feature(fi, (GenomeFeature*) gn2, env);

  /* set the Range for the diagram */
  dr1.start = 400UL;
  dr1.end   = 900UL;

  /* create a config object */
  if (!had_err)
    cfg = config_new(false, env);

  /* create a diagram object and test it */
  if (!had_err)
    dia = diagram_new(fi, dr1, "test1", cfg, env);

  ensure(had_err, dia->config);
  ensure(had_err, dia->range.start == 400UL);
  ensure(had_err, dia->range.end == 900UL);
  if (!had_err &&
      !config_cstr_in_list(dia->config,"collapse","to_parent","gene", env)) {
    track_key = track_key_new("generated", gft_gene, env);
    ensure(had_err, hashtable_get(dia->tracks, str_get(track_key)));
    str_delete(track_key, env);
  }

  if (!had_err &&
      !config_cstr_in_list(dia->config,"collapse","to_parent","exon", env)) {
    track_key = track_key_new("generated", gft_exon, env);
    ensure(had_err, hashtable_get(dia->tracks, str_get(track_key)));
    str_delete(track_key, env);
  }
  ensure(had_err, range_compare(diagram_get_range(dia),dr1) == 0);

  /* create a diagram object and test it */
  if (!had_err) {
    dia2 = diagram_new(fi,dr1,"test2",cfg,env);
    ensure(had_err, dia->range.start == 400UL);
    ensure(had_err, dia->range.end == 900UL);
  }

  if (!had_err &&
      !config_cstr_in_list(dia2->config,"collapse","to_parent","gene", env)) {
    track_key = track_key_new("generated", gft_gene, env);
    ensure(had_err, hashtable_get(dia2->tracks, str_get(track_key)));
    str_delete(track_key, env);
  }

  if (!had_err &&
      !config_cstr_in_list(dia2->config,"collapse","to_parent","exon", env)) {
    track_key = track_key_new("generated", gft_exon, env);
    ensure(had_err, hashtable_get(dia2->tracks, str_get(track_key)));
    str_delete(track_key, env);
  }

  if (!had_err &&
      !config_cstr_in_list(dia2->config,"collapse","to_parent","CDS", env)) {
    track_key = track_key_new("generated", gft_CDS, env);
    ensure(had_err, hashtable_get(dia2->tracks, str_get(track_key)));
    str_delete(track_key, env);
  }
  ensure(had_err, range_compare(diagram_get_range(dia),dr1) == 0);

  /* delete all generated objects */
  config_delete(cfg, env);
  diagram_delete(dia,env);
  diagram_delete(dia2,env);
  feature_index_delete(fi, env);
  genome_node_rec_delete(gn1, env);
  genome_node_rec_delete(gn2, env);
  genome_node_rec_delete((GenomeNode*) sr1, env);
  genome_node_rec_delete((GenomeNode*) sr2, env);
  str_delete(seqid1, env);
  str_delete(seqid2, env);

  return had_err;
}

void diagram_delete(Diagram *diagram, Env *env)
{
  if (!diagram) return;
  (void) hashtable_foreach(diagram->tracks, diagram_track_delete, NULL, env);
  hashtable_delete(diagram->tracks, env);
  hashtable_delete(diagram->nodeinfo, env);
  env_ma_free(diagram, env);
}
