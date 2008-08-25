/*
  Copyright (c) 2006-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2008 Center for Bioinformatics, University of Hamburg

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

#include <assert.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "libgtcore/cstr.h"
#include "libgtcore/hashmap.h"
#include "libgtcore/ma.h"
#include "libgtcore/parseutils.h"
#include "libgtcore/splitter.h"
#include "libgtcore/undef.h"
#include "libgtcore/unused.h"
#include "libgtcore/warning.h"
#include "libgtext/comment.h"
#include "libgtext/feature_info.h"
#include "libgtext/genome_feature.h"
#include "libgtext/genome_node.h"
#include "libgtext/genome_node_iterator.h"
#include "libgtext/gff3_escaping.h"
#include "libgtext/gff3_parser.h"
#include "libgtext/mapping.h"
#include "libgtext/sequence_region.h"

struct GFF3Parser {
  FeatureInfo *feature_info;
  Hashmap *seqid_to_ssr_mapping, /* maps seqids to simple sequence regions */
    *source_to_str_mapping,
    *undefined_sequence_regions; /* contains all (automatically created)
                                  * sequence regions */
  bool incomplete_node, /* at least on node is potentially incomplete */
       checkids,
       tidy,
       fasta_parsing; /* parser is in FASTA parsing mode */
  long offset;
  Mapping *offset_mapping;
  FeatureTypeFactory *feature_type_factory;
  unsigned int last_terminator; /* line number of the last terminator */
};

typedef struct {
  GenomeNode *sequence_region; /* the automatically created sequence region */
  Array *genome_features; /* the genome features which belong to this reagion */
} AutomaticSequenceRegion;

static AutomaticSequenceRegion* automatic_sequence_region_new(void)
{
  AutomaticSequenceRegion *auto_sr;
  auto_sr = ma_malloc(sizeof (AutomaticSequenceRegion));
  auto_sr->genome_features = array_new(sizeof (GenomeFeature*));
  return auto_sr;
}

static void automatic_sequence_region_delete(AutomaticSequenceRegion *auto_sr)
{
  unsigned long i;
  if (!auto_sr) return;
  genome_node_delete(auto_sr->sequence_region);
  for (i = 0; i < array_size(auto_sr->genome_features); i++)
    genome_node_delete(*(GenomeNode**) array_get(auto_sr->genome_features, i));
  array_delete(auto_sr->genome_features);
  ma_free(auto_sr);
}

typedef struct {
  Str *seqid_str;
  Range range;
  unsigned int line_number;
} SimpleSequenceRegion;

static SimpleSequenceRegion* simple_sequence_region_new(const char *seqid,
                                                        Range range,
                                                        unsigned int
                                                        line_number)
{
  SimpleSequenceRegion *ssr = ma_malloc(sizeof *ssr);
  ssr->seqid_str = str_new_cstr(seqid);
  ssr->range = range;
  ssr->line_number = line_number;
  return ssr;
}

static void simple_sequence_region_delete(SimpleSequenceRegion *ssr)
{
  if (!ssr) return;
  str_delete(ssr->seqid_str);
  ma_free(ssr);
}

GFF3Parser* gff3parser_new(bool checkids,
                           FeatureTypeFactory *feature_type_factory)
{
  GFF3Parser *parser;
  assert(feature_type_factory);
  parser = ma_malloc(sizeof (GFF3Parser));
  parser->feature_info = feature_info_new();
  parser->seqid_to_ssr_mapping = hashmap_new(
    HASH_STRING, NULL, (FreeFunc) simple_sequence_region_delete);
  parser->source_to_str_mapping = hashmap_new(HASH_STRING, NULL,
                                              (FreeFunc) str_delete);
  parser->undefined_sequence_regions = hashmap_new(
    HASH_STRING, NULL, (FreeFunc) automatic_sequence_region_delete);
  parser->incomplete_node = false;
  parser->checkids = checkids;
  parser->tidy = false;
  parser->fasta_parsing = false;
  parser->offset = UNDEF_LONG;
  parser->offset_mapping = NULL;
  parser->feature_type_factory = feature_type_factory;
  parser->last_terminator = 0;
  return parser;
}

void gff3parser_set_offset(GFF3Parser *parser, long offset)
{
  assert(parser);
  assert(!parser->offset_mapping);
  parser->offset = offset;
}

int gff3parser_set_offsetfile(GFF3Parser *parser, Str *offsetfile, Error *err)
{
  error_check(err);
  assert(parser);
  assert(parser->offset == UNDEF_LONG);
  parser->offset_mapping = mapping_new(offsetfile, "offsets",
                                       MAPPINGTYPE_INTEGER, err);
  if (parser->offset_mapping)
    return 0;
  return -1;
}

void gff3parser_enable_tidy_mode(GFF3Parser *parser)
{
  assert(parser);
  parser->tidy = true;
}

static int add_offset_if_necessary(Range *range, GFF3Parser *parser,
                                   const char *seqid, Error *err)
{
  long offset;
  int had_err = 0;
  error_check(err);
  if (parser->offset != UNDEF_LONG)
    *range = range_offset(*range, parser->offset);
  else if (parser->offset_mapping) {
    had_err = mapping_map_integer(parser->offset_mapping, &offset, seqid, err);
    if (!had_err)
      *range = range_offset(*range, offset);
  }
  return had_err;
}

static int parse_target_attribute(const char *value, Str *target_id,
                                  Range *target_range, Strand *target_strand,
                                  const char *filename,
                                  unsigned int line_number, Error *err)
{
  unsigned long num_of_tokens;
  Str *unescaped_target;
  char *escaped_target;
  Strand parsed_strand;
  Splitter *splitter;
  Range parsed_range;
  int had_err = 0;
  error_check(err);
  assert(value && filename);
  splitter = splitter_new();
  unescaped_target = str_new();
  escaped_target = cstr_dup(value);
  splitter_split(splitter, escaped_target, strlen(escaped_target), ' ');
  num_of_tokens = splitter_size(splitter);
  if (!(num_of_tokens == 3 || num_of_tokens == 4)) {
    error_set(err, "Target attribute value '%s' on line %u in file \"%s\" "
              "must have 3 or 4 blank separated entries", value, line_number,
              filename);
    had_err = -1;
  }
  /* parse target id */
  if (!had_err) {
    had_err = gff3_unescape(unescaped_target, splitter_get_token(splitter, 0),
                            strlen(splitter_get_token(splitter, 0)), err);
  }
  if (!had_err && target_id) str_append_str(target_id, unescaped_target);
  /* parse target range */
  if (!had_err) {
    had_err = parse_range(&parsed_range, splitter_get_token(splitter, 1),
                          splitter_get_token(splitter, 2), line_number,
                          filename, err);
  }
  if (!had_err && target_range) *target_range = parsed_range;
  /* parse target strand (if given) */
  if (!had_err) {
    if (splitter_size(splitter) == 4) {
      had_err = parse_strand(&parsed_strand, splitter_get_token(splitter, 3),
                             line_number, filename, err);
      if (!had_err && target_strand)
        *target_strand = parsed_strand;
    }
    else if (target_strand)
      *target_strand = NUM_OF_STRAND_TYPES; /* undefined */
  }
  ma_free(escaped_target);
  str_delete(unescaped_target);
  splitter_delete(splitter);
  return had_err;
}

int gff3parser_parse_target_attributes(const char *values,
                                       unsigned long *num_of_targets,
                                       Str *first_target_id,
                                       Range *first_target_range,
                                       Strand *first_target_strand,
                                       const char *filename,
                                       unsigned int line_number, Error *err)
{
  Splitter *splitter;
  unsigned long i;
  char *targets;
  int had_err = 0;
  error_check(err);
  assert(values && filename);
  targets = cstr_dup(values);
  splitter = splitter_new();
  splitter_split(splitter, targets, strlen(targets), ',');
  if (num_of_targets)
    *num_of_targets = splitter_size(splitter);
  for (i = 0; !had_err && i < splitter_size(splitter); i++) {
    had_err = parse_target_attribute(splitter_get_token(splitter, i),
                                     i ? NULL : first_target_id,
                                     i ? NULL : first_target_range,
                                     i ? NULL : first_target_strand, filename,
                                     line_number, err);
  }
  ma_free(targets);
  splitter_delete(splitter);
  return had_err;
}

static int set_seqid(GenomeNode *genome_feature, const char *seqid, Range range,
                     AutomaticSequenceRegion **auto_sr, GFF3Parser *parser,
                     const char *filename, unsigned int line_number, Error *err)
{
  bool seqid_str_created = false;
  SimpleSequenceRegion *ssr;
  Str *seqid_str = NULL;
  int had_err = 0;

  error_check(err);
  assert(genome_feature);

  ssr = hashmap_get(parser->seqid_to_ssr_mapping, seqid);
  if (!ssr) {
    /* sequence region has not been previously introduced -> check if one has
       already been created automatically */
    *auto_sr = hashmap_get(parser->undefined_sequence_regions, seqid);
    if (!*auto_sr) {
      /* sequence region has not been createad automatically -> do it now */
      warning("seqid \"%s\" on line %u in file \"%s\" has not been "
              "previously introduced with a \"%s\" line, create such a line "
              "automatically", seqid, line_number, filename,
              GFF_SEQUENCE_REGION);
      *auto_sr = automatic_sequence_region_new();
      seqid_str = str_new_cstr(seqid);
      seqid_str_created = true;
      (*auto_sr)->sequence_region = sequence_region_new(seqid_str, range);
      hashmap_add(parser->undefined_sequence_regions, str_get(seqid_str),
                  *auto_sr);
    }
    else {
      /* get seqid string */
      seqid_str = genome_node_get_seqid((*auto_sr)->sequence_region);
      /* update the range of the sequence region */
      genome_node_set_range((*auto_sr)->sequence_region,
                            range_join(range,
                                       genome_node_get_range((*auto_sr)
                                                         ->sequence_region)));
    }
  }
  else {
    seqid_str = ssr->seqid_str;
    /* perform range check */
    if (!range_contains(ssr->range, range)) {
      error_set(err, "range (%lu,%lu) of feature on line %u in file \"%s\" "
                "is not contained in range (%lu,%lu) of corresponding "
                "sequence region on line %u", range.start, range.end,
                line_number, filename, ssr->range.start, ssr->range.end,
                ssr->line_number);
      had_err = -1;
    }
  }
  if (!had_err) {
    assert(seqid_str);
    genome_node_set_seqid(genome_feature, seqid_str);
  }
  if (seqid_str_created)
    str_delete(seqid_str);

  return had_err;
}

typedef struct {
  GenomeNode *node_to_replace,
             *replacing_node;
} ReplaceInfo;

static int replace_func(void **elem, void *info, UNUSED Error *err)
{
  ReplaceInfo *replace_info = info;
  GenomeNode **node = (GenomeNode**) elem;
  error_check(err);
  assert(node && replace_info);
  if (*node == replace_info->node_to_replace) {
    *node = replace_info->replacing_node;
    return 1;
  }
  return 0;
}

static void replace_node(GenomeNode *node_to_replace,
                         GenomeNode *replacing_node, Queue *genome_nodes,
                         AutomaticSequenceRegion *auto_sr)
{
  ReplaceInfo replace_info;
  int rval;
  assert(node_to_replace && replacing_node && genome_nodes);
  replace_info.node_to_replace = node_to_replace;
  replace_info.replacing_node = replacing_node;
  if (auto_sr) {
    rval = array_iterate(auto_sr->genome_features,
                         (ArrayProcessor) replace_func, &replace_info, NULL);
    assert(rval == 1);
  }
  else {
    rval = queue_iterate(genome_nodes, replace_func, &replace_info, NULL);
    assert(rval == 1);
  }
}

static void remove_node(GenomeNode *genome_node, Queue *genome_nodes,
                        AutomaticSequenceRegion *auto_sr)
{
  unsigned long i;
  assert(genome_node && genome_nodes);
  if (auto_sr) {
    for (i = 0; i < array_size(auto_sr->genome_features); i++) {
      if (genome_node == *(GenomeNode**) array_get(auto_sr->genome_features, i))
        break;
    }
    assert(i < array_size(auto_sr->genome_features));
    array_rem(auto_sr->genome_features, i);
  }
  else {
    /* XXX: use proper queue_rem() method */
    Queue *tmp_queue = queue_new();
    while (queue_size(genome_nodes)) {
      GenomeNode *gn = queue_get(genome_nodes);
      if (gn != genome_node)
        queue_add(tmp_queue, gn);
    }
    assert(!queue_size(genome_nodes));
    while (queue_size(tmp_queue))
      queue_add(genome_nodes, queue_get(tmp_queue));
    queue_delete(tmp_queue);
  }
}

static void update_pseudo_node_range(GenomeNode *pseudo_node,
                                     GenomeNode *genome_feature)
{
  assert(pseudo_node && genome_feature);
  assert(genome_feature_is_pseudo((GenomeFeature*) pseudo_node));
  assert(!genome_feature_is_pseudo((GenomeFeature*) genome_feature));
  genome_node_set_range(pseudo_node,
                        range_join(genome_node_get_range(pseudo_node),
                                   genome_node_get_range(genome_feature)));
}

static void genome_node_is_part_of_pseudo_node(GenomeNode *pseudo_node,
                                               GenomeNode *child,
                                               FeatureInfo *feature_info)
{
  const char *id;
  assert(pseudo_node && genome_feature_is_pseudo((GenomeFeature*) pseudo_node));
  assert(child && !genome_feature_is_pseudo((GenomeFeature*) child));
  assert(feature_info);
  genome_node_is_part_of_genome_node(pseudo_node, child);
  id = genome_feature_get_attribute(child, ID_STRING);
  assert(id);
  feature_info_add_pseudo_parent(feature_info, id, pseudo_node);
}

static int store_id(const char *id, GenomeNode *genome_feature, bool *is_child,
                    GFF3Parser *parser, Queue *genome_nodes,
                    AutomaticSequenceRegion *auto_sr, const char *filename,
                    unsigned int line_number, Error *err)
{
  GenomeNode *gn;
  int had_err = 0;

  error_check(err);
  assert(id && genome_feature && parser);

  if ((gn = feature_info_get(parser->feature_info, id))) {
    /* this id has been used already -> try to make this a multi-feature */
    if (genome_node_get_line_number(gn) < parser->last_terminator) {
      error_set(err, "the multi-feature with %s \"%s\" on line %u in file "
                "\"%s\" is separated from its counterpart on line %u by "
                "terminator %s on line %u", ID_STRING, id, line_number,
                filename, genome_node_get_line_number(gn), GFF_TERMINATOR,
                parser->last_terminator);
      had_err = -1;
    }
    if (!had_err) {
      GenomeNode *pseudo_parent;
      bool has_parent ;
      has_parent = genome_feature_get_attribute(gn, PARENT_STRING)
                   ? true : false;
      assert(!genome_feature_is_pseudo((GenomeFeature*) gn));
      pseudo_parent = feature_info_get_pseudo_parent(parser->feature_info, id);
      if (pseudo_parent || !genome_feature_is_multi((GenomeFeature*) gn)) {
        if (!pseudo_parent) {
          genome_feature_make_multi_representative((GenomeFeature*) gn);
          if (!has_parent) { /* create pseudo node */
            GenomeNode *pseudo_node = genome_feature_new_pseudo((GenomeFeature*)
                                                                gn);
            genome_node_is_part_of_pseudo_node(pseudo_node, gn,
                                               parser->feature_info);
            replace_node(gn, pseudo_node, genome_nodes, auto_sr);
            genome_node_is_part_of_genome_node(pseudo_node, genome_feature);
            *is_child = true;
          }
        }
        else {
          assert(pseudo_parent);
          update_pseudo_node_range(pseudo_parent, genome_feature);
          genome_node_is_part_of_genome_node(pseudo_parent, genome_feature);
          *is_child = true;
        }
      }
      else {
        assert(has_parent);
        assert(genome_feature_get_multi_representative((GenomeFeature*) gn) ==
               (GenomeFeature*) gn);
      }
      genome_feature_set_multi_representative((GenomeFeature*) genome_feature,
                                              (GenomeFeature*) gn);
    }
  }
  else
    feature_info_add(parser->feature_info, id, genome_feature);

  if (!had_err)
    parser->incomplete_node = true;

  return had_err;
}

static Array* find_roots(StrArray *parents, FeatureInfo *feature_info)
{
  Array *roots;
  unsigned long i;
  assert(parents);
  roots = array_new(sizeof (GenomeNode*));
  for (i = 0; i < strarray_size(parents); i++) {
    GenomeNode *root = feature_info_find_root(feature_info,
                                              strarray_get(parents, i));
    array_add(roots, root);
  }
  return roots;
}

static bool roots_differ(Array *roots)
{
  GenomeNode *first_root;
  unsigned long i;
  assert(roots);
  first_root = *(GenomeNode**) array_get(roots, 0);
  for (i = 1; i < array_size(roots); i++) {
    if (first_root != *(GenomeNode**) array_get(roots, i))
      return true;
  }
  return false;
}

static GenomeNode* merge_pseudo_roots(GenomeNode *pseudo_a,
                                      GenomeNode *pseudo_b,
                                      FeatureInfo *feature_info,
                                      Queue *genome_nodes,
                                      AutomaticSequenceRegion *auto_sr)
{
  GenomeNodeIterator *gni;
  GenomeNode *child;
  assert(pseudo_a && genome_feature_is_pseudo((GenomeFeature*) pseudo_a));
  assert(pseudo_b && genome_feature_is_pseudo((GenomeFeature*) pseudo_b));
  assert(feature_info && genome_nodes);
  /* add children of pseudo node b to pseudo node a */
  gni = genome_node_iterator_new_direct(pseudo_b);
  while ((child = genome_node_iterator_next(gni))) {
    genome_node_is_part_of_genome_node(pseudo_a, child);
    feature_info_replace_pseudo_parent(feature_info, child, pseudo_a);
  }
  genome_node_iterator_delete(gni);
  /* remove pseudo node b from buffer */
  remove_node(pseudo_b, genome_nodes, auto_sr);
  genome_node_delete(pseudo_b);
  return pseudo_a;
}

static GenomeNode* add_node_to_pseudo_node(GenomeNode *pseudo_node,
                                           GenomeNode *normal_node,
                                           FeatureInfo *feature_info,
                                           Queue *genome_nodes,
                                           AutomaticSequenceRegion *auto_sr)
{
  assert(pseudo_node && genome_feature_is_pseudo((GenomeFeature*) pseudo_node));
  assert(normal_node &&
         !genome_feature_is_pseudo((GenomeFeature*) normal_node));
  assert(feature_info && genome_nodes);
  genome_node_is_part_of_pseudo_node(pseudo_node, normal_node, feature_info);
  remove_node(normal_node, genome_nodes, auto_sr);
  return pseudo_node;
}

static GenomeNode* create_pseudo_node(GenomeNode *node_a, GenomeNode *node_b,
                                      FeatureInfo *feature_info,
                                      Queue *genome_nodes,
                                      AutomaticSequenceRegion *auto_sr)
{
  GenomeNode *pseudo_node;
  assert(node_a && !genome_feature_is_pseudo((GenomeFeature*) node_a));
  assert(node_b && !genome_feature_is_pseudo((GenomeFeature*) node_b));
  assert(feature_info && genome_nodes);
  pseudo_node = genome_feature_new_pseudo((GenomeFeature*) node_a);
  genome_node_is_part_of_pseudo_node(pseudo_node, node_a, feature_info);
  genome_node_is_part_of_pseudo_node(pseudo_node, node_b, feature_info);
  replace_node(node_a, pseudo_node, genome_nodes, auto_sr);
  remove_node(node_b, genome_nodes, auto_sr);
  return pseudo_node;
}

static GenomeNode* join_root_pair(GenomeNode *root_a, GenomeNode *root_b,
                                  FeatureInfo *feature_info,
                                  Queue *genome_nodes,
                                  AutomaticSequenceRegion *auto_sr)
{
  bool root_a_is_pseudo, root_b_is_pseudo;
  GenomeNode *master_root;
  assert(root_a && root_b && feature_info && genome_nodes);
  root_a_is_pseudo = genome_feature_is_pseudo((GenomeFeature*) root_a);
  root_b_is_pseudo = genome_feature_is_pseudo((GenomeFeature*) root_b);
  if (root_a_is_pseudo && root_b_is_pseudo) {
    master_root = merge_pseudo_roots(root_a, root_b, feature_info, genome_nodes,
                                     auto_sr);
  }
  else if (root_a_is_pseudo) { /* !root_b_is_pseudo */
    master_root = add_node_to_pseudo_node(root_a, root_b, feature_info,
                                          genome_nodes, auto_sr);
  }
  else if (root_b_is_pseudo) { /* !root_a_is_pseudo */
    master_root = add_node_to_pseudo_node(root_b, root_a, feature_info,
                                          genome_nodes, auto_sr);
  }
  else { /* !root_a_is_pseudo && !root_b_is_pseudo */
    master_root =  create_pseudo_node(root_a, root_b, feature_info,
                                      genome_nodes, auto_sr);
  }
  return master_root;
}

static void join_roots(Array *roots, FeatureInfo *feature_info,
                       Queue *genome_nodes, AutomaticSequenceRegion *auto_sr)
{
  GenomeNode *master_root;
  unsigned long i;
  assert(roots && feature_info && genome_nodes);
  master_root = *(GenomeNode**) array_get(roots, 0);
  for (i = 1; i < array_size(roots); i++) {
    master_root = join_root_pair(master_root,
                                 *(GenomeNode**) array_get(roots, i),
                                 feature_info, genome_nodes, auto_sr);
  }
}

static int process_parent_attr(char *parent_attr, GenomeNode *genome_feature,
                               bool *is_child, GFF3Parser *parser,
                               Queue *genome_nodes,
                               AutomaticSequenceRegion *auto_sr,
                               const char *filename, unsigned int line_number,
                               Error *err)
{
  Splitter *parent_splitter;
  StrArray *valid_parents;
  unsigned long i;
  int had_err = 0;

  error_check(err);
  assert(parent_attr);

  valid_parents = strarray_new();
  parent_splitter = splitter_new();
  splitter_split(parent_splitter, parent_attr, strlen(parent_attr), ',');
  assert(splitter_size(parent_splitter));

  for (i = 0; i < splitter_size(parent_splitter); i++) {
    GenomeNode* parent_gf;
    const char *parent = splitter_get_token(parent_splitter, i);
    parent_gf = feature_info_get(parser->feature_info, parent);
    if (!parent_gf) {
      if (!parser->tidy) {
        error_set(err, "%s \"%s\" on line %u in file \"%s\" has not been "
                  "previously defined (via \"%s=\")", PARENT_STRING, parent,
                  line_number, filename, ID_STRING);
        had_err = -1;
      }
      else {
        warning("%s \"%s\" on line %u in file \"%s\" has not been "
                "previously defined (via \"%s=\")", PARENT_STRING, parent,
                line_number, filename, ID_STRING);
      }
    }
    else if (str_cmp(genome_node_get_seqid(parent_gf),
                     genome_node_get_seqid(genome_feature))) {
      error_set(err, "child on line %u in file \"%s\" has different "
                "sequence id than its parent on line %u ('%s' vs. '%s')",
                genome_node_get_line_number(genome_feature), filename,
                genome_node_get_line_number(parent_gf),
                str_get(genome_node_get_seqid(genome_feature)),
                str_get(genome_node_get_seqid(parent_gf)));
      had_err = -1;
    }
    else {
      assert(parser->incomplete_node);
      genome_node_is_part_of_genome_node(parent_gf, genome_feature);
      *is_child = true;
      strarray_add_cstr(valid_parents, parent);
    }
  }

  if (!had_err && !parser->tidy) {
    assert(splitter_size(parent_splitter) == strarray_size(valid_parents));
  }

  splitter_delete(parent_splitter);

  /* make sure all (valid) parents have the same (pseudo-)root */
  if (!had_err && strarray_size(valid_parents) >= 2) {
    Array *roots = find_roots(valid_parents, parser->feature_info);
    if (roots_differ(roots))
        join_roots(roots, parser->feature_info, genome_nodes, auto_sr);
    array_delete(roots);
  }

  strarray_delete(valid_parents);

  return had_err;
}

static bool is_blank_attribute(const char *attribute)
{
  while (*attribute != '\0') {
    if (*attribute != ' ')
      return false;
    attribute++;
  }
  return true;
}

static int check_multi_feature_constrains(GenomeNode *new_gf,
                                          GenomeNode *old_gf, const char *id,
                                          const char *filename,
                                          unsigned int line_number, Error *err)
{
  int had_err = 0;
  error_check(err);
  assert(new_gf && old_gf);
  assert(!genome_feature_is_pseudo((GenomeFeature*) new_gf));
  assert(!genome_feature_is_pseudo((GenomeFeature*) old_gf));
  /* check seqid */
  if (str_cmp(genome_node_get_seqid(new_gf), genome_node_get_seqid(old_gf))) {
    error_set(err, "the multi-feature with %s \"%s\" on line %u in file \"%s\" "
              "has a different sequence id than its counterpart on line %u",
              ID_STRING, id, line_number, filename,
              genome_node_get_line_number(old_gf));
    had_err = -1;
  }
  /* check source */
  if (!had_err && strcmp(genome_feature_get_source((GenomeFeature*) new_gf),
                         genome_feature_get_source((GenomeFeature*) old_gf))) {
    error_set(err, "the multi-feature with %s \"%s\" on line %u in file \"%s\" "
              "has a different source than its counterpart on line %u",
              ID_STRING, id, line_number, filename,
              genome_node_get_line_number(old_gf));
    had_err = -1;
  }
  /* check type */
  if (!had_err && genome_feature_get_type((GenomeFeature*) new_gf) !=
                  genome_feature_get_type((GenomeFeature*) old_gf)) {
    error_set(err, "the multi-feature with %s \"%s\" on line %u in file \"%s\" "
              "has a different type than its counterpart on line %u",
              ID_STRING, id, line_number, filename,
              genome_node_get_line_number(old_gf));
    had_err = -1;
  }
  /* check strand */
  if (!had_err && genome_feature_get_strand((GenomeFeature*) new_gf) !=
                  genome_feature_get_strand((GenomeFeature*) old_gf)) {
    error_set(err, "the multi-feature with %s \"%s\" on line %u in file \"%s\" "
              "has a different strand than its counterpart on line %u",
              ID_STRING, id, line_number, filename,
              genome_node_get_line_number(old_gf));
    had_err = -1;
  }
  /* check attributes (for target attribute only the name) */
  /* XXX */
  return had_err;
}

static int parse_attributes(char *attributes, GenomeNode *genome_feature,
                            bool *is_child, GFF3Parser *parser,
                            Queue *genome_nodes,
                            AutomaticSequenceRegion *auto_sr,
                            const char *filename, unsigned int line_number,
                            Error *err)
{
  Splitter *attribute_splitter, *tmp_splitter, *parent_splitter;
  unsigned long i;
  int had_err = 0;

  error_check(err);
  assert(attributes);

  attribute_splitter = splitter_new();
  tmp_splitter = splitter_new();
  parent_splitter = splitter_new();
  splitter_split(attribute_splitter, attributes, strlen(attributes), ';');

  for (i = 0; !had_err && i < splitter_size(attribute_splitter); i++) {
    const char *attr_tag = NULL;
    char *attr_value = NULL, *token = splitter_get_token(attribute_splitter, i);
    if (strncmp(token, ".", 1) == 0) {
      if (splitter_size(attribute_splitter) > 1) {
        error_set(err, "more than one attribute token defined on line %u in "
                  "file \"%s\", altough the first one is '.'", line_number,
                  filename);
        had_err = -1;
      }
      else
        break; /* no attributes to parse */
    }
    else if (is_blank_attribute(token))
      continue;
    else {
      splitter_reset(tmp_splitter);
      splitter_split(tmp_splitter, token, strlen(token), '=');
      if (splitter_size(tmp_splitter) != 2) {
        error_set(err, "token \"%s\" on line %u in file \"%s\" does not "
                  "contain exactly one '='", token, line_number, filename);
        had_err = -1;
        break;
      }
      else {
        attr_tag = splitter_get_token(tmp_splitter, 0);
        /* Skip leading blanks of attribute tag.
           Iit is not mentioned in the GFF3 spec that attribute tags cannot
           start with blanks, but if a Parent or ID attribute is prepended by a
           blank (e.g. '; Parent=' instead of ';Parent=') the parent-child
           relations do not get reconstructed correctly, because then '; Parent'
           would be treated as an attribute without special meaning.
           Therefore we decided to skip leading blanks and do _not_ consider
           them as part of the attribute but rather as an artefact of the GFF3
           construction. */
        while (attr_tag[0] == ' ')
          attr_tag++;
        attr_value = splitter_get_token(tmp_splitter, 1);
      }
    }
    if (!had_err && !strlen(attr_tag)) {
      error_set(err, "attribute \"=%s\" on line %u in file \"%s\" has no tag",
                attr_value, line_number, filename);
      had_err = -1;
    }
    if (!had_err && !strlen(attr_value)) {
      error_set(err, "attribute \"%s=\" on line %u in file \"%s\" has no value",
                 attr_tag, line_number, filename);
      had_err = -1;
    }
    /* check for duplicate attributes */
    if (!had_err && genome_feature_get_attribute(genome_feature, attr_tag)) {
      error_set(err, "more then one %s attribute on line %u in file \"%s\"",
                attr_tag, line_number, filename);
      had_err = -1;
    }
    /* save all attributes, although the Parent and ID attribute is newly
       created in GFF3 output */
    if (!had_err) {
      genome_feature_add_attribute((GenomeFeature*) genome_feature, attr_tag,
                                   attr_value);
    }
    /* some attributes require special care */
    if (!had_err) {
      if (!strcmp(attr_tag, ID_STRING)) {
        had_err = store_id(attr_value, genome_feature, is_child, parser,
                           genome_nodes, auto_sr, filename, line_number, err);
      }
      else if (!strcmp(attr_tag, PARENT_STRING)) {
        had_err = process_parent_attr(attr_value, genome_feature, is_child,
                                      parser, genome_nodes, auto_sr, filename,
                                      line_number, err);
      }
      else if (!strcmp(attr_tag, "Target")) {
        /* the value of ``Target'' attributes have a special syntax which is
           checked here */
        had_err = gff3parser_parse_target_attributes(attr_value, NULL, NULL,
                                                     NULL, NULL, filename,
                                                     line_number, err);
      }
    }
  }

  if (!had_err && genome_feature_is_multi((GenomeFeature*) genome_feature)) {
    had_err = check_multi_feature_constrains(genome_feature,
                        (GenomeNode*)
                        genome_feature_get_multi_representative((GenomeFeature*)
                                                                genome_feature),
                        genome_feature_get_attribute(genome_feature, ID_STRING),
                        filename, line_number, err);
  }

  splitter_delete(parent_splitter);
  splitter_delete(tmp_splitter);
  splitter_delete(attribute_splitter);

  return had_err;
}

static void set_source(GenomeNode *genome_feature, const char *source,
                       Hashmap *source_to_str_mapping)
{
  Str *source_str;
  assert(genome_feature && source && source_to_str_mapping);
  source_str = hashmap_get(source_to_str_mapping, source);
  if (!source_str) {
    source_str = str_new_cstr(source);
    hashmap_add(source_to_str_mapping, str_get(source_str), source_str);
  }
  assert(source_str);
  genome_feature_set_source(genome_feature, source_str);
}

static int parse_regular_gff3_line(GFF3Parser *parser, Queue *genome_nodes,
                                   char *line, size_t line_length,
                                   Str *filenamestr, unsigned int line_number,
                                   Error *err)
{
  GenomeNode *gn = NULL, *genome_feature = NULL;
  GenomeFeatureType *gft = NULL;
  Splitter *splitter;
  AutomaticSequenceRegion *auto_sr = NULL;
  Str *changed_seqid = NULL;
  Strand strand_value;
  float score_value;
  Phase phase_value;
  Range range;
  char *seqid = NULL, *source = NULL, *type = NULL, *start = NULL,
       *end = NULL, *score = NULL, *strand = NULL, *phase = NULL,
       *attributes = NULL, **tokens;
  const char *filename;
  bool score_is_defined, is_child = false;
  int had_err = 0;

  error_check(err);

  filename = str_get(filenamestr);

  /* create splitter */
  splitter = splitter_new();

  /* parse */
  splitter_split(splitter, line, line_length, '\t');
  if (splitter_size(splitter) != 9UL) {
    error_set(err, "line %u in file \"%s\" does not contain 9 tab (\\t) "
                   "separated fields", line_number, filename);
    had_err = -1;
  }
  if (!had_err) {
    tokens = splitter_get_tokens(splitter);
    seqid      = tokens[0];
    source     = tokens[1];
    type       = tokens[2];
    start      = tokens[3];
    end        = tokens[4];
    score      = tokens[5];
    strand     = tokens[6];
    phase      = tokens[7];
    attributes = tokens[8];
  }

  /* parse the feature type */
  if (!had_err &&
      !(gft = feature_type_factory_create_gft(parser->feature_type_factory,
                                              type))) {
    error_set(err, "type \"%s\" on line %u in file \"%s\" is not a valid one",
              type, line_number, filename);
    had_err = -1;
  }

  /* parse the range */
  if (!had_err)
    had_err = parse_range(&range, start, end, line_number, filename, err);
  if (!had_err && range.start == 0) {
      error_set(err, "illegal feature start 0 on line %u in file \"%s\" "
                   "(GFF3 files are 1-based)", line_number, filename);
      had_err = -1;
  }
  if (!had_err)
    had_err = add_offset_if_necessary(&range, parser, seqid, err);

  /* parse the score */
  if (!had_err) {
    had_err = parse_score(&score_is_defined, &score_value, score, line_number,
                          filename, err);
  }

  /* parse the strand */
  if (!had_err)
    had_err = parse_strand(&strand_value, strand, line_number, filename, err);

  /* parse the phase */
  if (!had_err)
    had_err = parse_phase(&phase_value, phase, line_number, filename, err);

  /* create the feature */
  if (!had_err) {
    genome_feature = genome_feature_new(gft, range, strand_value);
    genome_node_set_origin(genome_feature, filenamestr, line_number);
  }

  /* set seqid */
  if (!had_err) {
    had_err = set_seqid(genome_feature, seqid, range, &auto_sr, parser,
                        filename, line_number, err);
  }

  /* set source */
  if (!had_err)
    set_source(genome_feature, source, parser->source_to_str_mapping);

  /* parse the attributes */
  if (!had_err) {
    had_err = parse_attributes(attributes, genome_feature, &is_child, parser,
                               genome_nodes, auto_sr, filename, line_number,
                               err);
  }

  if (!had_err && score_is_defined)
    genome_feature_set_score((GenomeFeature*) genome_feature, score_value);
  if (!had_err && phase_value != PHASE_UNDEFINED)
    genome_feature_set_phase(genome_feature, phase_value);

  if (!had_err) {
    gn = (is_child || auto_sr) ? NULL : genome_feature;
    if (auto_sr && !is_child)
      array_add(auto_sr->genome_features, genome_feature);
  }
  else if (!is_child)
    genome_node_delete(genome_feature);

  if (!had_err && gn)
    queue_add(genome_nodes, gn);

  /* free */
  str_delete(changed_seqid);
  splitter_delete(splitter);

  return had_err;
}

static int parse_first_gff3_line(const char *line, const char *filename,
                                 Error *err)
{
  int version, had_err = 0;
  error_check(err);
  assert(line && filename);
  if (strncmp(line, GFF_VERSION_PREFIX, strlen(GFF_VERSION_PREFIX))) {
    error_set(err, "line 1 in file \"%s\" does not begin with \"%s\"", filename,
              GFF_VERSION_PREFIX);
    had_err = -1;
  }
  if (!had_err) {
    line += strlen(GFF_VERSION_PREFIX);
    /* skip blanks */
    while (line[0] == ' ')
      line++;
    had_err = parse_int_line(&version, line, 1, filename, err);
  }
  if (!had_err && version != GFF_VERSION) {
    error_set(err, "GFF version %d does not equal required version %u ",
              version, GFF_VERSION);
    had_err = -1;
  }
  return had_err;
}

static int parse_fasta_entry(Queue *genome_nodes, const char *line,
                             Str *filename, unsigned int line_number,
                             GenFile *fpin, Error *err)
{
  int had_err = 0;
  error_check(err);
  assert(line && line_number && fpin);
  if (line[0] != '>') {
    error_set(err, "line %d does not start with '>' as expected", line_number);
    had_err = -1;
  }
  if (!had_err) {
    GenomeNode *sequence_node;
    Str *sequence = str_new();
    int cc;
    while ((cc = genfile_xfgetc(fpin)) != EOF) {
      if (cc == '>') {
        genfile_unget_char(fpin, cc);
        break;
      }
      if (cc != '\n' && cc != '\r' && cc != ' ')
        str_append_char(sequence, cc);
    }
    sequence_node = sequence_node_new(line+1, sequence);
    genome_node_set_origin(sequence_node, filename, line_number);
    queue_add(genome_nodes, sequence_node);
  }
  return had_err;
}

static int add_auto_sr_to_queue(UNUSED void *key, void *value, void *data,
                                UNUSED Error *err)
{
  AutomaticSequenceRegion *auto_sr = value;
  Queue *genome_nodes = data;
  GenomeNode *gf;
  unsigned int i;
  error_check(err);
  assert(key && value && data);
  if (array_size(auto_sr->genome_features)) {
    queue_add(genome_nodes, auto_sr->sequence_region);
    auto_sr->sequence_region = NULL;
    for (i = 0; i < array_size(auto_sr->genome_features); i++) {
      gf = *(GenomeNode**) array_get(auto_sr->genome_features, i);
      queue_add(genome_nodes, gf);
    }
    array_reset(auto_sr->genome_features);
  }
  return 0;
}

static int parse_meta_gff3_line(GFF3Parser *parser, Queue *genome_nodes,
                                char *line, size_t line_length,
                                Str *filenamestr, unsigned int line_number,
                                Error *err)
{
  char *tmpline, *tmplineend, *seqstart, *seqid = NULL;
  GenomeNode *gn;
  Str *changed_seqid = NULL;
  SimpleSequenceRegion *ssr = NULL;
  Range range;
  const char *filename;
  int had_err = 0;

  error_check(err);
  assert(line[0] == '#');

  filename = str_get(filenamestr);

  if (line_length == 1 || line[1] != '#') {
    /* storing comment */
    gn = comment_new(line+1);
    genome_node_set_origin(gn, filenamestr, line_number);
    queue_add(genome_nodes, gn);
  }
  else if (strcmp(line, GFF_FASTA_DIRECTIVE) == 0) {
    if (!parser->fasta_parsing) {
      parser->fasta_parsing = true;
      had_err = hashmap_foreach(parser->undefined_sequence_regions,
                                add_auto_sr_to_queue, genome_nodes, NULL);
      assert(!had_err); /* add_auto_sr_to_queue() is sane */
    }
  }
  else if ((strncmp(line, GFF_SEQUENCE_REGION,
                    strlen(GFF_SEQUENCE_REGION)) == 0)) {
    /* we are in a line starting with "##sequence-region" */
    tmpline = line + strlen(GFF_SEQUENCE_REGION);
    tmplineend = line + line_length - 1;

    /* skip blanks */
    while (tmpline[0] == ' ')
      tmpline++;
    if (tmpline > tmplineend) {
      error_set(err, "missing sequence region name on line %u in file \"%s\"",
                line_number, filename);
      had_err = -1;
    }
    if (!had_err) {
      seqid = tmpline; /* save seqid */
      /* skip non-blanks */
      while (tmpline < tmplineend && !(tmpline[0] == ' '))
        tmpline++;
      /* terminate seqid */
      *tmpline++ = '\0';
      /* skip blanks */
      while (tmpline < tmplineend && tmpline[0] == ' ')
        tmpline++;
      if (tmpline > tmplineend) {
        error_set(err, "missing sequence region start on line %u in file "
                  "\"%s\"", line_number, filename);
        had_err = -1;
      }
      else
        seqstart = tmpline;
    }

    if (!had_err) {
      /* skip non-blanks */
      while (tmpline <= tmplineend && !(tmpline[0] == ' '))
        tmpline++;
      /* terminate seqstart */
      *tmpline++ = '\0';
      /* skip blanks */
      while (tmpline < tmplineend && tmpline[0] == ' ')
        tmpline++;
      if (tmpline > tmplineend) {
        error_set(err, "missing sequence region end on line %u in file \"%s\"",
                  line_number, filename);
        had_err = -1;
      }
    }
    if (!had_err) {
      had_err = parse_range(&range, seqstart, tmpline, line_number, filename,
                            err);
    }
    if (!had_err  && range.start == 0) {
      error_set(err, "illegal region start 0 on line %u in file \"%s\" "
                "(GFF3 files are 1-based)", line_number, filename);
      had_err = -1;
    }
    if (!had_err)
      had_err = add_offset_if_necessary(&range, parser, seqid, err);
    if (!had_err) {
      if (hashmap_get(parser->undefined_sequence_regions, seqid)) {
        error_set(err, "genome feature with id \"%s\" has been defined before "
                  "the corresponding \"%s\" definition on line %u in file "
                  "\"%s\"", seqid, GFF_SEQUENCE_REGION, line_number, filename);
        had_err = -1;
      }
    }
    if (!had_err) {
      /* now we can create a sequence region node */
      assert(seqid);
      ssr = hashmap_get(parser->seqid_to_ssr_mapping, seqid);
      if (ssr) {
        error_set(err, "the sequence region \"%s\" on line %u in file \"%s\" "
                  "has already been defined", str_get(ssr->seqid_str),
                  line_number, filename);
        had_err = -1;
      }
      else {
        ssr = simple_sequence_region_new(seqid, range, line_number);
        hashmap_add(parser->seqid_to_ssr_mapping, str_get(ssr->seqid_str), ssr);
      }
    }
    if (!had_err) {
      assert(ssr);
      gn = sequence_region_new(ssr->seqid_str, range);
      genome_node_set_origin(gn, filenamestr, line_number);
      queue_add(genome_nodes, gn);
    }
  }
  else if (strcmp(line, GFF_TERMINATOR) == 0) { /* terminator */
    /* now all nodes are complete */
    parser->incomplete_node = false;
    if (!parser->checkids)
      feature_info_reset(parser->feature_info);
    parser->last_terminator = line_number;
  }
  else {
    warning("skipping unknown meta line %u in file \"%s\": %s", line_number,
            filename, line);
  }
  str_delete(changed_seqid);
  return had_err;
}

int gff3parser_parse_genome_nodes(int *status_code, GFF3Parser *parser,
                                  Queue *genome_nodes, Str *filenamestr,
                                  unsigned long long *line_number,
                                  GenFile *fpin, Error *err)
{
  size_t line_length;
  Str *line_buffer;
  char *line;
  const char *filename;
  int rval, had_err = 0;

  error_check(err);

  filename = str_get(filenamestr);

  /* init */
  line_buffer = str_new();

  while ((rval = str_read_next_line_generic(line_buffer, fpin)) != EOF) {
    line = str_get(line_buffer);
    line_length = str_length(line_buffer);
    (*line_number)++;

    if (*line_number == 1) {
      if ((had_err = parse_first_gff3_line(line, filename, err)))
        break;
    }
    else if (line_length == 0) {
      warning("skipping blank line %llu in file \"%s\"", *line_number,
              filename);
    }
    else if (parser->fasta_parsing || line[0] == '>') {
      if (!parser->fasta_parsing) {
        parser->fasta_parsing = true;
        had_err = hashmap_foreach(parser->undefined_sequence_regions,
                                  add_auto_sr_to_queue, genome_nodes, NULL);
        assert(!had_err); /* add_auto_sr_to_queue() is sane */
      }
      had_err = parse_fasta_entry(genome_nodes, line, filenamestr, *line_number,
                                  fpin, err);
      break;
    }
    else if (line[0] == '#') {
      had_err = parse_meta_gff3_line(parser, genome_nodes, line, line_length,
                                     filenamestr, *line_number, err);
      if (had_err ||
          (!parser->incomplete_node && queue_size(genome_nodes))) {
        break;
      }
    }
    else {
      had_err = parse_regular_gff3_line(parser, genome_nodes, line, line_length,
                                        filenamestr, *line_number, err);
      if (had_err || (!parser->incomplete_node && queue_size(genome_nodes)))
        break;
    }
    str_reset(line_buffer);
  }

  if (had_err) {
    while (queue_size(genome_nodes))
      genome_node_rec_delete(queue_get(genome_nodes));
  }
  else if (rval == EOF) {
    /* the file has been parsed completely, add automatically created sequence
       regions to queue */
    had_err = hashmap_foreach(parser->undefined_sequence_regions,
                              add_auto_sr_to_queue, genome_nodes, NULL);
    assert(!had_err); /* add_auto_sr_to_queue() is sane */
  }

  str_delete(line_buffer);
  if (queue_size(genome_nodes))
    *status_code = 0; /* at least one node was created */
  else
    *status_code = EOF;
  return had_err;
}

void gff3parser_reset(GFF3Parser *parser)
{
  assert(parser);
  parser->fasta_parsing = false;
  feature_info_reset(parser->feature_info);
  hashmap_reset(parser->seqid_to_ssr_mapping);
  hashmap_reset(parser->source_to_str_mapping);
  hashmap_reset(parser->undefined_sequence_regions);
  parser->last_terminator = 0;
}

void gff3parser_delete(GFF3Parser *parser)
{
  if (!parser) return;
  feature_info_delete(parser->feature_info);
  hashmap_delete(parser->seqid_to_ssr_mapping);
  hashmap_delete(parser->source_to_str_mapping);
  hashmap_delete(parser->undefined_sequence_regions);
  mapping_delete(parser->offset_mapping);
  ma_free(parser);
}
