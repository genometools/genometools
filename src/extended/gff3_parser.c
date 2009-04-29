/*
  Copyright (c) 2006-2009 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "core/assert_api.h"
#include "core/cstr.h"
#include "core/hashmap.h"
#include "core/ma.h"
#include "core/parseutils.h"
#include "core/splitter.h"
#include "core/undef.h"
#include "core/unused_api.h"
#include "core/warning_api.h"
#include "extended/comment_node_api.h"
#include "extended/feature_info.h"
#include "extended/feature_node.h"
#include "extended/feature_node_iterator_api.h"
#include "extended/genome_node.h"
#include "extended/gff3_escaping.h"
#include "extended/gff3_parser.h"
#include "extended/mapping.h"
#include "extended/region_node.h"

struct GtGFF3Parser {
  GtFeatureInfo *feature_info;
  GtHashmap *seqid_to_ssr_mapping, /* maps seqids to simple sequence regions */
            *source_to_str_mapping,
            *undefined_sequence_regions; /* contains all (automatically created)
                                            sequence regions */
  bool incomplete_node, /* at least on node is potentially incomplete */
       checkids,
       tidy,
       fasta_parsing; /* parser is in FASTA parsing mode */
  long offset;
  GtMapping *offset_mapping;
  GtTypeChecker *type_checker;
  unsigned int last_terminator; /* line number of the last terminator */
};

typedef struct {
  GtGenomeNode *sequence_region; /* the automatically created sequence region */
  GtArray *feature_nodes; /* the features nodes which belong to this region */
} AutomaticSequenceRegion;

static AutomaticSequenceRegion* automatic_sequence_region_new(void)
{
  AutomaticSequenceRegion *auto_sr;
  auto_sr = gt_malloc(sizeof (AutomaticSequenceRegion));
  auto_sr->feature_nodes = gt_array_new(sizeof (GtFeatureNode*));
  return auto_sr;
}

static void automatic_sequence_region_delete(AutomaticSequenceRegion *auto_sr)
{
  unsigned long i;
  if (!auto_sr) return;
  gt_genome_node_delete(auto_sr->sequence_region);
  for (i = 0; i < gt_array_size(auto_sr->feature_nodes); i++) {
    gt_genome_node_delete(*(GtGenomeNode**)
                              gt_array_get(auto_sr->feature_nodes, i));
  }
  gt_array_delete(auto_sr->feature_nodes);
  gt_free(auto_sr);
}

typedef struct {
  GtStr *seqid_str;
  GtRange range;
  unsigned int line_number;
} SimpleSequenceRegion;

static SimpleSequenceRegion* simple_sequence_region_new(const char *seqid,
                                                        GtRange range,
                                                        unsigned int
                                                        line_number)
{
  SimpleSequenceRegion *ssr = gt_malloc(sizeof *ssr);
  ssr->seqid_str = gt_str_new_cstr(seqid);
  ssr->range = range;
  ssr->line_number = line_number;
  return ssr;
}

static void simple_sequence_region_delete(SimpleSequenceRegion *ssr)
{
  if (!ssr) return;
  gt_str_delete(ssr->seqid_str);
  gt_free(ssr);
}

GtGFF3Parser* gt_gff3_parser_new(GtTypeChecker *type_checker)
{
  GtGFF3Parser *parser;
  parser = gt_malloc(sizeof *parser);
  parser->feature_info = gt_feature_info_new();
  parser->seqid_to_ssr_mapping = gt_hashmap_new(
    HASH_STRING, NULL, (GtFree) simple_sequence_region_delete);
  parser->source_to_str_mapping = gt_hashmap_new(HASH_STRING, NULL,
                                              (GtFree) gt_str_delete);
  parser->undefined_sequence_regions = gt_hashmap_new(HASH_STRING, NULL,
                                     (GtFree) automatic_sequence_region_delete);
  parser->incomplete_node = false;
  parser->checkids = false;
  parser->tidy = false;
  parser->fasta_parsing = false;
  parser->offset = UNDEF_LONG;
  parser->offset_mapping = NULL;
  parser->type_checker = type_checker ? gt_type_checker_ref(type_checker)
                                      : NULL;
  parser->last_terminator = 0;
  return parser;
}

void gt_gff3_parser_check_id_attributes(GtGFF3Parser *parser)
{
  gt_assert(parser);
  parser->checkids = true;
}

void gt_gff3_parser_set_offset(GtGFF3Parser *parser, long offset)
{
  gt_assert(parser);
  gt_assert(!parser->offset_mapping);
  parser->offset = offset;
}

int gt_gff3_parser_set_offsetfile(GtGFF3Parser *parser, GtStr *offsetfile,
                                  GtError *err)
{
  gt_error_check(err);
  gt_assert(parser);
  gt_assert(parser->offset == UNDEF_LONG);
  parser->offset_mapping = gt_mapping_new(offsetfile, "offsets",
                                       MAPPINGTYPE_INTEGER, err);
  if (parser->offset_mapping)
    return 0;
  return -1;
}

void gt_gff3_parser_enable_tidy_mode(GtGFF3Parser *parser)
{
  gt_assert(parser);
  parser->tidy = true;
}

static int add_offset_if_necessary(GtRange *range, GtGFF3Parser *parser,
                                   const char *seqid, GtError *err)
{
  long offset;
  int had_err = 0;
  gt_error_check(err);
  if (parser->offset != UNDEF_LONG)
    *range = gt_range_offset(range, parser->offset);
  else if (parser->offset_mapping) {
    had_err = gt_mapping_map_integer(parser->offset_mapping, &offset, seqid,
                                     err);
    if (!had_err)
      *range = gt_range_offset(range, offset);
  }
  return had_err;
}

static int parse_target_attribute(const char *value, GtStr *target_id,
                                  GtRange *target_range,
                                  GtStrand *target_strand,
                                  const char *filename,
                                  unsigned int line_number, GtError *err)
{
  unsigned long num_of_tokens;
  GtStr *unescaped_target;
  char *escaped_target;
  GtStrand parsed_strand;
  GtSplitter *splitter;
  GtRange parsed_range;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(value && filename);
  splitter = gt_splitter_new();
  unescaped_target = gt_str_new();
  escaped_target = gt_cstr_dup(value);
  gt_splitter_split(splitter, escaped_target, strlen(escaped_target), ' ');
  num_of_tokens = gt_splitter_size(splitter);
  if (!(num_of_tokens == 3 || num_of_tokens == 4)) {
    gt_error_set(err, "Target attribute value '%s' on line %u in file \"%s\" "
              "must have 3 or 4 blank separated entries", value, line_number,
              filename);
    had_err = -1;
  }
  /* parse target id */
  if (!had_err) {
    had_err = gt_gff3_unescape(unescaped_target,
                               gt_splitter_get_token(splitter, 0),
                               strlen(gt_splitter_get_token(splitter, 0)), err);
  }
  if (!had_err && target_id) gt_str_append_str(target_id, unescaped_target);
  /* parse target range */
  if (!had_err) {
    had_err = gt_parse_range(&parsed_range, gt_splitter_get_token(splitter, 1),
                          gt_splitter_get_token(splitter, 2), line_number,
                          filename, err);
  }
  if (!had_err && target_range) *target_range = parsed_range;
  /* parse target strand (if given) */
  if (!had_err) {
    if (gt_splitter_size(splitter) == 4) {
      had_err = gt_parse_strand(&parsed_strand, gt_splitter_get_token(splitter,
                                                                      3),
                             line_number, filename, err);
      if (!had_err && target_strand)
        *target_strand = parsed_strand;
    }
    else if (target_strand)
      *target_strand = GT_NUM_OF_STRAND_TYPES; /* undefined */
  }
  gt_free(escaped_target);
  gt_str_delete(unescaped_target);
  gt_splitter_delete(splitter);
  return had_err;
}

int gt_gff3_parser_parse_target_attributes(const char *values,
                                           unsigned long *num_of_targets,
                                           GtStr *first_target_id,
                                           GtRange *first_target_range,
                                           GtStrand *first_target_strand,
                                           const char *filename,
                                           unsigned int line_number,
                                           GtError *err)
{
  GtSplitter *splitter;
  unsigned long i;
  char *targets;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(values && filename);
  targets = gt_cstr_dup(values);
  splitter = gt_splitter_new();
  gt_splitter_split(splitter, targets, strlen(targets), ',');
  if (num_of_targets)
    *num_of_targets = gt_splitter_size(splitter);
  for (i = 0; !had_err && i < gt_splitter_size(splitter); i++) {
    had_err = parse_target_attribute(gt_splitter_get_token(splitter, i),
                                     i ? NULL : first_target_id,
                                     i ? NULL : first_target_range,
                                     i ? NULL : first_target_strand, filename,
                                     line_number, err);
  }
  gt_free(targets);
  gt_splitter_delete(splitter);
  return had_err;
}

static int get_seqid_str(GtStr **seqid_str, const char *seqid, GtRange range,
                         AutomaticSequenceRegion **auto_sr,
                         GtGFF3Parser *parser, const char *filename,
                         unsigned int line_number, GtError *err)
{
  SimpleSequenceRegion *ssr;
  int had_err = 0;

  gt_error_check(err);

  ssr = gt_hashmap_get(parser->seqid_to_ssr_mapping, seqid);
  if (!ssr) {
    /* sequence region has not been previously introduced -> check if one has
       already been created automatically */
    *auto_sr = gt_hashmap_get(parser->undefined_sequence_regions, seqid);
    if (!*auto_sr) {
      /* sequence region has not been createad automatically -> do it now */
      gt_warning("seqid \"%s\" on line %u in file \"%s\" has not been "
                 "previously introduced with a \"%s\" line, create such a line "
                 "automatically", seqid, line_number, filename,
                 GFF_SEQUENCE_REGION);
      *auto_sr = automatic_sequence_region_new();
      *seqid_str = gt_str_new_cstr(seqid);
      (*auto_sr)->sequence_region = gt_region_node_new(*seqid_str, range.start,
                                                                   range.end);
      gt_hashmap_add(parser->undefined_sequence_regions, gt_str_get(*seqid_str),
                     *auto_sr);
    }
    else {
      GtRange joined_range,
              sr_range = gt_genome_node_get_range((*auto_sr)->sequence_region);
      /* get seqid string */
      *seqid_str =
        gt_str_ref(gt_genome_node_get_seqid((*auto_sr)->sequence_region));
      /* update the range of the sequence region */
      joined_range = gt_range_join(&range, &sr_range);
      gt_genome_node_set_range((*auto_sr)->sequence_region, &joined_range);
    }
  }
  else {
    /* perform range check */
    if (!gt_range_contains(&ssr->range, &range)) {
      gt_error_set(err, "range (%lu,%lu) of feature on line %u in file \"%s\" "
                "is not contained in range (%lu,%lu) of corresponding "
                "sequence region on line %u", range.start, range.end,
                line_number, filename, ssr->range.start, ssr->range.end,
                ssr->line_number);
      had_err = -1;
    }
    else
      *seqid_str = gt_str_ref(ssr->seqid_str);
  }

  return had_err;
}

typedef struct {
  GtFeatureNode *node_to_replace,
                *replacing_node;
} ReplaceInfo;

static int replace_func(void **elem, void *info, GT_UNUSED GtError *err)
{
  ReplaceInfo *replace_info = info;
  GtGenomeNode **node = (GtGenomeNode**) elem;
  gt_error_check(err);
  gt_assert(node && replace_info);
  if (*node == (GtGenomeNode*) replace_info->node_to_replace) {
    *node = (GtGenomeNode*) replace_info->replacing_node;
    return 1;
  }
  return 0;
}

static void replace_node(GtFeatureNode *node_to_replace,
                         GtFeatureNode *replacing_node, GtQueue *genome_nodes,
                         AutomaticSequenceRegion *auto_sr)
{
  ReplaceInfo replace_info;
  int rval;
  gt_assert(node_to_replace && replacing_node && genome_nodes);
  replace_info.node_to_replace = node_to_replace;
  replace_info.replacing_node = replacing_node;
  /* we go backwards in both cases, because we expect that the <node_to_replace>
     is near the end of the queue */
  if (auto_sr) {
    rval = gt_array_iterate_reverse(auto_sr->feature_nodes,
                                    (GtArrayProcessor) replace_func,
                                    &replace_info, NULL);
    gt_assert(rval == 1);
  }
  else {
    rval = gt_queue_iterate_reverse(genome_nodes, replace_func, &replace_info,
                                    NULL);
    gt_assert(rval == 1);
  }
}

static void remove_node(GtGenomeNode *genome_node, GtQueue *genome_nodes,
                        AutomaticSequenceRegion *auto_sr)
{
  gt_assert(genome_node && genome_nodes);
  if (auto_sr) {
    long i;
    /* we go backwards, because we expect that the <genome_node> is near the end
       of the array */
    for (i = gt_array_size(auto_sr->feature_nodes) - 1; i >= 0; i--) {
      if (genome_node ==
          *(GtGenomeNode**) gt_array_get(auto_sr->feature_nodes, i)) {
        break;
      }
    }
    gt_assert(i < gt_array_size(auto_sr->feature_nodes));
    gt_array_rem(auto_sr->feature_nodes, i);
  }
  else
    gt_queue_remove(genome_nodes, genome_node); /* traverses in reverse order */
}

static void update_pseudo_node_range(GtFeatureNode *pseudo_node,
                                     GtFeatureNode *feature_node)
{
  GtRange pseudo_range, feature_range, joined_range;
  gt_assert(pseudo_node && feature_node);
  gt_assert(gt_feature_node_is_pseudo(pseudo_node));
  gt_assert(!gt_feature_node_is_pseudo(feature_node));
  pseudo_range = gt_genome_node_get_range((GtGenomeNode*) pseudo_node);
  feature_range = gt_genome_node_get_range((GtGenomeNode*) feature_node);
  joined_range = gt_range_join(&pseudo_range, &feature_range);
  gt_genome_node_set_range((GtGenomeNode*) pseudo_node, &joined_range);
}

static void feature_node_is_part_of_pseudo_node(GtFeatureNode *pseudo_node,
                                                GtFeatureNode *child,
                                                GtFeatureInfo *feature_info)
{
  const char *id;
  gt_assert(pseudo_node &&
         gt_feature_node_is_pseudo((GtFeatureNode*) pseudo_node));
  gt_assert(child && !gt_feature_node_is_pseudo((GtFeatureNode*) child));
  gt_assert(feature_info);
  gt_feature_node_add_child(pseudo_node, child);
  id = gt_feature_node_get_attribute(child, ID_STRING);
  gt_assert(id);
  gt_feature_info_add_pseudo_parent(feature_info, id, pseudo_node);
}

static int store_id(const char *id, GtFeatureNode *feature_node,
                    bool *is_child, GtGFF3Parser *parser,
                    GtQueue *genome_nodes, AutomaticSequenceRegion *auto_sr,
                    const char *filename, unsigned int line_number,
                    GtError *err)
{
  GtFeatureNode *fn;
  int had_err = 0;

  gt_error_check(err);
  gt_assert(id && feature_node && parser);

  if ((fn = gt_feature_info_get(parser->feature_info, id))) {
    /* this id has been used already -> try to make this a multi-feature */
    if (gt_genome_node_get_line_number((GtGenomeNode*) fn) <
        parser->last_terminator) {
      gt_error_set(err, "the multi-feature with %s \"%s\" on line %u in file "
                   "\"%s\" is separated from its counterpart on line %u by "
                   "terminator %s on line %u", ID_STRING, id, line_number,
                   filename, gt_genome_node_get_line_number((GtGenomeNode*) fn),
                   GFF_TERMINATOR, parser->last_terminator);
      had_err = -1;
    }
    /* check seqid */
    if (!had_err &&
        gt_str_cmp(gt_genome_node_get_seqid((GtGenomeNode*) feature_node),
                   gt_genome_node_get_seqid((GtGenomeNode*) fn))) {
      gt_error_set(err, "the multi-feature with %s \"%s\" on line %u in file "
                   "\"%s\" has a different sequence id than its counterpart on "
                   "line %u", ID_STRING, id, line_number, filename,
                   gt_genome_node_get_line_number((GtGenomeNode*) fn));
      had_err = -1;
    }
    if (!had_err) {
      GtFeatureNode *pseudo_parent;
      bool has_parent ;
      has_parent = gt_feature_node_get_attribute(fn, PARENT_STRING)
                   ? true : false;
      gt_assert(!gt_feature_node_is_pseudo(fn));
      pseudo_parent = gt_feature_info_get_pseudo_parent(parser->feature_info,
                                                        id);
      if (pseudo_parent || !gt_feature_node_is_multi(fn)) {
        if (!pseudo_parent) {
          gt_feature_node_make_multi_representative(fn);
          if (!has_parent) { /* create pseudo node */
            GtFeatureNode *pseudo_node = (GtFeatureNode*)
                                         gt_feature_node_new_pseudo(fn);
            feature_node_is_part_of_pseudo_node(pseudo_node, fn,
                                                parser->feature_info);
            replace_node(fn, pseudo_node, genome_nodes, auto_sr);
            gt_feature_node_add_child(pseudo_node, feature_node);
            *is_child = true;
          }
        }
        else {
          gt_assert(pseudo_parent);
          update_pseudo_node_range(pseudo_parent, feature_node);
          gt_feature_node_add_child(pseudo_parent, feature_node);
          *is_child = true;
        }
      }
      else {
        gt_assert(has_parent);
        gt_assert(gt_feature_node_get_multi_representative(fn) == fn);
      }
      gt_feature_node_set_multi_representative(feature_node, fn);
    }
  }
  else
    gt_feature_info_add(parser->feature_info, id, feature_node);

  if (!had_err)
    parser->incomplete_node = true;

  return had_err;
}

static GtArray* find_roots(GtStrArray *parents, GtFeatureInfo *feature_info)
{
  GtArray *roots;
  unsigned long i;
  gt_assert(parents);
  roots = gt_array_new(sizeof (GtGenomeNode*));
  for (i = 0; i < gt_str_array_size(parents); i++) {
    GtFeatureNode *root = gt_feature_info_find_root(feature_info,
                                                 gt_str_array_get(parents, i));
    gt_array_add(roots, root);
  }
  return roots;
}

static bool roots_differ(GtArray *roots)
{
  GtGenomeNode *first_root;
  unsigned long i;
  gt_assert(roots);
  first_root = *(GtGenomeNode**) gt_array_get(roots, 0);
  for (i = 1; i < gt_array_size(roots); i++) {
    if (first_root != *(GtGenomeNode**) gt_array_get(roots, i))
      return true;
  }
  return false;
}

static GtFeatureNode* merge_pseudo_roots(GtFeatureNode *pseudo_a,
                                         GtFeatureNode *pseudo_b,
                                         GtFeatureInfo *feature_info,
                                         GtQueue *genome_nodes,
                                         AutomaticSequenceRegion *auto_sr)
{
  GtFeatureNodeIterator *fni;
  GtFeatureNode *child;
  gt_assert(pseudo_a && gt_feature_node_is_pseudo(pseudo_a));
  gt_assert(pseudo_b && gt_feature_node_is_pseudo(pseudo_b));
  gt_assert(feature_info && genome_nodes);
  /* add children of pseudo node b to pseudo node a */
  fni = gt_feature_node_iterator_new_direct(pseudo_b);
  while ((child = gt_feature_node_iterator_next(fni))) {
    gt_feature_node_add_child(pseudo_a, child);
    gt_feature_info_replace_pseudo_parent(feature_info, child, pseudo_a);
  }
  gt_feature_node_iterator_delete(fni);
  /* remove pseudo node b from buffer */
  remove_node((GtGenomeNode*) pseudo_b, genome_nodes, auto_sr);
  gt_genome_node_delete((GtGenomeNode*) pseudo_b);
  return pseudo_a;
}

static GtFeatureNode* add_node_to_pseudo_node(GtFeatureNode *pseudo_node,
                                              GtFeatureNode *normal_node,
                                              GtFeatureInfo *feature_info,
                                              GtQueue *genome_nodes,
                                              AutomaticSequenceRegion
                                                 *auto_sr)
{
  gt_assert(pseudo_node &&
         gt_feature_node_is_pseudo((GtFeatureNode*) pseudo_node));
  gt_assert(normal_node &&
         !gt_feature_node_is_pseudo((GtFeatureNode*) normal_node));
  gt_assert(feature_info && genome_nodes);
  feature_node_is_part_of_pseudo_node(pseudo_node, normal_node, feature_info);
  remove_node((GtGenomeNode*) normal_node, genome_nodes, auto_sr);
  return pseudo_node;
}

static GtFeatureNode* create_pseudo_node(GtFeatureNode *node_a,
                                         GtFeatureNode *node_b,
                                         GtFeatureInfo *feature_info,
                                         GtQueue *genome_nodes,
                                         AutomaticSequenceRegion *auto_sr)
{
  GtFeatureNode *pseudo_node;
  gt_assert(node_a && !gt_feature_node_is_pseudo((GtFeatureNode*) node_a));
  gt_assert(node_b && !gt_feature_node_is_pseudo((GtFeatureNode*) node_b));
  gt_assert(feature_info && genome_nodes);
  pseudo_node = (GtFeatureNode*)
                gt_feature_node_new_pseudo((GtFeatureNode*) node_a);
  feature_node_is_part_of_pseudo_node(pseudo_node, node_a, feature_info);
  feature_node_is_part_of_pseudo_node(pseudo_node, node_b, feature_info);
  replace_node(node_a, pseudo_node, genome_nodes, auto_sr);
  remove_node((GtGenomeNode*) node_b, genome_nodes, auto_sr);
  return pseudo_node;
}

static GtFeatureNode* join_root_pair(GtFeatureNode *root_a,
                                     GtFeatureNode *root_b,
                                     GtFeatureInfo *feature_info,
                                     GtQueue *genome_nodes,
                                     AutomaticSequenceRegion *auto_sr)
{
  bool root_a_is_pseudo, root_b_is_pseudo;
  GtFeatureNode *master_root;
  gt_assert(root_a && root_b && feature_info && genome_nodes);
  root_a_is_pseudo = gt_feature_node_is_pseudo((GtFeatureNode*) root_a);
  root_b_is_pseudo = gt_feature_node_is_pseudo((GtFeatureNode*) root_b);
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

static void join_roots(GtArray *roots, GtFeatureInfo *feature_info,
                       GtQueue *genome_nodes,
                       AutomaticSequenceRegion *auto_sr)
{
  GtFeatureNode *master_root;
  unsigned long i;
  gt_assert(roots && feature_info && genome_nodes);
  master_root = *(GtFeatureNode**) gt_array_get(roots, 0);
  for (i = 1; i < gt_array_size(roots); i++) {
    master_root = join_root_pair(master_root,
                                 *(GtFeatureNode**) gt_array_get(roots, i),
                                 feature_info, genome_nodes, auto_sr);
  }
}

static int process_parent_attr(char *parent_attr, GtGenomeNode *feature_node,
                               bool *is_child, GtGFF3Parser *parser,
                               GtQueue *genome_nodes,
                               AutomaticSequenceRegion *auto_sr,
                               const char *filename, unsigned int line_number,
                               GtError *err)
{
  GtSplitter *parent_splitter;
  GtStrArray *valid_parents;
  unsigned long i;
  int had_err = 0;

  gt_error_check(err);
  gt_assert(parent_attr);

  valid_parents = gt_str_array_new();
  parent_splitter = gt_splitter_new();
  gt_splitter_split(parent_splitter, parent_attr, strlen(parent_attr), ',');
  gt_assert(gt_splitter_size(parent_splitter));

  for (i = 0; i < gt_splitter_size(parent_splitter); i++) {
    GtGenomeNode* parent_gf;
    const char *parent = gt_splitter_get_token(parent_splitter, i);
    parent_gf = (GtGenomeNode*) gt_feature_info_get(parser->feature_info,
                                                    parent);
    if (!parent_gf) {
      if (!parser->tidy) {
        gt_error_set(err, "%s \"%s\" on line %u in file \"%s\" has not been "
                     "previously defined (via \"%s=\")", PARENT_STRING, parent,
                     line_number, filename, ID_STRING);
        had_err = -1;
      }
      else {
        gt_warning("%s \"%s\" on line %u in file \"%s\" has not been "
                   "previously defined (via \"%s=\")", PARENT_STRING, parent,
                   line_number, filename, ID_STRING);
      }
    }
    else if (gt_str_cmp(gt_genome_node_get_seqid(parent_gf),
                        gt_genome_node_get_seqid(feature_node))) {
      gt_error_set(err, "child on line %u in file \"%s\" has different "
                   "sequence id than its parent on line %u ('%s' vs. '%s')",
                   gt_genome_node_get_line_number(feature_node), filename,
                   gt_genome_node_get_line_number(parent_gf),
                   gt_str_get(gt_genome_node_get_seqid(feature_node)),
                   gt_str_get(gt_genome_node_get_seqid(parent_gf)));
      had_err = -1;
    }
    else {
      gt_assert(parser->incomplete_node);
      if (i)
        feature_node = gt_genome_node_ref(feature_node);
      gt_feature_node_add_child((GtFeatureNode*) parent_gf,
                                (GtFeatureNode*) feature_node);
      *is_child = true;
      gt_str_array_add_cstr(valid_parents, parent);
    }
  }

  if (!had_err && !parser->tidy) {
    gt_assert(gt_splitter_size(parent_splitter) ==
              gt_str_array_size(valid_parents));
  }

  gt_splitter_delete(parent_splitter);

  /* make sure all (valid) parents have the same (pseudo-)root */
  if (!had_err && gt_str_array_size(valid_parents) >= 2) {
    GtArray *roots = find_roots(valid_parents, parser->feature_info);
    if (roots_differ(roots))
        join_roots(roots, parser->feature_info, genome_nodes, auto_sr);
    gt_array_delete(roots);
  }

  gt_str_array_delete(valid_parents);

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

static int check_missing_attributes(GtGenomeNode *this_feature,
                                    GtStrArray *this_attributes,
                                    GtFeatureNode *other_feature,
                                    const char *id, const char *filename,
                                    GtError *err)
{
  unsigned long i;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(this_feature && this_attributes && other_feature);
  for (i = 0; !had_err && i < gt_str_array_size(this_attributes); i++) {
    if (!gt_feature_node_get_attribute(other_feature,
                                      gt_str_array_get(this_attributes, i))) {
      gt_error_set(err, "the multi-feature with %s \"%s\" on line %u in file "
                "\"%s\" does not have a '%s' attribute which is present in its "
                "counterpart on line %u", ID_STRING, id,
                gt_genome_node_get_line_number((GtGenomeNode*) other_feature),
                filename, gt_str_array_get(this_attributes, i),
                gt_genome_node_get_line_number(this_feature));
      had_err = -1;
      break;
    }
  }
  return had_err;
}

static int compare_target_attribute(GtFeatureNode *new_gf,
                                    GtFeatureNode *old_gf,
                                    const char *id, GtError *err)
{
  unsigned long new_target_num, old_target_num;
  GtStr *new_target_str, *old_target_str;
  const char *new_target, *old_target;
  int had_err;
  gt_error_check(err);
  gt_assert(new_gf && old_gf);
  new_target = gt_feature_node_get_attribute(new_gf, TARGET_STRING);
  old_target = gt_feature_node_get_attribute(old_gf, TARGET_STRING);
  new_target_str = gt_str_new();
  old_target_str = gt_str_new();
  had_err = gt_gff3_parser_parse_target_attributes(new_target, &new_target_num,
                                                   new_target_str, NULL, NULL,
                                                   "", 0, NULL);
  gt_assert(!had_err); /* has been parsed already */
  had_err = gt_gff3_parser_parse_target_attributes(old_target, &old_target_num,
                                                   old_target_str, NULL, NULL,
                                                   "", 0, NULL);
  gt_assert(!had_err); /* has been parsed already */
  if (gt_str_cmp(new_target_str, old_target_str)) {
    gt_error_set(err, "the multi-feature with %s \"%s\" on line %u in file "
                 "\"%s\" has a different %s name than its counterpart on line "
                 "%u", ID_STRING, id,
                 gt_genome_node_get_line_number((GtGenomeNode*) new_gf),
                 gt_genome_node_get_filename((GtGenomeNode*) new_gf),
                 TARGET_STRING,
                 gt_genome_node_get_line_number((GtGenomeNode*) old_gf));
    had_err = -1;
  }
  gt_str_delete(old_target_str);
  gt_str_delete(new_target_str);
  return had_err;
}

static int compare_other_attribute(const char *attr_name, GtFeatureNode *new_gf,
                                   GtFeatureNode *old_gf, const char *id,
                                   GtError *err)
{
  gt_error_check(err);
  gt_assert(attr_name && new_gf && old_gf);
  if (strcmp(gt_feature_node_get_attribute(new_gf, attr_name),
             gt_feature_node_get_attribute(old_gf, attr_name))) {
    gt_error_set(err, "the multi-feature with %s \"%s\" on line %u in file "
                 "\"%s\" has a different attribute '%s' than its counterpart "
                 "on line %u ('%s' vs. '%s')",
                 ID_STRING, id,
                 gt_genome_node_get_line_number((GtGenomeNode*) new_gf),
                 gt_genome_node_get_filename((GtGenomeNode*) new_gf), attr_name,
                 gt_genome_node_get_line_number((GtGenomeNode*) old_gf),
                 gt_feature_node_get_attribute(new_gf, attr_name),
                 gt_feature_node_get_attribute(old_gf, attr_name));
    return -1;
  }
  return 0;
}

static int check_multi_feature_constrains(GtGenomeNode *new_gf,
                                          GtGenomeNode *old_gf, const char *id,
                                          const char *filename,
                                          unsigned int line_number,
                                          GtError *err)
{
  int had_err = 0;
  gt_error_check(err);
  gt_assert(new_gf && old_gf);
  gt_assert(!gt_feature_node_is_pseudo((GtFeatureNode*) new_gf));
  gt_assert(!gt_feature_node_is_pseudo((GtFeatureNode*) old_gf));
  /* the seqids have been compared already */
  gt_assert(!gt_str_cmp(gt_genome_node_get_seqid(new_gf),
                        gt_genome_node_get_seqid(old_gf)));
  /* check source */
  if (!had_err &&
      strcmp(gt_feature_node_get_source((GtFeatureNode*) new_gf),
             gt_feature_node_get_source((GtFeatureNode*) old_gf))) {
    gt_error_set(err, "the multi-feature with %s \"%s\" on line %u in file "
                 "\"%s\" has a different source than its counterpart on line "
                 "%u", ID_STRING, id, line_number, filename,
                 gt_genome_node_get_line_number(old_gf));
    had_err = -1;
  }
  /* check type */
  if (!had_err && gt_feature_node_get_type((GtFeatureNode*) new_gf) !=
                  gt_feature_node_get_type((GtFeatureNode*) old_gf)) {
    gt_error_set(err, "the multi-feature with %s \"%s\" on line %u in file "
                 "\"%s\" has a different type than its counterpart on line %u",
                 ID_STRING, id, line_number, filename,
                 gt_genome_node_get_line_number(old_gf));
    had_err = -1;
  }
  /* check strand */
  if (!had_err && gt_feature_node_get_strand((GtFeatureNode*) new_gf) !=
                  gt_feature_node_get_strand((GtFeatureNode*) old_gf)) {
    gt_error_set(err, "the multi-feature with %s \"%s\" on line %u in file "
                 "\"%s\" has a different strand than its counterpart on line "
                 "%u", ID_STRING, id, line_number, filename,
                 gt_genome_node_get_line_number(old_gf));
    had_err = -1;
  }
  /* check attributes (for target attribute only the name) */
  if (!had_err) {
    GtStrArray *new_attributes, *old_attributes;
    new_attributes =
      gt_feature_node_get_attribute_list((GtFeatureNode*) new_gf);
    old_attributes =
      gt_feature_node_get_attribute_list((GtFeatureNode*) old_gf);
    had_err = check_missing_attributes(new_gf, new_attributes,
                                       (GtFeatureNode*) old_gf, id,
                                       filename, err);
    if (!had_err) {
      had_err = check_missing_attributes(old_gf, old_attributes,
                                         (GtFeatureNode*) new_gf, id, filename,
                                         err);
    }
    if (!had_err) {
      unsigned long i;
      gt_assert(gt_str_array_size(new_attributes) ==
             gt_str_array_size(old_attributes));
      for (i = 0; !had_err && i < gt_str_array_size(new_attributes); i++) {
        const char *attr_name = gt_str_array_get(new_attributes, i);
        if (!strcmp(attr_name, TARGET_STRING)) {
          had_err = compare_target_attribute((GtFeatureNode*) new_gf,
                                             (GtFeatureNode*) old_gf, id, err);
        }
        else {
          had_err = compare_other_attribute(attr_name, (GtFeatureNode*) new_gf,
                                            (GtFeatureNode*) old_gf, id, err);
        }
      }
    }
    gt_str_array_delete(new_attributes);
    gt_str_array_delete(old_attributes);
  }
  return had_err;
}

static int parse_attributes(char *attributes, GtGenomeNode *feature_node,
                            bool *is_child, GtGFF3Parser *parser,
                            GtQueue *genome_nodes,
                            AutomaticSequenceRegion *auto_sr,
                            const char *filename, unsigned int line_number,
                            GtError *err)
{
  GtSplitter *attribute_splitter, *tmp_splitter, *parent_splitter;
  unsigned long i;
  int had_err = 0;

  gt_error_check(err);
  gt_assert(attributes);

  attribute_splitter = gt_splitter_new();
  tmp_splitter = gt_splitter_new();
  parent_splitter = gt_splitter_new();
  gt_splitter_split(attribute_splitter, attributes, strlen(attributes), ';');

  for (i = 0; !had_err && i < gt_splitter_size(attribute_splitter); i++) {
    const char *attr_tag = NULL;
    bool attr_valid = true;
    char *attr_value = NULL,
         *token = gt_splitter_get_token(attribute_splitter, i);
    if (strncmp(token, ".", 1) == 0) {
      if (gt_splitter_size(attribute_splitter) > 1) {
        gt_error_set(err, "more than one attribute token defined on line %u in "
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
      gt_splitter_reset(tmp_splitter);
      gt_splitter_split(tmp_splitter, token, strlen(token), '=');
      if (gt_splitter_size(tmp_splitter) != 2) {
        gt_error_set(err, "token \"%s\" on line %u in file \"%s\" does not "
                     "contain exactly one '='", token, line_number, filename);
        had_err = -1;
        break;
      }
      else {
        attr_tag = gt_splitter_get_token(tmp_splitter, 0);
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
        attr_value = gt_splitter_get_token(tmp_splitter, 1);
      }
    }
    if (!had_err && !strlen(attr_tag)) {
      attr_valid = false;
      if (parser->tidy) {
        gt_warning("attribute \"=%s\" on line %u in file \"%s\" has no tag; "
                   "skip it", attr_value, line_number, filename);
      }
      else {
        gt_error_set(err, "attribute \"=%s\" on line %u in file \"%s\" has no "
                     "tag", attr_value, line_number, filename);
        had_err = -1;
      }
    }
    if (!had_err && !strlen(attr_value)) {
      attr_valid = false;
      if (parser->tidy) {
        gt_warning("attribute \"%s=\" on line %u in file \"%s\" has no value; "
                   "skip it", attr_tag, line_number, filename);
      }
      else {
        gt_error_set(err, "attribute \"%s=\" on line %u in file \"%s\" has no "
                     "value", attr_tag, line_number, filename);
        had_err = -1;
      }
    }
    /* check for duplicate attributes */
    if (!had_err && attr_valid &&
         gt_feature_node_get_attribute((GtFeatureNode*) feature_node,
                                       attr_tag)) {
      gt_error_set(err, "more then one %s attribute on line %u in file \"%s\"",
                attr_tag, line_number, filename);
      had_err = -1;
    }
    /* save all attributes, although the Parent and ID attribute is newly
       created in GFF3 output */
    if (!had_err && attr_valid) {
      gt_feature_node_add_attribute((GtFeatureNode*) feature_node,
                                      attr_tag, attr_value);
    }
    /* some attributes require special care */
    if (!had_err && attr_valid) {
      if (!strcmp(attr_tag, ID_STRING)) {
        had_err = store_id(attr_value, (GtFeatureNode*) feature_node,
                           is_child, parser, genome_nodes, auto_sr, filename,
                           line_number, err);
      }
      else if (!strcmp(attr_tag, PARENT_STRING)) {
        had_err = process_parent_attr(attr_value, feature_node, is_child,
                                      parser, genome_nodes, auto_sr, filename,
                                      line_number, err);
      }
      else if (!strcmp(attr_tag, TARGET_STRING)) {
        /* the value of ``Target'' attributes have a special syntax which is
           checked here */
        had_err = gt_gff3_parser_parse_target_attributes(attr_value, NULL, NULL,
                                                         NULL, NULL, filename,
                                                         line_number, err);
      }
    }
  }

  if (!had_err &&
      gt_feature_node_is_multi((GtFeatureNode*) feature_node)) {
    had_err =
      check_multi_feature_constrains(feature_node, (GtGenomeNode*)
                  gt_feature_node_get_multi_representative((GtFeatureNode*)
                                                             feature_node),
                  gt_feature_node_get_attribute((GtFeatureNode*) feature_node,
                                                ID_STRING), filename,
                                                line_number, err);
  }

  gt_splitter_delete(parent_splitter);
  gt_splitter_delete(tmp_splitter);
  gt_splitter_delete(attribute_splitter);

  return had_err;
}

static void set_source(GtFeatureNode *feature_node, const char *source,
                       GtHashmap *source_to_str_mapping)
{
  GtStr *source_str;
  gt_assert(feature_node && source && source_to_str_mapping);
  source_str = gt_hashmap_get(source_to_str_mapping, source);
  if (!source_str) {
    source_str = gt_str_new_cstr(source);
    gt_hashmap_add(source_to_str_mapping, gt_str_get(source_str), source_str);
  }
  gt_assert(source_str);
  gt_feature_node_set_source(feature_node, source_str);
}

static int parse_regular_gff3_line(GtGFF3Parser *parser,
                                   GtQueue *genome_nodes,
                                   GtCstrTable *used_types, char *line,
                                   size_t line_length, GtStr *filenamestr,
                                   unsigned int line_number, GtError *err)
{
  GtGenomeNode *gn = NULL, *feature_node = NULL;
  GtSplitter *splitter;
  AutomaticSequenceRegion *auto_sr = NULL;
  GtStr *seqid_str = NULL;
  GtStrand gt_strand_value;
  float score_value;
  GtPhase phase_value;
  GtRange range;
  char *seqid = NULL, *source = NULL, *type = NULL, *start = NULL,
       *end = NULL, *score = NULL, *strand = NULL, *phase = NULL,
       *attributes = NULL, **tokens;
  const char *filename;
  bool score_is_defined, is_child = false;
  int had_err = 0;

  gt_error_check(err);

  filename = gt_str_get(filenamestr);

  /* create splitter */
  splitter = gt_splitter_new();

  /* parse */
  gt_splitter_split(splitter, line, line_length, '\t');
  if (gt_splitter_size(splitter) != 9UL) {
    gt_error_set(err, "line %u in file \"%s\" does not contain 9 tab (\\t) "
                   "separated fields", line_number, filename);
    had_err = -1;
  }
  if (!had_err) {
    tokens = gt_splitter_get_tokens(splitter);
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
  if (!had_err) {
    if (parser->type_checker &&
        !gt_type_checker_is_valid(parser->type_checker, type)) {
      gt_error_set(err, "type \"%s\" on line %u in file \"%s\" is not a valid "
                   "one", type, line_number, filename);
      had_err = -1;
    }
    else if (!gt_cstr_table_get(used_types, type))
      gt_cstr_table_add(used_types, type);
  }

  /* parse the range */
  if (!had_err) {
    if (parser->tidy) {
      had_err = gt_parse_range_tidy(&range, start, end, line_number, filename,
                                    err);
    }
    else
      had_err = gt_parse_range(&range, start, end, line_number, filename, err);
  }
  if (!had_err && range.start == 0) {
      gt_error_set(err, "illegal feature start 0 on line %u in file \"%s\" "
                   "(GFF3 files are 1-based)", line_number, filename);
      had_err = -1;
  }
  if (!had_err)
    had_err = add_offset_if_necessary(&range, parser, seqid, err);

  /* parse the score */
  if (!had_err) {
    had_err = gt_parse_score(&score_is_defined, &score_value, score,
                             line_number, filename, err);
  }

  /* parse the strand */
  if (!had_err) {
    had_err = gt_parse_strand(&gt_strand_value, strand, line_number, filename,
                              err);
  }

  /* parse the phase */
  if (!had_err)
    had_err = gt_parse_phase(&phase_value, phase, line_number, filename, err);

  /* get seqid */
  if (!had_err) {
    had_err = get_seqid_str(&seqid_str, seqid, range, &auto_sr, parser,
                            filename, line_number, err);
  }

  /* create the feature */
  if (!had_err) {
    feature_node = gt_feature_node_new(seqid_str, type, range.start, range.end,
                                       gt_strand_value);
    gt_genome_node_set_origin(feature_node, filenamestr, line_number);
  }

  /* set source */
  if (!had_err) {
    set_source((GtFeatureNode*) feature_node, source,
               parser->source_to_str_mapping);
  }

  /* parse the attributes */
  if (!had_err) {
    had_err = parse_attributes(attributes, feature_node, &is_child, parser,
                               genome_nodes, auto_sr, filename, line_number,
                               err);
  }

  if (!had_err && score_is_defined)
    gt_feature_node_set_score((GtFeatureNode*) feature_node, score_value);
  if (!had_err && phase_value != GT_PHASE_UNDEFINED)
    gt_feature_node_set_phase((GtFeatureNode*) feature_node, phase_value);

  if (!had_err) {
    gn = (is_child || auto_sr) ? NULL : feature_node;
    if (auto_sr && !is_child)
      gt_array_add(auto_sr->feature_nodes, feature_node);
  }
  else if (!is_child)
    gt_genome_node_delete(feature_node);

  if (!had_err && gn)
    gt_queue_add(genome_nodes, gn);

  /* free */
  gt_str_delete(seqid_str);
  gt_splitter_delete(splitter);

  return had_err;
}

static int parse_first_gff3_line(const char *line, const char *filename,
                                 bool tidy, GtError *err)
{
  int version, had_err = 0;
  gt_error_check(err);
  gt_assert(line && filename);
  if (strncmp(line, GFF_VERSION_PREFIX, strlen(GFF_VERSION_PREFIX))) {
    if (tidy) {
      gt_warning("line 1 in file \"%s\" does not begin with \"%s\", create "
                 "\"%s %d\" line automaticallly", filename, GFF_VERSION_PREFIX,
                 GFF_VERSION_PREFIX, GFF_VERSION);
      return 0;
    }
    else {
      gt_error_set(err, "line 1 in file \"%s\" does not begin with \"%s\"",
                   filename, GFF_VERSION_PREFIX);
      had_err = -1;
    }
  }
  if (!had_err) {
    line += strlen(GFF_VERSION_PREFIX);
    /* skip blanks */
    while (line[0] == ' ')
      line++;
    had_err = gt_parse_int_line(&version, line, 1, filename, err);
  }
  if (!had_err && version != GFF_VERSION) {
    gt_error_set(err, "GFF version %d does not equal required version %u ",
              version, GFF_VERSION);
    had_err = -1;
  }
  if (!had_err)
    return 1;
  return had_err;
}

static int gff3_parser_parse_fasta_entry(GtQueue *genome_nodes,
                                         const char *line, GtStr *filename,
                                         unsigned int line_number,
                                         GtGenFile *fpin, GtError *err)
{
  int had_err = 0;
  gt_error_check(err);
  gt_assert(line && line_number && fpin);
  if (line[0] != '>') {
    gt_error_set(err, "line %d does not start with '>' as expected",
                 line_number);
    had_err = -1;
  }
  if (!had_err) {
    GtGenomeNode *sequence_node;
    GtStr *sequence = gt_str_new();
    int cc;
    while ((cc = gt_genfile_xfgetc(fpin)) != EOF) {
      if (cc == '>') {
        gt_genfile_unget_char(fpin, cc);
        break;
      }
      if (cc != '\n' && cc != '\r' && cc != ' ')
        gt_str_append_char(sequence, cc);
    }
    sequence_node = gt_sequence_node_new(line+1, sequence);
    gt_genome_node_set_origin(sequence_node, filename, line_number);
    gt_queue_add(genome_nodes, sequence_node);
  }
  return had_err;
}

static int add_auto_sr_to_queue(GT_UNUSED void *key, void *value, void *data,
                                GT_UNUSED GtError *err)
{
  AutomaticSequenceRegion *auto_sr = value;
  GtQueue *genome_nodes = data;
  GtGenomeNode *gf;
  unsigned int i;
  gt_error_check(err);
  gt_assert(key && value && data);
  if (gt_array_size(auto_sr->feature_nodes)) {
    gt_queue_add(genome_nodes, auto_sr->sequence_region);
    auto_sr->sequence_region = NULL;
    for (i = 0; i < gt_array_size(auto_sr->feature_nodes); i++) {
      gf = *(GtGenomeNode**) gt_array_get(auto_sr->feature_nodes, i);
      gt_queue_add(genome_nodes, gf);
    }
    gt_array_reset(auto_sr->feature_nodes);
  }
  return 0;
}

static void process_undefined_sequence_regions(GtHashmap
                                               *undefined_sequence_regions,
                                               GtQueue *genome_nodes)
{
  int had_err;
  gt_assert(undefined_sequence_regions && genome_nodes);
  had_err = gt_hashmap_foreach(undefined_sequence_regions,
                               add_auto_sr_to_queue, genome_nodes, NULL);
  gt_assert(!had_err); /* add_auto_sr_to_queue() is sane */
  gt_hashmap_reset(undefined_sequence_regions);
}

static int parse_meta_gff3_line(GtGFF3Parser *parser, GtQueue *genome_nodes,
                                char *line, size_t line_length,
                                GtStr *filenamestr, unsigned int line_number,
                                GtError *err)
{
  char *tmpline, *tmplineend, *seqstart, *seqid = NULL;
  GtGenomeNode *gn;
  GtStr *changed_seqid = NULL;
  SimpleSequenceRegion *ssr = NULL;
  GtRange range;
  const char *filename;
  int had_err = 0;

  gt_error_check(err);
  gt_assert(line[0] == '#');

  filename = gt_str_get(filenamestr);

  if (line_length == 1 || line[1] != '#') {
    /* storing comment */
    gn = gt_comment_node_new(line+1);
    gt_genome_node_set_origin(gn, filenamestr, line_number);
    gt_queue_add(genome_nodes, gn);
  }
  else if (strcmp(line, GFF_FASTA_DIRECTIVE) == 0) {
    if (!parser->fasta_parsing) {
      parser->fasta_parsing = true;
      process_undefined_sequence_regions(parser->undefined_sequence_regions,
                                         genome_nodes);
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
      gt_error_set(err, "missing sequence region name on line %u in file "
                   "\"%s\"", line_number, filename);
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
        gt_error_set(err, "missing sequence region start on line %u in file "
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
        gt_error_set(err, "missing sequence region end on line %u in file "
                     "\"%s\"", line_number, filename);
        had_err = -1;
      }
    }
    if (!had_err) {
      had_err = gt_parse_range(&range, seqstart, tmpline, line_number, filename,
                               err);
    }
    if (!had_err  && range.start == 0) {
      gt_error_set(err, "illegal region start 0 on line %u in file \"%s\" "
                "(GFF3 files are 1-based)", line_number, filename);
      had_err = -1;
    }
    if (!had_err)
      had_err = add_offset_if_necessary(&range, parser, seqid, err);
    if (!had_err) {
      if (gt_hashmap_get(parser->undefined_sequence_regions, seqid)) {
        gt_error_set(err, "genome feature with id \"%s\" has been defined "
                     "before the corresponding \"%s\" definition on line %u in "
                     "file \"%s\"", seqid, GFF_SEQUENCE_REGION, line_number,
                     filename);
        had_err = -1;
      }
    }
    if (!had_err) {
      /* now we can create a sequence region node */
      gt_assert(seqid);
      ssr = gt_hashmap_get(parser->seqid_to_ssr_mapping, seqid);
      if (ssr) {
        gt_error_set(err, "the sequence region \"%s\" on line %u in file "
                     "\"%s\" has already been defined",
                     gt_str_get(ssr->seqid_str), line_number, filename);
        had_err = -1;
      }
      else {
        ssr = simple_sequence_region_new(seqid, range, line_number);
        gt_hashmap_add(parser->seqid_to_ssr_mapping, gt_str_get(ssr->seqid_str),
                    ssr);
      }
    }
    if (!had_err) {
      gt_assert(ssr);
      gn = gt_region_node_new(ssr->seqid_str, range.start, range.end);
      gt_genome_node_set_origin(gn, filenamestr, line_number);
      gt_queue_add(genome_nodes, gn);
    }
  }
  else if (strcmp(line, GFF_TERMINATOR) == 0) { /* terminator */
    /* now all nodes are complete */
    parser->incomplete_node = false;
    if (!parser->checkids)
      gt_feature_info_reset(parser->feature_info);
    parser->last_terminator = line_number;
  }
  else {
    gt_warning("skipping unknown meta line %u in file \"%s\": %s", line_number,
               filename, line);
  }
  gt_str_delete(changed_seqid);
  return had_err;
}

int gt_gff3_parser_parse_genome_nodes(GtGFF3Parser *parser, int *status_code,
                                      GtQueue *genome_nodes,
                                      GtCstrTable *used_types,
                                      GtStr *filenamestr,
                                      unsigned long long *line_number,
                                      GtGenFile *fpin, GtError *err)
{
  size_t line_length;
  GtStr *line_buffer;
  char *line;
  const char *filename;
  int rval, had_err = 0;

  gt_error_check(err);

  filename = gt_str_get(filenamestr);

  /* init */
  line_buffer = gt_str_new();

  while ((rval = gt_str_read_next_line_generic(line_buffer, fpin)) != EOF) {
    line = gt_str_get(line_buffer);
    line_length = gt_str_length(line_buffer);
    (*line_number)++;

    if (*line_number == 1) {
      had_err = parse_first_gff3_line(line, filename, parser->tidy, err);
      if (had_err == -1) /* error */
        break;
      if (had_err == 1) { /* line processed */
        gt_str_reset(line_buffer);
        had_err = 0;
        continue;
      }
      gt_assert(had_err == 0); /* line not processed */
    }
    if (line_length == 0) {
      gt_warning("skipping blank line %llu in file \"%s\"", *line_number,
                 filename);
    }
    else if (parser->fasta_parsing || line[0] == '>') {
      if (!parser->fasta_parsing) {
        parser->fasta_parsing = true;
        process_undefined_sequence_regions(parser->undefined_sequence_regions,
                                           genome_nodes);
      }
      had_err = gff3_parser_parse_fasta_entry(genome_nodes, line, filenamestr,
                                              *line_number, fpin, err);
      break;
    }
    else if (line[0] == '#') {
      had_err = parse_meta_gff3_line(parser, genome_nodes, line, line_length,
                                     filenamestr, *line_number, err);
      if (had_err ||
          (!parser->incomplete_node && gt_queue_size(genome_nodes))) {
        break;
      }
    }
    else {
      had_err = parse_regular_gff3_line(parser, genome_nodes, used_types, line,
                                        line_length, filenamestr, *line_number,
                                        err);
      if (had_err || (!parser->incomplete_node && gt_queue_size(genome_nodes)))
        break;
    }
    gt_str_reset(line_buffer);
  }

  if (had_err) {
    while (gt_queue_size(genome_nodes))
      gt_genome_node_delete(gt_queue_get(genome_nodes));
  }
  else if (rval == EOF) {
    /* the file has been parsed completely, add automatically created sequence
       regions to queue */
    process_undefined_sequence_regions(parser->undefined_sequence_regions,
                                       genome_nodes);
  }

  gt_str_delete(line_buffer);
  if (gt_queue_size(genome_nodes))
    *status_code = 0; /* at least one node was created */
  else
    *status_code = EOF;
  return had_err;
}

void gt_gff3_parser_reset(GtGFF3Parser *parser)
{
  gt_assert(parser);
  parser->fasta_parsing = false;
  gt_feature_info_reset(parser->feature_info);
  gt_hashmap_reset(parser->seqid_to_ssr_mapping);
  gt_hashmap_reset(parser->source_to_str_mapping);
  gt_hashmap_reset(parser->undefined_sequence_regions);
  parser->last_terminator = 0;
}

void gt_gff3_parser_delete(GtGFF3Parser *parser)
{
  if (!parser) return;
  gt_feature_info_delete(parser->feature_info);
  gt_hashmap_delete(parser->seqid_to_ssr_mapping);
  gt_hashmap_delete(parser->source_to_str_mapping);
  gt_hashmap_delete(parser->undefined_sequence_regions);
  gt_mapping_delete(parser->offset_mapping);
  gt_type_checker_delete(parser->type_checker);
  gt_free(parser);
}
