/*
  Copyright (c) 2006-2013 Gordon Gremme <gordon@gremme.org>
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
#include "core/array.h"
#include "core/assert_api.h"
#include "core/compat.h"
#include "core/cstr_api.h"
#include "core/hashmap.h"
#include "core/ma.h"
#include "core/md5_seqid.h"
#include "core/parseutils.h"
#include "core/queue.h"
#include "core/splitter.h"
#include "core/symbol_api.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "core/warning_api.h"
#include "extended/comment_node_api.h"
#include "extended/feature_info.h"
#include "extended/feature_node.h"
#include "extended/feature_node_iterator_api.h"
#include "extended/feature_type.h"
#include "extended/gap_str.h"
#include "extended/genome_node.h"
#include "extended/gff3_defines.h"
#include "extended/gff3_escaping.h"
#include "extended/gff3_parser.h"
#include "extended/mapping.h"
#include "extended/orphanage.h"
#include "extended/region_node.h"
#include "extended/xrf_checker_api.h"

struct GtGFF3Parser {
  GtFeatureInfo *feature_info;
  GtHashmap *seqid_to_ssr_mapping, /* maps seqids to simple sequence regions */
            *source_to_str_mapping;
  bool incomplete_node, /* at least one node is potentially incomplete */
       checkids,
       checkregions,
       strict,
       tidy,
       fasta_parsing, /* parser is in FASTA parsing mode */
       eof_emitted,
       gvf_mode;
  GtGenomeNode *gff3_pragma;
  GtWord offset;
  GtMapping *offset_mapping;
  GtOrphanage *orphanage;
  GtTypeChecker *type_checker;
  GtXRFChecker *xrf_checker;
  unsigned int last_terminator; /* line number of the last terminator */
};

typedef struct {
  GtStr *seqid_str;
  GtRange range;
  unsigned int line_number;
  bool pseudo,
       is_circular;
} SimpleSequenceRegion;

static SimpleSequenceRegion* simple_sequence_region_new(const char *seqid,
                                                        GtRange range,
                                                        unsigned int
                                                        line_number)
{
  SimpleSequenceRegion *ssr = gt_calloc(1, sizeof *ssr);
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
  parser = gt_calloc(1, sizeof *parser);
  parser->feature_info = gt_feature_info_new();
  parser->seqid_to_ssr_mapping = gt_hashmap_new(GT_HASH_STRING, NULL,
                                        (GtFree) simple_sequence_region_delete);
  parser->source_to_str_mapping = gt_hashmap_new(GT_HASH_STRING, NULL,
                                                 (GtFree) gt_str_delete);
  parser->offset = GT_UNDEF_WORD;
  parser->orphanage = gt_orphanage_new();
  parser->type_checker = type_checker ? gt_type_checker_ref(type_checker)
                                      : NULL;
  parser->xrf_checker = NULL;
  return parser;
}

void gt_gff3_parser_set_xrf_checker(GtGFF3Parser *parser,
                                    GtXRFChecker *xrf_checker)
{
  gt_assert(parser && xrf_checker);
  gt_xrf_checker_delete(parser->xrf_checker);
  parser->xrf_checker = gt_xrf_checker_ref(xrf_checker);
}

void gt_gff3_parser_check_id_attributes(GtGFF3Parser *parser)
{
  gt_assert(parser);
  parser->checkids = true;
}

void gt_gff3_parser_check_region_boundaries(GtGFF3Parser *parser)
{
  gt_assert(parser);
  parser->checkregions = true;
}

void gt_gff3_parser_do_not_check_region_boundaries(GtGFF3Parser *parser)
{
  gt_assert(parser);
  parser->checkregions = false;
}

void gt_gff3_parser_set_offset(GtGFF3Parser *parser, GtWord offset)
{
  gt_assert(parser);
  gt_assert(!parser->offset_mapping);
  parser->offset = offset;
}

void gt_gff3_parser_set_type_checker(GtGFF3Parser *parser,
                                     GtTypeChecker *type_checker)
{
  gt_assert(parser && type_checker);
  gt_type_checker_delete(parser->type_checker);
  parser->type_checker = gt_type_checker_ref(type_checker);
}

int gt_gff3_parser_set_offsetfile(GtGFF3Parser *parser, GtStr *offsetfile,
                                  GtError *err)
{
  gt_error_check(err);
  gt_assert(parser);
  gt_assert(parser->offset == GT_UNDEF_WORD);
  parser->offset_mapping = gt_mapping_new(offsetfile, "offsets",
                                          GT_MAPPINGTYPE_INTEGER, err);
  if (parser->offset_mapping)
    return 0;
  return -1;
}

void gt_gff3_parser_enable_strict_mode(GtGFF3Parser *parser)
{
  gt_assert(parser && !parser->tidy);
  parser->strict = true;
}

void gt_gff3_parser_enable_tidy_mode(GtGFF3Parser *parser)
{
  gt_assert(parser && !parser->strict);
  parser->tidy = true;
}

static int offset_possible(const GtRange *range, GtWord offset,
                           const char *filename, unsigned int line_number,
                           GtError *err)
{
  int had_err = 0;
  gt_error_check(err);
  if (offset < 0) {
    /* check for underflow */
    GtUword result = range->start + offset;
    if (result == 0) {
      gt_error_set(err, "==0");
      gt_error_set(err, "adding offset "GT_WD" to node on line %u in file "
                   "\"%s\" leads to start 0 (GFF3 files are 1-based)",
                   offset, line_number, filename);
      had_err = -1;
    }
    else if (result > range->start) {
      gt_error_set(err, "adding offset "GT_WD" to node on line %u in file "
                   "\"%s\" leads to underflow",
                   offset, line_number, filename);
      had_err = -1;
    }
  }
  return had_err;
}

static int add_offset_if_necessary(GtRange *range, GtGFF3Parser *parser,
                                   const char *seqid, const char *filename,
                                   GtUword line_number, GtError *err)
{
  int had_err = 0;
  gt_error_check(err);
  if (parser->offset != GT_UNDEF_WORD) {
    had_err = offset_possible(range, parser->offset, filename, line_number,
                              err);
    if (!had_err)
      *range = gt_range_offset(range, parser->offset);
  }
  else if (parser->offset_mapping) {
    GtWord offset;
    had_err = gt_mapping_map_integer(parser->offset_mapping, &offset, seqid,
                                     err);
    if (!had_err)
      had_err = offset_possible(range, offset, filename, line_number, err);
    if (!had_err)
      *range = gt_range_offset(range, offset);
  }
  return had_err;
}

static int verify_seqid(const GtStr *seqid, const char *filename,
                        unsigned int line_number, GtError *err)
{
  int had_err = 0;
  gt_error_check(err);
  gt_assert(seqid);
  if (gt_md5_seqid_has_prefix(gt_str_get(seqid))) {
    if (gt_str_length(seqid) < GT_MD5_SEQID_PREFIX_LEN +
                               GT_MD5_SEQID_HASH_LEN) {
      gt_error_set(err, "MD5 sequence ID '%s' on line %u in file \"%s\" is too "
                   "short", gt_str_get(seqid), line_number, filename);
      had_err = -1;
    }
    if (!had_err && gt_str_length(seqid) >= GT_MD5_SEQID_TOTAL_LEN) {
      const char *sid = gt_str_get(seqid);
      if (sid[GT_MD5_SEQID_TOTAL_LEN-1] != GT_MD5_SEQID_SEPARATOR) {
        gt_error_set(err, "MD5 sequence ID '%s' on line %u in file \"%s\" has "
                     "wrong separator '%c' (must be '%c')",
                     gt_str_get(seqid), line_number, filename,
                     sid[GT_MD5_SEQID_TOTAL_LEN-1], GT_MD5_SEQID_SEPARATOR);
        had_err = -1;
      }
    }
    if (!had_err && gt_str_length(seqid) == GT_MD5_SEQID_TOTAL_LEN) {
      gt_error_set(err, "MD5 sequence ID '%s' on line %u in file \"%s\" has "
                   "missing sequence ID after separator '%c'",
                   gt_str_get(seqid), line_number, filename,
                   GT_MD5_SEQID_SEPARATOR);
      had_err = -1;
    }
  }
  return had_err;
}

static int parse_target_attribute(const char *value, bool tidy,
                                  GtStr *target_id, GtRange *target_range,
                                  GtStrand *target_strand,
                                  GtStrArray *target_ids,
                                  GtArray *target_ranges,
                                  GtArray *target_strands,
                                  const char *filename,
                                  unsigned int line_number, GtError *err)
{
  GtUword num_of_tokens;
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
  if (!had_err)
    had_err = verify_seqid(unescaped_target, filename, line_number, err);
  if (!had_err && target_id)
    gt_str_append_str(target_id, unescaped_target);
  if (!had_err && target_ids)
    gt_str_array_add(target_ids, unescaped_target);
  /* parse target range */
  if (!had_err) {
    if (tidy) {
      had_err = gt_parse_range_tidy(&parsed_range,
                                    gt_splitter_get_token(splitter, 1),
                                    gt_splitter_get_token(splitter, 2),
                                    line_number, filename, err);
    }
    else {
      had_err = gt_parse_range(&parsed_range,
                               gt_splitter_get_token(splitter, 1),
                               gt_splitter_get_token(splitter, 2), line_number,
                               filename, err);
    }
  }
  if (!had_err && target_range)
    *target_range = parsed_range;
  if (!had_err && target_ranges)
    gt_array_add(target_ranges, parsed_range);
  /* parse target strand (if given) */
  if (!had_err) {
    if (gt_splitter_size(splitter) == 4) {
      had_err = gt_parse_strand(&parsed_strand,
                                gt_splitter_get_token(splitter, 3),
                                line_number, filename, err);
      if (!had_err && target_strand)
        *target_strand = parsed_strand;
      if (!had_err && target_strands)
        gt_array_add(target_strands, parsed_strand);
    }
    else {
      if (target_strand)
        *target_strand = GT_NUM_OF_STRAND_TYPES; /* undefined */
      if (target_strands) {
        parsed_strand = GT_NUM_OF_STRAND_TYPES; /* undefined */
        gt_array_add(target_strands, parsed_strand);
      }
    }
  }
  gt_free(escaped_target);
  gt_str_delete(unescaped_target);
  gt_splitter_delete(splitter);
  return had_err;
}

static int parse_target_attributes(const char *values, bool tidy,
                                   GtUword *num_of_targets,
                                   GtStr *first_target_id,
                                   GtRange *first_target_range,
                                   GtStrand *first_target_strand,
                                   GtStrArray *target_ids,
                                   GtArray *target_ranges,
                                   GtArray *target_strands,
                                   const char *filename,
                                   unsigned int line_number,
                                   GtError *err)
{
  GtSplitter *splitter;
  GtUword i;
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
    had_err = parse_target_attribute(gt_splitter_get_token(splitter, i), tidy,
                                     i ? NULL : first_target_id,
                                     i ? NULL : first_target_range,
                                     i ? NULL : first_target_strand,
                                     target_ids, target_ranges, target_strands,
                                     filename, line_number, err);
  }
  gt_free(targets);
  gt_splitter_delete(splitter);
  return had_err;
}

int gt_gff3_parser_parse_target_attributes(const char *values,
                                           GtUword *num_of_targets,
                                           GtStr *first_target_id,
                                           GtRange *first_target_range,
                                           GtStrand *first_target_strand,
                                           const char *filename,
                                           unsigned int line_number,
                                           GtError *err)
{
  return parse_target_attributes(values, false, num_of_targets, first_target_id,
                                 first_target_range, first_target_strand, NULL,
                                 NULL, NULL, filename, line_number, err);
}

int gt_gff3_parser_parse_all_target_attributes(const char *values, bool tidy,
                                               GtStrArray *target_ids,
                                               GtArray *target_ranges,
                                               GtArray *target_strands,
                                               const char *filename,
                                               unsigned int line_number,
                                               GtError *err)
{
  return parse_target_attributes(values, tidy, NULL, NULL, NULL, NULL,
                                 target_ids, target_ranges, target_strands,
                                 filename, line_number, err);
}

static int get_seqid_str(GtStr **seqid_str, const char *seqid, GtRange range,
                         GtGFF3Parser *parser, const char *filename,
                         unsigned int line_number, GtError *err)
{
  SimpleSequenceRegion *ssr;
  int had_err = 0;

  gt_error_check(err);

  ssr = gt_hashmap_get(parser->seqid_to_ssr_mapping, seqid);
  if (!ssr) {
    GtRange range;
    range.start = 0;
    range.end = ULONG_MAX;
    ssr = simple_sequence_region_new(seqid, range, line_number);
    ssr->pseudo = true;
    gt_hashmap_add(parser->seqid_to_ssr_mapping, gt_str_get(ssr->seqid_str),
                   ssr);
  }
  else if (parser->checkregions && !ssr->is_circular &&
           !gt_range_contains(&ssr->range, &range) && parser->checkregions) {
    gt_error_set(err, "range ("GT_WU","GT_WU") of feature on line %u in file "
                 "\"%s\" is not contained in range ("GT_WU","GT_WU") of "
                 "corresponding sequence region on line %u",
                 range.start, range.end, line_number, filename,
                 ssr->range.start, ssr->range.end, ssr->line_number);
    had_err = -1;
  }

  if (!had_err)
    *seqid_str = gt_str_ref(ssr->seqid_str);

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
                         GtFeatureNode *replacing_node, GtQueue *genome_nodes)
{
  ReplaceInfo replace_info;
  GT_UNUSED int rval;
  gt_assert(node_to_replace && replacing_node && genome_nodes);
  replace_info.node_to_replace = node_to_replace;
  replace_info.replacing_node = replacing_node;
  /* we go backwards in both cases, because we expect that the <node_to_replace>
     is near the end of the queue */
  rval = gt_queue_iterate_reverse(genome_nodes, replace_func, &replace_info,
                                  NULL);
  gt_assert(rval == 1);
}

static void remove_node(GtGenomeNode *genome_node, GtQueue *genome_nodes)
{
  gt_assert(genome_node && genome_nodes);
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
  id = gt_feature_node_get_attribute(child, GT_GFF_ID);
  gt_assert(id);
  gt_feature_info_add_pseudo_parent(feature_info, id, pseudo_node);
}

static int store_id(const char *id, GtFeatureNode *feature_node,
                    bool *is_child, GtGFF3Parser *parser, GtQueue *genome_nodes,
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
                   "terminator %s on line %u", GT_GFF_ID, id, line_number,
                   filename, gt_genome_node_get_line_number((GtGenomeNode*) fn),
                   GT_GFF_TERMINATOR, parser->last_terminator);
      had_err = -1;
    }
    /* check seqid */
    if (!had_err &&
        gt_str_cmp(gt_genome_node_get_seqid((GtGenomeNode*) feature_node),
                   gt_genome_node_get_seqid((GtGenomeNode*) fn))) {
      gt_error_set(err, "the multi-feature with %s \"%s\" on line %u in file "
                   "\"%s\" has a different sequence id than its counterpart on "
                   "line %u", GT_GFF_ID, id, line_number, filename,
                   gt_genome_node_get_line_number((GtGenomeNode*) fn));
      had_err = -1;
    }
    if (!had_err) {
      GtFeatureNode *pseudo_parent;
      bool has_parent = gt_feature_node_get_attribute(fn, GT_GFF_PARENT)
                        ? true : false;
      gt_assert(!gt_feature_node_is_pseudo(fn));
      pseudo_parent = gt_feature_info_get_pseudo_parent(parser->feature_info,
                                                        id);
      if (pseudo_parent || !gt_feature_node_is_multi(fn)) {
        if (!pseudo_parent) {
          gt_feature_node_make_multi_representative(fn);
          if (!has_parent) { /* create pseudo node */
            GtFeatureNode *pseudo_node = (GtFeatureNode*)
                                        gt_feature_node_new_pseudo_template(fn);
            feature_node_is_part_of_pseudo_node(pseudo_node, fn,
                                                parser->feature_info);
            replace_node(fn, pseudo_node, genome_nodes);
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
  else {
    gt_feature_info_add(parser->feature_info, id, feature_node);
    if (!parser->strict)
      gt_orphanage_reg_parent(parser->orphanage, id);
  }

  if (!had_err)
    parser->incomplete_node = true;

  return had_err;
}

static GtArray* find_roots(GtStrArray *parents, GtFeatureInfo *feature_info)
{
  GtArray *roots;
  GtUword i;
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
  GtUword i;
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
                                         GtQueue *genome_nodes)
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
  remove_node((GtGenomeNode*) pseudo_b, genome_nodes);
  gt_genome_node_delete((GtGenomeNode*) pseudo_b);
  return pseudo_a;
}

static GtFeatureNode* add_node_to_pseudo_node(GtFeatureNode *pseudo_node,
                                              GtFeatureNode *normal_node,
                                              GtFeatureInfo *feature_info,
                                              GtQueue *genome_nodes)
{
  gt_assert(pseudo_node &&
         gt_feature_node_is_pseudo((GtFeatureNode*) pseudo_node));
  gt_assert(normal_node &&
         !gt_feature_node_is_pseudo((GtFeatureNode*) normal_node));
  gt_assert(feature_info && genome_nodes);
  feature_node_is_part_of_pseudo_node(pseudo_node, normal_node, feature_info);
  remove_node((GtGenomeNode*) normal_node, genome_nodes);
  return pseudo_node;
}

static GtFeatureNode* create_pseudo_node(GtFeatureNode *node_a,
                                         GtFeatureNode *node_b,
                                         GtFeatureInfo *feature_info,
                                         GtQueue *genome_nodes)
{
  GtFeatureNode *pseudo_node;
  gt_assert(node_a && !gt_feature_node_is_pseudo((GtFeatureNode*) node_a));
  gt_assert(node_b && !gt_feature_node_is_pseudo((GtFeatureNode*) node_b));
  gt_assert(feature_info && genome_nodes);
  pseudo_node = (GtFeatureNode*)
                gt_feature_node_new_pseudo_template((GtFeatureNode*) node_a);
  feature_node_is_part_of_pseudo_node(pseudo_node, node_a, feature_info);
  feature_node_is_part_of_pseudo_node(pseudo_node, node_b, feature_info);
  replace_node(node_a, pseudo_node, genome_nodes);
  remove_node((GtGenomeNode*) node_b, genome_nodes);
  return pseudo_node;
}

static GtFeatureNode* join_root_pair(GtFeatureNode *root_a,
                                     GtFeatureNode *root_b,
                                     GtFeatureInfo *feature_info,
                                     GtQueue *genome_nodes)
{
  bool root_a_is_pseudo, root_b_is_pseudo;
  GtFeatureNode *master_root;
  gt_assert(root_a && root_b && feature_info && genome_nodes);
  root_a_is_pseudo = gt_feature_node_is_pseudo((GtFeatureNode*) root_a);
  root_b_is_pseudo = gt_feature_node_is_pseudo((GtFeatureNode*) root_b);
  if (root_a_is_pseudo && root_b_is_pseudo) {
    master_root = merge_pseudo_roots(root_a, root_b, feature_info,
                                     genome_nodes);
  }
  else if (root_a_is_pseudo) { /* !root_b_is_pseudo */
    master_root = add_node_to_pseudo_node(root_a, root_b, feature_info,
                                          genome_nodes);
  }
  else if (root_b_is_pseudo) { /* !root_a_is_pseudo */
    master_root = add_node_to_pseudo_node(root_b, root_a, feature_info,
                                          genome_nodes);
  }
  else { /* !root_a_is_pseudo && !root_b_is_pseudo */
    master_root = create_pseudo_node(root_a, root_b, feature_info,
                                     genome_nodes);
  }
  return master_root;
}

static void join_roots(GtArray *roots, GtFeatureInfo *feature_info,
                       GtQueue *genome_nodes)
{
  GtFeatureNode *master_root;
  GtUword i;
  gt_assert(roots && feature_info && genome_nodes);
  master_root = *(GtFeatureNode**) gt_array_get(roots, 0);
  for (i = 1; i < gt_array_size(roots); i++) {
    master_root = join_root_pair(master_root,
                                 *(GtFeatureNode**) gt_array_get(roots, i),
                                 feature_info, genome_nodes);
  }
}

static int process_child(GtGenomeNode *child, GtSplitter *parent_splitter,
                         GtFeatureInfo *feature_info, bool strict,
                         unsigned int last_terminator,
                         GtTypeChecker *type_checker, GtQueue *genome_nodes,
                         GtError *err)
{
  GtStrArray *valid_parents;
  GtGenomeNode* parent_gf;
  GtUword i;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(child && parent_splitter && feature_info && genome_nodes);
  valid_parents = gt_str_array_new();
  for (i = 0; !had_err && i < gt_splitter_size(parent_splitter); i++) {
    const char *parent = gt_splitter_get_token(parent_splitter, i);
    parent_gf = (GtGenomeNode*) gt_feature_info_get(feature_info,
                                                    parent);
    gt_assert(parent_gf);
    if (gt_genome_node_get_line_number(parent_gf) < last_terminator) {
      gt_error_set(err, "the child with %s \"%s\" on line %u in file "
                   "\"%s\" is separated from its corresponding %s on line %u "
                   "by terminator %s on line %u", GT_GFF_PARENT, parent,
                   gt_genome_node_get_line_number(child),
                   gt_genome_node_get_filename(child), GT_GFF_PARENT,
                   gt_genome_node_get_line_number(parent_gf), GT_GFF_TERMINATOR,
                   last_terminator);
      gt_genome_node_delete(child);
      had_err = -1;
    }
    if (!had_err) {
      if (i)
        child = gt_genome_node_ref(child);
      /* check for cycles (in strict mode, no cycles can occur) */
      if (!strict) {
        GtFeatureNodeIterator *fni;
        GtFeatureNode *node;
        fni = gt_feature_node_iterator_new((GtFeatureNode*) child);
        while (!had_err && (node = gt_feature_node_iterator_next(fni))) {
          if (node == (GtFeatureNode*) parent_gf) {
            gt_error_set(err, "linking the feature on line %u in file \"%s\" "
                         "to its %s with %s \"%s\" would cause a cycle",
                         gt_genome_node_get_line_number(child),
                         gt_genome_node_get_filename(child), GT_GFF_PARENT,
                         GT_GFF_ID, parent);
            gt_genome_node_delete(child);
            had_err = -1;
          }
        }
        gt_feature_node_iterator_delete(fni);
      }
    }
    if (!had_err && type_checker) {
      const char *parent_type, *child_type;
      /* check partof relationships */
      parent_type = gt_feature_node_get_type((GtFeatureNode*) parent_gf);
      child_type = gt_feature_node_get_type((GtFeatureNode*) child);
      if (!gt_type_checker_is_partof(type_checker, parent_type, child_type)) {
        gt_error_set(err, "the child feature with type '%s' on line %u in file "
                     "\"%s\" is not part-of parent feature with type '%s' "
                     "given on line %u (according to type checker '%s')",
                     child_type, gt_genome_node_get_line_number(child),
                     gt_genome_node_get_filename(child), parent_type,
                     gt_genome_node_get_line_number(parent_gf),
                     gt_type_checker_description(type_checker));
        gt_genome_node_delete(child);
        had_err = -1;
      }
    }
    if (!had_err) {
      gt_feature_node_add_child((GtFeatureNode*) parent_gf,
                                (GtFeatureNode*) child);
      gt_str_array_add_cstr(valid_parents, parent);
    }
  }
  if (!had_err) {
    gt_assert(gt_splitter_size(parent_splitter) ==
              gt_str_array_size(valid_parents));
    /* make sure all (valid) parents have the same (pseudo-)root */
    if (gt_str_array_size(valid_parents) >= 2) {
      GtArray *roots = find_roots(valid_parents, feature_info);
      if (roots_differ(roots))
        join_roots(roots, feature_info, genome_nodes);
      gt_array_delete(roots);
    }
  }
  gt_str_array_delete(valid_parents);
  return had_err;
}

static int process_parent_attr(char *parent_attr, GtGenomeNode *feature_node,
                               const char *id, bool *is_child,
                               GtGFF3Parser *parser, GtQueue *genome_nodes,
                               const char *filename, unsigned int line_number,
                               GtError *err)
{
  GtSplitter *parent_splitter;
  GtStrArray *missing_parents = NULL;
  bool orphaned_parent = false;
  GtGenomeNode* parent_gf;
  const char *parent;
  GtUword i;
  int had_err = 0;

  gt_error_check(err);
  gt_assert(parent_attr);

  parent_splitter = gt_splitter_new();
  gt_splitter_split(parent_splitter, parent_attr, strlen(parent_attr), ',');
  gt_assert(gt_splitter_size(parent_splitter));

  for (i = 0; i < gt_splitter_size(parent_splitter); i++) {
    parent = gt_splitter_get_token(parent_splitter, i);
    parent_gf = (GtGenomeNode*) gt_feature_info_get(parser->feature_info,
                                                    parent);
    if (!parent_gf) {
      if (parser->strict) {
        gt_error_set(err, "%s \"%s\" on line %u in file \"%s\" was not "
                     "previously defined (via \"%s=\")", GT_GFF_PARENT, parent,
                     line_number, filename, GT_GFF_ID);
        had_err = -1;
      }
      else {
        if (!missing_parents)
          missing_parents = gt_str_array_new();
        gt_str_array_add_cstr(missing_parents, parent);
      }
    }
    else if (!parser->strict &&
             gt_orphanage_is_orphan(parser->orphanage, parent)) {
      /* children of orphaned parends are orphans themself */
      orphaned_parent = true;
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
    else if (parent_gf == feature_node) {
      gt_error_set(err, "feature on line %u in file \"%s\" is "
                   "self-referential (%s and %s are the same)",
                   gt_genome_node_get_line_number(feature_node), filename,
                   GT_GFF_PARENT, GT_GFF_ID);
      had_err = -1;
    }
  }

  if (!had_err) {
    if (!missing_parents && !orphaned_parent) {
      had_err = process_child(feature_node, parent_splitter,
                              parser->feature_info, parser->strict,
                              parser->last_terminator, parser->type_checker,
                              genome_nodes, err);
    }
    else {
      gt_assert(!parser->strict);
      gt_orphanage_add(parser->orphanage, feature_node, id, missing_parents);
      parser->incomplete_node = true;
    }
    *is_child = true;
  }

  gt_splitter_delete(parent_splitter);
  gt_str_array_delete(missing_parents);

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
  GtUword i;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(this_feature && this_attributes && other_feature);
  for (i = 0; !had_err && i < gt_str_array_size(this_attributes); i++) {
    const char *attrkey = gt_str_array_get(this_attributes, i);
    if (strcmp(attrkey, "ID") != 0 && strcmp(attrkey, "Parent"))
      continue;

    if (!gt_feature_node_get_attribute(other_feature, attrkey)) {
      gt_error_set(err, "the multi-feature with %s \"%s\" on line %u in file "
                   "\"%s\" does not have a '%s' attribute which is present in "
                   "its counterpart on line %u", GT_GFF_ID, id,
                   gt_genome_node_get_line_number((GtGenomeNode*)
                                                  other_feature),
                   filename, attrkey,
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
  GtUword new_target_num, old_target_num;
  GtStr *new_target_str, *old_target_str;
  const char *new_target, *old_target;
  int had_err;
  gt_error_check(err);
  gt_assert(new_gf && old_gf);
  new_target = gt_feature_node_get_attribute(new_gf, GT_GFF_TARGET);
  old_target = gt_feature_node_get_attribute(old_gf, GT_GFF_TARGET);
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
                 "%u", GT_GFF_ID, id,
                 gt_genome_node_get_line_number((GtGenomeNode*) new_gf),
                 gt_genome_node_get_filename((GtGenomeNode*) new_gf),
                 GT_GFF_TARGET,
                 gt_genome_node_get_line_number((GtGenomeNode*) old_gf));
    had_err = -1;
  }
  gt_str_delete(old_target_str);
  gt_str_delete(new_target_str);
  return had_err;
}

static void tidy_multi_feature_with_different_parent(GtFeatureNode *new_gf,
                                                     GtFeatureNode *old_gf,
                                                     const char *id,
                                                     GtFeatureInfo
                                                     *feature_info)
{
  GtFeatureNode *parent;
  gt_assert(new_gf && old_gf && id && feature_info);
  gt_warning("the multi-feature with %s \"%s\" on line %u in file \"%s\" has a "
             "different attribute '%s' than its counterpart on line %u "
             "('%s' vs. '%s') -> tidy this as normal feature", GT_GFF_ID, id,
             gt_genome_node_get_line_number((GtGenomeNode*) new_gf),
             gt_genome_node_get_filename((GtGenomeNode*) new_gf), GT_GFF_PARENT,
             gt_genome_node_get_line_number((GtGenomeNode*) old_gf),
             gt_feature_node_get_attribute(new_gf, GT_GFF_PARENT),
             gt_feature_node_get_attribute(old_gf, GT_GFF_PARENT));
  /* the new feature is not a multi-feature anymore */
  gt_feature_node_unset_multi(new_gf);
  /* unset multi-feature status of the old feature, if it is the only child of
     its parent (with the given type) */
  parent = gt_feature_info_get(feature_info,
                               gt_feature_node_get_attribute(old_gf,
                                                             GT_GFF_PARENT));
  if (gt_feature_node_number_of_children_of_type(parent, old_gf) == 1)
    gt_feature_node_unset_multi(old_gf);
}

static int compare_other_attribute(const char *attr_name, GtFeatureNode *new_gf,
                                   GtFeatureNode *old_gf, const char *id,
                                   GtGFF3Parser *parser, GtError *err)
{
  gt_error_check(err);
  gt_assert(attr_name && new_gf && old_gf && parser);
  if (strcmp(gt_feature_node_get_attribute(new_gf, attr_name),
             gt_feature_node_get_attribute(old_gf, attr_name))) {
    if (parser->tidy && !strcmp(attr_name, GT_GFF_PARENT)) {
      tidy_multi_feature_with_different_parent(new_gf, old_gf, id,
                                               parser->feature_info);
    }
    else {
      gt_error_set(err, "the multi-feature with %s \"%s\" on line %u in file "
                   "\"%s\" has a different attribute '%s' than its counterpart "
                   "on line %u ('%s' vs. '%s')",
                   GT_GFF_ID, id,
                   gt_genome_node_get_line_number((GtGenomeNode*) new_gf),
                   gt_genome_node_get_filename((GtGenomeNode*) new_gf),
                   attr_name,
                   gt_genome_node_get_line_number((GtGenomeNode*) old_gf),
                   gt_feature_node_get_attribute(new_gf, attr_name),
                   gt_feature_node_get_attribute(old_gf, attr_name));
      return -1;
    }
  }
  return 0;
}

static int check_multi_feature_constrains(GtGenomeNode *new_gf,
                                          GtGenomeNode *old_gf, const char *id,
                                          GtGFF3Parser *parser,
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
                 "%u", GT_GFF_ID, id, line_number, filename,
                 gt_genome_node_get_line_number(old_gf));
    had_err = -1;
  }
  /* check type */
  if (!had_err && gt_feature_node_get_type((GtFeatureNode*) new_gf) !=
                  gt_feature_node_get_type((GtFeatureNode*) old_gf)) {
    gt_error_set(err, "the multi-feature with %s \"%s\" on line %u in file "
                 "\"%s\" has a different type than its counterpart on line %u",
                 GT_GFF_ID, id, line_number, filename,
                 gt_genome_node_get_line_number(old_gf));
    had_err = -1;
  }
  /* check strand */
  if (!had_err && gt_feature_node_get_strand((GtFeatureNode*) new_gf) !=
                  gt_feature_node_get_strand((GtFeatureNode*) old_gf)) {
    gt_error_set(err, "the multi-feature with %s \"%s\" on line %u in file "
                 "\"%s\" has a different strand than its counterpart on line "
                 "%u", GT_GFF_ID, id, line_number, filename,
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
                                       (GtFeatureNode*) old_gf, id, filename,
                                       err);
    if (!had_err) {
      had_err = check_missing_attributes(old_gf, old_attributes,
                                         (GtFeatureNode*) new_gf, id, filename,
                                         err);
    }
    if (!had_err) {
      GtUword i;
      for (i = 0; !had_err && i < gt_str_array_size(new_attributes); i++) {
        const char *attr_name = gt_str_array_get(new_attributes, i);
        if (!strcmp(attr_name, GT_GFF_TARGET)) {
          had_err = compare_target_attribute((GtFeatureNode*) new_gf,
                                             (GtFeatureNode*) old_gf, id, err);
        }
        else if (!strcmp(attr_name, GT_GFF_PARENT) ||
                 !strcmp(attr_name, GT_GFF_NAME)) {
          /* we relaxed the check and now only require the 'Parent' and the
             'Name' attribute to match*/
          had_err = compare_other_attribute(attr_name, (GtFeatureNode*) new_gf,
                                            (GtFeatureNode*) old_gf, id, parser,
                                            err);
        }
      }
    }
    gt_str_array_delete(new_attributes);
    gt_str_array_delete(old_attributes);
  }
  return had_err;
}

void gt_gff3_parser_build_target_str(GtStr *target, GtStrArray *target_ids,
                                     GtArray *target_ranges,
                                     GtArray *target_strands)
{
  GtUword i;
  gt_assert(target && target_ids && target_ranges && target_strands);
  for (i = 0; i < gt_str_array_size(target_ids); i++) {
    GtRange *range;
    GtStrand *strand;
    range = gt_array_get(target_ranges, i);
    strand = gt_array_get(target_strands, i);
    if (i)
      gt_str_append_char(target, ',');
    gt_str_append_cstr(target, gt_str_array_get(target_ids, i));
    gt_str_append_char(target, ' ');
    gt_str_append_ulong(target, range->start);
    gt_str_append_char(target, ' ');
    gt_str_append_ulong(target, range->end);
    if (*strand != GT_NUM_OF_STRAND_TYPES) {
      gt_str_append_char(target, ' ');
      gt_str_append_char(target, GT_STRAND_CHARS[*strand]);
    }
  }
}

static bool invalid_uppercase_gff3_attribute(const char *attr_tag)
{
  return (strcmp(attr_tag, GT_GFF_ID) &&
          strcmp(attr_tag, GT_GFF_NAME) &&
          strcmp(attr_tag, GT_GFF_ALIAS) &&
          strcmp(attr_tag, GT_GFF_PARENT) &&
          strcmp(attr_tag, GT_GFF_TARGET) &&
          strcmp(attr_tag, GT_GFF_GAP) &&
          strcmp(attr_tag, GT_GFF_DERIVES_FROM) &&
          strcmp(attr_tag, GT_GFF_NOTE) &&
          strcmp(attr_tag, GT_GFF_DBXREF) &&
          strcmp(attr_tag, GT_GFF_ONTOLOGY_TERM) &&
          strcmp(attr_tag, GT_GFF_START_RANGE) &&
          strcmp(attr_tag, GT_GFF_END_RANGE) &&
          strcmp(attr_tag, GT_GFF_IS_CIRCULAR));
}

static bool invalid_uppercase_gvf_attribute(const char *attr_tag)
{
  return (strcmp(attr_tag, GT_GVF_GENOTYPE) &&
          strcmp(attr_tag, GT_GVF_REFERENCE_SEQ) &&
          strcmp(attr_tag, GT_GVF_VARIANT_SEQ) &&
          strcmp(attr_tag, GT_GVF_VARIANT_FREQ) &&
          strcmp(attr_tag, GT_GVF_VARIANT_EFFECT) &&
          strcmp(attr_tag, GT_GVF_VARIANT_READS) &&
          strcmp(attr_tag, GT_GVF_TOTAL_READS) &&
          strcmp(attr_tag, GT_GVF_PHASED) &&
          strcmp(attr_tag, GT_GVF_START_RANGE) &&
          strcmp(attr_tag, GT_GVF_END_RANGE) &&
          strcmp(attr_tag, GT_GVF_INDIVIDUAL) &&
          strcmp(attr_tag, GT_GVF_REFERENCE_CODON) &&
          strcmp(attr_tag, GT_GVF_VARIANT_CODON) &&
          strcmp(attr_tag, GT_GVF_REFERENCE_AA) &&
          strcmp(attr_tag, GT_GVF_VARIANT_AA) &&
          strcmp(attr_tag, GT_GVF_BREAKPOINT_DETAIL) &&
          strcmp(attr_tag, GT_GVF_SEQUENCE_CONTEXT) &&
          strcmp(attr_tag, GT_GVF_ZYGOSITY));
}

static int parse_attributes(char *attributes, GtGenomeNode *feature_node,
                            bool *is_child, GtGFF3Parser *parser,
                            const char *seqid, GtQueue *genome_nodes,
                            const char *filename, unsigned int line_number,
                            GtError *err)
{
  GtSplitter *attribute_splitter, *tmp_splitter, *parent_splitter;
  char *id_value = NULL, *parent_value = NULL;
  GtUword i;
  int had_err = 0;

  gt_error_check(err);
  gt_assert(attributes);

  attribute_splitter = gt_splitter_new();
  tmp_splitter = gt_splitter_new();
  parent_splitter = gt_splitter_new();
  gt_splitter_split(attribute_splitter, attributes, strlen(attributes), ';');

  for (i = 0; !had_err && i < gt_splitter_size(attribute_splitter); i++) {
    const char *old_value;
    bool attr_valid = true;
    char *attr_tag = NULL,
         *attr_value = NULL,
         *token = gt_splitter_get_token(attribute_splitter, i);
    if (strncmp(token, ".", 1) == 0) {
      if (gt_splitter_size(attribute_splitter) > 1) {
        gt_error_set(err, "more than one attribute token defined on line %u in "
                     "file \"%s\", although the first one is '.'", line_number,
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
        if (parser->tidy && gt_splitter_size(tmp_splitter) == 1) {
          gt_warning("token \"%s\" on line %u in file \"%s\" does not "
                     "contain exactly one '='", token, line_number, filename);
          continue;
        }
        else {
          gt_error_set(err, "token \"%s\" on line %u in file \"%s\" does not "
                       "contain exactly one '='", token, line_number, filename);
          had_err = -1;
          break;
        }
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
    if (!had_err && attr_valid && isupper(attr_tag[0])) {
      /* check if uppercase attributes are the predefined ones */
      bool invalid;
      if (parser->gvf_mode)
        invalid = invalid_uppercase_gff3_attribute(attr_tag)
                    && invalid_uppercase_gvf_attribute(attr_tag);
      else
        invalid = invalid_uppercase_gff3_attribute(attr_tag);
      if (invalid) {
        if (parser->tidy) {
          gt_warning("illegal uppercase attribute \"%s\" on line %u in file "
                     "\"%s\"; change to lowercase", attr_tag, line_number,
                     filename);
          attr_tag[0] = tolower(attr_tag[0]);
        }
        else {
          gt_error_set(err, "illegal uppercase attribute \"%s\" on line %u in "
                            "file \"%s\" (uppercase attributes are reserved)",
                            attr_tag, line_number, filename);
          had_err = -1;
        }
      }
    }
    /* save all attributes, although the Parent and ID attributes are newly
       created in GFF3 output */
    if (!had_err && attr_valid) {
      if ((old_value = gt_feature_node_get_attribute((GtFeatureNode*)
                                                     feature_node, attr_tag))) {
        /* handle duplicate attribute */
        if (parser->tidy) {
          GtStr *combined_value;
          gt_warning("more than one %s attribute on line %u in file \"%s\"; "
                     "join them", attr_tag, line_number, filename);
          combined_value = gt_str_new_cstr(old_value);
          gt_str_append_char(combined_value, ',');
          gt_str_append_cstr(combined_value, attr_value);
          gt_feature_node_set_attribute((GtFeatureNode*) feature_node,
                                        attr_tag, gt_str_get(combined_value));
          gt_str_delete(combined_value);
        }
        else {
          gt_error_set(err, "more than one %s attribute on line %u in file "
                            "\"%s\"", attr_tag, line_number, filename);
          had_err = -1;
        }
      }
      else {
        gt_feature_node_add_attribute((GtFeatureNode*) feature_node, attr_tag,
                                      attr_value);
      }
    }
    /* some attributes require special care */
    if (!had_err && attr_valid) {
      if (!strcmp(attr_tag, GT_GFF_ID))
        id_value = attr_value; /* process later */
      else if (!strcmp(attr_tag, GT_GFF_PARENT))
        parent_value = attr_value; /* process later */
      else if (!strcmp(attr_tag, GT_GFF_IS_CIRCULAR)) {
        SimpleSequenceRegion *ssr;
        if (strcmp(attr_value, "true")) {
          gt_error_set(err, "value \"%s\" of %s attribute on line %u in file "
                       "\"%s\" does not equal \"true\"", attr_value,
                       GT_GFF_IS_CIRCULAR, line_number, filename);
          had_err = -1;
        }
        ssr = gt_hashmap_get(parser->seqid_to_ssr_mapping, seqid);
        gt_assert(ssr); /* XXX */
        gt_assert(!ssr->is_circular); /* XXX */
        ssr->is_circular = true;
      }
      else if (!strcmp(attr_tag, GT_GFF_TARGET)) {
        /* the value of ``Target'' attributes have a special syntax which is
           checked here */
        had_err = gt_gff3_parser_parse_target_attributes(attr_value, NULL, NULL,
                                                         NULL, NULL, filename,
                                                         line_number, err);
        if (had_err && parser->tidy) {
          GtStrArray *target_ids;
          GtArray *target_ranges, *target_strands;
          /* try to tidy up the ``Target'' attributes */
          gt_error_unset(err);
          target_ids = gt_str_array_new();
          target_ranges = gt_array_new(sizeof (GtRange));
          target_strands = gt_array_new(sizeof (GtStrand));
          had_err = gt_gff3_parser_parse_all_target_attributes(attr_value, true,
                                                               target_ids,
                                                               target_ranges,
                                                               target_strands,
                                                               filename,
                                                               line_number,
                                                               err);
          if (!had_err) {
            GtStr *new_target = gt_str_new();
            gt_gff3_parser_build_target_str(new_target, target_ids,
                                            target_ranges, target_strands);
            gt_feature_node_set_attribute((GtFeatureNode*) feature_node,
                                          GT_GFF_TARGET,
                                          gt_str_get(new_target));
            gt_str_delete(new_target);
          }
          gt_array_delete(target_strands);
          gt_array_delete(target_ranges);
          gt_str_array_delete(target_ids);
        }
      }
      else if (!strcmp(attr_tag, GT_GFF_DBXREF)
                 || !strcmp(attr_tag, GT_GFF_ONTOLOGY_TERM)) {
        if (parser->xrf_checker) {
          if (!gt_xrf_checker_is_valid(parser->xrf_checker, attr_value, err)) {
            had_err = -1;
          }
        }
      }
      else if (parser->type_checker && !strcmp(attr_tag, GT_GFF_GAP)) {
        GtGapStr *gs = NULL;
        GtRange rng = gt_genome_node_get_range(feature_node);
        if (gt_type_checker_is_a(parser->type_checker,
                                 gt_symbol("protein_match"),
                                 gt_feature_node_get_type((GtFeatureNode*)
                                                          feature_node))) {
          gs = gt_gap_str_new_protein(attr_value, err);
        } else {
          gs = gt_gap_str_new_nucleotide(attr_value, err);
        }
        if (!gs) {
          gt_assert(gt_error_is_set(err));
          had_err = -1;
        }
        if (!had_err) {
          if (gt_range_length(&rng) != gt_gap_str_length_reference(gs)) {
            gt_error_set(err, "length of aligned reference in %s attribute on "
                              "line %u in file \"%s\" (" GT_WU ") does not "
                              "match the length of its %s feature (" GT_WU ")",
                         GT_GFF_GAP, line_number, filename,
                         gt_gap_str_length_reference(gs),
                         gt_feature_node_get_type((GtFeatureNode*)
                                                  feature_node),
                         gt_range_length(&rng));
            had_err = -1;
          }
        }
        gt_gap_str_delete(gs);
      }
    }
  }

  if (!had_err && id_value) {
    had_err = store_id(id_value, (GtFeatureNode*) feature_node, is_child,
                       parser, genome_nodes, filename, line_number, err);
  }
  if (!had_err && parent_value) {
    had_err = process_parent_attr(parent_value, feature_node, id_value,
                                  is_child, parser, genome_nodes, filename,
                                  line_number, err);
  }

  if (!had_err && gt_feature_node_is_multi((GtFeatureNode*) feature_node)) {
    had_err =
      check_multi_feature_constrains(feature_node, (GtGenomeNode*)
                    gt_feature_node_get_multi_representative((GtFeatureNode*)
                                                             feature_node),
                    gt_feature_node_get_attribute((GtFeatureNode*) feature_node,
                                                  GT_GFF_ID),
                                     parser, filename, line_number, err);
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

static int parse_gff3_feature_line(GtGFF3Parser *parser,
                                   GtQueue *genome_nodes,
                                   GtCstrTable *used_types, char *line,
                                   size_t line_length, GtStr *filenamestr,
                                   unsigned int line_number, GtError *err)
{
  GtGenomeNode *gn = NULL, *feature_node = NULL;
  GtSplitter *splitter;
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
  if (gt_splitter_size(splitter) != 9) {
    if (parser->tidy && gt_splitter_size(splitter) == 10) {
      gt_warning("line %u in file \"%s\" does not contain 9 tab (\\t) "
                 "separated fields, dropping 10th field",
                 line_number, filename);
    }
    else {
      gt_error_set(err, "line %u in file \"%s\" does not contain 9 tab (\\t) "
                        "separated fields", line_number, filename);
      had_err = -1;
    }
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

  if (!had_err && parser->tidy && (start[0] == '.' || end[0] == '.')) {
    gt_warning("feature \"%s\" on line %u in file \"%s\" has undefined "
               "range, discarding feature", type, line_number, filename);
    gt_splitter_delete(splitter);
    return 0;
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
    if (parser->strict)
      had_err = gt_parse_range(&range, start, end, line_number, filename, err);
    else if (parser->tidy) {
      had_err = gt_parse_range_tidy(&range, start, end, line_number, filename,
                                    err);
    }
    else {
      had_err = gt_parse_range_correct_neg(&range, start, end, line_number,
                                           filename, err);
    }
  }
  if (!had_err && range.start == 0) {
      gt_error_set(err, "illegal feature start 0 on line %u in file \"%s\" "
                   "(GFF3 files are 1-based)", line_number, filename);
      had_err = -1;
  }
  if (!had_err) {
    had_err = add_offset_if_necessary(&range, parser, seqid, filename,
                                      line_number, err);
  }

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
    had_err = get_seqid_str(&seqid_str, seqid, range, parser, filename,
                            line_number, err);
  }
  /* verify seqid */
  if (!had_err)
    had_err = verify_seqid(seqid_str, filename, line_number, err);

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
                               seqid, genome_nodes, filename, line_number, err);
  }

  if (!had_err && score_is_defined)
    gt_feature_node_set_score((GtFeatureNode*) feature_node, score_value);
  if (!had_err && phase_value != GT_PHASE_UNDEFINED)
    gt_feature_node_set_phase((GtFeatureNode*) feature_node, phase_value);

  if (!had_err)
    gn = is_child ? NULL : feature_node;
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
                                 GtQueue *genome_nodes, GtStr *filenamestr,
                                 GtUint64 *line_number,
                                 bool *gvf_mode,
                                 bool tidy, GtError *err)
{
  int had_err = 0;
  const char *directive;
  char *data;
  gt_error_check(err);
  gt_assert(line && filename);

  if (strncmp(line, GT_GFF_VERSION_PREFIX, strlen(GT_GFF_VERSION_PREFIX))) {
    if (strncmp(line, GT_GVF_VERSION_PREFIX, strlen(GT_GVF_VERSION_PREFIX))) {
      if (tidy) {
        gt_warning("line 1 in file \"%s\" does not begin with \"%s\" or "
                   "\"%s\", create \"%s %d\" line automatically",
                   filename,
                   GT_GFF_VERSION_PREFIX, GT_GVF_VERSION_PREFIX,
                   GT_GFF_VERSION_PREFIX, GT_GFF_VERSION);
        return 0;
      }
      else {
        gt_error_set(err, "line 1 in file \"%s\" does not begin with \"%s\" "
                     "or \"%s\"",
                     filename, GT_GFF_VERSION_PREFIX, GT_GFF_VERSION_PREFIX);
        had_err = -1;
      }
    } else {
      *gvf_mode = true;
    }
  }
  directive = line + 2;
  data = strchr(directive, ' ');
  if (!data)
    data = strchr(line+2, '\t');
  if (data) {
    data[0] = '\0';
    data++;
  }
  else {
    gt_error_set(err, "version pragma encountered in line %u in file "
                      "\"%s\" does not have a version number",
                      (unsigned int) *line_number, filename);
    had_err = -1;
  }
  if (!had_err) {
    line += *gvf_mode
              ? strlen(GT_GVF_VERSION_PREFIX)
              : strlen(GT_GFF_VERSION_PREFIX);
    /* skip blanks */
    while (line[0] == ' ')
      line++;
  }
  if (!had_err && !(*gvf_mode)) {
    int version;
    had_err = gt_parse_int_line(&version, data, (unsigned int) *line_number,
                                filename, err);
    if (!had_err && version != GT_GFF_VERSION) {
      gt_error_set(err, "GFF version %s does not equal required version %s ",
                   data, GT_GFF_VERSION_STRING);
      had_err = -1;
    }
  }
  if (!had_err && *gvf_mode) {
    GtGenomeNode *gn;
    gn = gt_meta_node_new(directive, data);
    gt_genome_node_set_origin(gn, filenamestr, (unsigned int) *line_number);
    gt_queue_add(genome_nodes, gn);
  }
  if (!had_err)
    return 1;
  return had_err;
}

static int gff3_parser_parse_fasta_entry(GtQueue *genome_nodes,
                                         const char *line, GtStr *filename,
                                         unsigned int line_number,
                                         GtFile *fpin, GtError *err)
{
  int had_err = 0;
  gt_error_check(err);
  gt_assert(line && line_number);
  if (line[0] != '>') {
    gt_error_set(err, "line %d does not start with '>' as expected",
                 line_number);
    had_err = -1;
  }
  if (!had_err) {
    GtGenomeNode *sequence_node;
    GtStr *sequence = gt_str_new();
    int cc;
    while ((cc = gt_file_xfgetc(fpin)) != EOF) {
      if (cc == '>') {
        gt_file_unget_char(fpin, cc);
        break;
      }
      if (cc != '\n' && cc != '\r' && cc != ' ')
        gt_str_append_char(sequence, cc);
    }
    sequence_node = gt_sequence_node_new(line+1, sequence);
    gt_genome_node_set_origin(sequence_node, filename, line_number);
    gt_queue_add(genome_nodes, sequence_node);
    gt_str_delete(sequence);
  }
  return had_err;
}

static int process_orphans(GtOrphanage *orphanage, GtFeatureInfo *feature_info,
                           bool strict, unsigned int last_terminator,
                           GtTypeChecker *type_checker, GtQueue *genome_nodes,
                           GtError *err)
{
  GtGenomeNode *orphan;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(orphanage && feature_info && genome_nodes);
  while (!had_err && (orphan = gt_orphanage_get_orphan(orphanage))) {
    const char *parent_attr;
    char *parent_attr_dup;
    GtSplitter *splitter;
    GtUword i;
    parent_attr = gt_feature_node_get_attribute((GtFeatureNode*) orphan,
                                                GT_GFF_PARENT);
    gt_assert(parent_attr);
    parent_attr_dup = gt_cstr_dup(parent_attr);
    splitter = gt_splitter_new();
    gt_splitter_split(splitter, parent_attr_dup, strlen(parent_attr_dup), ',');
    gt_assert(gt_splitter_size(splitter));
    /* check for missing parents */
    for (i = 0; !had_err && i < gt_splitter_size(splitter); i++) {
      const char *parent = gt_splitter_get_token(splitter, i);
      if (gt_orphanage_parent_is_missing(orphanage, parent)) {
        gt_error_set(err, "%s \"%s\" on line %u in file \"%s\" was not defined "
                     "(via \"%s=\")", GT_GFF_PARENT, parent,
                     gt_genome_node_get_line_number(orphan),
                     gt_genome_node_get_filename(orphan) , GT_GFF_ID);
        gt_genome_node_delete(orphan);
        had_err = -1;
      }
    }
    if (!had_err) {
      had_err = process_child(orphan, splitter, feature_info, strict,
                              last_terminator, type_checker, genome_nodes, err);
    }
    gt_splitter_delete(splitter);
    gt_free(parent_attr_dup);
  }
  return had_err;
}

static bool invalid_gvf_pragma(const char *line)
{
  return (strncmp(line, GT_GVF_REFERENCE_FASTA, strlen(GT_GVF_REFERENCE_FASTA))
            && strncmp(line, GT_GVF_FEATURE_GFF3, strlen(GT_GVF_FEATURE_GFF3))
            && strncmp(line, GT_GVF_FILE_VERSION, strlen(GT_GVF_FILE_VERSION))
            && strncmp(line, GT_GVF_FILE_DATE, strlen(GT_GVF_FILE_DATE))
            && strncmp(line, GT_GVF_INDIVIDUAL_ID, strlen(GT_GVF_INDIVIDUAL_ID))
            && strncmp(line, GT_GVF_POPULATION, strlen(GT_GVF_POPULATION))
            && strncmp(line, GT_GVF_SEX, strlen(GT_GVF_SEX))
            && strncmp(line, GT_GVF_TECHNOLOGY_PLATFORM,
                       strlen(GT_GVF_TECHNOLOGY_PLATFORM))
            && strncmp(line, GT_GVF_TECHNOLOGY_PLATFORM_CLASS,
                       strlen(GT_GVF_TECHNOLOGY_PLATFORM_CLASS))
            && strncmp(line, GT_GVF_TECHNOLOGY_PLATFORM_NAME,
                       strlen(GT_GVF_TECHNOLOGY_PLATFORM_NAME))
            && strncmp(line, GT_GVF_TECHNOLOGY_PLATFORM_VERSION,
                       strlen(GT_GVF_TECHNOLOGY_PLATFORM_VERSION))
            && strncmp(line, GT_GVF_TECHNOLOGY_PLATFORM_MACHINE_ID,
                       strlen(GT_GVF_TECHNOLOGY_PLATFORM_MACHINE_ID))
            && strncmp(line, GT_GVF_TECHNOLOGY_PLATFORM_READ_LENGTH,
                       strlen(GT_GVF_TECHNOLOGY_PLATFORM_READ_LENGTH))
            && strncmp(line, GT_GVF_TECHNOLOGY_PLATFORM_READ_TYPE,
                       strlen(GT_GVF_TECHNOLOGY_PLATFORM_READ_TYPE))
            && strncmp(line, GT_GVF_TECHNOLOGY_PLATFORM_READ_PAIR_SPAN,
                       strlen(GT_GVF_TECHNOLOGY_PLATFORM_READ_PAIR_SPAN))
            && strncmp(line, GT_GVF_TECHNOLOGY_PLATFORM_AVERAGE_COVERAGE,
                       strlen(GT_GVF_TECHNOLOGY_PLATFORM_AVERAGE_COVERAGE))
            && strncmp(line, GT_GVF_SEQUENCING_SCOPE,
                       strlen(GT_GVF_SEQUENCING_SCOPE))
            && strncmp(line, GT_GVF_CAPTURE_METHOD,
                       strlen(GT_GVF_CAPTURE_METHOD))
            && strncmp(line, GT_GVF_CAPTURE_REGIONS,
                       strlen(GT_GVF_CAPTURE_REGIONS))
            && strncmp(line, GT_GVF_SEQUENCE_ALIGNMENT,
                       strlen(GT_GVF_SEQUENCE_ALIGNMENT))
            && strncmp(line, GT_GVF_VARIANT_CALLING,
                       strlen(GT_GVF_VARIANT_CALLING))
            && strncmp(line, GT_GVF_SAMPLE_DESCRIPTION,
                       strlen(GT_GVF_SAMPLE_DESCRIPTION))
            && strncmp(line, GT_GVF_GENOMIC_SOURCE,
                       strlen(GT_GVF_GENOMIC_SOURCE))
            && strncmp(line, GT_GVF_MULTI_INDIVIDUAL,
                       strlen(GT_GVF_MULTI_INDIVIDUAL))
            && strncmp(line, GT_GVF_DATA_SOURCE, strlen(GT_GVF_DATA_SOURCE))
            && strncmp(line, GT_GVF_SCORE_METHOD, strlen(GT_GVF_SCORE_METHOD))
            && strncmp(line, GT_GVF_SOURCE_METHOD, strlen(GT_GVF_SOURCE_METHOD))
            && strncmp(line, GT_GVF_ATTRIBUTE_METHOD,
                       strlen(GT_GVF_ATTRIBUTE_METHOD))
            && strncmp(line, GT_GVF_PHENOTYPE_DESCRIPTION,
                       strlen(GT_GVF_PHENOTYPE_DESCRIPTION))
            && strncmp(line, GT_GVF_PHASED_GENOTYPES,
                       strlen(GT_GVF_PHASED_GENOTYPES)));
}

static bool invalid_gff3_pragma(const char *line)
{
  return (strncmp(line, GT_GFF_SPECIES, strlen(GT_GFF_SPECIES))
            && strncmp(line, GT_GFF_FEATURE_ONTOLOGY,
                       strlen(GT_GFF_FEATURE_ONTOLOGY))
            && strncmp(line, GT_GFF_ATTRIBUTE_ONTOLOGY,
                       strlen(GT_GFF_ATTRIBUTE_ONTOLOGY))
            && strncmp(line, GT_GFF_SOURCE_ONTOLOGY,
                       strlen(GT_GFF_SOURCE_ONTOLOGY))
            && strncmp(line, GT_GFF_GENOME_BUILD,
                       strlen(GT_GFF_GENOME_BUILD)));
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
  else if (strcmp(line, GT_GFF_FASTA_DIRECTIVE) == 0)
    parser->fasta_parsing = true;
  else if ((strncmp(line, GT_GFF_SEQUENCE_REGION,
                    strlen(GT_GFF_SEQUENCE_REGION)) == 0)) {
    /* we are in a line starting with "##sequence-region" */
    tmpline = line + strlen(GT_GFF_SEQUENCE_REGION);
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
      while (tmpline < tmplineend && !(tmpline[0] == ' ' || tmpline[0] == '\t'))
        tmpline++;
      /* terminate seqid */
      *tmpline++ = '\0';
      /* skip blanks */
      while (tmpline < tmplineend && (tmpline[0] == ' ' || tmpline[0] == '\t'))
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
      while (tmpline <= tmplineend &&
             !(tmpline[0] == ' ' || tmpline[0] == '\t')) {
        tmpline++;
      }
      /* terminate seqstart */
      *tmpline++ = '\0';
      /* skip blanks */
      while (tmpline < tmplineend && (tmpline[0] == ' ' || tmpline[0] == '\t'))
        tmpline++;
      if (tmpline > tmplineend) {
        gt_error_set(err, "missing sequence region end on line %u in file "
                     "\"%s\"", line_number, filename);
        had_err = -1;
      }
    }
    if (!had_err) {
      if (parser->strict) {
        had_err = gt_parse_range(&range, seqstart, tmpline, line_number,
                                 filename, err);
      }
      else if (parser->tidy) {
        had_err = gt_parse_range_tidy(&range, seqstart, tmpline, line_number,
                                      filename, err);
      }
      else {
        had_err = gt_parse_range_correct_neg(&range, seqstart, tmpline,
                                             line_number, filename, err);
      }
    }
    if (!had_err && range.start == 0) {
      gt_error_set(err, "illegal region start 0 on line %u in file \"%s\" "
                   "(GFF3 files are 1-based)", line_number, filename);
      had_err = -1;
    }
    if (!had_err) {
      had_err = add_offset_if_necessary(&range, parser, seqid, filename,
                                        line_number, err);
    }
    if (!had_err) {
      /* now we can create a sequence region node */
      gt_assert(seqid);
      ssr = gt_hashmap_get(parser->seqid_to_ssr_mapping, seqid);
      if (ssr) {
        if (!ssr->pseudo) {
          gt_error_set(err, "the sequence region \"%s\" on line %u in file "
                       "\"%s\" has already been defined",
                       gt_str_get(ssr->seqid_str), line_number, filename);
          had_err = -1;
        }
        else {
          ssr->range = range;
          ssr->line_number = line_number;
          ssr->pseudo = false;
        }
      }
      else {
        ssr = simple_sequence_region_new(seqid, range, line_number);
        gt_hashmap_add(parser->seqid_to_ssr_mapping, gt_str_get(ssr->seqid_str),
                       ssr);
      }
    }
    if (!had_err)
      had_err = verify_seqid(ssr->seqid_str, filename, line_number, err);
    if (!had_err) {
      gt_assert(ssr);
      gn = gt_region_node_new(ssr->seqid_str, range.start, range.end);
      gt_genome_node_set_origin(gn, filenamestr, line_number);
      gt_queue_add(genome_nodes, gn);
    }
  }
  else if (strncmp(line, GT_GFF_TERMINATOR,
                   strlen(GT_GFF_TERMINATOR)) == 0) { /* terminator */
    if (line_length > strlen(GT_GFF_TERMINATOR)) {
      gt_warning("superfluous information after terminator in line %u of file "
                 "\"%s\": %s", line_number, filename, line);
    }
    /* now all nodes are complete */
    if (!parser->strict) {
      had_err = process_orphans(parser->orphanage, parser->feature_info,
                                parser->strict, parser->last_terminator,
                                parser->type_checker, genome_nodes, err);
    }
    parser->incomplete_node = false;
    if (!parser->checkids)
      gt_feature_info_reset(parser->feature_info);
    parser->last_terminator = line_number;
  }
  else if (strncmp(line, GT_GFF_VERSION_PREFIX,
                   strlen(GT_GFF_VERSION_PREFIX)) == 0) {
    if (parser->tidy) {
      gt_warning("skipping illegal GFF version pragma in line %u of file "
                 "\"%s\": %s", line_number, filename, line);
    }
    else {
      gt_error_set(err, "illegal GFF version pragma in line %u of file \"%s\": "
                   "%s", line_number, filename, line);
      had_err = -1;
    }
  }
  else if (strncmp(line, GT_GVF_VERSION_PREFIX,
                   strlen(GT_GVF_VERSION_PREFIX)) == 0) {
    char *data;
    bool make_node = false;

    if (parser->gvf_mode) {
      if (parser->tidy) {
        gt_warning("skipping illegal GVF version pragma in line %u of file "
                   "\"%s\": %s", line_number, filename, line);
      }
      else {
        gt_error_set(err, "illegal GVF version pragma in line %u of file "
                          "\"%s\": %s", line_number, filename, line);
        had_err = -1;
      }
    } else make_node = true;

    if (make_node) {
      if (!had_err) {
        data = strchr(line+2, ' ');
        if (!data)
          data = strchr(line+2, '\t');
        if (data) {
          data[0] = '\0';
          data++;
        } else {
          gt_error_set(err, "meta-directive encountered in line %u in file "
                            "\"%s\" does not have data: %s", line_number,
                            filename, line);
          had_err = -1;
        }
      }
      if (!had_err) {
        parser->gvf_mode = true;
        gn = gt_meta_node_new(line+2, data);
        gt_genome_node_set_origin(gn, filenamestr, line_number);
        gt_queue_add(genome_nodes, gn);
      }
    }
  }
  else {
    char *data;
    bool invalid;
    if (parser->gvf_mode)
      invalid = invalid_gff3_pragma(line) && invalid_gvf_pragma(line);
    else
      invalid = invalid_gff3_pragma(line);
    if (invalid) {
      gt_warning("unknown meta-directive encountered in line %u in file "
                 "\"%s\", keep as comment: %s", line_number, filename, line);
    }
    data = strchr(line+2, ' ');
    if (data) {
      data[0] = '\0';
      data++;
    }
    else {
      gt_error_set(err, "meta-directive encountered in line %u in file "
                 "\"%s\" does not have data: %s", line_number, filename, line);
      had_err = -1;
    }
    if (!had_err) {
      /* storing meta node */
      gn = gt_meta_node_new(line+2, data);
      gt_genome_node_set_origin(gn, filenamestr, line_number);
      gt_queue_add(genome_nodes, gn);
    }
  }
  gt_str_delete(changed_seqid);
  return had_err;
}

int gt_gff3_parser_parse_genome_nodes(GtGFF3Parser *parser, int *status_code,
                                      GtQueue *genome_nodes,
                                      GtCstrTable *used_types,
                                      GtStr *filenamestr,
                                      GtUint64 *line_number,
                                      GtFile *fpin, GtError *err)
{
  size_t line_length;
  GtStr *line_buffer;
  char *line;
  const char *filename;
  int rval, had_err = 0;

  gt_error_check(err);
  gt_assert(status_code && genome_nodes && used_types);

  filename = gt_str_get(filenamestr);

  /* init */
  line_buffer = gt_str_new();

  while ((rval = gt_str_read_next_line_generic(line_buffer, fpin)) != EOF) {
    line = gt_str_get(line_buffer);
    line_length = gt_str_length(line_buffer);
    (*line_number)++;

    if (*line_number == 1) {
      had_err = parse_first_gff3_line(line, filename, genome_nodes, filenamestr,
                                      line_number, &parser->gvf_mode,
                                      parser->tidy, err);
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
      gt_warning("skipping blank line "GT_LLU" in file \"%s\"", *line_number,
                 filename);
    }
    else if (parser->fasta_parsing || line[0] == '>') {
      parser->fasta_parsing = true;
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
      had_err = parse_gff3_feature_line(parser, genome_nodes, used_types, line,
                                        line_length, filenamestr, *line_number,
                                        err);
      if (had_err || (!parser->incomplete_node && gt_queue_size(genome_nodes)))
        break;
    }
    gt_str_reset(line_buffer);
  }

  if (!had_err && rval == EOF && *line_number == 0) {
    if (parser->tidy) {
      gt_warning("GFF3 file \"%s\" is empty", gt_str_get(filenamestr));
    } else {
      gt_error_set(err, "GFF3 file \"%s\" is empty", gt_str_get(filenamestr));
      had_err = -1;
    }
  }

  if (!had_err && !parser->strict) {
    had_err = process_orphans(parser->orphanage, parser->feature_info,
                              parser->strict, parser->last_terminator,
                              parser->type_checker, genome_nodes, err);
  }

  if (had_err) {
    while (gt_queue_size(genome_nodes))
      gt_genome_node_delete(gt_queue_get(genome_nodes));
  }
  else if (rval == EOF && !parser->eof_emitted) {
    GtGenomeNode *eofn = gt_eof_node_new();
    gt_genome_node_set_origin(eofn, filenamestr, *line_number+1);
    gt_queue_add(genome_nodes, eofn);
    parser->eof_emitted = true;
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
  parser->eof_emitted = false;
  gt_feature_info_reset(parser->feature_info);
  gt_hashmap_reset(parser->seqid_to_ssr_mapping);
  gt_hashmap_reset(parser->source_to_str_mapping);
  gt_orphanage_reset(parser->orphanage);
  parser->last_terminator = 0;
}

void gt_gff3_parser_delete(GtGFF3Parser *parser)
{
  if (!parser) return;
  gt_feature_info_delete(parser->feature_info);
  gt_hashmap_delete(parser->seqid_to_ssr_mapping);
  gt_hashmap_delete(parser->source_to_str_mapping);
  gt_mapping_delete(parser->offset_mapping);
  gt_orphanage_delete(parser->orphanage);
  gt_type_checker_delete(parser->type_checker);
  gt_xrf_checker_delete(parser->xrf_checker);
  gt_free(parser);
}
