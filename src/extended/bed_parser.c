/*
  Copyright (c) 2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>

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
#include "core/io.h"
#include "core/ma.h"
#include "core/parseutils.h"
#include "core/splitter.h"
#include "core/str.h"
#include "extended/bed_parser.h"
#include "extended/feature_node.h"
#include "extended/genome_node.h"
#include "extended/gff3_defines.h"
#include "extended/region_node_builder.h"

#define BROWSER_KEYWORD  "browser"
#define TRACK_KEYWORD    "track"
#define OFFSET_KEYWORD   "offset"

#define BLANK_CHAR       ' '
#define COMMENT_CHAR     '#'
#define TABULATOR_CHAR   '\t'
#define PAIR_SEPARATOR   '='
#define QUOTE_CHAR       '"'

struct GtBEDParser {
  GtRegionNodeBuilder *region_node_builder;
  GtQueue *feature_nodes;
  GtHashmap *seqid_to_str_mapping;
  GtStr *word,
        *another_word;
  char *feature_type,
       *thick_feature_type,
       *block_type;
  long offset;
};

GtBEDParser* gt_bed_parser_new(void)
{
  GtBEDParser *bed_parser = gt_calloc(1, sizeof *bed_parser);
  bed_parser->region_node_builder = gt_region_node_builder_new();
  bed_parser->feature_nodes = gt_queue_new();
  bed_parser->seqid_to_str_mapping = gt_hashmap_new(GT_HASH_STRING, NULL,
                                                    (GtFree) gt_str_delete);
  bed_parser->word = gt_str_new();
  bed_parser->another_word = gt_str_new();
  return bed_parser;
}

void gt_bed_parser_delete(GtBEDParser *bed_parser)
{
  if (!bed_parser) return;
  gt_free(bed_parser->block_type);
  gt_free(bed_parser->thick_feature_type);
  gt_free(bed_parser->feature_type);
  gt_str_delete(bed_parser->another_word);
  gt_str_delete(bed_parser->word);
  gt_hashmap_delete(bed_parser->seqid_to_str_mapping);
  while (gt_queue_size(bed_parser->feature_nodes))
    gt_genome_node_delete(gt_queue_get(bed_parser->feature_nodes));
  gt_queue_delete(bed_parser->feature_nodes);
  gt_region_node_builder_delete(bed_parser->region_node_builder);
  gt_free(bed_parser);
}

static int bed_parser_blank_line(GtIO *bed_file, GtError *err)
{
  int had_err = 0;
  gt_error_check(err);
  had_err = gt_io_expect(bed_file, BLANK_CHAR, err);
  while (!had_err) {
    char cc = gt_io_peek(bed_file);
    if (cc == GT_CARRIAGE_RETURN) {
      gt_io_next(bed_file);
      if (gt_io_peek(bed_file) == GT_END_OF_LINE)
        gt_io_next(bed_file);
      break;
    }
    else if ((cc == GT_END_OF_LINE) || (cc == GT_END_OF_FILE)) {
      gt_io_next(bed_file);
      break;
    }
    else
      had_err = gt_io_expect(bed_file, BLANK_CHAR, err);
  }
  return had_err;
}

static void rest_line(GtIO *bed_file)
{
  for (;;) {
    switch (gt_io_peek(bed_file)) {
      case GT_CARRIAGE_RETURN:
        gt_io_next(bed_file);
        if (gt_io_peek(bed_file) == GT_END_OF_LINE)
          gt_io_next(bed_file);
        return;
      case GT_END_OF_LINE:
        gt_io_next(bed_file);
        /*@fallthrough@*/
      case GT_END_OF_FILE:
        return;
      default:
        gt_io_next(bed_file);
    }
  }
}

static int bed_parser_comment_line(GtIO *bed_file, GtError *err)
{
  int had_err;
  gt_error_check(err);
  had_err = gt_io_expect(bed_file, COMMENT_CHAR, err);
  if (!had_err)
    rest_line(bed_file);
  return had_err;
}

static void word(GtStr *word, GtIO *bed_file)
{
  gt_str_reset(word);
  for (;;) {
    switch (gt_io_peek(bed_file)) {
      case BLANK_CHAR:
      case TABULATOR_CHAR:
      case PAIR_SEPARATOR:
      case GT_CARRIAGE_RETURN:
      case GT_END_OF_LINE:
      case GT_END_OF_FILE:
        return;
      default:
        gt_str_append_char(word, gt_io_next(bed_file));
    }
  }
}

static int quoted_word(GtStr *word, GtIO *bed_file, GtError *err)
{
  bool break_while = false;
  int had_err;
  gt_error_check(err);
  gt_str_reset(word);
  had_err = gt_io_expect(bed_file, QUOTE_CHAR, err);
  while (!had_err) {
    switch (gt_io_peek(bed_file)) {
      case QUOTE_CHAR:
      case GT_CARRIAGE_RETURN:
      case GT_END_OF_LINE:
      case GT_END_OF_FILE:
        break_while = true;
        break;
      default:
        gt_str_append_char(word, gt_io_next(bed_file));
    }
    if (break_while)
      break;
  }
  if (!had_err)
    had_err = gt_io_expect(bed_file, QUOTE_CHAR, err);
  return had_err;
}

static GtStr* get_seqid(GtBEDParser *bed_parser)
{
  GtStr *seqid = gt_hashmap_get(bed_parser->seqid_to_str_mapping,
                                gt_str_get(bed_parser->word));
  if (!seqid) {
    seqid = gt_str_new_cstr(gt_str_get(bed_parser->word));
    gt_hashmap_add(bed_parser->seqid_to_str_mapping, gt_str_get(seqid), seqid);
  }
  gt_assert(seqid);
  return seqid;
}

static bool bed_separator(GtIO *bed_file)
{
  char cc = gt_io_peek(bed_file);
  if (cc == BLANK_CHAR || cc == TABULATOR_CHAR)
    return true;
  return false;
}

static int skip_blanks(GtIO *bed_file, GtError *err)
{
  gt_error_check(err);
  if (!bed_separator(bed_file)) {
    gt_error_set(err, "file \"%s\": line %lu: expected blank or tabulator, got "
                      "'%c'", gt_io_get_filename(bed_file),
                      gt_io_get_line_number(bed_file), gt_io_peek(bed_file));
    return -1;
  }
  while (bed_separator(bed_file))
    gt_io_next(bed_file);
  return 0;
}

static int track_rest(GtBEDParser *bed_parser, GtIO *bed_file, GtError *err)
{
  char cc;
  int had_err = 0;
  gt_error_check(err);
  bed_parser->offset = 0; /* reset offset for new track line */
  if (bed_separator(bed_file)) /* skip to first attribute=value pair */
    had_err = skip_blanks(bed_file, err);
  while (!had_err &&
         (cc = gt_io_peek(bed_file)) != GT_END_OF_LINE &&
         cc != GT_CARRIAGE_RETURN) {
    /* parse attribute */
    word(bed_parser->word, bed_file);
    had_err = gt_io_expect(bed_file, PAIR_SEPARATOR, err);
    /* parse value */
    if (!had_err) {
      if (gt_io_peek(bed_file) == QUOTE_CHAR)
        had_err = quoted_word(bed_parser->another_word, bed_file, err);
      else
        word(bed_parser->another_word, bed_file);
    }
    /* process offset if necessary */
    if (!had_err && !strcmp(gt_str_get(bed_parser->word), OFFSET_KEYWORD)) {
      if (gt_parse_long(&bed_parser->offset,
                         gt_str_get(bed_parser->another_word))) {
        gt_error_set(err, "file \"%s\": line %lu: could not parse offset value "
                     "'%s'", gt_io_get_filename(bed_file),
                     gt_io_get_line_number(bed_file),
                     gt_str_get(bed_parser->another_word));
        had_err = -1;
      }
    }
    /* skip blanks up to next attribute or end-of-line */
    if (!had_err && bed_separator(bed_file))
      had_err = skip_blanks(bed_file, err);
  }
  /* the end of the line should now be reached */
  if (!had_err)
    had_err = gt_io_expect(bed_file, GT_END_OF_LINE, err);
  return had_err;
}

static int parse_bed_range(GtRange *range, GtStr *start, GtStr *end,
                           long offset, GtIO *bed_file, bool thick,
                           GtError *err)
{
  int had_err;
  gt_error_check(err);
  had_err = gt_parse_range(range, gt_str_get(start), gt_str_get(end),
                           gt_io_get_line_number(bed_file),
                           gt_io_get_filename(bed_file), err);
  /* BED has a weird numbering scheme: positions are 0-based, but the end
     position is not part of the feature. Transform to 1-based coordinates. */
  range->start++;
  /* Ranges defining a 'thick' region sometimes come with length 0 to
     designate that there are no thick regions. So do not fail here and
     handle that case later. */
  if (!thick) {
    if (!had_err && range->start > range->end) {
      gt_error_set(err, "file \"%s\": line %lu: BED feature has length 0",
                   gt_io_get_filename(bed_file),
                   gt_io_get_line_number(bed_file));
      had_err = -1;
    }
  }
  if (offset)
    *range = gt_range_offset(range, offset);
  return had_err;
}

static void construct_thick_feature(GtBEDParser *bed_parser, GtFeatureNode *fn,
                                    GtRange range)
{
  GtGenomeNode *thick_feature;
  const char *name;
  gt_assert(fn);
  thick_feature = gt_feature_node_new(gt_genome_node_get_seqid((GtGenomeNode*)
                                                               fn),
                                      bed_parser->thick_feature_type
                                      ? bed_parser->thick_feature_type
                                      : BED_THICK_FEATURE_TYPE,
                                      range.start, range.end,
                                      gt_feature_node_get_strand(fn));
  if ((name = gt_feature_node_get_attribute(fn, GT_GFF_NAME))) {
    gt_feature_node_add_attribute((GtFeatureNode*) thick_feature, GT_GFF_NAME,
                                  name);
  }
  gt_feature_node_set_score((GtFeatureNode*) thick_feature,
                            gt_feature_node_get_score(fn));
  gt_feature_node_set_strand((GtFeatureNode*) thick_feature,
                             gt_feature_node_get_strand(fn));
  gt_feature_node_add_child(fn, (GtFeatureNode*) thick_feature);
}

static int create_block_features(GtBEDParser *bed_parser, GtFeatureNode *fn,
                                 unsigned long block_count,
                                 GtSplitter *size_splitter,
                                 GtSplitter *start_splitter, GtIO *bed_file,
                                 GtError *err)
{
  unsigned long i;
  int had_err = 0;
  gt_assert(fn && block_count && size_splitter && start_splitter);
  gt_assert(gt_splitter_size(size_splitter) == block_count);
  gt_assert(gt_splitter_size(start_splitter) == block_count);
  for (i = 0; !had_err && i < block_count; i++) {
    unsigned long block_size, block_start, start, end;
    GtGenomeNode *block;
    const char *name;
    if (gt_parse_ulong(&block_size, gt_splitter_get_token(size_splitter, i))) {
      gt_error_set(err, "file \"%s\": line %lu: could not parse blockSize '%s'",
                   gt_io_get_filename(bed_file),
                   gt_io_get_line_number(bed_file),
                   gt_splitter_get_token(size_splitter, i));
      had_err = -1;
    }
    if (!had_err && gt_parse_ulong(&block_start,
                                   gt_splitter_get_token(start_splitter, i))) {
      gt_error_set(err, "file \"%s\": line %lu: could not parse blockStart "
                   "'%s'", gt_io_get_filename(bed_file),
                   gt_io_get_line_number(bed_file),
                   gt_splitter_get_token(start_splitter, i));
      had_err = -1;
    }
    if (!had_err) {
      start = gt_genome_node_get_start((GtGenomeNode*) fn) + block_start;
      end = start + block_size - 1;
      block = gt_feature_node_new(gt_genome_node_get_seqid((GtGenomeNode*) fn),
                                  bed_parser->block_type
                                  ? bed_parser->block_type
                                  : BED_BLOCK_TYPE,
                                  start, end, gt_feature_node_get_strand(fn));
      if ((name = gt_feature_node_get_attribute(fn, GT_GFF_NAME))) {
        gt_feature_node_add_attribute((GtFeatureNode*) block, GT_GFF_NAME,
                                      name);
      }
      gt_feature_node_set_score((GtFeatureNode*) block,
                                gt_feature_node_get_score(fn));
      gt_feature_node_set_strand((GtFeatureNode*) block,
                                 gt_feature_node_get_strand(fn));
      gt_feature_node_add_child(fn, (GtFeatureNode*) block);
    }
  }
  return had_err;
}

static void remove_terminal_comma(GtStr *str)
{
  gt_assert(str && gt_str_length(str));
  if (gt_str_get(str)[gt_str_length(str)-1] == ',')
    gt_str_set_length(str, gt_str_length(str)-1);
}

static int process_blocks(GtBEDParser *bed_parser, GtFeatureNode *fn,
                          unsigned long block_count, GtStr *block_sizes,
                          GtStr *block_starts, GtIO *bed_file, GtError *err)
{
  GtSplitter *size_splitter = NULL , *start_splitter = NULL;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(fn && block_count && block_sizes && block_starts);
  if (!gt_str_length(block_sizes)) {
    gt_error_set(err,
                 "file \"%s\": line %lu: blockCount given without blockSizes",
                 gt_io_get_filename(bed_file),
                 gt_io_get_line_number(bed_file));
    had_err = -1;
  }
  if (!had_err && !gt_str_length(block_starts)) {
    gt_error_set(err,
                 "file \"%s\": line %lu: blockCount given without blockStarts",
                 gt_io_get_filename(bed_file),
                 gt_io_get_line_number(bed_file));
    had_err = -1;
  }
  if (!had_err) {
    /* remove terminal commas found in real-world BED files */
    remove_terminal_comma(block_sizes);
    remove_terminal_comma(block_starts);
  }
  if (!had_err) {
    size_splitter = gt_splitter_new();
    gt_splitter_split(size_splitter, gt_str_get(block_sizes),
                      gt_str_length(block_sizes), ',');
    if (gt_splitter_size(size_splitter) != block_count) {
      gt_error_set(err, "file \"%s\": line %lu: blockSizes column does not "
                        "have blockCount=%lu many comma separated fields",
                   gt_io_get_filename(bed_file),
                   gt_io_get_line_number(bed_file), block_count);
      had_err = -1;
    }
  }
  if (!had_err) {
    start_splitter = gt_splitter_new();
    gt_splitter_split(start_splitter, gt_str_get(block_starts),
                      gt_str_length(block_starts), ',');
    if (gt_splitter_size(start_splitter) != block_count) {
      gt_error_set(err, "file \"%s\": line %lu: blockStarts column does not "
                        "have " "blockCount=%lu many comma separated fields",
                   gt_io_get_filename(bed_file),
                   gt_io_get_line_number(bed_file), block_count);
      had_err = -1;
    }
  }
  if (!had_err) {
    had_err = create_block_features(bed_parser, fn, block_count, size_splitter,
                                    start_splitter, bed_file, err);
  }
  gt_splitter_delete(start_splitter);
  gt_splitter_delete(size_splitter);
  return had_err;
}

static int bed_rest(GtBEDParser *bed_parser, GtIO *bed_file, GtError *err)
{
  unsigned long block_count = 0;
  GtGenomeNode *gn = NULL;
  GtRange range;
  GtStr *seqid;
  int had_err;
  gt_error_check(err);
  /* column 1.: chrom */
  seqid = get_seqid(bed_parser);
  had_err = skip_blanks(bed_file, err);
  /* column 2.: chromStart */
  if (!had_err) {
    word(bed_parser->word, bed_file);
    had_err = skip_blanks(bed_file, err);
  }
  /* column 3.: chromEnd */
  if (!had_err) {
    word(bed_parser->another_word, bed_file);
    had_err = parse_bed_range(&range, bed_parser->word,
                              bed_parser->another_word, bed_parser->offset,
                              bed_file, false, err);
  }
  if (!had_err) {
    /* add region */
    gt_region_node_builder_add_region(bed_parser->region_node_builder,
                                      gt_str_get(seqid), range);
    /* create feature */
    gn = gt_feature_node_new(seqid,
                             bed_parser->feature_type
                             ? bed_parser->feature_type
                             : BED_FEATURE_TYPE,
                             range.start, range.end, GT_STRAND_BOTH);
    gt_queue_add(bed_parser->feature_nodes, gn);
    if (bed_separator(bed_file))
      had_err = skip_blanks(bed_file, err);
  }
  /* optional column 4.: name */
  if (!had_err) {
    word(bed_parser->word, bed_file);
    if (gt_str_length(bed_parser->word)) {
      gt_feature_node_add_attribute((GtFeatureNode*) gn, GT_GFF_NAME,
                                    gt_str_get(bed_parser->word));
    }
    if (bed_separator(bed_file))
      had_err = skip_blanks(bed_file, err);
  }
  /* optional column 5.: score */
  if (!had_err) {
    word(bed_parser->word, bed_file);
    if (gt_str_length(bed_parser->word)) {
      bool score_is_defined;
      float score_value;
      had_err = gt_parse_score(&score_is_defined, &score_value,
                               gt_str_get(bed_parser->word),
                               gt_io_get_line_number(bed_file),
                               gt_io_get_filename(bed_file), err);
      if (!had_err && score_is_defined)
        gt_feature_node_set_score((GtFeatureNode*) gn, score_value);
    }
  }
  if (!had_err && bed_separator(bed_file))
    had_err = skip_blanks(bed_file, err);
  /* optional column 6.: strand */
  if (!had_err) {
    word(bed_parser->word, bed_file);
    if (gt_str_length(bed_parser->word)) {
      GtStrand strand;
      had_err = gt_parse_strand(&strand, gt_str_get(bed_parser->word),
                                gt_io_get_line_number(bed_file),
                                gt_io_get_filename(bed_file), err);
      if (!had_err)
        gt_feature_node_set_strand((GtFeatureNode*) gn, strand);
    }
  }
  if (!had_err && bed_separator(bed_file))
    had_err = skip_blanks(bed_file, err);
  /* optional column 7.: thickStart */
  if (!had_err) {
    word(bed_parser->word, bed_file);
    if (bed_separator(bed_file))
      had_err = skip_blanks(bed_file, err);
  }
  /* optional column 8.: thickEnd */
  if (!had_err) {
    word(bed_parser->another_word, bed_file);
    if (gt_str_length(bed_parser->another_word)) {
      gt_assert(gt_str_length(bed_parser->word));
      /* got a thickStart and a thickEnd -> construct corresponding feature */
      had_err = parse_bed_range(&range, bed_parser->word,
                                bed_parser->another_word, bed_parser->offset,
                                bed_file, true, err);
      if (!had_err && range.start <= range.end)
        construct_thick_feature(bed_parser, (GtFeatureNode*) gn, range);
    }
  }
  if (!had_err && bed_separator(bed_file))
    had_err = skip_blanks(bed_file, err);
  /* optional column 9.: itemRgb */
  if (!had_err) {
    word(bed_parser->word, bed_file);
    /* we do not use the RGB values */
    if (bed_separator(bed_file))
      had_err = skip_blanks(bed_file, err);
  }
  /* optional column 10.: blockCount */
  if (!had_err) {
    word(bed_parser->word, bed_file);
    if (gt_str_length(bed_parser->word)) {
      if (gt_parse_ulong(&block_count, gt_str_get(bed_parser->word))) {
        gt_error_set(err, "file \"%s\": line %lu: could not parse blockCount",
                     gt_io_get_filename(bed_file),
                     gt_io_get_line_number(bed_file));
        had_err = -1;
      }
      else {
        /* reset to parse/process blockSizes and blockStarts properly */
        gt_str_reset(bed_parser->word);
        gt_str_reset(bed_parser->another_word);
      }
    }
  }
  if (!had_err && bed_separator(bed_file))
    had_err = skip_blanks(bed_file, err);
  /* optional column 11.: blockSizes */
  if (!had_err) {
    word(bed_parser->word, bed_file);
    if (bed_separator(bed_file))
      had_err = skip_blanks(bed_file, err);
  }
  /* optional column 12.: blockStarts */
  if (!had_err) {
    word(bed_parser->another_word, bed_file);
    if (bed_separator(bed_file))
      had_err = skip_blanks(bed_file, err);
  }
  /* process blocks if necessary */
  if (!had_err && block_count) {
    had_err = process_blocks(bed_parser, (GtFeatureNode*) gn, block_count,
                             bed_parser->word, bed_parser->another_word,
                             bed_file, err);
  }
  /* the end of the line should now be reached */
  if (!had_err)
    had_err = gt_io_expect(bed_file, GT_END_OF_LINE, err);
  return had_err;
}

static int bed_line(GtBEDParser *bed_parser, GtIO *bed_file, GtError *err)
{
  int had_err = 0;
  gt_error_check(err);
  word(bed_parser->word, bed_file);
  if (!strcmp(gt_str_get(bed_parser->word), BROWSER_KEYWORD))
    rest_line(bed_file); /* ignore browser lines completely */
  else if (!strcmp(gt_str_get(bed_parser->word), TRACK_KEYWORD))
    had_err = track_rest(bed_parser, bed_file, err);
  else
    had_err = bed_rest(bed_parser, bed_file, err);
  return had_err;
}

static int parse_bed_file(GtBEDParser *bed_parser, GtIO *bed_file, GtError *err)
{
  int had_err = 0;
  gt_error_check(err);
  gt_assert(bed_file);
  while (!had_err && gt_io_has_char(bed_file)) {
    switch (gt_io_peek(bed_file)) {
      case BLANK_CHAR:
        had_err = bed_parser_blank_line(bed_file, err);
        break;
      case COMMENT_CHAR:
        had_err = bed_parser_comment_line(bed_file, err);
        break;
      case GT_CARRIAGE_RETURN:
        gt_io_next(bed_file);
        if (gt_io_peek(bed_file) == GT_END_OF_LINE)
          gt_io_next(bed_file);
        break;
      case GT_END_OF_LINE:
        gt_io_next(bed_file);
        break;
      default:
        had_err = bed_line(bed_parser, bed_file, err);
    }
  }
  if (!had_err)
    had_err = gt_io_expect(bed_file, GT_END_OF_FILE, err);
  return had_err;
}

int gt_bed_parser_parse(GtBEDParser *bed_parser, GtQueue *genome_nodes,
                        const char *filename, GtError *err)
{
  GtIO *bed_file;
  int had_err;
  gt_error_check(err);
  gt_assert(bed_parser && genome_nodes);
  bed_file = gt_io_new(filename, "r");
  /* parse BED file */
  had_err = parse_bed_file(bed_parser, bed_file, err);
  /* process created region and feature nodes */
  gt_region_node_builder_build(bed_parser->region_node_builder, genome_nodes);
  gt_region_node_builder_reset(bed_parser->region_node_builder);
  while (gt_queue_size(bed_parser->feature_nodes))
    gt_queue_add(genome_nodes, gt_queue_get(bed_parser->feature_nodes));
  gt_io_delete(bed_file);
  return had_err;
}

void gt_bed_parser_set_feature_type(GtBEDParser *bed_parser, const char *type)
{
  gt_assert(bed_parser && type);
  gt_free(bed_parser->feature_type);
  bed_parser->feature_type = gt_cstr_dup(type);
}

void gt_bed_parser_set_thick_feature_type(GtBEDParser *bed_parser,
                                          const char *type)
{
  gt_assert(bed_parser && type);
  gt_free(bed_parser->thick_feature_type);
  bed_parser->thick_feature_type = gt_cstr_dup(type);
}

void gt_bed_parser_set_block_type(GtBEDParser *bed_parser, const char *type)
{
  gt_assert(bed_parser && type);
  gt_free(bed_parser->block_type);
  bed_parser->block_type = gt_cstr_dup(type);
}
