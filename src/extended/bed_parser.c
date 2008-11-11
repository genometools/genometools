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
#include "core/hashmap.h"
#include "core/io.h"
#include "core/ma.h"
#include "core/parseutils.h"
#include "core/str.h"
#include "extended/bed_parser.h"
#include "extended/feature_node.h"

#define BROWSER_KEYWORD  "browser"
#define TRACK_KEYWORD    "track"

#define BED_FEATURE_TYPE "BED_feature"

#define BLANK_CHAR       ' '
#define COMMENT_CHAR     '#'
#define TABULATOR_CHAR   '\t'

struct GtBEDParser {
  GtHashmap *seqid_to_str_mapping;
  GtStr *word,
        *another_word;
};

GtBEDParser* gt_bed_parser_new(void)
{
  GtBEDParser *bed_parser = gt_malloc(sizeof *bed_parser);
  bed_parser->seqid_to_str_mapping = gt_hashmap_new(HASH_STRING, NULL,
                                                    (GtFree) gt_str_delete);
  bed_parser->word = gt_str_new();
  bed_parser->another_word = gt_str_new();
  return bed_parser;
}

void gt_bed_parser_delete(GtBEDParser *bed_parser)
{
  if (!bed_parser) return;
  gt_str_delete(bed_parser->another_word);
  gt_str_delete(bed_parser->word);
  gt_hashmap_delete(bed_parser->seqid_to_str_mapping);
  gt_free(bed_parser);
}

static int blank_line(GtIO *bed_file, GtError *err)
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

static int comment_line(GtIO *bed_file, GtError *err)
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
      case GT_CARRIAGE_RETURN:
      case GT_END_OF_LINE:
      case GT_END_OF_FILE:
        return;
      default:
        gt_str_append_char(word, gt_io_next(bed_file));
    }
  }
}

static GtStr* get_seqid(GtBEDParser *bed_parser)
{
  GtStr *seqid = gt_hashmap_get(bed_parser->seqid_to_str_mapping,
                                bed_parser->word);
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

static int bed_rest(GtBEDParser *bed_parser, GtQueue *genome_nodes,
                    GtIO *bed_file, GtError *err)
{
  GtGenomeNode *gn;
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
    had_err = gt_parse_range(&range, gt_str_get(bed_parser->word),
                             gt_str_get(bed_parser->another_word),
                             gt_io_get_line_number(bed_file),
                             gt_io_get_filename(bed_file), err);
  }
  if (!had_err && range.start == range.end) {
    gt_error_set(err, "file \"%s\": line %lu: BED feature has length 0",
                 gt_io_get_filename(bed_file), gt_io_get_line_number(bed_file));
    had_err = -1;
  }
  if (!had_err) {
    /* BED has a weird numbering scheme: positions are 0-based, but the end
       position is not part of the feature. Transform to 1-based coordinates. */
    range.start++;
    /* create feature */
    gn = gt_feature_node_new(seqid, BED_FEATURE_TYPE, range.start, range.end,
                             GT_STRAND_BOTH);
    gt_queue_add(genome_nodes, gn);
    if (bed_separator(bed_file))
      had_err = skip_blanks(bed_file, err);
  }
  /* optional column 4.: name */
  if (!had_err) {
    word(bed_parser->word, bed_file);
    if (gt_str_length(bed_parser->word)) {
      gt_feature_node_add_attribute((GtFeatureNode*) gn, "Name",
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
        gt_feature_node_set_strand(gn, strand);
    }
  }
  if (!had_err && bed_separator(bed_file))
    had_err = skip_blanks(bed_file, err);
  /* optional column 7.: thickStart */
  if (!had_err) {
    word(bed_parser->word, bed_file);
    if (gt_str_length(bed_parser->word)) {
      /* XXX */
    }
    if (bed_separator(bed_file))
      had_err = skip_blanks(bed_file, err);
  }
  /* optional column 8.: thickEnd */
  if (!had_err) {
    word(bed_parser->word, bed_file);
    if (gt_str_length(bed_parser->word)) {
      /* XXX */
    }
    if (bed_separator(bed_file))
      had_err = skip_blanks(bed_file, err);
  }
  /* optional column 9.: itemRgb */
  if (!had_err) {
    word(bed_parser->word, bed_file);
    if (gt_str_length(bed_parser->word)) {
      /* XXX */
    }
    if (bed_separator(bed_file))
      had_err = skip_blanks(bed_file, err);
  }
  /* optional column 10.: blockCount */
  if (!had_err) {
    word(bed_parser->word, bed_file);
    if (gt_str_length(bed_parser->word)) {
      /* XXX */
    }
    if (bed_separator(bed_file))
      had_err = skip_blanks(bed_file, err);
  }
  /* optional column 11.: blockSizes */
  if (!had_err) {
    word(bed_parser->word, bed_file);
    if (gt_str_length(bed_parser->word)) {
      /* XXX */
    }
    if (bed_separator(bed_file))
      had_err = skip_blanks(bed_file, err);
  }
  /* optional column 12.: blockStarts */
  if (!had_err) {
    word(bed_parser->word, bed_file);
    if (gt_str_length(bed_parser->word)) {
      /* XXX */
    }
    if (bed_separator(bed_file))
      had_err = skip_blanks(bed_file, err);
  }
  /* the end of the line should now be reached */
  if (!had_err)
    had_err = gt_io_expect(bed_file, GT_END_OF_LINE, err);
  return had_err;
}

static int bed_line(GtBEDParser *bed_parser, GtQueue *genome_nodes,
                    GtIO *bed_file, GtError *err)
{
  int had_err = 0;
  gt_error_check(err);
  word(bed_parser->word, bed_file);
  if (!strcmp(gt_str_get(bed_parser->word), BROWSER_KEYWORD))
    rest_line(bed_file); /* ignore browser lines completely */
  else if (!strcmp(gt_str_get(bed_parser->word), TRACK_KEYWORD))
    rest_line(bed_file); /* ignore track lines completely */
  else
    had_err = bed_rest(bed_parser, genome_nodes, bed_file, err);
  return had_err;
}

static int parse_bed_file(GtBEDParser *bed_parser, GtQueue *genome_nodes,
                          GtIO *bed_file, GtError *err)
{
  int had_err = 0;
  gt_error_check(err);
  gt_assert(bed_file);
  while (!had_err && gt_io_has_char(bed_file)) {
    switch (gt_io_peek(bed_file)) {
      case BLANK_CHAR:
        had_err = blank_line(bed_file, err);
        break;
      case COMMENT_CHAR:
        had_err = comment_line(bed_file, err);
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
        had_err = bed_line(bed_parser, genome_nodes, bed_file, err);
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
  had_err = parse_bed_file(bed_parser, genome_nodes, bed_file, err);
  gt_io_delete(bed_file);
  return had_err;
}
