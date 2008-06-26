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
#include "libgtcore/hashtable.h"
#include "libgtcore/ma.h"
#include "libgtcore/parseutils.h"
#include "libgtcore/splitter.h"
#include "libgtcore/undef.h"
#include "libgtcore/unused.h"
#include "libgtcore/warning.h"
#include "libgtext/comment.h"
#include "libgtext/genome_feature.h"
#include "libgtext/genome_node.h"
#include "libgtext/gff3_escaping.h"
#include "libgtext/gff3_parser.h"
#include "libgtext/mapping.h"
#include "libgtext/sequence_region.h"

struct GFF3Parser {
  Hashtable *id_to_genome_node_mapping,
            *seqid_to_ssr_mapping, /* maps seqids to simple sequence regions */
            *source_to_str_mapping,
            *undefined_sequence_regions; /* contains all (automatically created)
                                            sequence regions */
  bool incomplete_node, /* at least on node is potentially incomplete */
       checkids;
  long offset;
  Mapping *offset_mapping;
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
  unsigned long line;
} SimpleSequenceRegion;

static SimpleSequenceRegion* simple_sequence_region_new(const char *seqid,
                                                        Range range,
                                                        unsigned long line)
{
  SimpleSequenceRegion *ssr = ma_malloc(sizeof *ssr);
  ssr->seqid_str = str_new_cstr(seqid);
  ssr->range = range;
  ssr->line = line;
  return ssr;
}

static void simple_sequence_region_delete(SimpleSequenceRegion *ssr)
{
  if (!ssr) return;
  str_delete(ssr->seqid_str);
  ma_free(ssr);
}

GFF3Parser* gff3parser_new(bool checkids)
{
  GFF3Parser *gff3_parser = ma_malloc(sizeof (GFF3Parser));
  gff3_parser->id_to_genome_node_mapping = hashtable_new(HASH_STRING,
                                                         ma_free_func,
                                                         (FreeFunc)
                                                         genome_node_delete);
  gff3_parser->seqid_to_ssr_mapping = hashtable_new(HASH_STRING, NULL,
                                                    (FreeFunc)
                                                 simple_sequence_region_delete);
  gff3_parser->source_to_str_mapping = hashtable_new(HASH_STRING, NULL,
                                                     (FreeFunc) str_delete);
  gff3_parser->undefined_sequence_regions = hashtable_new(HASH_STRING, NULL,
                              (FreeFunc) automatic_sequence_region_delete);
  gff3_parser->incomplete_node = false;
  gff3_parser->checkids = checkids;
  gff3_parser->offset = UNDEF_LONG;
  gff3_parser->offset_mapping = NULL;
  return gff3_parser;
}

void gff3parser_set_offset(GFF3Parser *gff3_parser, long offset)
{
  assert(gff3_parser);
  assert(!gff3_parser->offset_mapping);
  gff3_parser->offset = offset;
}

int gff3parser_set_offsetfile(GFF3Parser *gff3_parser, Str *offsetfile,
                              Error *err)
{
  error_check(err);
  assert(gff3_parser);
  assert(gff3_parser->offset == UNDEF_LONG);
  gff3_parser->offset_mapping = mapping_new(offsetfile, "offsets",
                                            MAPPINGTYPE_INTEGER, err);
  if (gff3_parser->offset_mapping)
    return 0;
  return -1;

}

static int add_offset_if_necessary(Range *range, GFF3Parser *gff3_parser,
                                   const char *seqid, Error *err)
{
  long offset;
  int had_err = 0;
  error_check(err);
  if (gff3_parser->offset != UNDEF_LONG)
    *range = range_offset(*range, gff3_parser->offset);
  else if (gff3_parser->offset_mapping) {
    had_err = mapping_map_integer(gff3_parser->offset_mapping, &offset, seqid,
                                  err);
    if (!had_err)
      *range = range_offset(*range, offset);
  }
  return had_err;
}

static int parse_target_attribute(const char *value, Str *target_id,
                                  Range *target_range, Strand *target_strand,
                                  const char *filename,
                                  unsigned long line_number, Error *err)
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
    error_set(err, "Target attribute value '%s' on line %lu in file \"%s\" "
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
                                       unsigned long line_number, Error *err)
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

static int parse_regular_gff3_line(GFF3Parser *gff3_parser, Queue *genome_nodes,
                                   char *line, size_t line_length,
                                   Str *filenamestr, unsigned long line_number,
                                   Error *err)
{
  GenomeNode *gn = NULL, *genome_feature = NULL, *parent_gf;
  GenomeFeatureType gft;
  Splitter *splitter, *attribute_splitter, *tmp_splitter, *parents_splitter;
  AutomaticSequenceRegion *auto_sr = NULL;
  Str *seqid_str = NULL, *source_str, *changed_seqid = NULL;
  SimpleSequenceRegion *ssr;
  Strand strand_value;
  double score_value;
  Phase phase_value;
  Range range;
  char *id = NULL, *seqid = NULL, *source = NULL, *type = NULL, *start = NULL,
       *end = NULL, *score = NULL, *strand = NULL, *phase = NULL,
       *attributes = NULL, *token, *tmp_token, **tokens;
  const char *filename;
  bool is_child = false, seqid_str_created = false;
  unsigned long i;
  int had_err = 0;

  error_check(err);

  filename = str_get(filenamestr);

  /* create splitter */
  splitter = splitter_new();
  attribute_splitter = splitter_new();
  tmp_splitter = splitter_new();
  parents_splitter = splitter_new();

  /* parse */
  splitter_split(splitter, line, line_length, '\t');
  if (splitter_size(splitter) != 9UL) {
    error_set(err, "line %lu in file \"%s\" does not contain 9 tab (\\t) "
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
  if (!had_err && genome_feature_type_get(&gft, type) == -1) {
    error_set(err, "type \"%s\" on line %lu in file \"%s\" is not a valid one",
              type, line_number, filename);
    had_err = -1;
  }

  /* parse the range */
  if (!had_err)
    had_err = parse_range(&range, start, end, line_number, filename, err);
  if (!had_err && range.start == 0) {
      error_set(err, "illegal feature start 0 on line %lu in file \"%s\" "
                   "(GFF3 files are 1-based)", line_number, filename);
      had_err = -1;
  }
  if (!had_err)
    had_err = add_offset_if_necessary(&range, gff3_parser, seqid, err);

  /* parse the score */
  if (!had_err)
    had_err = parse_score(&score_value, score, line_number, filename, err);

  /* parse the strand */
  if (!had_err)
    had_err = parse_strand(&strand_value, strand, line_number, filename, err);

  /* parse the phase */
  if (!had_err)
    had_err = parse_phase(&phase_value, phase, line_number, filename, err);

  /* create the feature */
  if (!had_err) {
    genome_feature = genome_feature_new(gft, range, strand_value, filenamestr,
                                        line_number);
  }

  /* parse the unique id and the parents */
  if (!had_err) {
    splitter_split(attribute_splitter, attributes, strlen(attributes), ';');
    for (i = 0; i < splitter_size(attribute_splitter); i++) {
      const char *attr_tag, *attr_value;
      token = splitter_get_token(attribute_splitter, i);
      if (strncmp(token, ".", 1) == 0) {
        if (splitter_size(attribute_splitter) > 1) {
          error_set(err, "more than one attribute token defined on line %lu in "
                    "file \"%s\", altough the first one is '.'", line_number,
                    filename);
          had_err = -1;
        }
        if (!had_err)
          break; /* no attributes to parse */
      }
      else if (strncmp(token, ID_STRING, strlen(ID_STRING)) == 0) {
        if (id) {
          error_set(err, "more then one %s token on line %lu in file \"%s\"",
                    ID_STRING, line_number, filename);
          had_err = -1;
          break;
        }
        splitter_reset(tmp_splitter);
        splitter_split(tmp_splitter, token, strlen(token), '=');
        if (splitter_size(tmp_splitter) != 2) {
          error_set(err, "token \"%s\" on line %lu in file \"%s\" does not "
                    "contain exactly one '='", token, line_number, filename);
          had_err = -1;
          break;
        }
        id = cstr_dup(splitter_get_token(tmp_splitter, 1));
      }
      else if (strncmp(token, PARENT_STRING, strlen(PARENT_STRING)) == 0) {
        splitter_reset(tmp_splitter);
        splitter_split(tmp_splitter, token, strlen(token), '=');
        if (splitter_size(tmp_splitter) != 2) {
          error_set(err, "token \"%s\" on line %lu in file \"%s\" does not "
                    "contain exactly one '='", token, line_number, filename);
          had_err = -1;
          break;
        }
        tmp_token = splitter_get_token(tmp_splitter, 1);
        splitter_split(parents_splitter, tmp_token, strlen(tmp_token), ',');
        assert(splitter_size(parents_splitter));
      }
      else {
        /* add other attributes here */
        splitter_reset(tmp_splitter);
        splitter_split(tmp_splitter, token, strlen(token), '=');
        if (splitter_size(tmp_splitter) != 2) {
          error_set(err, "token \"%s\" on line %lu in file \"%s\" does not "
                    "contain exactly one '='", token, line_number, filename);
          had_err = -1;
          break;
        }
      }
      if (!had_err) {
        attr_tag = splitter_get_token(tmp_splitter, 0);
        attr_value = splitter_get_token(tmp_splitter, 1);
        if (!strcmp(attr_tag, "Target")) {
          /* the value of ``Target'' attributes have a special syntax which is
             checked here */
          had_err = gff3parser_parse_target_attributes(attr_value, NULL, NULL,
                                                       NULL, NULL, filename,
                                                       line_number, err);
        }
      }
      if (!had_err) {
        /* save all attributes, although the ID and Parent attributes are
           generated */
        genome_feature_add_attribute((GenomeFeature*) genome_feature, attr_tag,
                                     attr_value);
      }
    }
  }

  /* set seqid */
  if (!had_err) {
    ssr = hashtable_get(gff3_parser->seqid_to_ssr_mapping, seqid);
    if (!ssr) {
      /* sequence region has not been previously introduced -> check if one has
         already been created automatically */
      auto_sr = hashtable_get(gff3_parser->undefined_sequence_regions, seqid);
      if (!auto_sr) {
        /* sequence region has not been createad automatically -> do it now */
        warning("seqid \"%s\" on line %lu in file \"%s\" has not been "
                "previously introduced with a \"%s\" line, create such a line "
                "automatically", seqid, line_number, filename,
                GFF_SEQUENCE_REGION);
        auto_sr = automatic_sequence_region_new();
        seqid_str = str_new_cstr(seqid);
        seqid_str_created = true;
        auto_sr->sequence_region = sequence_region_new(seqid_str, range, NULL,
                                                       0);
        hashtable_add(gff3_parser->undefined_sequence_regions,
                      str_get(seqid_str), auto_sr);
      }
      else {
        /* get seqid string */
        seqid_str = genome_node_get_seqid(auto_sr->sequence_region);
        /* update the range of the sequence region */
        genome_node_set_range(auto_sr->sequence_region,
                              range_join(range,
                                         genome_node_get_range(auto_sr
                                                           ->sequence_region)));
      }
    }
    else {
      seqid_str = ssr->seqid_str;
      /* perform range check */
      if (!range_contains(ssr->range, range)) {
        error_set(err, "range (%lu,%lu) of feature on line %lu in file \"%s\" "
                  "is not contained in range (%lu,%lu) of corresponding "
                  "sequence region on line %lu", range.start, range.end,
                  line_number, filename, ssr->range.start, ssr->range.end,
                  ssr->line);
        had_err = -1;
      }
    }
  }
  if (!had_err) {
    assert(seqid_str);
    genome_node_set_seqid(genome_feature, seqid_str);
  }
  if (seqid_str_created)
    str_delete(seqid_str);

  /* set source */
  if (!had_err) {
    source_str = hashtable_get(gff3_parser->source_to_str_mapping, source);
    if (!source_str) {
      source_str = str_new_cstr(source);
      hashtable_add(gff3_parser->source_to_str_mapping, str_get(source_str),
                    source_str);
    }
    assert(source_str);
    genome_feature_set_source(genome_feature, source_str);
  }

  if (!had_err && score_value != UNDEF_DOUBLE)
    genome_feature_set_score((GenomeFeature*) genome_feature, score_value);
  if (!had_err && phase_value != PHASE_UNDEFINED)
    genome_feature_set_phase(genome_feature, phase_value);

  /* store id */
  if (!had_err && id) {
    if ((gn = hashtable_get(gff3_parser->id_to_genome_node_mapping, id))) {
      /* this id has been used already -> raise error */
      error_set(err, "the %s \"%s\" on line %lu in file \"%s\" has been used "
                "already for the feature defined on line %lu", ID_STRING, id,
                line_number, filename, genome_node_get_line_number(gn));
      had_err = -1;
      ma_free(id);
    }
    else {
      gff3_parser->incomplete_node = true;
      hashtable_add(gff3_parser->id_to_genome_node_mapping, id,
                    genome_node_ref(genome_feature));
    }
  }
  else
    ma_free(id);

  if (!had_err) {
    for (i = 0; i < splitter_size(parents_splitter); i++) {
      parent_gf = hashtable_get(gff3_parser->id_to_genome_node_mapping,
                                splitter_get_token(parents_splitter, i));
      if (!parent_gf) {
        error_set(err, "%s \"%s\" on line %lu in file \"%s\" has not been "
                  "previously defined (via \"%s=\")", PARENT_STRING,
                  splitter_get_token(parents_splitter, i), line_number,
                  filename, ID_STRING);
        had_err = -1;
      }
      else if (str_cmp(genome_node_get_seqid(parent_gf),
                       genome_node_get_seqid(genome_feature))) {
        error_set(err, "child on line %lu in file \"%s\" has different "
                  "sequence id than its parent on line %lu ('%s' vs. '%s')",
                  genome_node_get_line_number(genome_feature), filename,
                  genome_node_get_line_number(parent_gf),
                  str_get(genome_node_get_seqid(genome_feature)),
                  str_get(genome_node_get_seqid(parent_gf)));
        had_err = -1;
      }
      else {
        assert(gff3_parser->incomplete_node);
        genome_node_is_part_of_genome_node(parent_gf, genome_feature);
        is_child = true;
      }
    }
  }

  if (!had_err) {
    gn = (is_child || auto_sr) ? NULL : genome_feature;
    if (auto_sr && !is_child)
      array_add(auto_sr->genome_features, genome_feature);
  }
  else
    genome_node_delete(genome_feature);

  if (!had_err && gn)
    queue_add(genome_nodes, gn);

  /* free */
  str_delete(changed_seqid);
  splitter_delete(splitter);
  splitter_delete(attribute_splitter);
  splitter_delete(tmp_splitter);
  splitter_delete(parents_splitter);

  return had_err;
}

static int parse_meta_gff3_line(GFF3Parser *gff3_parser, Queue *genome_nodes,
                                char *line, size_t line_length,
                                Str *filenamestr, unsigned long line_number,
                                Error *err)
{
  char *tmpline, *tmplineend, *seqid = NULL;
  GenomeNode *gn;
  Str *changed_seqid = NULL;
  SimpleSequenceRegion *ssr = NULL;
  Range range;
  const char *filename;
  int rval, had_err = 0;

  error_check(err);
  assert(line[0] == '#');

  filename = str_get(filenamestr);

  if (line_length == 1 || line[1] != '#') {
    /* storing comment */
    gn = comment_new(line+1, filenamestr, line_number);
    queue_add(genome_nodes, gn);
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
      error_set(err, "missing sequence region name on line %lu in file \"%s\"",
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
        error_set(err, "missing sequence region start on line %lu in file "
                  "\"%s\"", line_number, filename);
        had_err = -1;
      }
    }
    if (!had_err) {
      if ((rval = sscanf(tmpline, "%lu", &range.start)) != 1) {
        error_set(err, "could not parse region start on line %lu in file "
                  "\"%s\"", line_number, filename);
        had_err = -1;
      }
    }
    if (!had_err  && range.start == 0) {
      error_set(err, "illegal region start 0 on line %lu in file \"%s\" "
                "(GFF3 files are 1-based)", line_number, filename);
      had_err = -1;
    }
    if (!had_err) {
      /* skip non-blanks */
      while (tmpline <= tmplineend && !(tmpline[0] == ' '))
        tmpline++;
      /* skip blanks */
      while (tmpline < tmplineend && tmpline[0] == ' ')
        tmpline++;
      if (tmpline > tmplineend) {
        error_set(err, "missing sequence region end on line %lu in file \"%s\"",
                  line_number, filename);
        had_err = -1;
      }
    }
    if (!had_err && (rval = sscanf(tmpline, "%lu", &range.end)) != 1) {
      error_set(err, "could not parse region end on line %lu in file \"%s\"",
                line_number, filename);
      had_err = -1;
    }
    if (!had_err && (range.start > range.end)) {
      error_set(err, "region start %lu is larger then region end %lu on line "
                "%lu in file \"%s\"", range.start, range.end, line_number,
                filename);
      had_err = -1;
    }
    if (!had_err)
      had_err = add_offset_if_necessary(&range, gff3_parser, seqid, err);
   if (!had_err) {
      if (hashtable_get(gff3_parser->undefined_sequence_regions, seqid)) {
        error_set(err, "genome feature with id \"%s\" has been defined before "
                  "the corresponding \"%s\" definition on line %lu in file "
                  "\"%s\"", seqid, GFF_SEQUENCE_REGION, line_number, filename);
        had_err = -1;
      }
    }
    if (!had_err) {
      /* now we can create a sequence region node */
      assert(seqid);
      ssr = hashtable_get(gff3_parser->seqid_to_ssr_mapping, seqid);
      if (ssr) {
        error_set(err, "the sequence region \"%s\" on line %lu in file \"%s\" "
                  "has already been defined", str_get(ssr->seqid_str),
                  line_number, filename);
        had_err = -1;
      }
      else {
        ssr = simple_sequence_region_new(seqid, range, line_number);
        hashtable_add(gff3_parser->seqid_to_ssr_mapping,
                      str_get(ssr->seqid_str), ssr);
      }
    }
    if (!had_err) {
      assert(ssr);
      gn = sequence_region_new(ssr->seqid_str, range, filenamestr, line_number);
      queue_add(genome_nodes, gn);
    }
  }
  else if (strcmp(line, GFF_TERMINATOR) == 0) { /* terminator */
    /* now all nodes are complete */
    gff3_parser->incomplete_node = false;
    if (!gff3_parser->checkids)
      hashtable_reset(gff3_parser->id_to_genome_node_mapping);
  }
  else {
    warning("skipping unknown meta line %lu in file \"%s\": %s", line_number,
            filename, line);
  }
  str_delete(changed_seqid);
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

int gff3parser_parse_genome_nodes(int *status_code, GFF3Parser *gff3_parser,
                                  Queue *genome_nodes, Str *filenamestr,
                                  unsigned long long *line_number,
                                  GenFile *fpin, Error *err)
{
  size_t line_length;
  Str *line_buffer;
  char *line;
  const char *filename;
  int rval, version, had_err = 0;

  error_check(err);

  filename = str_get(filenamestr);

  /* init */
  line_buffer = str_new();

  while ((rval = str_read_next_line_generic(line_buffer, fpin)) != EOF) {
    line = str_get(line_buffer);
    line_length = str_length(line_buffer);
    (*line_number)++;

    if (*line_number == 1) {
      if (strncmp(line, GFF_VERSION_PREFIX, strlen(GFF_VERSION_PREFIX))) {
        error_set(err, "line %llu in file \"%s\" does not begin with \"%s\"",
                  *line_number, filename, GFF_VERSION_PREFIX);
        had_err = -1;
        break;
      }
      if (!had_err) {
        line += strlen(GFF_VERSION_PREFIX);
        /* skip blanks */
        while (line[0] == ' ')
          line++;
        if (parse_int_line(&version, line, *line_number, filename, err)) {
          had_err = -1;
          break;
        }
      }
      if (!had_err) {
        if (version != GFF_VERSION) {
          error_set(err, "GFF version %d does not equal required version %u ",
                    version, GFF_VERSION);
          had_err = -1;
          break;
        }
      }
    }
    else if (line_length == 0) {
      warning("skipping blank line %llu in file \"%s\"", *line_number,
              filename);
    }
    else if (line[0] == '#') {
      had_err = parse_meta_gff3_line(gff3_parser, genome_nodes, line,
                                     line_length, filenamestr, *line_number,
                                     err);
      if (had_err ||
          (!gff3_parser->incomplete_node && queue_size(genome_nodes))) {
        break;
      }
    }
    else {
      had_err = parse_regular_gff3_line(gff3_parser, genome_nodes, line,
                                        line_length, filenamestr, *line_number,
                                        err);
      if (had_err ||
          (!gff3_parser->incomplete_node && queue_size(genome_nodes))) {
        break;
      }
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
    had_err = hashtable_foreach(gff3_parser->undefined_sequence_regions,
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

void gff3parser_reset(GFF3Parser *gff3_parser)
{
  assert(gff3_parser && gff3_parser->id_to_genome_node_mapping);
  hashtable_reset(gff3_parser->id_to_genome_node_mapping);
  hashtable_reset(gff3_parser->seqid_to_ssr_mapping);
  hashtable_reset(gff3_parser->source_to_str_mapping);
  hashtable_reset(gff3_parser->undefined_sequence_regions);
}

void gff3parser_delete(GFF3Parser *gff3_parser)
{
  if (!gff3_parser) return;
  assert(gff3_parser->id_to_genome_node_mapping);
  hashtable_delete(gff3_parser->id_to_genome_node_mapping);
  hashtable_delete(gff3_parser->seqid_to_ssr_mapping);
  hashtable_delete(gff3_parser->source_to_str_mapping);
  hashtable_delete(gff3_parser->undefined_sequence_regions);
  mapping_delete(gff3_parser->offset_mapping);
  ma_free(gff3_parser);
}
