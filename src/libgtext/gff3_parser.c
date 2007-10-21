/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg

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
#include "libgtcore/parseutils.h"
#include "libgtcore/splitter.h"
#include "libgtcore/undef.h"
#include "libgtcore/warning.h"
#include "libgtext/comment.h"
#include "libgtext/genome_feature.h"
#include "libgtext/genome_node.h"
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

static AutomaticSequenceRegion* automatic_sequence_region_new(Env *env)
{
  AutomaticSequenceRegion *auto_sr;
  env_error_check(env);
  auto_sr = env_ma_malloc(env, sizeof (AutomaticSequenceRegion));
  auto_sr->genome_features = array_new(sizeof (GenomeFeature*), env);
  return auto_sr;
}

static void automatic_sequence_region_delete(AutomaticSequenceRegion *auto_sr,
                                             Env *env)
{
  unsigned long i;
  if (!auto_sr) return;
  genome_node_delete(auto_sr->sequence_region, env);
  for (i = 0; i < array_size(auto_sr->genome_features); i++) {
    genome_node_delete(*(GenomeNode**) array_get(auto_sr->genome_features, i),
                       env);
  }
  array_delete(auto_sr->genome_features, env);
  env_ma_free(auto_sr, env);
}

typedef struct {
  Str *seqid_str;
  Range range;
  unsigned long line;
} SimpleSequenceRegion;

static SimpleSequenceRegion* simple_sequence_region_new(const char *seqid,
                                                        Range range,
                                                        unsigned long line,
                                                        Env *env)
{
  SimpleSequenceRegion *ssr = env_ma_malloc(env, sizeof *ssr);
  ssr->seqid_str = str_new_cstr(seqid, env);
  ssr->range = range;
  ssr->line = line;
  return ssr;
}

static void simple_sequence_region_delete(SimpleSequenceRegion *ssr, Env *env)
{
  if (!ssr) return;
  str_delete(ssr->seqid_str, env);
  env_ma_free(ssr, env);
}

GFF3Parser* gff3parser_new(bool checkids, Env *env)
{
  GFF3Parser *gff3_parser = env_ma_malloc(env, sizeof (GFF3Parser));
  gff3_parser->id_to_genome_node_mapping = hashtable_new(HASH_STRING,
                                                         env_ma_free_func,
                                                         (FreeFunc)
                                                         genome_node_delete,
                                                         env);
  gff3_parser->seqid_to_ssr_mapping = hashtable_new(HASH_STRING, NULL,
                                                    (FreeFunc)
                                                  simple_sequence_region_delete,
                                                    env);
  gff3_parser->source_to_str_mapping = hashtable_new(HASH_STRING, NULL,
                                                     (FreeFunc) str_delete,
                                                     env);
  gff3_parser->undefined_sequence_regions = hashtable_new(HASH_STRING, NULL,
                              (FreeFunc) automatic_sequence_region_delete, env);
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
                              Env *env)
{
  env_error_check(env);
  assert(gff3_parser);
  assert(gff3_parser->offset == UNDEF_LONG);
  gff3_parser->offset_mapping = mapping_new(offsetfile, "offsets",
                                            MAPPINGTYPE_INTEGER, env);
  if (gff3_parser->offset_mapping)
    return 0;
  return -1;

}

static int add_offset_if_necessary(Range *range, GFF3Parser *gff3_parser,
                                   const char *seqid, Env *env)
{
  long offset;
  int had_err = 0;
  env_error_check(env);
  if (gff3_parser->offset != UNDEF_LONG)
    *range = range_offset(*range, gff3_parser->offset);
  else if (gff3_parser->offset_mapping) {
    had_err = mapping_map_integer(gff3_parser->offset_mapping, &offset, seqid,
                                  env);
    if (!had_err)
      *range = range_offset(*range, offset);
  }
  return had_err;
}

static int parse_regular_gff3_line(GFF3Parser *gff3_parser, Queue *genome_nodes,
                                   char *line, size_t line_length,
                                   Str *filenamestr, unsigned long line_number,
                                   Env *env)
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
  char *id = NULL, *seqid = NULL, *source = NULL, *type, *start, *end, *score,
    *strand = NULL, *phase = NULL, *attributes = NULL, *token, *tmp_token,
    **tokens;
  const char *filename;
  bool is_child = false, seqid_str_created = false;
  unsigned long i;
  int had_err = 0;

  env_error_check(env);

  filename = str_get(filenamestr);

  /* create splitter */
  splitter = splitter_new(env);
  attribute_splitter = splitter_new(env);
  tmp_splitter = splitter_new(env);
  parents_splitter = splitter_new(env);

  /* parse */
  splitter_split(splitter, line, line_length, '\t', env);
  if (splitter_size(splitter) != 9UL) {
    env_error_set(env, "line %lu in file \"%s\" does not contain 9 tab (\\t) "
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
    env_error_set(env, "type \"%s\" on line %lu in file \"%s\" is not a valid "
                       "one", type, line_number, filename);
    had_err = -1;
  }

  /* parse the range */
  if (!had_err)
    had_err = parse_range(&range, start, end, line_number, filename, env);
  if (!had_err && range.start == 0) {
      env_error_set(env, "illegal feature start 0 on line %lu in file \"%s\" "
                    "(GFF3 files are 1-based)", line_number, filename);
      had_err = -1;
  }
  if (!had_err)
    had_err = add_offset_if_necessary(&range, gff3_parser, seqid, env);

  /* parse the score */
  if (!had_err)
    had_err = parse_score(&score_value, score, line_number, filename, env);

  /* parse the strand */
  if (!had_err)
    had_err = parse_strand(&strand_value, strand, line_number, filename, env);

  /* parse the phase */
  if (!had_err)
    had_err = parse_phase(&phase_value, phase, line_number, filename, env);

  /* create the feature */
  if (!had_err) {
    genome_feature = genome_feature_new(gft, range, strand_value, filenamestr,
                                        line_number, env);
  }

  /* parse the unique id and the parents */
  if (!had_err) {
    splitter_split(attribute_splitter, attributes, strlen(attributes), ';',
                   env);
    for (i = 0; i < splitter_size(attribute_splitter); i++) {
      token = splitter_get_token(attribute_splitter, i);
      if (strncmp(token, ".", 1) == 0) {
        if (splitter_size(attribute_splitter) > 1) {
          env_error_set(env, "more than one attribute token defined on line "
                        "%lu in file \"%s\", altough the first one is '.'",
                        line_number, filename);
          had_err = -1;
        }
        if (!had_err)
          break; /* no attributes to parse */
      }
      else if (strncmp(token, ID_STRING, strlen(ID_STRING)) == 0) {
        if (id) {
          env_error_set(env, "more then one %s token on line %lu in file "
                        "\"%s\"", ID_STRING, line_number, filename);
          had_err = -1;
          break;
        }
        splitter_reset(tmp_splitter);
        splitter_split(tmp_splitter, token, strlen(token), '=', env);
        if (splitter_size(tmp_splitter) != 2) {
          env_error_set(env, "token \"%s\" on line %lu in file \"%s\" does not "
                             "contain exactly one '='", token, line_number,
                        filename);
          had_err = -1;
          break;
        }
        id = cstr_dup(splitter_get_token(tmp_splitter, 1), env);
      }
      else if (strncmp(token, PARENT_STRING, strlen(PARENT_STRING)) == 0) {
        splitter_reset(tmp_splitter);
        splitter_split(tmp_splitter, token, strlen(token), '=', env);
        if (splitter_size(tmp_splitter) != 2) {
          env_error_set(env, "token \"%s\" on line %lu in file \"%s\" does not "
                             "contain exactly one '='", token, line_number,
                        filename);
          had_err = -1;
          break;
        }
        tmp_token = splitter_get_token(tmp_splitter, 1);
        splitter_split(parents_splitter, tmp_token, strlen(tmp_token), ',',
                       env);
        assert(splitter_size(parents_splitter));
      }
      else {
        /* add other attributes here */
        splitter_reset(tmp_splitter);
        splitter_split(tmp_splitter, token, strlen(token), '=', env);
        if (splitter_size(tmp_splitter) != 2) {
          env_error_set(env, "token \"%s\" on line %lu in file \"%s\" does not "
                             "contain exactly one '='", token, line_number,
                        filename);
          had_err = -1;
          break;
        }
      }
      if (!had_err) {
        /* save all attributes, although the ID and Parent attributes are
           generated */
        genome_feature_add_attribute((GenomeFeature*) genome_feature,
                                     splitter_get_token(tmp_splitter, 0),
                                     splitter_get_token(tmp_splitter, 1), env);
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
        auto_sr = automatic_sequence_region_new(env);
        seqid_str = str_new_cstr(seqid, env);
        seqid_str_created = true;
        auto_sr->sequence_region = sequence_region_new(seqid_str, range, NULL,
                                                       0, env);
        hashtable_add(gff3_parser->undefined_sequence_regions,
                      str_get(seqid_str), auto_sr, env);
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
        env_error_set(env, "range (%lu,%lu) of feature on line %lu in file "
                      "\"%s\" is not contained in range (%lu,%lu) of "
                      "corresponding sequence region on line %lu",
                      range.start, range.end, line_number, filename,
                      ssr->range.start, ssr->range.end, ssr->line);
        had_err = -1;
      }
    }
  }
  if (!had_err) {
    assert(seqid_str);
    genome_node_set_seqid(genome_feature, seqid_str, env);
  }
  if (seqid_str_created)
    str_delete(seqid_str, env);

  /* set source */
  if (!had_err) {
    source_str = hashtable_get(gff3_parser->source_to_str_mapping, source);
    if (!source_str) {
      source_str = str_new_cstr(source, env);
      hashtable_add(gff3_parser->source_to_str_mapping, str_get(source_str),
                    source_str, env);
    }
    assert(source_str);
    genome_node_set_source(genome_feature, source_str);
  }

  if (!had_err && score_value != UNDEF_DOUBLE)
    genome_feature_set_score((GenomeFeature*) genome_feature, score_value);
  if (!had_err && phase_value != PHASE_UNDEFINED)
    genome_node_set_phase(genome_feature, phase_value);

  /* store id */
  if (!had_err && id) {
    if ((gn = hashtable_get(gff3_parser->id_to_genome_node_mapping, id))) {
      /* this id has been used already -> raise error */
      env_error_set(env, "the %s \"%s\" on line %lu in file \"%s\" has been "
                    "used already for the feature defined on line %lu",
                    ID_STRING, id, line_number, filename,
                    genome_node_get_line_number(gn));
      had_err = -1;
      env_ma_free(id, env);
    }
    else {
      gff3_parser->incomplete_node = true;
      hashtable_add(gff3_parser->id_to_genome_node_mapping, id,
                    genome_node_ref(genome_feature), env);
    }
  }
  else
    env_ma_free(id, env);

  if (!had_err) {
    for (i = 0; i < splitter_size(parents_splitter); i++) {
      parent_gf = hashtable_get(gff3_parser->id_to_genome_node_mapping,
                                splitter_get_token(parents_splitter, i));
      if (!parent_gf) {
        env_error_set(env, "%s \"%s\" on line %lu in file \"%s\" has not been "
                           "previously defined (via \"%s=\")", PARENT_STRING,
                      splitter_get_token(parents_splitter, i), line_number,
                      filename, ID_STRING);
        had_err = -1;
      }
      else if (str_cmp(genome_node_get_seqid(parent_gf),
                       genome_node_get_seqid(genome_feature))) {
        env_error_set(env, "child on line %lu in file \"%s\" has different "
                      "sequence id than its parent on line %lu ('%s' vs. '%s')",
                      genome_node_get_line_number(genome_feature), filename,
                      genome_node_get_line_number(parent_gf),
                      str_get(genome_node_get_seqid(genome_feature)),
                      str_get(genome_node_get_seqid(parent_gf)));
        had_err = -1;
      }
      else {
        assert(gff3_parser->incomplete_node);
        genome_node_is_part_of_genome_node(parent_gf, genome_feature, env);
        is_child = true;
      }
    }
  }

  if (!had_err) {
    gn = (is_child || auto_sr) ? NULL : genome_feature;
    if (auto_sr && !is_child)
      array_add(auto_sr->genome_features, genome_feature, env);
  }
  else
    genome_node_delete(genome_feature, env);

  if (!had_err && gn)
    queue_add(genome_nodes, gn, env);

  /* free */
  str_delete(changed_seqid, env);
  splitter_delete(splitter, env);
  splitter_delete(attribute_splitter, env);
  splitter_delete(tmp_splitter, env);
  splitter_delete(parents_splitter, env);

  return had_err;
}

static int parse_meta_gff3_line(GFF3Parser *gff3_parser, Queue *genome_nodes,
                                char *line, size_t line_length,
                                Str *filenamestr, unsigned long line_number,
                                Env *env)
{
  char *tmpline, *tmplineend, *seqid = NULL;
  GenomeNode *gn;
  Str *changed_seqid = NULL;
  SimpleSequenceRegion *ssr;
  Range range;
  const char *filename;
  int rval, had_err = 0;

  env_error_check(env);
  assert(line[0] == '#');

  filename = str_get(filenamestr);

  if (line_length == 1 || line[1] != '#') {
    /* storing comment */
    gn = comment_new(line+1, filenamestr, line_number, env);
    queue_add(genome_nodes, gn, env);
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
      env_error_set(env, "missing sequence region name on line %lu in file "
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
        env_error_set(env, "missing sequence region start on line %lu in file "
                           "\"%s\"", line_number, filename);
        had_err = -1;
      }
    }
    if (!had_err) {
      if ((rval = sscanf(tmpline, "%lu", &range.start)) != 1) {
        env_error_set(env, "could not parse region start on line %lu in file "
                           "\"%s\"", line_number, filename);
        had_err = -1;
      }
    }
    if (!had_err  && range.start == 0) {
      env_error_set(env, "illegal region start 0 on line %lu in file \"%s\" "
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
        env_error_set(env, "missing sequence region end on line %lu in file "
                           "\"%s\"", line_number, filename);
        had_err = -1;
      }
    }
    if (!had_err && (rval = sscanf(tmpline, "%lu", &range.end)) != 1) {
      env_error_set(env, "could not parse region end on line %lu in file "
                         "\"%s\"", line_number, filename);
      had_err = -1;
    }
    if (!had_err && (range.start > range.end)) {
      env_error_set(env, "region start %lu is larger then region end %lu on "
                         "line %lu in file \"%s\"", range.start, range.end,
                    line_number, filename);
      had_err = -1;
    }
    if (!had_err)
      had_err = add_offset_if_necessary(&range, gff3_parser, seqid, env);
   if (!had_err) {
      if (hashtable_get(gff3_parser->undefined_sequence_regions, seqid)) {
        env_error_set(env, "genome feature with id \"%s\" has been defined "
                      "before the corresponding \"%s\" definition on line %lu "
                      "in file \"%s\"", seqid, GFF_SEQUENCE_REGION, line_number,
                      filename);
        had_err = -1;
      }
    }
    if (!had_err) {
      /* now we can create a sequence region node */
      assert(seqid);
      ssr = hashtable_get(gff3_parser->seqid_to_ssr_mapping, seqid);
      if (ssr) {
        env_error_set(env, "the sequence region \"%s\" on line %lu in file "
                      "\"%s\" has already been defined",
                      str_get(ssr->seqid_str), line_number, filename);
        had_err = -1;
      }
      else {
        ssr = simple_sequence_region_new(seqid, range, line_number, env);
        hashtable_add(gff3_parser->seqid_to_ssr_mapping,
                      str_get(ssr->seqid_str), ssr, env);
      }
    }
    if (!had_err) {
      assert(ssr);
      gn = sequence_region_new(ssr->seqid_str, range, filenamestr, line_number,
                               env);
      queue_add(genome_nodes, gn, env);
    }
  }
  else if (strcmp(line, GFF_TERMINATOR) == 0) { /* terminator */
    /* now all nodes are complete */
    gff3_parser->incomplete_node = false;
    if (!gff3_parser->checkids)
      hashtable_reset(gff3_parser->id_to_genome_node_mapping, env);
  }
  else {
    warning("skipping unknown meta line %lu in file \"%s\": %s", line_number,
            filename, line);
  }
  str_delete(changed_seqid, env);
  return had_err;
}

static int add_auto_sr_to_queue(void *key, void *value, void *data, Env *env)
{
  AutomaticSequenceRegion *auto_sr = value;
  Queue *genome_nodes = data;
  GenomeNode *gf;
  unsigned int i;
  env_error_check(env);
  assert(key && value && data);
  if (array_size(auto_sr->genome_features)) {
    queue_add(genome_nodes, auto_sr->sequence_region, env);
    auto_sr->sequence_region = NULL;
    for (i = 0; i < array_size(auto_sr->genome_features); i++) {
      gf = *(GenomeNode**) array_get(auto_sr->genome_features, i);
      queue_add(genome_nodes, gf, env);
    }
    array_reset(auto_sr->genome_features);
  }
  return 0;
}

int gff3parser_parse_genome_nodes(int *status_code, GFF3Parser *gff3_parser,
                                  Queue *genome_nodes, Str *filenamestr,
                                  unsigned long long *line_number,
                                  GenFile *fpin, Env *env)
{
  size_t line_length;
  Str *line_buffer;
  char *line;
  const char *filename;
  int rval, version, had_err = 0;

  env_error_check(env);

  filename = str_get(filenamestr);

  /* init */
  line_buffer = str_new(env);

  while ((rval = str_read_next_line_generic(line_buffer, fpin, env)) != EOF) {
    line = str_get(line_buffer);
    line_length = str_length(line_buffer);
    (*line_number)++;

    if (*line_number == 1) {
      if (strncmp(line, GFF_VERSION_PREFIX, strlen(GFF_VERSION_PREFIX))) {
        env_error_set(env,
                      "line %llu in file \"%s\" does not begin with \"%s\"",
                      *line_number, filename, GFF_VERSION_PREFIX);
        had_err = -1;
        break;
      }
      if (!had_err) {
        line += strlen(GFF_VERSION_PREFIX);
        /* skip blanks */
        while (line[0] == ' ')
          line++;
        had_err = parse_int(&version, line, *line_number, filename, env);
      }
      if (!had_err) {
        if (version != GFF_VERSION) {
          env_error_set(env, "GFF version %d does not equal required version "
                             "%u ", version, GFF_VERSION);
          had_err = -1;
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
                                     env);
      if (had_err ||
          (!gff3_parser->incomplete_node && queue_size(genome_nodes))) {
        break;
      }
    }
    else {
      had_err = parse_regular_gff3_line(gff3_parser, genome_nodes, line,
                                        line_length, filenamestr, *line_number,
                                        env);
      if (had_err ||
          (!gff3_parser->incomplete_node && queue_size(genome_nodes))) {
        break;
      }
    }
    str_reset(line_buffer);
  }

  if (had_err) {
    while (queue_size(genome_nodes))
      genome_node_rec_delete(queue_get(genome_nodes, env), env);
  }
  else if (rval == EOF) {
    /* the file has been parsed completely, add automatically created sequence
       regions to queue */
    had_err = hashtable_foreach(gff3_parser->undefined_sequence_regions,
                                add_auto_sr_to_queue, genome_nodes, env);
    assert(!had_err); /* add_auto_sr_to_queue() is sane */
  }

  str_delete(line_buffer, env);
  if (queue_size(genome_nodes))
    *status_code = 0; /* at least one node was created */
  else
    *status_code = EOF;
  return had_err;
}

void gff3parser_reset(GFF3Parser *gff3_parser, Env *env)
{
  assert(gff3_parser && gff3_parser->id_to_genome_node_mapping);
  hashtable_reset(gff3_parser->id_to_genome_node_mapping, env);
  hashtable_reset(gff3_parser->seqid_to_ssr_mapping, env);
  hashtable_reset(gff3_parser->source_to_str_mapping, env);
  hashtable_reset(gff3_parser->undefined_sequence_regions, env);
}

void gff3parser_delete(GFF3Parser *gff3_parser, Env *env)
{
  if (!gff3_parser) return;
  assert(gff3_parser->id_to_genome_node_mapping);
  hashtable_delete(gff3_parser->id_to_genome_node_mapping, env);
  hashtable_delete(gff3_parser->seqid_to_ssr_mapping, env);
  hashtable_delete(gff3_parser->source_to_str_mapping, env);
  hashtable_delete(gff3_parser->undefined_sequence_regions, env);
  mapping_delete(gff3_parser->offset_mapping, env);
  env_ma_free(gff3_parser, env);
}
