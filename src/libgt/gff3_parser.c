/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "array.h"
#include "comment.h"
#include "cstr.h"
#include "error.h"
#include "fptr.h"
#include "genome_feature.h"
#include "genome_node.h"
#include "gff3_parser.h"
#include "hashtable.h"
#include "parseutils.h"
#include "phase.h"
#include "range.h"
#include "sequence_region.h"
#include "str.h"
#include "undef.h"
#include "splitter.h"
#include "warning.h"
#include "xansi.h"

struct GFF3Parser {
  Hashtable *id_to_genome_node_mapping,
            *seqid_to_str_mapping,
            *source_to_str_mapping;
  long offset;
};

GFF3Parser* gff3parser_new(Env *env)
{
  GFF3Parser *gff3_parser = env_ma_malloc(env, sizeof (GFF3Parser));
  gff3_parser->id_to_genome_node_mapping = hashtable_new(HASH_STRING,
                                                         env_ma_free_func, NULL,
                                                         env);
  gff3_parser->seqid_to_str_mapping = hashtable_new(HASH_STRING, NULL,
                                                    (FreeFunc) str_delete, env);
  gff3_parser->source_to_str_mapping = hashtable_new(HASH_STRING, NULL,
                                                     (FreeFunc) str_delete,
                                                     env);
  gff3_parser->offset = UNDEFLONG;
  return gff3_parser;
}

void gff3parser_set_offset(GFF3Parser *gff3_parser, long offset)
{
  assert(gff3_parser);
  gff3_parser->offset = offset;
}

static int parse_regular_gff3_line(GFF3Parser *gff3_parser,
                                   Queue *genome_nodes, char *line,
                                   size_t line_length, const char *filename,
                                   unsigned long line_number,
                                   bool *break_loop, Env *env)
{
  GenomeNode *gn = NULL, *genome_feature = NULL, *parent_gf;
  GenomeFeatureType gft;
  Splitter *splitter, *attribute_splitter, *tmp_splitter, *parents_splitter;
  Str *seqid_str, *source_str;
  Strand strand_value;
  double score_value;
  Phase phase_value;
  Range range;
  char *id = NULL, *seqid, *source, *type, *start, *end, *score, *strand,
       *phase, *attributes, *token, *tmp_token, **tokens;
  bool out_node_complete = true, is_child = false;
  unsigned long i;
  int has_err = 0;

  env_error_check(env);

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
    has_err = -1;
  }
  if (!has_err) {
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
  if (!has_err && genome_feature_type_get(&gft, type) == -1) {
    env_error_set(env, "type \"%s\" on line %lu in file \"%s\" is not a valid "
                  "one", type, line_number, filename);
    has_err = -1;
  }

  /* parse the range */
  if (!has_err)
    has_err = parse_range(&range, start, end, line_number, filename, env);

  if (!has_err && gff3_parser->offset != UNDEFLONG)
    range = range_offset(range, gff3_parser->offset, line_number);

  /* parse the score */
  if (!has_err)
    has_err = parse_score(&score_value, score, line_number, filename, env);

  /* parse the strand */
  if (!has_err)
    has_err = parse_strand(&strand_value, strand, line_number, filename, env);

  /* parse the phase */
  if (!has_err)
    has_err = parse_phase(&phase_value, phase, line_number, filename, env);

  /* create the feature */
  if (!has_err) {
    genome_feature = genome_feature_new(gft, range, strand_value, filename,
                                        line_number, env);
  }

  /* parse the unique id and the parents */
  if (!has_err) {
    splitter_split(attribute_splitter, attributes, strlen(attributes), ';',
                   env);
    for (i = 0; i < splitter_size(attribute_splitter); i++) {
      token = splitter_get_token(attribute_splitter, i);
      if (strncmp(token, ".", 1) == 0) {
        if (splitter_size(attribute_splitter) > 1) {
          env_error_set(env, "more than one attribute token defined on line "
                        "%lu in file \"%s\", altough the first one is '.'",
                        line_number, filename);
          has_err = -1;
        }
        if (!has_err)
          break; /* no attributes to parse */
      }
      else if (strncmp(token, ID_STRING, strlen(ID_STRING)) == 0) {
        if (id) {
          env_error_set(env, "more then one %s token on line %lu in file "
                        "\"%s\"", ID_STRING, line_number, filename);
          has_err = -1;
          break;
        }
        splitter_reset(tmp_splitter);
        splitter_split(tmp_splitter, token, strlen(token), '=', env);
        if (splitter_size(tmp_splitter) != 2) {
          env_error_set(env, "token \"%s\" on line %lu in file \"%s\" does not "
                    "contain exactly one '='", token, line_number, filename);
          has_err = -1;
          break;
        }
        id = cstr_dup(splitter_get_token(tmp_splitter, 1), env);
      }
      else if (strncmp(token, PARENT_STRING, strlen(PARENT_STRING)) == 0) {
        splitter_reset(tmp_splitter);
        splitter_split(tmp_splitter, token, strlen(token), '=', env);
        if (splitter_size(tmp_splitter) != 2) {
          env_error_set(env, "token \"%s\" on line %lu in file \"%s\" does not "
                    "contain exactly one '='", token, line_number, filename);
          has_err = -1;
          break;
        }
        tmp_token = splitter_get_token(tmp_splitter, 1);
        splitter_split(parents_splitter, tmp_token, strlen(tmp_token), ',',
                       env);
        assert(splitter_size(parents_splitter)); /* XXX: should be an error */
      }
      else {
        /* add other attributes here */
        splitter_reset(tmp_splitter);
        splitter_split(tmp_splitter, token, strlen(token), '=', env);
        if (splitter_size(tmp_splitter) != 2) {
          env_error_set(env, "token \"%s\" on line %lu in file \"%s\" does not "
                    "contain exactly one '='", token, line_number, filename);
          has_err = -1;
          break;
        }
        genome_feature_add_attribute((GenomeFeature*) genome_feature,
                                     splitter_get_token(tmp_splitter, 0),
                                     splitter_get_token(tmp_splitter, 1), env);
      }
    }
  }

  /* set seqid */
  if (!has_err) {
    seqid_str = hashtable_get(gff3_parser->seqid_to_str_mapping, seqid);
    if (!seqid_str) {
      env_error_set(env, "seqid \"%s\" on line %lu in file \"%s\" has not been "
                "previously introduced with a \"%s\" line", seqid, line_number,
                filename, GFF_SEQUENCE_REGION);
      has_err = -1;
    }
    else
      genome_node_set_seqid(genome_feature, seqid_str);
  }

  /* set source */
  if (!has_err) {
    source_str = hashtable_get(gff3_parser->source_to_str_mapping, source);
    if (!source_str) {
      source_str = str_new_cstr(source, env);
      hashtable_add(gff3_parser->source_to_str_mapping, str_get(source_str),
                    source_str, env);
    }
    assert(source_str);
    genome_node_set_source(genome_feature, source_str);
  }

  if (!has_err && score_value != UNDEFDOUBLE)
    genome_feature_set_score((GenomeFeature*) genome_feature, score_value);
  if (!has_err && phase_value != PHASE_UNDEFINED)
    genome_node_set_phase(genome_feature, phase_value);

  /* store id */
  if (!has_err && id) {
    if ((gn = hashtable_get(gff3_parser->id_to_genome_node_mapping, id))) {
      /* this id has been used already -> raise error */
      env_error_set(env, "the %s \"%s\" on line %lu in file \"%s\" has been "
                    "used already for the feature defined on line %lu",
                    ID_STRING, id, line_number, filename,
                    genome_node_get_line_number(gn));
      has_err = -1;
      env_ma_free(id, env);
    }
    else {
      out_node_complete = false;
      hashtable_add(gff3_parser->id_to_genome_node_mapping, id, genome_feature,
                    env);
    }
  }
  else
    env_ma_free(id, env);

  if (!has_err) {
    for (i = 0; i < splitter_size(parents_splitter); i++) {
      parent_gf = hashtable_get(gff3_parser->id_to_genome_node_mapping,
                                splitter_get_token(parents_splitter, i));
      if (!parent_gf) {
        env_error_set(env, "%s \"%s\" on line %lu in file \"%s\" has not been "
                  "previously defined (via \"%s=\")", PARENT_STRING,
                  splitter_get_token(parents_splitter, i), line_number,
                  filename, ID_STRING);
        has_err = -1;
      }
      else {
        genome_node_is_part_of_genome_node(parent_gf, genome_feature, env);
        out_node_complete = false;
        is_child = true;
      }
    }
  }

  if (!has_err)
    gn = is_child ? NULL : genome_feature;
  else
    genome_node_delete(genome_feature, env);

  if (!has_err && gn)
    queue_add(genome_nodes, gn, env);
  if (out_node_complete)
    *break_loop = true;

  /* free */
  splitter_delete(splitter, env);
  splitter_delete(attribute_splitter, env);
  splitter_delete(tmp_splitter, env);
  splitter_delete(parents_splitter, env);

  return has_err;
}

static int parse_meta_gff3_line(GFF3Parser *gff3_parser, Queue *genome_nodes,
                                char *line, size_t line_length,
                                const char *filename,
                                unsigned long line_number,
                                bool *break_loop, Env *env)
{
  char *tmpline, *tmplineend, *seqid = NULL;
  GenomeNode *gn;
  Str *seqid_str;
  Range range;
  int rval, has_err = 0;

  env_error_check(env);
  assert(line[0] == '#');

  if (line_length == 1 || line[1] != '#') {
    /* storing comment */
    gn = comment_new(line+1, filename, line_number, env);
    queue_add(genome_nodes, gn, env);
    *break_loop = true;
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
      env_error_set(env, "missing sequence region on line %lu in file \"%s\"",
                line_number, filename);
      has_err = -1;
    }
    if (!has_err) {
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
        has_err = -1;
      }
    }
    if (!has_err) {
      if ((rval = sscanf(tmpline, "%lu", &range.start)) != 1) {
        env_error_set(env, "could not parse region start on line %lu in file "
                  "\"%s\"", line_number, filename);
        has_err = -1;
      }
    }
    if (!has_err) {
      /* skip non-blanks */
      while (tmpline <= tmplineend && !(tmpline[0] == ' '))
        tmpline++;
      /* skip blanks */
      while (tmpline < tmplineend && tmpline[0] == ' ')
        tmpline++;
      if (tmpline > tmplineend) {
        env_error_set(env, "missing sequence region end on line %lu in file "
                      "\"%s\"", line_number, filename);
        has_err = -1;
      }
    }
    if (!has_err && (rval = sscanf(tmpline, "%lu", &range.end)) != 1) {
      env_error_set(env, "could not parse region end on line %lu in file "
                    "\"%s\"", line_number, filename);
      has_err = -1;
    }
    if (!has_err && (range.start > range.end)) {
      env_error_set(env, "region start %lu is larger then region end %lu on "
                    "line %lu in file \"%s\"", range.start, range.end,
                    line_number, filename);
      has_err = -1;
    }
    if (!has_err) {
      if (gff3_parser->offset != UNDEFLONG)
        range = range_offset(range, gff3_parser->offset, line_number);
      /* now we can create a sequence region node */
      assert(seqid);
      seqid_str = hashtable_get(gff3_parser->seqid_to_str_mapping, seqid);
      if (seqid_str) {
        env_error_set(env, "the sequence region \"%s\" on line %lu in file "
                  "\"%s\" is already defined", str_get(seqid_str), line_number,
                  filename);
        has_err = -1;
      }
      else {
        seqid_str = str_new_cstr(seqid, env);
        hashtable_add(gff3_parser->seqid_to_str_mapping, str_get(seqid_str),
                      seqid_str, env);
      }
    }
    if (!has_err) {
      assert(seqid_str);
      gn = sequence_region_new(seqid_str, range, filename, line_number, env);
      queue_add(genome_nodes, gn, env);
    }
    *break_loop = true;
  }
  else if (strcmp(line, GFF_TERMINATOR) == 0) { /* terminator */
    *break_loop = true;
  }
  else {
    warning("skipping unknown meta line %lu in file \"%s\": %s", line_number,
            filename, line);
  }
  return has_err;
}

int gff3parser_parse_genome_nodes(int *status_code, GFF3Parser *gff3_parser,
                                  Queue *genome_nodes,
                                  const char *filename,
                                  unsigned long *line_number,
                                  FILE *fpin, Env *env)
{
  bool break_loop = false;
  size_t line_length;
  Str *line_buffer;
  char *line;
  int version, has_err = 0;

  env_error_check(env);

  /* init */
  line_buffer = str_new(env);

  /* the given (buffer) queue is empty */
  assert(!queue_size(genome_nodes));

  while (str_read_next_line(line_buffer, fpin, env) != EOF) {
    line = str_get(line_buffer);
    line_length = str_length(line_buffer);
    (*line_number)++;

    if (*line_number == 1) {
      if (strncmp(line, GFF_VERSION_PREFIX, strlen(GFF_VERSION_PREFIX))) {
        env_error_set(env, "line %lu in file \"%s\" does begin with \"%s\"",
                  *line_number, filename, GFF_VERSION_PREFIX);
        has_err = -1;
        break;
      }
      if (!has_err) {
        line += strlen(GFF_VERSION_PREFIX);
        /* skip blanks */
        while (line[0] == ' ')
          line++;
        has_err = parse_int(&version, line, *line_number, filename, env);
      }
      if (!has_err) {
        if (version != GFF_VERSION) {
          env_error_set(env, "GFF version %d does not equal required version "
                        "%u ", version, GFF_VERSION);
          has_err = -1;
        }
      }
    }
    else if (line_length == 0)
      warning("skipping blank line %lu in file \"%s\"", *line_number, filename);
    else if (line[0] == '#') {
      has_err = parse_meta_gff3_line(gff3_parser, genome_nodes, line,
                                     line_length, filename, *line_number,
                                     &break_loop, env);
      if (has_err || break_loop)
        break;
    }
    else {
      has_err = parse_regular_gff3_line(gff3_parser, genome_nodes, line,
                                        line_length, filename, *line_number,
                                        &break_loop, env);
      if (has_err || break_loop)
        break;
    }
    str_reset(line_buffer);
  }

  if (has_err) {
    while (queue_size(genome_nodes))
      genome_node_rec_delete(*(GenomeNode**) queue_get(genome_nodes), env);  
  }

  str_delete(line_buffer, env);
  if (queue_size(genome_nodes))
    *status_code = 0; /* at least one node was created */
  else
    *status_code = EOF;
  return has_err;
}

void gff3parser_reset(GFF3Parser *gff3_parser, Env *env)
{
  assert(gff3_parser && gff3_parser->id_to_genome_node_mapping);
  hashtable_reset(gff3_parser->id_to_genome_node_mapping, env);
  hashtable_reset(gff3_parser->seqid_to_str_mapping, env);
  hashtable_reset(gff3_parser->source_to_str_mapping, env);
}

void gff3parser_delete(GFF3Parser *gff3_parser, Env *env)
{
  if (!gff3_parser) return;
  assert(gff3_parser->id_to_genome_node_mapping);
  hashtable_delete(gff3_parser->id_to_genome_node_mapping, env);
  hashtable_delete(gff3_parser->seqid_to_str_mapping, env);
  hashtable_delete(gff3_parser->source_to_str_mapping, env);
  env_ma_free(gff3_parser, env);
}
