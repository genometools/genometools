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
#include <string.h>
#include "libgtcore/cstr.h"
#include "libgtcore/hashtable.h"
#include "libgtcore/ma.h"
#include "libgtcore/parseutils.h"
#include "libgtcore/splitter.h"
#include "libgtcore/undef.h"
#include "libgtcore/unused.h"
#include "libgtcore/warning.h"
#include "libgtext/compare.h"
#include "libgtext/genome_node.h"
#include "libgtext/gtf_parser.h"

#define GENE_ID_ATTRIBUTE       "gene_id"
#define TRANSCRIPT_ID_ATTRIBUTE "transcript_id"

struct GTF_parser {
  Hashtable *sequence_region_to_range, /* map from sequence regions to ranges */
            *gene_id_hash, /* map from gene_id to transcript_id hash */
            *seqid_to_str_mapping,
            *source_to_str_mapping;
};

typedef enum {
  GTF_CDS,
  GTF_exon,
  GTF_stop_codon
} GTF_feature_type;

static const char *GTF_feature_type_strings[] = { "CDS",
                                                  "exon",
                                                  "stop_codon" };

static int GTF_feature_type_get(GTF_feature_type *type, char *feature_string)
{
  void *result;

  assert(type && feature_string);

  result = bsearch(&feature_string,
                   GTF_feature_type_strings,
                   sizeof (GTF_feature_type_strings) /
                   sizeof (GTF_feature_type_strings[0]),
                   sizeof (char*),
                   compare);

  if (result) {
    *type = (GTF_feature_type)
            ((char**) result - (char**) GTF_feature_type_strings);
    return 0;
  }
  /* else type not found */
  return -1;
}

GTF_parser* gtf_parser_new(void)
{
  GTF_parser *parser = ma_malloc(sizeof (GTF_parser));
  parser->sequence_region_to_range = hashtable_new(HASH_STRING,
                                                   ma_free_func,
                                                   ma_free_func);
  parser->gene_id_hash = hashtable_new(HASH_STRING, ma_free_func,
                                       (FreeFunc) hashtable_delete);
  parser->seqid_to_str_mapping = hashtable_new(HASH_STRING, NULL,
                                               (FreeFunc) str_delete);
  parser->source_to_str_mapping = hashtable_new(HASH_STRING, NULL,
                                                (FreeFunc) str_delete);
  return parser;
}

static int construct_sequence_regions(void *key, void *value, void *data,
                                      UNUSED Error *err)
{
  Str *seqid;
  Range range;
  GenomeNode *gn;
  Queue *genome_nodes = (Queue*) data;
  error_check(err);
  assert(key && value && data);
  seqid = str_new_cstr(key);
  range = *(Range*) value;
  gn = sequence_region_new(seqid, range, NULL, 0);
  queue_add(genome_nodes, gn);
  str_delete(seqid);
  return 0;
}

static int construct_mRNAs(UNUSED void *key, void *value, void *data,
                           Error *err)
{
  Array *genome_node_array = (Array*) value,
        *mRNAs = (Array*) data;
  GenomeNode *mRNA_node, *first_node, *gn;
  Strand mRNA_strand;
  Range mRNA_range;
  Str *mRNA_seqid;
  unsigned long i;
  int had_err = 0;

  error_check(err);
  assert(key && value && data);
  assert(array_size(genome_node_array)); /* at least one node in array */

  /* determine the range and the strand of the mRNA */
  first_node = *(GenomeNode**) array_get(genome_node_array, 0);
  mRNA_range = genome_node_get_range(first_node);
  mRNA_strand = genome_feature_get_strand((GenomeFeature*) first_node);
  mRNA_seqid = genome_node_get_seqid(first_node);
  for (i = 1; i < array_size(genome_node_array); i++) {
    gn = *(GenomeNode**) array_get(genome_node_array, i);
    mRNA_range = range_join(mRNA_range, genome_node_get_range(gn));
    /* XXX: an error check is necessary here, otherwise strand_join() can cause
       a failed assertion */
    mRNA_strand = strand_join(mRNA_strand,
                              genome_feature_get_strand((GenomeFeature*) gn));
    if (str_cmp(mRNA_seqid, genome_node_get_seqid(gn))) {
      error_set(err, "The features on lines %lu and %lu refer to different "
                "genomic sequences (``seqname''), although they have the same "
                "gene IDs (``gene_id'') which must be globally unique",
                genome_node_get_line_number(first_node),
                genome_node_get_line_number(gn));
      had_err = -1;
      break;
    }
  }

  if (!had_err) {
    mRNA_node = genome_feature_new(gft_mRNA, mRNA_range, mRNA_strand, NULL, 0);
    genome_node_set_seqid(mRNA_node, mRNA_seqid);

    /* register children */
    for (i = 0; i < array_size(genome_node_array); i++) {
      gn = *(GenomeNode**) array_get(genome_node_array, i);
      genome_node_is_part_of_genome_node(mRNA_node, gn);
    }

    /* store the mRNA */
    array_add(mRNAs, mRNA_node);
  }

  return had_err;
}

static int construct_genes(UNUSED void *key, void *value, void *data,
                           Error *err)
{
  Hashtable *transcript_id_hash = (Hashtable*) value;
  Queue *genome_nodes = (Queue*) data;
  Array *mRNAs = array_new(sizeof (GenomeNode*));
  GenomeNode *gene_node, *gn;
  Strand gene_strand;
  Range gene_range;
  Str *gene_seqid;
  unsigned long i;
  int had_err = 0;

  error_check(err);
  assert(key && value && data);

  had_err = hashtable_foreach(transcript_id_hash, construct_mRNAs, mRNAs, err);
  if (!had_err) {
    assert(array_size(mRNAs)); /* at least one mRNA constructed */

    /* determine the range and the strand of the gene */
    gn = *(GenomeNode**) array_get(mRNAs, 0);
    gene_range = genome_node_get_range(gn);
    gene_strand = genome_feature_get_strand((GenomeFeature*) gn);
    gene_seqid = genome_node_get_seqid(gn);
    for (i = 1; i < array_size(mRNAs); i++) {
      gn = *(GenomeNode**) array_get(mRNAs, i);
      gene_range = range_join(gene_range, genome_node_get_range(gn));
      gene_strand = strand_join(gene_strand,
                                genome_feature_get_strand((GenomeFeature*) gn));
      assert(str_cmp(gene_seqid, genome_node_get_seqid(gn)) == 0);
    }

    gene_node = genome_feature_new(gft_gene, gene_range, gene_strand, NULL, 0);
    genome_node_set_seqid(gene_node, gene_seqid);

    /* register children */
    for (i = 0; i < array_size(mRNAs); i++) {
      gn = *(GenomeNode**) array_get(mRNAs, i);
      genome_node_is_part_of_genome_node(gene_node, gn);
    }

    /* store the gene */
    queue_add(genome_nodes, gene_node);

    /* free */
    array_delete(mRNAs);
  }

  return had_err;
}

int gtf_parser_parse(GTF_parser *parser, Queue *genome_nodes,
                     Str *filenamestr, FILE *fpin, unsigned int be_tolerant,
                     Error *err)
{
  Str *seqid_str, *source_str, *line_buffer;
  char *line;
  size_t line_length;
  unsigned long i, line_number = 0;
  GenomeNode *gn;
  Range range, *rangeptr;
  Phase phase_value;
  Strand strand_value;
  Splitter *splitter, *attribute_splitter;
  double score_value;
  char *seqname,
       *source,
       *feature,
       *start,
       *end,
       *score,
       *strand,
       *frame,
       *attributes,
       *token,
       *gene_id,
       *transcript_id,
       **tokens;
  Hashtable *transcript_id_hash; /* map from transcript id to array of genome
                                    nodes */
  Array *genome_node_array;
  GTF_feature_type gtf_feature_type;
  /* abuse gft_TF_binding_site as an undefined value */
  GenomeFeatureType gff_feature_type = gft_TF_binding_site;
  const char *filename;
  int had_err = 0;

  assert(parser && genome_nodes && fpin);
  error_check(err);

  filename = str_get(filenamestr);

  /* alloc */
  line_buffer = str_new();
  splitter = splitter_new(),
  attribute_splitter = splitter_new();

#define HANDLE_ERROR                                                \
        if (had_err) {                                              \
          if (be_tolerant) {                                        \
            fprintf(stderr, "skipping line: %s\n", error_get(err)); \
            error_unset(err);                                       \
            str_reset(line_buffer);                                 \
            had_err = 0;                                            \
            continue;                                               \
          }                                                         \
          else {                                                    \
            had_err = -1;                                           \
            break;                                                  \
          }                                                         \
        }

  while (str_read_next_line(line_buffer, fpin) != EOF) {
    line = str_get(line_buffer);
    line_length = str_length(line_buffer);
    line_number++;
    had_err = 0;

    if (line_length == 0)
      warning("skipping blank line %lu in file \"%s\"", line_number, filename);
    else if (line[0] == '#') {
      /* storing comment */
      gn = comment_new(line+1, filenamestr, line_number);
      queue_add(genome_nodes, gn);
    }
    else {
      /* process tab delimited GTF line */
      splitter_reset(splitter);
      splitter_split(splitter, line, line_length, '\t');
      if (splitter_size(splitter) != 9UL) {
        error_set(err, "line %lu in file \"%s\" contains %lu tab (\\t) "
                  "separated fields instead of 9", line_number, filename,
                  splitter_size(splitter));
        had_err = -1;
        break;
      }
      tokens = splitter_get_tokens(splitter);
      seqname    = tokens[0];
      source     = tokens[1];
      feature    = tokens[2];
      start      = tokens[3];
      end        = tokens[4];
      score      = tokens[5];
      strand     = tokens[6];
      frame      = tokens[7];
      attributes = tokens[8];

      /* parse feature */
      if (GTF_feature_type_get(&gtf_feature_type, feature) == -1) {
        /* we skip unknown features */
        fprintf(stderr, "skipping line %lu in file \"%s\": unknown feature: "
                        "\"%s\"\n", line_number, filename, feature);
        str_reset(line_buffer);
        continue;
      }

      /* translate into GFF3 feature type */
      switch (gtf_feature_type) {
        case GTF_CDS:
        case GTF_stop_codon:
          gff_feature_type = gft_CDS;
          break;
        case GTF_exon:
          gff_feature_type = gft_exon;
      }
      assert(gff_feature_type != gft_TF_binding_site);

      /* parse the range */
      had_err = parse_range(&range, start, end, line_number, filename, err);
      HANDLE_ERROR;

      /* process seqname (we have to do it here because we need the range) */
      if ((rangeptr = hashtable_get(parser->sequence_region_to_range,
                                    seqname))) {
        /* sequence region is already defined -> update range */
        *rangeptr = range_join(range, *rangeptr);
      }
      else {
        /* sequence region is not already defined -> define it */
        rangeptr = ma_malloc(sizeof (Range));
        *rangeptr = range;
        hashtable_add(parser->sequence_region_to_range, cstr_dup(seqname),
                      rangeptr);
      }

      /* parse the score */
      had_err = parse_score(&score_value, score, line_number, filename, err);
      HANDLE_ERROR;

      /* parse the strand */
      had_err = parse_strand(&strand_value, strand, line_number, filename, err);
      HANDLE_ERROR;

      /* parse the frame */
      had_err = parse_phase(&phase_value, frame, line_number, filename, err);
      HANDLE_ERROR;

      /* parse the attributes */
      splitter_reset(attribute_splitter);
      gene_id = NULL;
      transcript_id = NULL;
      splitter_split(attribute_splitter, attributes, strlen(attributes), ';');
      for (i = 0; i < splitter_size(attribute_splitter); i++) {
        token = splitter_get_token(attribute_splitter, i);
        /* skip blank newline, if necessary and possible */
        if (i && token[0])
          token++;
        /* look for the two mandatory attributes */
        if (strncmp(token, GENE_ID_ATTRIBUTE, strlen(GENE_ID_ATTRIBUTE)) == 0) {
          if (strlen(token) + 2 < strlen(GENE_ID_ATTRIBUTE)) {
            error_set(err, "missing value to attribute \"%s\" on line %lu in "
                      "file \"%s\"", GENE_ID_ATTRIBUTE, line_number, filename);
            had_err = -1;
          }
          HANDLE_ERROR;
          gene_id = token + strlen(GENE_ID_ATTRIBUTE) + 1;
        }
        else if (strncmp(token, TRANSCRIPT_ID_ATTRIBUTE,
                         strlen(TRANSCRIPT_ID_ATTRIBUTE)) == 0) {
          if (strlen(token) + 2 < strlen(TRANSCRIPT_ID_ATTRIBUTE)) {
            error_set(err, "missing value to attribute \"%s\" on line %lu in "
                      "file \"%s\"", TRANSCRIPT_ID_ATTRIBUTE, line_number,
                      filename);
            had_err = -1;
          }
          HANDLE_ERROR;
          transcript_id = token + strlen(TRANSCRIPT_ID_ATTRIBUTE) + 1;
        }
      }

      /* check for the manadatory attributes */
      if (!gene_id) {
        error_set(err, "missing attribute \"%s\" on line %lu in file \"%s\"",
                  GENE_ID_ATTRIBUTE, line_number, filename);
        had_err = -1;
      }
      HANDLE_ERROR;
      if (!transcript_id) {
        error_set(err, "missing attribute \"%s\" on line %lu in file \"%s\"",
                  TRANSCRIPT_ID_ATTRIBUTE, line_number, filename);
        had_err = -1;
      }
      HANDLE_ERROR;

      /* process the mandatory attributes */
      if (!(transcript_id_hash = hashtable_get(parser->gene_id_hash,
                                               gene_id))) {
        transcript_id_hash = hashtable_new(HASH_STRING, ma_free_func,
                                           (FreeFunc) array_delete);
        hashtable_add(parser->gene_id_hash, cstr_dup(gene_id),
                      transcript_id_hash);
      }
      assert(transcript_id_hash);

      if (!(genome_node_array = hashtable_get(transcript_id_hash,
                                              transcript_id))) {
        genome_node_array = array_new(sizeof (GenomeNode*));
        hashtable_add(transcript_id_hash, cstr_dup(transcript_id),
                      genome_node_array);
      }
      assert(genome_node_array);

      /* construct the new feature */
      gn = genome_feature_new(gff_feature_type, range, strand_value,
                              filenamestr, line_number);

      /* set seqid */
      seqid_str = hashtable_get(parser->seqid_to_str_mapping, seqname);
      if (!seqid_str) {
        seqid_str = str_new_cstr(seqname);
        hashtable_add(parser->seqid_to_str_mapping, str_get(seqid_str),
                      seqid_str);
      }
      assert(seqid_str);
      genome_node_set_seqid(gn, seqid_str);

      /* set source */
      source_str = hashtable_get(parser->source_to_str_mapping, source);
      if (!source_str) {
        source_str = str_new_cstr(source);
        hashtable_add(parser->source_to_str_mapping, str_get(source_str),
                      source_str);
      }
      assert(source_str);
      genome_feature_set_source(gn, source_str);

      if (score_value != UNDEF_DOUBLE)
        genome_feature_set_score((GenomeFeature*) gn, score_value);
      if (phase_value != PHASE_UNDEFINED)
        genome_feature_set_phase(gn, phase_value);
      array_add(genome_node_array, gn);
    }

    str_reset(line_buffer);
  }

  /* process all comments features */
  if (!had_err) {
    had_err = hashtable_foreach(parser->sequence_region_to_range,
                                construct_sequence_regions, genome_nodes, NULL);
    assert(!had_err); /* construct_sequence_regions() is sane */
  }

  /* process all genome_features */
  if (!had_err) {
    had_err = hashtable_foreach(parser->gene_id_hash, construct_genes,
                                genome_nodes, err);
  }

  /* free */
  splitter_delete(splitter);
  splitter_delete(attribute_splitter);
  str_delete(line_buffer);

  return had_err;
}

void gtf_parser_delete(GTF_parser *parser)
{
  if (!parser) return;
  hashtable_delete(parser->sequence_region_to_range);
  hashtable_delete(parser->gene_id_hash);
  hashtable_delete(parser->seqid_to_str_mapping);
  hashtable_delete(parser->source_to_str_mapping);
  ma_free(parser);
}
