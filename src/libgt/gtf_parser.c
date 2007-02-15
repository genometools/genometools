/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include "error.h"
#include "fptr.h"
#include "genome_node.h"
#include "gtf_parser.h"
#include "hashtable.h"
#include "parseutils.h"
#include "splitter.h"
#include "str.h"
#include "undef.h"
#include "warning.h"
#include "xansi.h"

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

static int compare (const void *a, const void *b)
{
  return strcmp(*((char**) a), *((char**) b));
}

static int GTF_feature_type_get(GTF_feature_type *type, char *feature_string)
{
  void *result;

  assert(type && feature_string);

  result = bsearch(&feature_string,
                   GTF_feature_type_strings,
                   sizeof(GTF_feature_type_strings) /
                   sizeof(GTF_feature_type_strings[0]),
                   sizeof(char*),
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
  GTF_parser *parser = xmalloc(sizeof(GTF_parser));

  parser->sequence_region_to_range = hashtable_new(HASH_STRING, free, free);
  parser->gene_id_hash = hashtable_new(HASH_STRING, free,
                                       (Free) hashtable_free);
  parser->seqid_to_str_mapping = hashtable_new(HASH_STRING, NULL,
                                               (Free) str_free);
  parser->source_to_str_mapping = hashtable_new(HASH_STRING, NULL,
                                                (Free) str_free);

  return parser;
}

static void construct_sequence_regions(void *key, void *value, void *data)
{
  Str *seqid;
  Range range;
  GenomeNode *gn;
  Queue *genome_nodes = (Queue*) data;

  assert(key && value && data);

  seqid = str_new_cstr(key);
  range = *(Range*) value;
  gn = sequence_region_new(seqid, range, NULL, 0);
  queue_add(genome_nodes, gn);

  str_free(seqid);
}

static void construct_mRNAs(void *key, void *value, void *data)
{
  Array *genome_node_array = (Array*) value,
        *mRNAs = (Array*) data;
  GenomeNode *mRNA_node, *first_node, *gn;
  Strand mRNA_strand;
  Range mRNA_range;
  Str *mRNA_seqid;
  unsigned long i;

  assert(key && value && data);
  assert(array_size(genome_node_array)); /* at least one node in array */

  /* determine the range and the strand of the mRNA */
  first_node = *(GenomeNode**) array_get(genome_node_array, 0);
  mRNA_range = genome_node_get_range(first_node);
  mRNA_strand = genome_feature_get_strand((Genome_feature*) first_node);
  mRNA_seqid = genome_node_get_seqid(first_node);
  for (i = 1; i < array_size(genome_node_array); i++) {
    gn = *(GenomeNode**) array_get(genome_node_array, i);
    mRNA_range = range_join(mRNA_range, genome_node_get_range(gn));
    /* XXX: an error check is necessary here, otherwise strand_join() can cause
       a failed assertion */
    mRNA_strand = strand_join(mRNA_strand,
                              genome_feature_get_strand((Genome_feature*) gn));
    if (str_cmp(mRNA_seqid, genome_node_get_seqid(gn))) {
      error("The features on lines %lu and %lu refer to different genomic "
            "sequences (``seqname''), although they have the same gene IDs "
            "(``gene_id'') which must be globally unique",
            genome_node_get_line_number(first_node),
            genome_node_get_line_number(gn));
    }
  }

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

static void construct_genes(void *key, void *value, void *data)
{
  Hashtable *transcript_id_hash = (Hashtable*) value;
  Queue *genome_nodes = (Queue*) data;
  Array *mRNAs = array_new(sizeof(GenomeNode*));
  GenomeNode *gene_node, *gn;
  Strand gene_strand;
  Range gene_range;
  Str *gene_seqid;
  unsigned long i;

  assert(key && value && data);

  hashtable_foreach(transcript_id_hash, construct_mRNAs, mRNAs);
  assert(array_size(mRNAs)); /* at least one mRNA constructed */

  /* determine the range and the strand of the gene */
  gn = *(GenomeNode**) array_get(mRNAs, 0);
  gene_range = genome_node_get_range(gn);
  gene_strand = genome_feature_get_strand((Genome_feature*) gn);
  gene_seqid = genome_node_get_seqid(gn);
  for (i = 1; i < array_size(mRNAs); i++) {
    gn = *(GenomeNode**) array_get(mRNAs, i);
    gene_range = range_join(gene_range, genome_node_get_range(gn));
    gene_strand = strand_join(gene_strand,
                              genome_feature_get_strand((Genome_feature*) gn));
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
  array_free(mRNAs);
}

int gtf_parser_parse(GTF_parser *parser, Queue *genome_nodes,
                     const char *filename, FILE *fpin,
                     unsigned int be_tolerant, Error *err)
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

  assert(parser && genome_nodes && filename && fpin);
  error_check(err);

  /* alloc */
  line_buffer = str_new();
  splitter = splitter_new(),
  attribute_splitter = splitter_new();

#define HANDLE_ERROR                                                \
        if (error_is_set(err)) {                                    \
          if (be_tolerant) {                                        \
            fprintf(stderr, "skipping line: %s\n", error_get(err)); \
            error_unset(err);                                       \
            str_reset(line_buffer);                                 \
            continue;                                               \
          }                                                         \
          else {                                                    \
            splitter_free(splitter);                                \
            splitter_free(attribute_splitter);                      \
            str_free(line_buffer);                                  \
            return -1;                                              \
          }                                                         \
        }

  while (str_read_next_line(line_buffer, fpin) != EOF) {
    line = str_get(line_buffer);
    line_length = str_length(line_buffer);
    line_number++;

    if (line_length == 0)
      warning("skipping blank line %lu in file \"%s\"", line_number, filename);
    else if (line[0] == '#') {
      /* storing comment */
      gn = comment_new(line+1, filename, line_number);
      queue_add(genome_nodes, gn);
    }
    else {
      /* process tab delimited GTF line */
      splitter_reset(splitter);
      splitter_split(splitter, line, line_length, '\t');
      if (splitter_size(splitter) != 9UL) {
        error("line %lu in file \"%s\" contains %lu tab (\\t) separated fields "
              "instead of 9", line_number, filename,
              splitter_size(splitter));
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
      /* XXX: process return value differently */
      (void) parse_range(&range, start, end, line_number, filename, err);
      HANDLE_ERROR;

      /* process seqname (we have to do it here because we need the range) */
      if ((rangeptr = hashtable_get(parser->sequence_region_to_range,
                                    seqname))) {
        /* sequence region is already defined -> update range */
        *rangeptr = range_join(range, *rangeptr);
      }
      else {
        /* sequence region is not already defined -> define it */
        rangeptr = xmalloc(sizeof(Range));
        *rangeptr = range;
        hashtable_add(parser->sequence_region_to_range, xstrdup(seqname),
                      rangeptr);
      }

      /* parse the score */
      /* XXX: process return value differently */
      (void) parse_score(&score_value, score, line_number, filename, err);
      HANDLE_ERROR;

      /* parse the strand */
      strand_value = parse_strand(strand, line_number, filename, err);
      HANDLE_ERROR;

      /* parse the frame */
      phase_value = parse_phase(frame, line_number, filename, err);
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
          }
          HANDLE_ERROR;
          transcript_id = token + strlen(TRANSCRIPT_ID_ATTRIBUTE) + 1;
        }
      }

      /* check for the manadatory attributes */
      if (!gene_id) {
        error_set(err, "missing attribute \"%s\" on line %lu in file \"%s\"",
                  GENE_ID_ATTRIBUTE, line_number, filename);
      }
      HANDLE_ERROR;
      if (!transcript_id) {
        error_set(err, "missing attribute \"%s\" on line %lu in file \"%s\"",
                  TRANSCRIPT_ID_ATTRIBUTE, line_number, filename);
      }
      HANDLE_ERROR;

      /* process the mandatory attributes */
      if (!(transcript_id_hash = hashtable_get(parser->gene_id_hash,
                                               gene_id))) {
        transcript_id_hash = hashtable_new(HASH_STRING, free,
                                           (Free) array_free);
        hashtable_add(parser->gene_id_hash, xstrdup(gene_id),
                      transcript_id_hash);
      }
      assert(transcript_id_hash);

      if (!(genome_node_array = hashtable_get(transcript_id_hash,
                                              transcript_id))) {
        genome_node_array = array_new(sizeof(GenomeNode*));
        hashtable_add(transcript_id_hash, xstrdup(transcript_id),
                      genome_node_array);
      }
      assert(genome_node_array);

      /* construct the new feature */
      gn = genome_feature_new(gff_feature_type, range, strand_value, filename,
                              line_number);

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
      genome_node_set_source(gn, source_str);

      if (score_value != UNDEFDOUBLE)
        genome_feature_set_score((Genome_feature*) gn, score_value);
      if (phase_value != PHASE_UNDEFINED)
        genome_node_set_phase(gn, phase_value);
      array_add(genome_node_array, gn);
    }

    str_reset(line_buffer);
  }

  /* process all comments features */
  hashtable_foreach(parser->sequence_region_to_range,
                    construct_sequence_regions, genome_nodes);

  /* process all genome_features */
  hashtable_foreach(parser->gene_id_hash, construct_genes, genome_nodes);

  /* free */
  splitter_free(splitter);
  splitter_free(attribute_splitter);
  str_free(line_buffer);

  return 0;
}

void gtf_parser_free(GTF_parser *parser)
{
  if (!parser) return;
  hashtable_free(parser->sequence_region_to_range);
  hashtable_free(parser->gene_id_hash);
  hashtable_free(parser->seqid_to_str_mapping);
  hashtable_free(parser->source_to_str_mapping);
  free(parser);
}
