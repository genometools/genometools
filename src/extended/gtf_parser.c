/*
  Copyright (c) 2006-2008, 2012 Gordon Gremme <gordon@gremme.org>
  Copyright (c) 2006-2008       Center for Bioinformatics, University of Hamburg

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
#include "core/assert_api.h"
#include "core/cstr_api.h"
#include "core/hashmap.h"
#include "core/ma.h"
#include "core/parseutils.h"
#include "core/splitter.h"
#include "core/strand.h"
#include "core/strcmp.h"
#include "core/symbol.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "core/warning_api.h"
#include "extended/feature_node.h"
#include "extended/genome_node.h"
#include "extended/gff3_defines.h"
#include "extended/gtf_parser.h"
#include "extended/region_node_builder.h"

#define GENE_ID_ATTRIBUTE          "gene_id"
#define GENE_NAME_ATTRIBUTE        "gene_name"
#define TRANSCRIPT_ID_ATTRIBUTE    "transcript_id"
#define TRANSCRIPT_NAME_ATTRIBUTE  "transcript_name"

#define GTF_PARSER_STOP_CODON_FLAG "stop_codon"

struct GtGTFParser {
  GtHashmap *gene_id_hash, /* map from gene_id to transcript_id hash */
            *seqid_to_str_mapping,
            *source_to_str_mapping,
            *gene_id_to_name_mapping,
            *transcript_id_to_name_mapping;
  GtRegionNodeBuilder *region_node_builder;
  GtTypeChecker *type_checker;
};

typedef struct {
  GtQueue *genome_nodes;
  GtArray *mRNAs;
  GtHashmap *gene_id_to_name_mapping,
            *transcript_id_to_name_mapping;
  bool tidy;
} ConstructionInfo;

typedef enum {
  GTF_CDS,
  GTF_exon,
  GTF_start_codon,
  GTF_stop_codon
} GTF_feature_type;

static const char *GTF_feature_type_strings[] = { "CDS",
                                                  "exon",
                                                  "start_codon",
                                                  "stop_codon" };

static int GTF_feature_type_get(GTF_feature_type *type, char *feature_string)
{
  void *result;

  gt_assert(type && feature_string);

  result = bsearch(&feature_string,
                   GTF_feature_type_strings,
                   sizeof (GTF_feature_type_strings) /
                   sizeof (GTF_feature_type_strings[0]),
                   sizeof (char*),
                   gt_strcmpptr);

  if (result) {
    *type = (GTF_feature_type)
            ((char**) result - (char**) GTF_feature_type_strings);
    return 0;
  }
  /* else type not found */
  return -1;
}

GtGTFParser* gt_gtf_parser_new(GtTypeChecker *type_checker)
{
  GtGTFParser *parser = gt_malloc(sizeof (GtGTFParser));
  parser->gene_id_hash = gt_hashmap_new(GT_HASH_STRING, gt_free_func,
                                        (GtFree) gt_hashmap_delete);
  parser->seqid_to_str_mapping = gt_hashmap_new(GT_HASH_STRING, NULL,
                                                (GtFree) gt_str_delete);
  parser->source_to_str_mapping = gt_hashmap_new(GT_HASH_STRING, NULL,
                                                 (GtFree) gt_str_delete);
  parser->gene_id_to_name_mapping = gt_hashmap_new(GT_HASH_STRING, gt_free_func,
                                                   gt_free_func);
  parser->transcript_id_to_name_mapping = gt_hashmap_new(GT_HASH_STRING,
                                                         gt_free_func,
                                                         gt_free_func);
  parser->region_node_builder = gt_region_node_builder_new();
  parser->type_checker = type_checker;
  return parser;
}

static int construct_mRNAs(GT_UNUSED void *key, void *value, void *data,
                           GtError *err)
{
  ConstructionInfo *cinfo = (ConstructionInfo*) data;
  GtArray *gt_genome_node_array = (GtArray*) value,
          *mRNAs = (GtArray*) cinfo->mRNAs;
  GtGenomeNode *mRNA_node, *first_node, *gn;
  const char *tname;
  GtStrand mRNA_strand;
  GtRange mRNA_range;
  GtStr *mRNA_seqid;
  GtUword i;
  int had_err = 0;

  gt_error_check(err);
  gt_assert(key && value && data);
   /* at least one node in array */
  gt_assert(gt_array_size(gt_genome_node_array));

  /* determine the range and the strand of the mRNA */
  first_node = *(GtGenomeNode**) gt_array_get(gt_genome_node_array, 0);
  mRNA_range = gt_genome_node_get_range(first_node);
  mRNA_strand = gt_feature_node_get_strand((GtFeatureNode*) first_node);
  mRNA_seqid = gt_genome_node_get_seqid(first_node);

  /* TODO: support discontinuous start/stop codons */
  for (i = 0; !had_err && i < gt_array_size(gt_genome_node_array); i++) {
    gn = *(GtGenomeNode**) gt_array_get(gt_genome_node_array, i);
    if (gt_feature_node_get_attribute((GtFeatureNode*) gn,
        GTF_PARSER_STOP_CODON_FLAG)) {
      GtUword j;
      GtRange stop_codon_rng = gt_genome_node_get_range(gn);
      bool found_cds = false;
      for (j = 0; !had_err && j < gt_array_size(gt_genome_node_array); j++) {
        GtGenomeNode* gn2;
        GtRange this_rng;
        const char *this_type;
        gn2 = *(GtGenomeNode**) gt_array_get(gt_genome_node_array, j);
        if (gn == gn2) continue;
        this_rng = gt_genome_node_get_range(gn2);
        this_type = gt_feature_node_get_type((GtFeatureNode*) gn2);
        if (this_type == gt_symbol(gt_ft_CDS)) {
          if (gt_range_contains(&this_rng, &stop_codon_rng)) {
            if (cinfo->tidy) {
              gt_warning("stop codon on line %u in file %s is contained in "
                         "CDS in line %u",
                         gt_genome_node_get_line_number(gn),
                         gt_genome_node_get_filename(gn),
                         gt_genome_node_get_line_number(gn2));
              found_cds = true;
            } else {
              gt_error_set(err, "stop codon on line %u in file %s is "
                                "contained in CDS in line %u",
                           gt_genome_node_get_line_number(gn),
                           gt_genome_node_get_filename(gn),
                           gt_genome_node_get_line_number(gn2));
              had_err = -1;
            }
            break;
          }
          if (this_rng.end + 1 == stop_codon_rng.start) {
            this_rng.end = stop_codon_rng.end;
            gt_genome_node_set_range(gn2, &this_rng);
            found_cds = true;
            break;
          }
          if (this_rng.start == stop_codon_rng.end + 1) {
            this_rng.start = stop_codon_rng.start;
            gt_genome_node_set_range(gn2, &this_rng);
            found_cds = true;
            break;
          }
        }
      }
      if (!found_cds) {
        if (!had_err) {
          if (cinfo->tidy) {
            gt_warning("found stop codon on line %u in file %s with no "
                       "flanking CDS, ignoring it",
                       gt_genome_node_get_line_number(gn),
                       gt_genome_node_get_filename(gn));
          } else {
            gt_error_set(err, "found stop codon on line %u in file %s with no "
                              "flanking CDS",
                         gt_genome_node_get_line_number(gn),
                         gt_genome_node_get_filename(gn));
            had_err = -1;
            break;
          }
        }
      } else {
        gt_array_rem(gt_genome_node_array, i);
        gt_genome_node_delete(gn);
      }
    }
  }

  for (i = 1; !had_err && i < gt_array_size(gt_genome_node_array); i++) {
    GtRange range;
    gn = *(GtGenomeNode**) gt_array_get(gt_genome_node_array, i);
    range = gt_genome_node_get_range(gn);
    mRNA_range = gt_range_join(&mRNA_range, &range);
    /* XXX: an error check is necessary here, otherwise gt_strand_join() can
       cause a failed assertion */
    mRNA_strand = gt_strand_join(mRNA_strand,
                          gt_feature_node_get_strand((GtFeatureNode*) gn));
    if (gt_str_cmp(mRNA_seqid, gt_genome_node_get_seqid(gn))) {
      gt_error_set(err, "The features on lines %u and %u refer to different "
                "genomic sequences (``seqname''), although they have the same "
                "gene IDs (``gene_id'') which must be globally unique",
                gt_genome_node_get_line_number(first_node),
                gt_genome_node_get_line_number(gn));
      had_err = -1;
      break;
    }
  }

  if (!had_err) {
    mRNA_node = gt_feature_node_new(mRNA_seqid, gt_ft_mRNA, mRNA_range.start,
                                    mRNA_range.end, mRNA_strand);

    if ((tname = gt_hashmap_get(cinfo->transcript_id_to_name_mapping,
                              (const char*) key)) && strlen(tname) > 0) {
      gt_feature_node_add_attribute((GtFeatureNode*) mRNA_node, GT_GFF_NAME,
                                      tname);
    }

    /* register children */
    for (i = 0; i < gt_array_size(gt_genome_node_array); i++) {
      gn = *(GtGenomeNode**) gt_array_get(gt_genome_node_array, i);
      gt_feature_node_add_child((GtFeatureNode*) mRNA_node,
                                (GtFeatureNode*) gt_genome_node_ref(gn));
    }

    /* store the mRNA */
    gt_array_add(mRNAs, mRNA_node);
  }

  return had_err;
}

static int delete_mRNAs(GT_UNUSED void *key, void *value,
                        GT_UNUSED void *data, GT_UNUSED GtError *err)
{
  GtArray *gt_genome_node_array = (GtArray*) value;
  GtUword i;

  gt_assert(key && value);
  for (i = 0; i < gt_array_size(gt_genome_node_array); i++) {
    GtGenomeNode *gn = *(GtGenomeNode**) gt_array_get(gt_genome_node_array, i);
    gt_genome_node_delete(gn);
  }

  return 0;
}

static int construct_genes(GT_UNUSED void *key, void *value, void *data,
                           GtError *err)
{
  GtHashmap *transcript_id_hash = (GtHashmap*) value;
  ConstructionInfo *cinfo = (ConstructionInfo*) data;
  GtQueue *genome_nodes = cinfo->genome_nodes;
  const char *gname;
  GtArray *mRNAs = gt_array_new(sizeof (GtGenomeNode*));
  GtGenomeNode *gene_node, *gn;
  GtStrand gene_strand;
  GtRange gene_range;
  GtStr *gene_seqid;
  GtUword i;
  int had_err = 0;

  gt_error_check(err);
  gt_assert(key && value && data);
  cinfo->mRNAs = mRNAs;
  had_err = gt_hashmap_foreach(transcript_id_hash, construct_mRNAs, cinfo, err);
  if (!had_err) {
    gt_assert(gt_array_size(mRNAs)); /* at least one mRNA constructed */

    /* determine the range and the strand of the gene */
    gn = *(GtGenomeNode**) gt_array_get(mRNAs, 0);
    gene_range = gt_genome_node_get_range(gn);
    gene_strand = gt_feature_node_get_strand((GtFeatureNode*) gn);
    gene_seqid = gt_genome_node_get_seqid(gn);
    for (i = 1; i < gt_array_size(mRNAs); i++) {
      GtRange range;
      gn = *(GtGenomeNode**) gt_array_get(mRNAs, i);
      range = gt_genome_node_get_range(gn);
      gene_range = gt_range_join(&gene_range, &range);
      gene_strand = gt_strand_join(gene_strand,
                          gt_feature_node_get_strand((GtFeatureNode*) gn));
      gt_assert(gt_str_cmp(gene_seqid, gt_genome_node_get_seqid(gn)) == 0);
    }

    gene_node = gt_feature_node_new(gene_seqid, gt_ft_gene, gene_range.start,
                                    gene_range.end, gene_strand);

    if ((gname = gt_hashmap_get(cinfo->gene_id_to_name_mapping,
                              (const char*) key)) && strlen(gname) > 0) {
      gt_feature_node_add_attribute((GtFeatureNode*) gene_node, GT_GFF_NAME,
                                      gname);
    }

    /* register children */
    for (i = 0; i < gt_array_size(mRNAs); i++) {
      gn = *(GtGenomeNode**) gt_array_get(mRNAs, i);
      gt_feature_node_add_child((GtFeatureNode*) gene_node,
                                (GtFeatureNode*) gn);
    }

    /* store the gene */
    gt_queue_add(genome_nodes, gene_node);
  }

  /* free */
  gt_array_delete(mRNAs);

  return had_err;
}

static int delete_genes(GT_UNUSED void *key, void *value, GT_UNUSED void *data,
                        GT_UNUSED GtError *err)
{
  GtHashmap *transcript_id_hash = (GtHashmap*) value;
  GT_UNUSED int had_err = 0;
  gt_assert(key && value);
  had_err = gt_hashmap_foreach(transcript_id_hash, delete_mRNAs, NULL, err);
  gt_assert(!had_err);
  return had_err;
}

int gt_gtf_parser_parse(GtGTFParser *parser, GtQueue *genome_nodes,
                        GtStr *filenamestr, GtFile *fpin, bool be_tolerant,
                        GtError *err)
{
  GtStr *seqid_str, *source_str, *line_buffer;
  char *line;
  size_t line_length;
  GtUword i, line_number = 0;
  GtGenomeNode *gn;
  GtRange range;
  GtPhase phase_value;
  GtStrand gt_strand_value;
  GtSplitter *splitter, *attribute_splitter;
  float score_value;
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
       *gene_name = NULL,
       *transcript_id,
       *transcript_name = NULL,
       **tokens;
  GtHashmap *transcript_id_hash; /* map from transcript id to array of genome
                                    nodes */
  GtArray *gt_genome_node_array;
  ConstructionInfo cinfo;
  GTF_feature_type gtf_feature_type;
  GT_UNUSED bool gff_type_is_valid = false;
  const char *type = NULL;
  const char *filename;
  bool score_is_defined;
  int had_err = 0;

  gt_assert(parser && genome_nodes);
  gt_error_check(err);

  filename = gt_str_get(filenamestr);

  /* alloc */
  line_buffer = gt_str_new();
  splitter = gt_splitter_new(),
  attribute_splitter = gt_splitter_new();

#define HANDLE_ERROR                                                \
        if (had_err) {                                              \
          if (be_tolerant) {                                        \
            fprintf(stderr, "skipping line: %s\n", gt_error_get(err)); \
            gt_error_unset(err);                                       \
            gt_str_reset(line_buffer);                                 \
            had_err = 0;                                            \
            continue;                                               \
          }                                                         \
          else {                                                    \
            had_err = -1;                                           \
            break;                                                  \
          }                                                         \
        }

  while (gt_str_read_next_line_generic(line_buffer, fpin) != EOF) {
    line = gt_str_get(line_buffer);
    line_length = gt_str_length(line_buffer);
    line_number++;
    gene_name = gene_id = transcript_id = transcript_name = NULL;
    had_err = 0;

    if (line_length == 0) {
      gt_warning("skipping blank line " GT_WU " in file \"%s\"", line_number,
                 filename);
    }
    else if (line[0] == '#') {
      /* storing comment */
      if (line_length >= 2 && line[1] == '#')
        gn = gt_comment_node_new(line+2); /* store '##' line as '#' line */
      else
        gn = gt_comment_node_new(line+1);
      gt_genome_node_set_origin(gn, filenamestr, line_number);
      gt_queue_add(genome_nodes, gn);
    }
    else {
      bool stop_codon = false;

      /* process tab delimited GTF line */
      gt_splitter_reset(splitter);
      gt_splitter_split(splitter, line, line_length, '\t');
      if (gt_splitter_size(splitter) != 9UL) {
        gt_error_set(err, "line " GT_WU " in file \"%s\" contains " GT_WU
                     " tab (\\t) " "separated fields instead of 9", line_number,
                     filename,
                  gt_splitter_size(splitter));
        had_err = -1;
        break;
      }
      tokens = gt_splitter_get_tokens(splitter);
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
        fprintf(stderr, "skipping line " GT_WU " in file \"%s\": unknown "
                "feature: \"%s\"\n", line_number, filename, feature);
        gt_str_reset(line_buffer);
        continue;
      }

      /* translate into GFF3 feature type */
      switch (gtf_feature_type) {
        case GTF_stop_codon:
          stop_codon = true;
        case GTF_CDS:
          gff_type_is_valid = gt_type_checker_is_valid(parser->type_checker,
                                                       gt_ft_CDS);
          type = gt_ft_CDS;
          break;
        case GTF_exon:
          gff_type_is_valid = gt_type_checker_is_valid(parser->type_checker,
                                                       gt_ft_exon);
          type = gt_ft_exon;
          break;
        case GTF_start_codon:
          /* we can skip the start codons, they are part of the CDS anyway */
          gt_str_reset(line_buffer);
          continue;
      }
      gt_assert(gff_type_is_valid);

      /* parse the range */
      had_err = gt_parse_range(&range, start, end, line_number, filename, err);
      HANDLE_ERROR;

      /* process seqname (we have to do it here because we need the range) */
      gt_region_node_builder_add_region(parser->region_node_builder, seqname,
                                        range);

      /* parse the score */
      had_err = gt_parse_score(&score_is_defined, &score_value, score,
                               line_number, filename, err);
      HANDLE_ERROR;

      /* parse the strand */
      had_err = gt_parse_strand(&gt_strand_value, strand, line_number, filename,
                               err);
      HANDLE_ERROR;

      /* parse the frame */
      had_err = gt_parse_phase(&phase_value, frame, line_number, filename, err);
      HANDLE_ERROR;

      /* parse the attributes */
      gt_splitter_reset(attribute_splitter);
      gene_id = NULL;
      transcript_id = NULL;
      gt_splitter_split(attribute_splitter, attributes, strlen(attributes),
                        ';');
      for (i = 0; i < gt_splitter_size(attribute_splitter); i++) {
        token = gt_splitter_get_token(attribute_splitter, i);
        /* skip leading blanks */
        while (*token == ' ')
          token++;
        /* look for the two mandatory attributes */
        if (strncmp(token, GENE_ID_ATTRIBUTE, strlen(GENE_ID_ATTRIBUTE)) == 0) {
          if (strlen(token) + 2 < strlen(GENE_ID_ATTRIBUTE)) {
            gt_error_set(err, "missing value to attribute \"%s\" on line "
                         GT_WU "in file \"%s\"", GENE_ID_ATTRIBUTE, line_number,
                         filename);
            had_err = -1;
          }
          HANDLE_ERROR;
          gene_id = token + strlen(GENE_ID_ATTRIBUTE) + 1;
        }
        else if (strncmp(token, TRANSCRIPT_ID_ATTRIBUTE,
                         strlen(TRANSCRIPT_ID_ATTRIBUTE)) == 0) {
          if (strlen(token) + 2 < strlen(TRANSCRIPT_ID_ATTRIBUTE)) {
            gt_error_set(err, "missing value to attribute \"%s\" on line "
                         GT_WU "in file \"%s\"", TRANSCRIPT_ID_ATTRIBUTE,
                         line_number, filename);
            had_err = -1;
          }
          HANDLE_ERROR;
          transcript_id = token + strlen(TRANSCRIPT_ID_ATTRIBUTE) + 1;
        }
        else if (strncmp(token, GENE_NAME_ATTRIBUTE,
                         strlen(GENE_NAME_ATTRIBUTE)) == 0) {
          if (strlen(token) + 2 < strlen(GENE_NAME_ATTRIBUTE)) {
            gt_error_set(err, "missing value to attribute \"%s\" on line "
                         GT_WU "in file \"%s\"", GENE_NAME_ATTRIBUTE,
                         line_number, filename);
            had_err = -1;
          }
          HANDLE_ERROR;
          gene_name = token + strlen(GENE_NAME_ATTRIBUTE) + 1;
          /* for output we want to strip quotes */
          if (*gene_name == '"')
            gene_name++;
          if (gene_name[strlen(gene_name)-1] == '"')
            gene_name[strlen(gene_name)-1] = '\0';
        }
        else if (strncmp(token, TRANSCRIPT_NAME_ATTRIBUTE,
                         strlen(TRANSCRIPT_NAME_ATTRIBUTE)) == 0) {
          if (strlen(token) + 2 < strlen(TRANSCRIPT_NAME_ATTRIBUTE)) {
            gt_error_set(err, "missing value to attribute \"%s\" on line "
                         GT_WU "in file \"%s\"", TRANSCRIPT_NAME_ATTRIBUTE,
                         line_number, filename);
            had_err = -1;
          }
          HANDLE_ERROR;
          transcript_name = token + strlen(TRANSCRIPT_NAME_ATTRIBUTE) + 1;
          /* for output we want to strip quotes */
          if (*transcript_name == '"')
            transcript_name++;
          if (transcript_name[strlen(transcript_name)-1] == '"')
            transcript_name[strlen(transcript_name)-1] = '\0';
        }
      }

      /* check for the mandatory attributes */
      if (!gene_id) {
        gt_error_set(err, "missing attribute \"%s\" on line " GT_WU
                     " in file \"%s\"", GENE_ID_ATTRIBUTE, line_number,
                     filename);
        had_err = -1;
      }
      HANDLE_ERROR;
      if (!transcript_id) {
        gt_error_set(err, "missing attribute \"%s\" on line " GT_WU
                     " in file \"%s\"", TRANSCRIPT_ID_ATTRIBUTE, line_number,
                     filename);
        had_err = -1;
      }
      HANDLE_ERROR;

      /* process the mandatory attributes */
      if (!(transcript_id_hash = gt_hashmap_get(parser->gene_id_hash,
                                             gene_id))) {
        transcript_id_hash = gt_hashmap_new(GT_HASH_STRING, gt_free_func,
                                            (GtFree) gt_array_delete);
        gt_hashmap_add(parser->gene_id_hash, gt_cstr_dup(gene_id),
                    transcript_id_hash);
      }
      gt_assert(transcript_id_hash);

      if (!(gt_genome_node_array = gt_hashmap_get(transcript_id_hash,
                                            transcript_id))) {
        gt_genome_node_array = gt_array_new(sizeof (GtGenomeNode*));
        gt_hashmap_add(transcript_id_hash, gt_cstr_dup(transcript_id),
                    gt_genome_node_array);
      }
      gt_assert(gt_genome_node_array);

      /* save optional gene_name and transcript_name attributes */
      if (transcript_name && strlen(transcript_name) > 0
            && !gt_hashmap_get(parser->transcript_id_to_name_mapping,
                             transcript_id)) {
        gt_hashmap_add(parser->transcript_id_to_name_mapping,
                    gt_cstr_dup(transcript_id),
                    gt_cstr_dup(transcript_name));
      }
      if (gene_name && strlen(gene_name) > 0
            && !gt_hashmap_get(parser->gene_id_to_name_mapping,
                                    gene_id)) {
        gt_hashmap_add(parser->gene_id_to_name_mapping,
                    gt_cstr_dup(gene_id),
                    gt_cstr_dup(gene_name));
      }

      /* get seqid */
      seqid_str = gt_hashmap_get(parser->seqid_to_str_mapping, seqname);
      if (!seqid_str) {
        seqid_str = gt_str_new_cstr(seqname);
        gt_hashmap_add(parser->seqid_to_str_mapping, gt_str_get(seqid_str),
                       seqid_str);
      }
      gt_assert(seqid_str);

      /* construct the new feature */
      gn = gt_feature_node_new(seqid_str, type, range.start, range.end,
                                 gt_strand_value);
      gt_genome_node_set_origin(gn, filenamestr, line_number);
      if (stop_codon) {
        gt_feature_node_add_attribute((GtFeatureNode*) gn,
                                      GTF_PARSER_STOP_CODON_FLAG, "true");
      }

      /* set source */
      source_str = gt_hashmap_get(parser->source_to_str_mapping, source);
      if (!source_str) {
        source_str = gt_str_new_cstr(source);
        gt_hashmap_add(parser->source_to_str_mapping, gt_str_get(source_str),
                    source_str);
      }
      gt_assert(source_str);
      gt_feature_node_set_source((GtFeatureNode*) gn, source_str);

      if (score_is_defined)
        gt_feature_node_set_score((GtFeatureNode*) gn, score_value);
      if (phase_value != GT_PHASE_UNDEFINED)
        gt_feature_node_set_phase((GtFeatureNode*) gn, phase_value);
      gt_array_add(gt_genome_node_array, gn);
    }

    gt_str_reset(line_buffer);
  }

  /* process all region nodes */
  if (!had_err)
    gt_region_node_builder_build(parser->region_node_builder, genome_nodes);

  /* process all feature nodes */
  cinfo.genome_nodes = genome_nodes;
  cinfo.tidy = be_tolerant;
  cinfo.gene_id_to_name_mapping = parser->gene_id_to_name_mapping;
  cinfo.transcript_id_to_name_mapping = parser->transcript_id_to_name_mapping;
  if (!had_err) {
    had_err = gt_hashmap_foreach(parser->gene_id_hash, construct_genes,
                                 &cinfo, err);
  }
  gt_hashmap_foreach(parser->gene_id_hash, delete_genes, NULL, err);

  /* free */
  gt_splitter_delete(splitter);
  gt_splitter_delete(attribute_splitter);
  gt_str_delete(line_buffer);

  return had_err;
}

void gt_gtf_parser_delete(GtGTFParser *parser)
{
  if (!parser) return;
  gt_region_node_builder_delete(parser->region_node_builder);
  gt_hashmap_delete(parser->gene_id_hash);
  gt_hashmap_delete(parser->seqid_to_str_mapping);
  gt_hashmap_delete(parser->source_to_str_mapping);
  gt_hashmap_delete(parser->transcript_id_to_name_mapping);
  gt_hashmap_delete(parser->gene_id_to_name_mapping);
  gt_free(parser);
}
