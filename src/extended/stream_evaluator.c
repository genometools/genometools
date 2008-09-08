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
#include "core/bsearch.h"
#include "core/cstr.h"
#include "core/hashmap.h"
#include "core/log.h"
#include "core/ma.h"
#include "core/unused_api.h"
#include "core/warning.h"
#include "core/xansi.h"
#include "extended/evaluator.h"
#include "extended/gff3_output.h"
#include "extended/stream_evaluator.h"
#include "extended/transcript_evaluators.h"
#include "extended/transcript_exons.h"
#include "extended/transcript_used_exons.h"

typedef struct {
  unsigned long TP, FP, FN;
} NucEval;

struct StreamEvaluator {
  GenomeStream *reality,
               *prediction;
  bool nuceval, evalLTR;
  unsigned long LTRdelta;
  Hashmap *slots; /* sequence id -> slot */
  Evaluator *gene_evaluator,
            *mRNA_evaluator,
            *LTR_evaluator;
  TranscriptEvaluators *mRNA_exon_evaluators,
                       *mRNA_exon_evaluators_collapsed,
                       *CDS_exon_evaluators,
                       *CDS_exon_evaluators_collapsed;
  unsigned long missing_genes,
                wrong_genes,
                missing_mRNAs,
                wrong_mRNAs,
                missing_LTRs,
                wrong_LTRs;
  NucEval mRNA_nucleotides,
          CDS_nucleotides;
};

typedef struct {
  GT_Array *genes_forward,
        *genes_reverse,
        *mRNAs_forward,
        *mRNAs_reverse,
        *LTRs;
  TranscriptExons *mRNA_exons_forward,
                  *mRNA_exons_reverse,
                  *CDS_exons_forward,
                  *CDS_exons_reverse;
  TranscriptCounts *mRNA_counts_forward,
                   *mRNA_counts_reverse,
                   *CDS_counts_forward,
                   *CDS_counts_reverse;
  GT_Range real_range;
  unsigned long FP_mRNA_nucleotides_forward,
                FP_mRNA_nucleotides_reverse,
                FP_CDS_nucleotides_forward,
                FP_CDS_nucleotides_reverse;
  GT_Bittab *real_mRNA_nucleotides_forward,
         *pred_mRNA_nucleotides_forward,
         *real_mRNA_nucleotides_reverse,
         *pred_mRNA_nucleotides_reverse,
         *real_CDS_nucleotides_forward,
         *pred_CDS_nucleotides_forward,
         *real_CDS_nucleotides_reverse,
         *pred_CDS_nucleotides_reverse,
         *true_genes_forward,
         *true_genes_reverse,
         *true_mRNAs_forward,
         *true_mRNAs_reverse,
         *true_LTRs,
         *overlapped_genes_forward,
         *overlapped_genes_reverse,
         *overlapped_mRNAs_forward,
         *overlapped_mRNAs_reverse,
         *overlapped_LTRs;
  TranscriptGT_Bittabs *mRNA_exon_bittabs_forward,
                    *mRNA_exon_bittabs_reverse,
                    *CDS_exon_bittabs_forward,
                    *CDS_exon_bittabs_reverse;
  TranscriptUsedExons *used_mRNA_exons_forward,
                      *used_mRNA_exons_reverse,
                      *used_CDS_exons_forward,
                      *used_CDS_exons_reverse;
} Slot;

typedef struct
{
  Slot *slot;
  bool nuceval,
       verbose;
} ProcessRealFeatureInfo;

typedef struct {
  Slot *slot;
  bool nuceval,
       verbose,
       exondiff;
  unsigned long LTRdelta;
  Evaluator *gene_evaluator,
            *mRNA_evaluator,
            *LTR_evaluator;
  TranscriptEvaluators *mRNA_exon_evaluators,
                       *mRNA_exon_evaluators_collapsed,
                       *CDS_exon_evaluators,
                       *CDS_exon_evaluators_collapsed;
  unsigned long *wrong_genes,
                *wrong_mRNAs,
                *wrong_LTRs;
} ProcessPredictedFeatureInfo;

static Slot* slot_new(bool nuceval, GT_Range range)
{
  unsigned long length;
  Slot *s = gt_calloc(1, sizeof (Slot));
  length = gt_range_length(range);
  s->genes_forward = gt_array_new(sizeof (GT_GenomeNode*));
  s->genes_reverse = gt_array_new(sizeof (GT_GenomeNode*));
  s->mRNAs_forward = gt_array_new(sizeof (GT_GenomeNode*));
  s->mRNAs_reverse = gt_array_new(sizeof (GT_GenomeNode*));
  s->LTRs          = gt_array_new(sizeof (GT_GenomeNode*));
  s->mRNA_exons_forward = transcript_exons_new();
  s->mRNA_exons_reverse = transcript_exons_new();
  s->CDS_exons_forward = transcript_exons_new();
  s->CDS_exons_reverse = transcript_exons_new();
  if (nuceval) {
    s->real_range = range;
    s->real_mRNA_nucleotides_forward = gt_bittab_new(length);
    s->pred_mRNA_nucleotides_forward = gt_bittab_new(length);
    s->real_mRNA_nucleotides_reverse = gt_bittab_new(length);
    s->pred_mRNA_nucleotides_reverse = gt_bittab_new(length);
    s->real_CDS_nucleotides_forward = gt_bittab_new(length);
    s->pred_CDS_nucleotides_forward = gt_bittab_new(length);
    s->real_CDS_nucleotides_reverse = gt_bittab_new(length);
    s->pred_CDS_nucleotides_reverse = gt_bittab_new(length);
  }
  s->used_mRNA_exons_forward = transcript_used_exons_new();
  s->used_mRNA_exons_reverse = transcript_used_exons_new();
  s->used_CDS_exons_forward = transcript_used_exons_new();
  s->used_CDS_exons_reverse = transcript_used_exons_new();
  return s;
}

static void slot_delete(Slot *s)
{
  unsigned long i;
  assert(s);
  for (i = 0; i < gt_array_size(s->genes_forward); i++)
    gt_genome_node_rec_delete(*(GT_GenomeNode**) gt_array_get(s->genes_forward, i));
  gt_array_delete(s->genes_forward);
  for (i = 0; i < gt_array_size(s->genes_reverse); i++)
    gt_genome_node_rec_delete(*(GT_GenomeNode**) gt_array_get(s->genes_reverse, i));
  gt_array_delete(s->genes_reverse);
  for (i = 0; i < gt_array_size(s->mRNAs_forward); i++)
    gt_genome_node_rec_delete(*(GT_GenomeNode**) gt_array_get(s->mRNAs_forward, i));
  gt_array_delete(s->mRNAs_forward);
  for (i = 0; i < gt_array_size(s->mRNAs_reverse); i++)
    gt_genome_node_rec_delete(*(GT_GenomeNode**) gt_array_get(s->mRNAs_reverse, i));
  gt_array_delete(s->mRNAs_reverse);
  for (i = 0; i < gt_array_size(s->LTRs); i++)
    gt_genome_node_rec_delete(*(GT_GenomeNode**) gt_array_get(s->LTRs, i));
  gt_array_delete(s->LTRs);
  transcript_exons_delete(s->mRNA_exons_forward);
  transcript_exons_delete(s->mRNA_exons_reverse);
  transcript_exons_delete(s->CDS_exons_forward);
  transcript_exons_delete(s->CDS_exons_reverse);
  transcript_counts_delete(s->mRNA_counts_forward);
  transcript_counts_delete(s->mRNA_counts_reverse);
  transcript_counts_delete(s->CDS_counts_forward);
  transcript_counts_delete(s->CDS_counts_reverse);
  gt_bittab_delete(s->real_mRNA_nucleotides_forward);
  gt_bittab_delete(s->pred_mRNA_nucleotides_forward);
  gt_bittab_delete(s->real_mRNA_nucleotides_reverse);
  gt_bittab_delete(s->pred_mRNA_nucleotides_reverse);
  gt_bittab_delete(s->real_CDS_nucleotides_forward);
  gt_bittab_delete(s->pred_CDS_nucleotides_forward);
  gt_bittab_delete(s->real_CDS_nucleotides_reverse);
  gt_bittab_delete(s->pred_CDS_nucleotides_reverse);
  gt_bittab_delete(s->true_genes_forward);
  gt_bittab_delete(s->true_genes_reverse);
  gt_bittab_delete(s->true_mRNAs_forward);
  gt_bittab_delete(s->true_mRNAs_reverse);
  gt_bittab_delete(s->true_LTRs);
  gt_bittab_delete(s->overlapped_genes_forward);
  gt_bittab_delete(s->overlapped_genes_reverse);
  gt_bittab_delete(s->overlapped_mRNAs_forward);
  gt_bittab_delete(s->overlapped_mRNAs_reverse);
  gt_bittab_delete(s->overlapped_LTRs);
  transcript_bittabs_delete(s->mRNA_exon_bittabs_forward);
  transcript_bittabs_delete(s->mRNA_exon_bittabs_reverse);
  transcript_bittabs_delete(s->CDS_exon_bittabs_forward);
  transcript_bittabs_delete(s->CDS_exon_bittabs_reverse);
  transcript_used_exons_delete(s->used_mRNA_exons_forward);
  transcript_used_exons_delete(s->used_mRNA_exons_reverse);
  transcript_used_exons_delete(s->used_CDS_exons_forward);
  transcript_used_exons_delete(s->used_CDS_exons_reverse);
  gt_free(s);
}

StreamEvaluator* stream_evaluator_new(GenomeStream *reality,
                                      GenomeStream *prediction, bool nuceval,
                                      bool evalLTR, unsigned long LTRdelta)
{
  StreamEvaluator *evaluator = gt_calloc(1, sizeof (StreamEvaluator));
  evaluator->reality = genome_stream_ref(reality);
  evaluator->prediction = genome_stream_ref(prediction);
  evaluator->nuceval = nuceval;
  evaluator->evalLTR = evalLTR;
  evaluator->LTRdelta = LTRdelta;
  evaluator->slots = hashmap_new(HASH_STRING, gt_free_func,
                                 (GT_FreeFunc) slot_delete);
  evaluator->gene_evaluator = evaluator_new();
  evaluator->mRNA_evaluator = evaluator_new();
  evaluator->LTR_evaluator = evaluator_new();
  evaluator->mRNA_exon_evaluators = transcript_evaluators_new();
  evaluator->mRNA_exon_evaluators_collapsed = transcript_evaluators_new();
  evaluator->CDS_exon_evaluators = transcript_evaluators_new();
  evaluator->CDS_exon_evaluators_collapsed = transcript_evaluators_new();
  return evaluator;
}

static int set_actuals_and_sort_them(GT_UNUSED void *key, void *value, void *data,
                                     GT_UNUSED GT_Error *err)
{
  StreamEvaluator *se = (StreamEvaluator*) data;
  Slot *s = (Slot*) value;

  gt_error_check(err);
  assert(key && value && data);

  /* set actual genes */
  evaluator_add_actual(se->gene_evaluator, gt_array_size(s->genes_forward));
  evaluator_add_actual(se->gene_evaluator, gt_array_size(s->genes_reverse));

  /* set actual mRNAs */
  evaluator_add_actual(se->mRNA_evaluator, gt_array_size(s->mRNAs_forward));
  evaluator_add_actual(se->mRNA_evaluator, gt_array_size(s->mRNAs_reverse));

  /* set actual LTRs */
  evaluator_add_actual(se->LTR_evaluator, gt_array_size(s->LTRs));

  /* set actual exons (before uniq!) */
  transcript_evaluators_add_actuals(se->mRNA_exon_evaluators,
                                    s->mRNA_exons_forward);
  transcript_evaluators_add_actuals(se->mRNA_exon_evaluators,
                                    s->mRNA_exons_reverse);
  transcript_evaluators_add_actuals(se->CDS_exon_evaluators,
                                    s->CDS_exons_forward);
  transcript_evaluators_add_actuals(se->CDS_exon_evaluators,
                                    s->CDS_exons_reverse);

  /* sort genes */
  genome_nodes_sort(s->genes_forward);
  genome_nodes_sort(s->genes_reverse);

  /* sort mRNAs */
  genome_nodes_sort(s->mRNAs_forward);
  genome_nodes_sort(s->mRNAs_reverse);

  /* sort LTRs */
  genome_nodes_sort(s->LTRs);

  /* sort exons */
  transcript_exons_sort(s->mRNA_exons_forward);
  transcript_exons_sort(s->mRNA_exons_reverse);
  transcript_exons_sort(s->CDS_exons_forward);
  transcript_exons_sort(s->CDS_exons_reverse);

  /* determine true exons */
  s->mRNA_counts_forward =
    transcript_exons_uniq_in_place_count(s->mRNA_exons_forward);
  s->mRNA_counts_reverse =
    transcript_exons_uniq_in_place_count(s->mRNA_exons_reverse);
  s->CDS_counts_forward =
    transcript_exons_uniq_in_place_count(s->CDS_exons_forward);
  s->CDS_counts_reverse =
    transcript_exons_uniq_in_place_count(s->CDS_exons_reverse);

  /* set actual exons for the collapsed case (after uniq!) */
  transcript_evaluators_add_actuals(se->mRNA_exon_evaluators_collapsed,
                                    s->mRNA_exons_forward);
  transcript_evaluators_add_actuals(se->mRNA_exon_evaluators_collapsed,
                                    s->mRNA_exons_reverse);
  transcript_evaluators_add_actuals(se->CDS_exon_evaluators_collapsed,
                                    s->CDS_exons_forward);
  transcript_evaluators_add_actuals(se->CDS_exon_evaluators_collapsed,
                                    s->CDS_exons_reverse);

  /* make sure that the genes are sorted */
  assert(genome_nodes_are_sorted(s->genes_forward));
  assert(genome_nodes_are_sorted(s->genes_reverse));

  /* make sure that the mRNAs are sorted */
  assert(genome_nodes_are_sorted(s->mRNAs_forward));
  assert(genome_nodes_are_sorted(s->mRNAs_reverse));

  /* make sure that the LTRs are sorted */
  assert(genome_nodes_are_sorted(s->LTRs));

  /* make sure that the exons are sorted */
  assert(transcript_exons_are_sorted(s->mRNA_exons_forward));
  assert(transcript_exons_are_sorted(s->mRNA_exons_reverse));
  assert(transcript_exons_are_sorted(s->CDS_exons_forward));
  assert(transcript_exons_are_sorted(s->CDS_exons_reverse));

  /* init true bittabs */
  s->true_genes_forward = gt_array_size(s->genes_forward)
                          ? gt_bittab_new(gt_array_size(s->genes_forward))
                          : NULL;
  s->true_genes_reverse = gt_array_size(s->genes_reverse)
                          ? gt_bittab_new(gt_array_size(s->genes_reverse))
                          : NULL;
  s->true_mRNAs_forward = gt_array_size(s->mRNAs_forward)
                          ? gt_bittab_new(gt_array_size(s->mRNAs_forward))
                          : NULL;
  s->true_mRNAs_reverse = gt_array_size(s->mRNAs_reverse)
                          ? gt_bittab_new(gt_array_size(s->mRNAs_reverse))
                          : NULL;
  s->true_LTRs          = gt_array_size(s->LTRs)
                          ? gt_bittab_new(gt_array_size(s->LTRs))
                          : NULL;

  /* init overlap bittabs */
  s->overlapped_genes_forward = gt_array_size(s->genes_forward)
                                ? gt_bittab_new(gt_array_size(s->genes_forward))
                                : NULL;
  s->overlapped_genes_reverse = gt_array_size(s->genes_reverse)
                                ? gt_bittab_new(gt_array_size(s->genes_reverse))
                                : NULL;
  s->overlapped_mRNAs_forward = gt_array_size(s->mRNAs_forward)
                                ? gt_bittab_new(gt_array_size(s->mRNAs_forward))
                                : NULL;
  s->overlapped_mRNAs_reverse = gt_array_size(s->mRNAs_reverse)
                                ? gt_bittab_new(gt_array_size(s->mRNAs_reverse))
                                : NULL;
  s->overlapped_LTRs          = gt_array_size(s->LTRs)
                                ? gt_bittab_new(gt_array_size(s->LTRs))
                                : NULL;

  /* init bittabs (for collapsed exons) */
  s->mRNA_exon_bittabs_forward =
    transcript_exons_create_bittabs(s->mRNA_exons_forward);
  s->mRNA_exon_bittabs_reverse =
    transcript_exons_create_bittabs(s->mRNA_exons_reverse);
  s->CDS_exon_bittabs_forward =
    transcript_exons_create_bittabs(s->CDS_exons_forward);
  s->CDS_exon_bittabs_reverse =
    transcript_exons_create_bittabs(s->CDS_exons_reverse);

  return 0;
}

static void add_real_exon(TranscriptExons *te, GT_Range range, GT_GenomeNode *gn)
{
  assert(te);
  gt_array_add(transcript_exons_get_all(te), range);
  switch (gt_genome_feature_get_transcriptfeaturetype((GT_GenomeFeature*) gn)) {
    case TRANSCRIPT_FEATURE_TYPE_UNDETERMINED:
      warning("type of feature (single, initial, internal, or terminal) given "
              "on line %u in file \"%s\" could not be determined, because the "
              "feature has no Parent attribute. Treating it as single.",
              gt_genome_node_get_line_number(gn), gt_genome_node_get_filename(gn));
      /*@fallthrough@*/
    case TRANSCRIPT_FEATURE_TYPE_SINGLE:
      gt_array_add(transcript_exons_get_single(te), range);
      break;
    case TRANSCRIPT_FEATURE_TYPE_INITIAL:
      gt_array_add(transcript_exons_get_initial(te), range);
      break;
    case TRANSCRIPT_FEATURE_TYPE_INTERNAL:
      gt_array_add(transcript_exons_get_internal(te), range);
      break;
    case TRANSCRIPT_FEATURE_TYPE_TERMINAL:
      gt_array_add(transcript_exons_get_terminal(te), range);
      break;
  }
}

static void add_nucleotide_exon(GT_Bittab *nucleotides, GT_Range range,
                                GT_Range real_range,
                                unsigned long *FP)
{
  unsigned long i;
  assert(nucleotides);
  for (i = range.start; i <= range.end; i++) {
    if (gt_range_within(real_range, i)) {
      assert(i >= real_range.start);
      gt_bittab_set_bit(nucleotides, i - real_range.start);
    }
    else {
      assert(FP);
      (*FP)++;
    }
  }
}

static int process_real_feature(GT_GenomeNode *gn, void *data, GT_UNUSED GT_Error *err)
{
  ProcessRealFeatureInfo *info = (ProcessRealFeatureInfo*) data;
  GT_GenomeNode *gn_ref;
  GT_GenomeFeature *gf;
  GT_Range range;

  gt_error_check(err);
  assert(gn && data);
  gf = (GT_GenomeFeature*) gn;

  if (gt_genome_feature_has_type(gf, gft_gene)) {
    switch (gt_genome_feature_get_strand(gf)) {
      case GT_STRAND_FORWARD:
        gn_ref = gt_genome_node_rec_ref(gn);
        gt_array_add(info->slot->genes_forward, gn_ref);
        break;
      case GT_STRAND_REVERSE:
        gn_ref = gt_genome_node_rec_ref(gn);
        gt_array_add(info->slot->genes_reverse, gn_ref);
        break;
      default:
        if (info->verbose) {
          fprintf(stderr, "skipping real gene with unknown orientation "
                  "(line %u)\n", gt_genome_node_get_line_number(gn));
        }
    }
  }
  else if (gt_genome_feature_has_type(gf, gft_mRNA)) {
    switch (gt_genome_feature_get_strand(gf)) {
      case GT_STRAND_FORWARD:
        gn_ref = gt_genome_node_rec_ref(gn);
        gt_array_add(info->slot->mRNAs_forward, gn_ref);
        break;
      case GT_STRAND_REVERSE:
        gn_ref = gt_genome_node_rec_ref(gn);
        gt_array_add(info->slot->mRNAs_reverse, gn_ref);
        break;
      default:
        if (info->verbose) {
          fprintf(stderr, "skipping real mRNA with unknown orientation "
                  "(line %u)\n", gt_genome_node_get_line_number(gn));
        }
    }
  }
  else if (gt_genome_feature_has_type(gf, gft_LTR_retrotransposon)) {
    gn_ref = gt_genome_node_rec_ref(gn);
    gt_array_add(info->slot->LTRs, gn_ref);
  }
  else if (gt_genome_feature_has_type(gf, gft_CDS)) {
    range = gt_genome_node_get_range(gn);
    switch (gt_genome_feature_get_strand(gf)) {
      case GT_STRAND_FORWARD:
        add_real_exon(info->slot->CDS_exons_forward, range, gn);
        /* nucleotide level */
        if (info->nuceval) {
          add_nucleotide_exon(info->slot->real_CDS_nucleotides_forward, range,
                              info->slot->real_range, NULL);
        }
        break;
      case GT_STRAND_REVERSE:
        add_real_exon(info->slot->CDS_exons_reverse, range, gn);
        /* nucleotide level */
        if (info->nuceval) {
          add_nucleotide_exon(info->slot->real_CDS_nucleotides_reverse, range,
                              info->slot->real_range, NULL);
        }
        break;
      default:
        if (info->verbose) {
          fprintf(stderr, "skipping real CDS exon with unknown orientation "
                  "(line %u)\n", gt_genome_node_get_line_number(gn));
        }
    }
  }
  else if (gt_genome_feature_has_type(gf, gft_exon)) {
    range = gt_genome_node_get_range(gn);
    switch (gt_genome_feature_get_strand(gf)) {
      case GT_STRAND_FORWARD:
        add_real_exon(info->slot->mRNA_exons_forward, range, gn);
        /* nucleotide level */
        if (info->nuceval) {
          add_nucleotide_exon(info->slot->real_mRNA_nucleotides_forward,
                              range, info->slot->real_range, NULL);
        }
        break;
      case GT_STRAND_REVERSE:
        add_real_exon(info->slot->mRNA_exons_reverse, range, gn);
        /* nucleotide level */
        if (info->nuceval) {
          add_nucleotide_exon(info->slot->real_mRNA_nucleotides_reverse,
                              range, info->slot->real_range, NULL);
        }
        break;
      default:
        if (info->verbose) {
          fprintf(stderr, "skipping real mRNA exon with unknown orientation "
                  "(line %u)\n", gt_genome_node_get_line_number(gn));
        }
    }
  }
  return 0;
}

static int store_exon(GT_GenomeNode *gn, void *data, GT_UNUSED GT_Error *err)
{
  GT_Array *exons = (GT_Array*) data;
  GT_Range range;
  GT_GenomeFeature *gf;
  gt_error_check(err);
  gf = gt_genome_node_cast(gt_genome_feature_class(), gn);
  assert(gf && exons);
  if (gt_genome_feature_has_type(gf, gft_exon)) {
    range = gt_genome_node_get_range(gn);
    gt_array_add(exons, range);
  }
  return 0;
}

static bool mRNAs_are_equal(GT_GenomeNode *gn_1, GT_GenomeNode *gn_2)
{
  GT_Array *exons_1, *exons_2;
  bool equal;
  int had_err;

  assert(gn_1 && gn_2);

  /* init */
  exons_1 = gt_array_new(sizeof (GT_Range));
  exons_2 = gt_array_new(sizeof (GT_Range));

  /* get exon ranges */
  had_err = gt_genome_node_traverse_children(gn_1, exons_1, store_exon, false,
                                          NULL);
  assert(!had_err); /* cannot happen, store_exon() is sane */
  had_err = gt_genome_node_traverse_children(gn_2, exons_2, store_exon, false,
                                          NULL);
  assert(!had_err); /* cannot happen, store_exon() is sane */

  /* sort exon ranges */
  ranges_sort(exons_1);
  ranges_sort(exons_2);

  /* compare exon ranges */
  equal = ranges_are_equal(exons_1, exons_2);

  /* free */
  gt_array_delete(exons_1);
  gt_array_delete(exons_2);

  return equal;
}

typedef struct {
  GT_Array *exons,
        *mRNAs;
} Store_gene_feature_info;

static int store_gene_feature(GT_GenomeNode *gn, void *data, GT_UNUSED GT_Error *err)
{
  GT_GenomeFeature *gf;
  Store_gene_feature_info *info = (Store_gene_feature_info*) data;
  GT_Range range;
  gt_error_check(err);
  gf = gt_genome_node_cast(gt_genome_feature_class(), gn);
  assert(gf && info);
  if (gt_genome_feature_has_type(gf, gft_mRNA)) {
    gt_array_add(info->mRNAs, gf);
  }
  else if (gt_genome_feature_has_type(gf, gft_exon)) {
    range = gt_genome_node_get_range(gn);
    gt_array_add(info->exons, range);
  }
  return 0;
}

static bool genes_are_equal(GT_GenomeNode *gn_1, GT_GenomeNode *gn_2)
{
  GT_Array *exons_1, *exons_2, *mRNAs_1, *mRNAs_2;
  Store_gene_feature_info info;
  unsigned long i;
  bool equal;
  int had_err;

  /* init */
  exons_1 = gt_array_new(sizeof (GT_Range));
  exons_2 = gt_array_new(sizeof (GT_Range));
  mRNAs_1 = gt_array_new(sizeof (GT_GenomeNode*));
  mRNAs_2 = gt_array_new(sizeof (GT_GenomeNode*));

  /* get (direct) gene features */
  info.exons = exons_1;
  info.mRNAs = mRNAs_1;
  had_err = gt_genome_node_traverse_direct_children(gn_1, &info,
                                                 store_gene_feature, NULL);
  assert(!had_err); /* cannot happen, store_gene_feature() is sane */
  info.exons = exons_2;
  info.mRNAs = mRNAs_2;
  had_err = gt_genome_node_traverse_direct_children(gn_2, &info,
                                                 store_gene_feature, NULL);
  assert(!had_err); /* cannot happen, store_gene_feature() is sane */

  /* sort exon ranges */
  ranges_sort(exons_1);
  ranges_sort(exons_2);

  /* compare exon ranges */
  equal = ranges_are_equal(exons_1, exons_2);

  /* compare mRNAs, if necessary */
  if (equal && gt_array_size(mRNAs_1) == gt_array_size(mRNAs_2)) {
    /* sort mRNAs */
    genome_nodes_sort(mRNAs_1);
    genome_nodes_sort(mRNAs_2);
    for (i = 0; i < gt_array_size(mRNAs_1); i++) {
      assert(equal);
      equal = mRNAs_are_equal(*(GT_GenomeNode**) gt_array_get(mRNAs_1, i),
                              *(GT_GenomeNode**) gt_array_get(mRNAs_2, i));
      if (!equal)
        break;
    }
  }

  /* free */
  gt_array_delete(exons_1);
  gt_array_delete(exons_2);
  gt_array_delete(mRNAs_1);
  gt_array_delete(mRNAs_2);

  return equal;
}

static void store_predicted_exon(TranscriptEvaluators *te, GT_GenomeNode *gn)
{
  assert(te && gn);
  evaluator_add_predicted(transcript_evaluators_get_all(te), 1);
  switch (gt_genome_feature_get_transcriptfeaturetype((GT_GenomeFeature*) gn)) {
    case TRANSCRIPT_FEATURE_TYPE_UNDETERMINED:
      warning("type of feature (single, initial, internal, or terminal) given "
              "on line %u in file \"%s\" could not be determined, because the "
              "feature has no Parent attribute. Treating it as single.",
              gt_genome_node_get_line_number(gn), gt_genome_node_get_filename(gn));
      /*@fallthrough@*/
    case TRANSCRIPT_FEATURE_TYPE_SINGLE:
      evaluator_add_predicted(transcript_evaluators_get_single(te), 1);
    break;
    case TRANSCRIPT_FEATURE_TYPE_INITIAL:
      evaluator_add_predicted(transcript_evaluators_get_initial(te), 1);
    break;
    case TRANSCRIPT_FEATURE_TYPE_INTERNAL:
      evaluator_add_predicted(transcript_evaluators_get_internal(te), 1);
    break;
    case TRANSCRIPT_FEATURE_TYPE_TERMINAL:
      evaluator_add_predicted(transcript_evaluators_get_terminal(te), 1);
    break;
  }
}

/* adds exon only if necessary */
static void add_predicted_collapsed(GT_Dlist *used_exons, GT_Range *predicted_range,
                                    Evaluator *exon_evaluator_collapsed)
{
  GT_Range *used_range;
  if (!gt_dlist_find(used_exons, predicted_range)) {
    used_range = gt_malloc(sizeof (GT_Range));
    used_range->start = predicted_range->start;
    used_range->end = predicted_range->end;
    gt_dlist_add(used_exons, used_range);
    evaluator_add_predicted(exon_evaluator_collapsed, 1);
  }
}

static void store_predicted_exon_collapsed(TranscriptUsedExons *used_exons,
                                           GT_Range *predicted_range,
                                           TranscriptEvaluators *te,
                                           GT_GenomeNode *gn)
{
  add_predicted_collapsed(transcript_used_exons_get_all(used_exons),
                          predicted_range, transcript_evaluators_get_all(te));
  switch (gt_genome_feature_get_transcriptfeaturetype((GT_GenomeFeature*) gn)) {
    case TRANSCRIPT_FEATURE_TYPE_UNDETERMINED:
      /* we do not show a warning here, because store_predicted_exon() has been
         called before and already shown one */
    case TRANSCRIPT_FEATURE_TYPE_SINGLE:
      add_predicted_collapsed(transcript_used_exons_get_single(used_exons),
                              predicted_range,
                              transcript_evaluators_get_single(te));
      break;
    case TRANSCRIPT_FEATURE_TYPE_INITIAL:
      add_predicted_collapsed(transcript_used_exons_get_initial(used_exons),
                              predicted_range,
                              transcript_evaluators_get_initial(te));
      break;
    case TRANSCRIPT_FEATURE_TYPE_INTERNAL:
      add_predicted_collapsed(transcript_used_exons_get_internal(used_exons),
                              predicted_range,
                              transcript_evaluators_get_internal(te));
      break;
    case TRANSCRIPT_FEATURE_TYPE_TERMINAL:
      add_predicted_collapsed(transcript_used_exons_get_terminal(used_exons),
                              predicted_range,
                              transcript_evaluators_get_terminal(te));
      break;
  }
}

static void mark_and_show_false_exon(GT_GenomeNode *gn, bool exondiff)
{
  gt_genome_node_mark(gn); /* mark false exons */
  if (exondiff) {
    gff3_output_leading((GT_GenomeFeature*) gn, NULL);
    printf(".\n");
  }
}

static void determine_true_exon(GT_GenomeNode *gn, GT_Strand predicted_strand,
                                bool exondiff, GT_Range *predicted_range,
                                GT_Array *exons_forward,
                                GT_Array *exons_reverse,
                                GT_Array *true_exons_forward,
                                GT_Array *true_exons_reverse,
                                GT_Bittab *true_exons_forward_collapsed,
                                GT_Bittab *true_exons_reverse_collapsed,
                                Evaluator *exon_evaluator,
                                Evaluator *exon_evaluator_collapsed)
{
  GT_Range *actual_range;
  unsigned long num, *ctr_ptr;

  if ((actual_range = bsearch(predicted_range,
                              predicted_strand == GT_STRAND_FORWARD
                              ? gt_array_get_space(exons_forward)
                              : gt_array_get_space(exons_reverse),
                              predicted_strand == GT_STRAND_FORWARD
                              ? gt_array_size(exons_forward)
                              : gt_array_size(exons_reverse), sizeof (GT_Range),
                              (GT_Compare) gt_range_compare_ptr))) {
    if (predicted_strand == GT_STRAND_FORWARD) {
      num = actual_range - (GT_Range*) gt_array_get_space(exons_forward);
      ctr_ptr = gt_array_get(true_exons_forward, num);
      if (*ctr_ptr) {
        (*ctr_ptr)--;
        evaluator_add_true(exon_evaluator);
      }
      else
        mark_and_show_false_exon(gn, exondiff);
      if (true_exons_forward_collapsed &&
          !gt_bittab_bit_is_set(true_exons_forward_collapsed, num)) {
        gt_bittab_set_bit(true_exons_forward_collapsed, num);
        evaluator_add_true(exon_evaluator_collapsed);
      }
    }
    else {
      num = actual_range - (GT_Range*) gt_array_get_space(exons_reverse);
      ctr_ptr = gt_array_get(true_exons_reverse, num);
      if (*ctr_ptr) {
        (*ctr_ptr)--;
        evaluator_add_true(exon_evaluator);
      }
      else
        mark_and_show_false_exon(gn, exondiff);
      if (true_exons_reverse_collapsed &&
          !gt_bittab_bit_is_set(true_exons_reverse_collapsed, num)) {
        gt_bittab_set_bit(true_exons_reverse_collapsed, num);
        evaluator_add_true(exon_evaluator_collapsed);
      }
    }
  }
  else
    mark_and_show_false_exon(gn, exondiff);
}

static void store_true_exon(GT_GenomeNode *gn, GT_Strand predicted_strand,
                            GT_Range *predicted_range, bool exondiff,
                            TranscriptExons *exons_forward,
                            TranscriptExons *exons_reverse,
                            TranscriptCounts *counts_forward,
                            TranscriptCounts *counts_reverse,
                            TranscriptGT_Bittabs *exon_bittabs_forward,
                            TranscriptGT_Bittabs *exon_bittabs_reverse,
                            TranscriptEvaluators *exon_evaluators,
                            TranscriptEvaluators *exon_evaluators_collapsed)
{
  assert(gn && predicted_range && exons_forward && exons_reverse);
  determine_true_exon(gn, predicted_strand, exondiff, predicted_range,
                      transcript_exons_get_all(exons_forward),
                      transcript_exons_get_all(exons_reverse),
                      transcript_counts_get_all(counts_forward),
                      transcript_counts_get_all(counts_reverse),
                      transcript_bittabs_get_all(exon_bittabs_forward),
                      transcript_bittabs_get_all(exon_bittabs_reverse),
                      transcript_evaluators_get_all(exon_evaluators),
                      transcript_evaluators_get_all(exon_evaluators_collapsed));
  switch (gt_genome_feature_get_transcriptfeaturetype((GT_GenomeFeature*) gn)) {
    case TRANSCRIPT_FEATURE_TYPE_UNDETERMINED:
    case TRANSCRIPT_FEATURE_TYPE_SINGLE:
      determine_true_exon(gn, predicted_strand, exondiff, predicted_range,
                          transcript_exons_get_single(exons_forward),
                          transcript_exons_get_single(exons_reverse),
                          transcript_counts_get_single(counts_forward),
                          transcript_counts_get_single(counts_reverse),
                          transcript_bittabs_get_single(exon_bittabs_forward),
                          transcript_bittabs_get_single(exon_bittabs_reverse),
                          transcript_evaluators_get_single(exon_evaluators),
                          transcript_evaluators_get_single(
                            exon_evaluators_collapsed));
      break;
    case TRANSCRIPT_FEATURE_TYPE_INITIAL:
      determine_true_exon(gn, predicted_strand, exondiff, predicted_range,
                          transcript_exons_get_initial(exons_forward),
                          transcript_exons_get_initial(exons_reverse),
                          transcript_counts_get_initial(counts_forward),
                          transcript_counts_get_initial(counts_reverse),
                          transcript_bittabs_get_initial(exon_bittabs_forward),
                          transcript_bittabs_get_initial(exon_bittabs_reverse),
                          transcript_evaluators_get_initial(exon_evaluators),
                          transcript_evaluators_get_initial(
                            exon_evaluators_collapsed));
      break;
    case TRANSCRIPT_FEATURE_TYPE_INTERNAL:
      determine_true_exon(gn, predicted_strand, exondiff, predicted_range,
                          transcript_exons_get_internal(exons_forward),
                          transcript_exons_get_internal(exons_reverse),
                          transcript_counts_get_internal(counts_forward),
                          transcript_counts_get_internal(counts_reverse),
                          transcript_bittabs_get_internal(exon_bittabs_forward),
                          transcript_bittabs_get_internal(exon_bittabs_reverse),
                          transcript_evaluators_get_internal(exon_evaluators),
                          transcript_evaluators_get_internal(
                            exon_evaluators_collapsed));
      break;
    case TRANSCRIPT_FEATURE_TYPE_TERMINAL:
      determine_true_exon(gn, predicted_strand, exondiff, predicted_range,
                          transcript_exons_get_terminal(exons_forward),
                          transcript_exons_get_terminal(exons_reverse),
                          transcript_counts_get_terminal(counts_forward),
                          transcript_counts_get_terminal(counts_reverse),
                          transcript_bittabs_get_terminal(exon_bittabs_forward),
                          transcript_bittabs_get_terminal(exon_bittabs_reverse),
                          transcript_evaluators_get_terminal(exon_evaluators),
                          transcript_evaluators_get_terminal(
                            exon_evaluators_collapsed));
      break;
  }
}

static int process_predicted_feature(GT_GenomeNode *gn, void *data,
                                     GT_UNUSED GT_Error *err)
{
  ProcessPredictedFeatureInfo *info = (ProcessPredictedFeatureInfo*) data;
  GT_Range predicted_range;
  unsigned long i, num;
  GT_Strand predicted_strand;
  GT_Array *real_genome_nodes;
  GT_GenomeNode **real_gn;

  gt_error_check(err);
  assert(gn && data);

  predicted_range = gt_genome_node_get_range(gn);
  predicted_strand = gt_genome_feature_get_strand((GT_GenomeFeature*) gn);
  real_genome_nodes = gt_array_new(sizeof (GT_GenomeNode**));

  if (gt_genome_feature_has_type((GT_GenomeFeature*) gn, gft_gene)) {
    /* store predicted gene */
    evaluator_add_predicted(info->gene_evaluator, 1);
    /* determine true gene */
    switch (predicted_strand) {
      case GT_STRAND_FORWARD:
      case GT_STRAND_REVERSE:
        bsearch_all_mark(real_genome_nodes, &gn,
                         predicted_strand == GT_STRAND_FORWARD
                         ? gt_array_get_space(info->slot->genes_forward)
                         : gt_array_get_space(info->slot->genes_reverse),
                         predicted_strand == GT_STRAND_FORWARD
                         ? gt_array_size(info->slot->genes_forward)
                         : gt_array_size(info->slot->genes_reverse),
                         sizeof (GT_GenomeNode*),
                         (GT_CompareWithData) gt_genome_node_compare_with_data,
                         NULL,
                         predicted_strand == GT_STRAND_FORWARD
                         ? info->slot->overlapped_genes_forward
                         : info->slot->overlapped_genes_reverse);
        if (gt_array_size(real_genome_nodes)) {
          /* gene(s) with the same range found -> check if they are equal */
          for (i = 0; i < gt_array_size(real_genome_nodes); i++) {
            real_gn = *(GT_GenomeNode***) gt_array_get(real_genome_nodes, i);
            if (genes_are_equal(gn, *real_gn)) {
              if (predicted_strand == GT_STRAND_FORWARD) {
                num = real_gn - (GT_GenomeNode**)
                      gt_array_get_space(info->slot->genes_forward);
                if (!gt_bittab_bit_is_set(info->slot->true_genes_forward, num)) {
                  gt_bittab_set_bit(info->slot->true_genes_forward, num);
                  evaluator_add_true(info->gene_evaluator);
                  /*@loopbreak@*/
                  break;
                }
              }
              else {
                num = real_gn - (GT_GenomeNode**)
                      gt_array_get_space(info->slot->genes_reverse);
                if (!gt_bittab_bit_is_set(info->slot->true_genes_reverse, num)) {
                  gt_bittab_set_bit(info->slot->true_genes_reverse, num);
                  evaluator_add_true(info->gene_evaluator);
                  /*@loopbreak@*/
                  break;
                }
              }
            }
          }
        }
        else {
          /* no gene with the same range found -> check if this is a wrong
             gene */
          if (!gt_genome_node_overlaps_nodes_mark(gn,
                                    predicted_strand == GT_STRAND_FORWARD
                                    ? info->slot->genes_forward
                                    : info->slot->genes_reverse,
                                    predicted_strand == GT_STRAND_FORWARD
                                    ? info->slot->overlapped_genes_forward
                                    : info->slot->overlapped_genes_reverse)) {
            (*info->wrong_genes)++;
          }
        }
        break;
      default:
        if (info->verbose) {
          fprintf(stderr, "skipping predicted gene with unknown orientation "
                  "(line %u)\n", gt_genome_node_get_line_number(gn));
        }
    }
  }
  else if (gt_genome_feature_has_type((GT_GenomeFeature*) gn, gft_mRNA)) {
    /* store predicted mRNA */
    evaluator_add_predicted(info->mRNA_evaluator, 1);
    /* determine true mRNA */
    switch (predicted_strand) {
      case GT_STRAND_FORWARD:
      case GT_STRAND_REVERSE:
        bsearch_all_mark(real_genome_nodes, &gn,
                         predicted_strand == GT_STRAND_FORWARD
                         ? gt_array_get_space(info->slot->mRNAs_forward)
                         : gt_array_get_space(info->slot->mRNAs_reverse),
                         predicted_strand == GT_STRAND_FORWARD
                         ? gt_array_size(info->slot->mRNAs_forward)
                         : gt_array_size(info->slot->mRNAs_reverse),
                         sizeof (GT_GenomeNode*),
                         (GT_CompareWithData) gt_genome_node_compare_with_data,
                         NULL,
                         predicted_strand == GT_STRAND_FORWARD
                         ? info->slot->overlapped_mRNAs_forward
                         : info->slot->overlapped_mRNAs_reverse);
        if (gt_array_size(real_genome_nodes)) {
          /* mRNA(s) with the same range found -> check if they are equal */
          for (i = 0; i < gt_array_size(real_genome_nodes); i++) {
            real_gn = *(GT_GenomeNode***) gt_array_get(real_genome_nodes, i);
            if (mRNAs_are_equal(gn, *real_gn)) {
              if (predicted_strand == GT_STRAND_FORWARD) {
                num = real_gn - (GT_GenomeNode**)
                      gt_array_get_space(info->slot->mRNAs_forward);
                if (!gt_bittab_bit_is_set(info->slot->true_mRNAs_forward, num)) {
                  gt_bittab_set_bit(info->slot->true_mRNAs_forward, num);
                  evaluator_add_true(info->mRNA_evaluator);
                  /*@loopbreak@*/
                  break;
                }
              }
              else {
                num = real_gn - (GT_GenomeNode**)
                      gt_array_get_space(info->slot->mRNAs_reverse);
                if (!gt_bittab_bit_is_set(info->slot->true_mRNAs_reverse, num)) {
                  gt_bittab_set_bit(info->slot->true_mRNAs_reverse, num);
                  evaluator_add_true(info->mRNA_evaluator);
                  /*@loopbreak@*/
                  break;
                }
              }
            }
          }
        }
        else {
          /* no mRNA with the same range found -> check if this is a wrong
             mRNA */
          if (!gt_genome_node_overlaps_nodes_mark(gn,
                                    predicted_strand == GT_STRAND_FORWARD
                                    ? info->slot->mRNAs_forward
                                    : info->slot->mRNAs_reverse,
                                    predicted_strand == GT_STRAND_FORWARD
                                    ? info->slot->overlapped_mRNAs_forward
                                    : info->slot->overlapped_mRNAs_reverse)) {
            (*info->wrong_mRNAs)++;
          }
        }
        break;
      default:
        if (info->verbose) {
          fprintf(stderr, "skipping predicted mRNA with unknown orientation "
                  "(line %u)\n", gt_genome_node_get_line_number(gn));
        }
    }
  }
  else if (gt_genome_feature_has_type((GT_GenomeFeature*) gn,
                                   gft_LTR_retrotransposon)) {
    /* store predicted LTR */
    evaluator_add_predicted(info->LTR_evaluator, 1);
    /* determine true LTR */
    bsearch_all_mark(real_genome_nodes, &gn,
                     gt_array_get_space(info->slot->LTRs),
                     gt_array_size(info->slot->LTRs), sizeof (GT_GenomeNode*),
                     (GT_CompareWithData) gt_genome_node_compare_delta,
                     &info->LTRdelta, info->slot->overlapped_LTRs);

    if (gt_array_size(real_genome_nodes)) {
      for (i = 0; i < gt_array_size(real_genome_nodes); i++) {
        real_gn = *(GT_GenomeNode***) gt_array_get(real_genome_nodes, i);
        num = real_gn - (GT_GenomeNode**) gt_array_get_space(info->slot->LTRs);
        if (!gt_bittab_bit_is_set(info->slot->true_LTRs, num)) {
          gt_bittab_set_bit(info->slot->true_LTRs, num);
          evaluator_add_true(info->LTR_evaluator);
          /*@loopbreak@*/
          break;
        }
      }
    }
    else {
      /* no LTR with the same range found -> check if this is a wrong LTR */
      if (!gt_genome_node_overlaps_nodes_mark(gn, info->slot->LTRs,
                                           info->slot->overlapped_LTRs)) {
        (*info->wrong_LTRs)++;
      }
    }
  }
  else if (gt_genome_feature_has_type((GT_GenomeFeature*) gn, gft_exon)) {
    /* store predicted exon (mRNA level)*/
    store_predicted_exon(info->mRNA_exon_evaluators, gn);

    /* store predicted exon (mRNA level, collapsed) */
    store_predicted_exon_collapsed(predicted_strand == GT_STRAND_FORWARD
                                   ? info->slot->used_mRNA_exons_forward
                                   : info->slot->used_mRNA_exons_reverse,
                                   &predicted_range,
                                   info->mRNA_exon_evaluators_collapsed, gn);

    /* determine true exon (mRNA level)*/
    switch (predicted_strand) {
      case GT_STRAND_FORWARD:
      case GT_STRAND_REVERSE:
        store_true_exon(gn, predicted_strand, &predicted_range,
                        info->exondiff,
                        info->slot->mRNA_exons_forward,
                        info->slot->mRNA_exons_reverse,
                        info->slot->mRNA_counts_forward,
                        info->slot->mRNA_counts_reverse,
                        info->slot->mRNA_exon_bittabs_forward,
                        info->slot->mRNA_exon_bittabs_reverse,
                        info->mRNA_exon_evaluators,
                        info->mRNA_exon_evaluators_collapsed);
        /* nucleotide level */
        if (info->nuceval) {
          add_nucleotide_exon(predicted_strand == GT_STRAND_FORWARD
                              ? info->slot->pred_mRNA_nucleotides_forward
                              : info->slot->pred_mRNA_nucleotides_reverse,
                              predicted_range, info->slot->real_range,
                              predicted_strand == GT_STRAND_FORWARD
                              ? &info->slot->FP_mRNA_nucleotides_forward
                              : &info->slot->FP_mRNA_nucleotides_reverse);
        }
        break;
      default:
        if (info->verbose) {
          fprintf(stderr, "skipping predicted exon with unknown orientation "
                  "(line %u)\n", gt_genome_node_get_line_number(gn));
        }
    }
  }
  else if (gt_genome_feature_has_type((GT_GenomeFeature*) gn, gft_CDS)) {
    /* store predicted exon (CDS level)*/
    store_predicted_exon(info->CDS_exon_evaluators, gn);

    /* store predicted exon (CDS level, collapsed) */
    store_predicted_exon_collapsed(predicted_strand == GT_STRAND_FORWARD
                                   ? info->slot->used_CDS_exons_forward
                                   : info->slot->used_CDS_exons_reverse,
                                   &predicted_range,
                                   info->CDS_exon_evaluators_collapsed, gn);

    /* determine true exon (CDS level) */
    switch (predicted_strand) {
      case GT_STRAND_FORWARD:
      case GT_STRAND_REVERSE:
        store_true_exon(gn, predicted_strand, &predicted_range,
                        info->exondiff,
                        info->slot->CDS_exons_forward,
                        info->slot->CDS_exons_reverse,
                        info->slot->CDS_counts_forward,
                        info->slot->CDS_counts_reverse,
                        info->slot->CDS_exon_bittabs_forward,
                        info->slot->CDS_exon_bittabs_reverse,
                        info->CDS_exon_evaluators,
                        info->CDS_exon_evaluators_collapsed);
        /* nucleotide level */
        if (info->nuceval) {
          add_nucleotide_exon(predicted_strand == GT_STRAND_FORWARD
                              ? info->slot->pred_CDS_nucleotides_forward
                              : info->slot->pred_CDS_nucleotides_reverse,
                              predicted_range, info->slot->real_range,
                              predicted_strand == GT_STRAND_FORWARD
                              ? &info->slot->FP_CDS_nucleotides_forward
                              : &info->slot->FP_CDS_nucleotides_reverse);
        }
        break;
      default:
        if (info->verbose) {
          fprintf(stderr, "skipping predicted exon with unknown orientation "
                  "(line %u)\n", gt_genome_node_get_line_number(gn));
        }
      }
  }
  gt_array_delete(real_genome_nodes);
  return 0;
}

int determine_missing_features(GT_UNUSED void *key, void *value, void *data,
                               GT_UNUSED GT_Error *err)
{
  StreamEvaluator *se = (StreamEvaluator*) data;
  Slot *slot = (Slot*) value;
  gt_error_check(err);
  assert(key && value && data);
  if (slot->overlapped_genes_forward) {
    se->missing_genes += gt_bittab_size(slot->overlapped_genes_forward) -
                         gt_bittab_count_set_bits(slot->overlapped_genes_forward);
  }
  if (slot->overlapped_genes_reverse) {
    se->missing_genes += gt_bittab_size(slot->overlapped_genes_reverse) -
                         gt_bittab_count_set_bits(slot->overlapped_genes_reverse);
  }
  if (slot->overlapped_mRNAs_forward) {
    se->missing_mRNAs += gt_bittab_size(slot->overlapped_mRNAs_forward) -
                         gt_bittab_count_set_bits(slot->overlapped_mRNAs_forward);
  }
  if (slot->overlapped_mRNAs_reverse) {
    se->missing_mRNAs += gt_bittab_size(slot->overlapped_mRNAs_reverse) -
                         gt_bittab_count_set_bits(slot->overlapped_mRNAs_reverse);
  }
  if (slot->overlapped_LTRs) {
    se->missing_LTRs  += gt_bittab_size(slot->overlapped_LTRs) -
                         gt_bittab_count_set_bits(slot->overlapped_LTRs);
  }
  return 0;
}

static void add_nucleotide_values(NucEval *nucleotides, GT_Bittab *real,
                                  GT_Bittab *pred, GT_Bittab *tmp, const char *level)
{
  assert(nucleotides && real && pred && tmp);
  if (log_enabled()) {
    log_log(level);
    log_log("reality:");
    gt_bittab_show(real, log_fp());
    log_log("prediction:");
    gt_bittab_show(pred, log_fp());
  }
  /* real & pred = TP */
  gt_bittab_and(tmp, real, pred);
  nucleotides->TP += gt_bittab_count_set_bits(tmp);
  /* ~real & pred = FP */;
  gt_bittab_complement(tmp, real);
  gt_bittab_and_equal(tmp, pred);
  nucleotides->FP += gt_bittab_count_set_bits(tmp);
  /* real & ~pred = FN */
  gt_bittab_complement(tmp, pred);
  gt_bittab_and_equal(tmp, real);
  nucleotides->FN += gt_bittab_count_set_bits(tmp);
}

int compute_nucleotides_values(GT_UNUSED void *key, void *value, void *data,
                               GT_UNUSED GT_Error *err)
{
  StreamEvaluator *se = (StreamEvaluator*) data;
  Slot *slot = (Slot*) value;
  GT_Bittab *tmp;
  gt_error_check(err);
  assert(key && value && data);
  /* add ``out of range'' FPs */
  se->mRNA_nucleotides.FP += slot->FP_mRNA_nucleotides_forward;
  se->mRNA_nucleotides.FP += slot->FP_mRNA_nucleotides_reverse;
  se->CDS_nucleotides.FP  += slot->FP_CDS_nucleotides_forward;
  se->CDS_nucleotides.FP  += slot->FP_CDS_nucleotides_reverse;
  /* add other values */
  tmp = gt_bittab_new(gt_range_length(slot->real_range));
  add_nucleotide_values(&se->mRNA_nucleotides,
                        slot->real_mRNA_nucleotides_forward,
                        slot->pred_mRNA_nucleotides_forward, tmp,
                        "mRNA forward");
  add_nucleotide_values(&se->mRNA_nucleotides,
                        slot->real_mRNA_nucleotides_reverse,
                        slot->pred_mRNA_nucleotides_reverse, tmp,
                        "mRNA reverse");
  add_nucleotide_values(&se->CDS_nucleotides,
                        slot->real_CDS_nucleotides_forward,
                        slot->pred_CDS_nucleotides_forward, tmp,
                        "CDS forward");
  add_nucleotide_values(&se->CDS_nucleotides,
                        slot->real_CDS_nucleotides_reverse,
                        slot->pred_CDS_nucleotides_reverse, tmp,
                        "CDS reverse");
  gt_bittab_delete(tmp);
  return 0;
}

int stream_evaluator_evaluate(StreamEvaluator *se, bool verbose, bool exondiff,
                              GenomeVisitor *gv, GT_Error *err)
{
  GT_GenomeNode *gn;
  GT_SequenceRegion *sr;
  GT_GenomeFeature *gf;
  Slot *slot;
  ProcessRealFeatureInfo real_info;
  ProcessPredictedFeatureInfo predicted_info;
  int had_err;

  gt_error_check(err);
  assert(se);

  /* init */
  real_info.nuceval = se->nuceval;
  real_info.verbose = verbose;
  predicted_info.nuceval = se->nuceval;
  predicted_info.verbose = verbose;
  predicted_info.exondiff = exondiff;
  predicted_info.LTRdelta = se->LTRdelta;
  predicted_info.gene_evaluator = se->gene_evaluator;
  predicted_info.mRNA_evaluator = se->mRNA_evaluator;
  predicted_info.LTR_evaluator  = se->LTR_evaluator;
  predicted_info.mRNA_exon_evaluators = se->mRNA_exon_evaluators;
  predicted_info.mRNA_exon_evaluators_collapsed =
    se->mRNA_exon_evaluators_collapsed;
  predicted_info.CDS_exon_evaluators = se->CDS_exon_evaluators;
  predicted_info.CDS_exon_evaluators_collapsed =
    se->CDS_exon_evaluators_collapsed;
  predicted_info.wrong_genes = &se->wrong_genes;
  predicted_info.wrong_mRNAs = &se->wrong_mRNAs;
  predicted_info.wrong_LTRs  = &se->wrong_LTRs;

  /* process the reality stream completely */
  while (!(had_err = genome_stream_next_tree(se->reality, &gn, err)) && gn) {
    sr = gt_genome_node_cast(gt_sequence_region_class(), gn);
    if (sr) {
      /* each sequence region gets its own ``slot'' */
      if (!(slot = hashmap_get(se->slots, gt_str_get(gt_genome_node_get_seqid(gn)))))
      {
        slot = slot_new(se->nuceval, gt_genome_node_get_range(gn));
        hashmap_add(se->slots,
                    cstr_dup(gt_str_get(gt_genome_node_get_seqid(gn))), slot);
      }
      assert(slot);
    }
    gf = gt_genome_node_cast(gt_genome_feature_class(), gn);
    /* we consider only genome features */
    if (gf) {
      /* each sequence must have its own ``slot'' at this point */
      slot = hashmap_get(se->slots, gt_str_get(gt_genome_node_get_seqid(gn)));
      assert(slot);
      /* store the exons */
      real_info.slot = slot;
      gt_genome_feature_determine_transcripttypes(gf);
      had_err = gt_genome_node_traverse_children(gn, &real_info,
                                              process_real_feature, false,
                                              NULL);
      assert(!had_err); /* cannot happen, process_real_feature() is sane */
    }
    if (gv)
      gt_genome_node_accept(gn, gv, err);
    gt_genome_node_rec_delete(gn);
  }

  /* set the actuals and sort them */
  if (!had_err) {
    had_err = hashmap_foreach(se->slots, set_actuals_and_sort_them, se, NULL);
    assert(!had_err); /* set_actuals_and_sort_them() is sane */
  }

  /* process the prediction stream */
  if (!had_err) {
    while (!(had_err = genome_stream_next_tree(se->prediction, &gn, err)) &&
           gn) {
      gf = gt_genome_node_cast(gt_genome_feature_class(), gn);
      /* we consider only genome features */
      if (gf) {
        /* get (real) slot */
        slot = hashmap_get(se->slots, gt_str_get(gt_genome_node_get_seqid(gn)));
        if (slot) {
          predicted_info.slot = slot;
          gt_genome_feature_determine_transcripttypes(gf);
          had_err = gt_genome_node_traverse_children(gn, &predicted_info,
                                                  process_predicted_feature,
                                                  false, NULL);
          assert(!had_err); /* cannot happen, process_predicted_feature() is
                               sane */
        }
        else {
          /* we got no (real) slot */
          warning("sequence id \"%s\" (with predictions) not given in "
                  "``reality''", gt_str_get(gt_genome_node_get_seqid(gn)));
        }
      }
      if (gv)
        had_err = gt_genome_node_accept(gn, gv, err);
      gt_genome_node_rec_delete(gn);
    }
  }

  /* determine the missing mRNAs */
  if (!had_err) {
    had_err = hashmap_foreach(se->slots, determine_missing_features, se, NULL);
    assert(!had_err); /* determine_missing_features() is sane */
  }

  /* compute the nucleotides values */
  if (!had_err && se->nuceval) {
    had_err = hashmap_foreach(se->slots, compute_nucleotides_values, se, NULL);
    assert(!had_err); /* compute_nucleotides_values() is sane */
  }

  return had_err;
}

static void show_transcript_values(TranscriptEvaluators *te, const char *level,
                                   const char *additional_info, FILE *outfp)
{
  assert(te);

  fprintf(outfp, "exon sensitivity (%s level, all%s): ", level,
          additional_info);
  evaluator_show_sensitivity(transcript_evaluators_get_all(te), outfp);
  xfputc('\n', outfp);

  fprintf(outfp, "exon specificity (%s level, all%s): ", level,
          additional_info);
  evaluator_show_specificity(transcript_evaluators_get_all(te), outfp);
  xfputc('\n', outfp);

  fprintf(outfp, "exon sensitivity (%s level, single%s): ", level,
          additional_info);
  evaluator_show_sensitivity(transcript_evaluators_get_single(te), outfp);
  xfputc('\n', outfp);

  fprintf(outfp, "exon specificity (%s level, single%s): ", level,
          additional_info );
  evaluator_show_specificity(transcript_evaluators_get_single(te), outfp);
  xfputc('\n', outfp);

  fprintf(outfp, "exon sensitivity (%s level, initial%s): ", level,
          additional_info);
  evaluator_show_sensitivity(transcript_evaluators_get_initial(te), outfp);
  xfputc('\n', outfp);

  fprintf(outfp, "exon specificity (%s level, initial%s): ", level,
          additional_info);
  evaluator_show_specificity(transcript_evaluators_get_initial(te), outfp);
  xfputc('\n', outfp);

  fprintf(outfp, "exon sensitivity (%s level, internal%s): ", level,
          additional_info);
  evaluator_show_sensitivity(transcript_evaluators_get_internal(te), outfp);
  xfputc('\n', outfp);

  fprintf(outfp, "exon specificity (%s level, internal%s): ", level,
          additional_info);
  evaluator_show_specificity(transcript_evaluators_get_internal(te), outfp);
  xfputc('\n', outfp);

  fprintf(outfp, "exon sensitivity (%s level, terminal%s): ", level,
          additional_info);
  evaluator_show_sensitivity(transcript_evaluators_get_terminal(te), outfp);
  xfputc('\n', outfp);

  fprintf(outfp, "exon specificity (%s level, terminal%s): ", level,
          additional_info);
  evaluator_show_specificity(transcript_evaluators_get_terminal(te), outfp);
  xfputc('\n', outfp);
}

static void show_nucleotide_values(NucEval *nucleotides, const char *level,
                                   FILE *outfp)
{
  double sensitivity = 1.0, specificity = 1.0;
  assert(nucleotides && level);
  if (nucleotides->TP || nucleotides->FN) {
    sensitivity = (double) nucleotides->TP /
                  (nucleotides->TP + nucleotides->FN);
  }
  if (nucleotides->TP || nucleotides->FP) {
    specificity = (double) nucleotides->TP /
                  (nucleotides->TP + nucleotides->FP);
  }
  fprintf(outfp, "nucleotide sensitivity (%s level): %6.2f%% "
          "(TP=%lu/(TP=%lu + FN=%lu))\n", level, sensitivity * 100.0,
          nucleotides->TP, nucleotides->TP, nucleotides->FN);
  fprintf(outfp, "nucleotide specificity (%s level): %6.2f%% "
          "(TP=%lu/(TP=%lu + FP=%lu))\n", level, specificity * 100.0,
          nucleotides->TP, nucleotides->TP, nucleotides->FP);
}

void stream_evaluator_show(StreamEvaluator *se, FILE *outfp)
{
  assert(se);

  if (!se->evalLTR) {
    /* gene level */
    fprintf(outfp, "gene sensitivity:              ");
    evaluator_show_sensitivity(se->gene_evaluator, outfp);
    fprintf(outfp, " (missing genes: %lu)\n", se->missing_genes);

    fprintf(outfp, "gene specificity:              ");
    evaluator_show_specificity(se->gene_evaluator, outfp);
    fprintf(outfp, " (wrong genes: %lu)\n", se->wrong_genes);

    /* mRNA level */
    fprintf(outfp, "mRNA sensitivity:              ");
    evaluator_show_sensitivity(se->mRNA_evaluator, outfp);
    fprintf(outfp, " (missing mRNAs: %lu)\n", se->missing_mRNAs);

    fprintf(outfp, "mRNA specificity:              ");
    evaluator_show_specificity(se->mRNA_evaluator, outfp);
    fprintf(outfp, " (wrong mRNAs: %lu)\n", se->wrong_mRNAs);

    /* mRNA exon level */
    show_transcript_values(se->mRNA_exon_evaluators, "mRNA", "", outfp);
    show_transcript_values(se->mRNA_exon_evaluators_collapsed, "mRNA",
                           ", collapsed", outfp);

    /* CDS exon level */
    show_transcript_values(se->CDS_exon_evaluators, "CDS", "", outfp);
    show_transcript_values(se->CDS_exon_evaluators_collapsed, "CDS",
                           ", collapsed", outfp);

    if (se->nuceval) {
      /* mRNA nucleotide level */
      show_nucleotide_values(&se->mRNA_nucleotides, "mRNA", outfp);
      /* CDS nucleotide level */
      show_nucleotide_values(&se->CDS_nucleotides, "CDS", outfp);
    }
  }
  else {
    /* LTR_retrotransposon prediction */
    fprintf(outfp, "LTR_retrotransposon sensitivity: ");
    evaluator_show_sensitivity(se->LTR_evaluator, outfp);
    fprintf(outfp, " (missing LTRs: %lu)\n", se->missing_LTRs);

    fprintf(outfp, "LTR_retrotransposon specificity: ");
    evaluator_show_specificity(se->LTR_evaluator, outfp);
    fprintf(outfp, " (wrong LTRs: %lu)\n", se->wrong_LTRs);
  }
}

void stream_evaluator_delete(StreamEvaluator *se)
{
  if (!se) return;
  genome_stream_delete(se->reality);
  genome_stream_delete(se->prediction);
  hashmap_delete(se->slots);
  evaluator_delete(se->gene_evaluator);
  evaluator_delete(se->mRNA_evaluator);
  evaluator_delete(se->LTR_evaluator);
  transcript_evaluators_delete(se->mRNA_exon_evaluators);
  transcript_evaluators_delete(se->mRNA_exon_evaluators_collapsed);
  transcript_evaluators_delete(se->CDS_exon_evaluators);
  transcript_evaluators_delete(se->CDS_exon_evaluators_collapsed);
  gt_free(se);
}
