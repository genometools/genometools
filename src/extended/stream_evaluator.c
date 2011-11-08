/*
  Copyright (c) 2006-2011 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#include "core/assert_api.h"
#include "core/bsearch.h"
#include "core/cstr_api.h"
#include "core/hashmap.h"
#include "core/log.h"
#include "core/ma.h"
#include "core/unused_api.h"
#include "core/warning_api.h"
#include "core/xansi_api.h"
#include "extended/evaluator.h"
#include "extended/gff3_output.h"
#include "extended/stream_evaluator.h"
#include "extended/transcript_evaluators.h"
#include "extended/transcript_exons.h"
#include "extended/transcript_used_exons.h"

typedef struct {
  unsigned long TP, FP, FN;
} NucEval;

struct GtStreamEvaluator {
  GtNodeStream *reference,
               *prediction;
  bool nuceval, evalLTR;
  unsigned long LTRdelta;
  GtHashmap *slots; /* sequence id -> slot */
  GtEvaluator *mRNA_gene_evaluator,
              *CDS_gene_evaluator,
              *mRNA_mRNA_evaluator,
              *CDS_mRNA_evaluator,
              *LTR_evaluator;
  GtTranscriptEvaluators *mRNA_exon_evaluators,
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
  GtArray *genes_forward,
          *genes_reverse,
          *mRNAs_forward,
          *mRNAs_reverse,
          *LTRs;
  GtTranscriptExons *mRNA_exons_forward,
                    *mRNA_exons_reverse,
                    *CDS_exons_forward,
                    *CDS_exons_reverse;
  GtTranscriptCounts *mRNA_counts_forward,
                     *mRNA_counts_reverse,
                     *CDS_counts_forward,
                     *CDS_counts_reverse;
  GtRange real_range;
  unsigned long FP_mRNA_nucleotides_forward,
                FP_mRNA_nucleotides_reverse,
                FP_CDS_nucleotides_forward,
                FP_CDS_nucleotides_reverse;
  GtBittab *real_mRNA_nucleotides_forward,
           *pred_mRNA_nucleotides_forward,
           *real_mRNA_nucleotides_reverse,
           *pred_mRNA_nucleotides_reverse,
           *real_CDS_nucleotides_forward,
           *pred_CDS_nucleotides_forward,
           *real_CDS_nucleotides_reverse,
           *pred_CDS_nucleotides_reverse,
           *true_mRNA_genes_forward,
           *true_mRNA_genes_reverse,
           *true_CDS_genes_forward,
           *true_CDS_genes_reverse,
           *true_mRNA_mRNAs_forward,
           *true_mRNA_mRNAs_reverse,
           *true_CDS_mRNAs_forward,
           *true_CDS_mRNAs_reverse,
           *true_LTRs,
           *overlapped_genes_forward,
           *overlapped_genes_reverse,
           *overlapped_mRNAs_forward,
           *overlapped_mRNAs_reverse,
           *overlapped_LTRs;
  GtTranscriptBittabs *mRNA_exon_bittabs_forward,
                      *mRNA_exon_bittabs_reverse,
                      *CDS_exon_bittabs_forward,
                      *CDS_exon_bittabs_reverse;
  GtTranscriptUsedExons *used_mRNA_exons_forward,
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
       exondiff,
       exondiffcollapsed;
  unsigned long LTRdelta;
  GtEvaluator *mRNA_gene_evaluator,
              *CDS_gene_evaluator,
              *mRNA_mRNA_evaluator,
              *CDS_mRNA_evaluator,
              *LTR_evaluator;
  GtTranscriptEvaluators *mRNA_exon_evaluators,
                         *mRNA_exon_evaluators_collapsed,
                         *CDS_exon_evaluators,
                         *CDS_exon_evaluators_collapsed;
  unsigned long *wrong_genes,
                *wrong_mRNAs,
                *wrong_LTRs;
} ProcessPredictedFeatureInfo;

static Slot* slot_new(bool nuceval, GtRange range)
{
  unsigned long length;
  Slot *s = gt_calloc(1, sizeof (Slot));
  length = gt_range_length(&range);
  s->genes_forward = gt_array_new(sizeof (GtGenomeNode*));
  s->genes_reverse = gt_array_new(sizeof (GtGenomeNode*));
  s->mRNAs_forward = gt_array_new(sizeof (GtGenomeNode*));
  s->mRNAs_reverse = gt_array_new(sizeof (GtGenomeNode*));
  s->LTRs          = gt_array_new(sizeof (GtGenomeNode*));
  s->mRNA_exons_forward = gt_transcript_exons_new();
  s->mRNA_exons_reverse = gt_transcript_exons_new();
  s->CDS_exons_forward = gt_transcript_exons_new();
  s->CDS_exons_reverse = gt_transcript_exons_new();
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
  s->used_mRNA_exons_forward = gt_transcript_used_exons_new();
  s->used_mRNA_exons_reverse = gt_transcript_used_exons_new();
  s->used_CDS_exons_forward = gt_transcript_used_exons_new();
  s->used_CDS_exons_reverse = gt_transcript_used_exons_new();
  return s;
}

static void slot_delete(Slot *s)
{
  unsigned long i;
  gt_assert(s);
  for (i = 0; i < gt_array_size(s->genes_forward); i++)
    gt_genome_node_delete(*(GtGenomeNode**) gt_array_get(s->genes_forward, i));
  gt_array_delete(s->genes_forward);
  for (i = 0; i < gt_array_size(s->genes_reverse); i++)
    gt_genome_node_delete(*(GtGenomeNode**) gt_array_get(s->genes_reverse, i));
  gt_array_delete(s->genes_reverse);
  for (i = 0; i < gt_array_size(s->mRNAs_forward); i++)
    gt_genome_node_delete(*(GtGenomeNode**) gt_array_get(s->mRNAs_forward, i));
  gt_array_delete(s->mRNAs_forward);
  for (i = 0; i < gt_array_size(s->mRNAs_reverse); i++)
    gt_genome_node_delete(*(GtGenomeNode**) gt_array_get(s->mRNAs_reverse, i));
  gt_array_delete(s->mRNAs_reverse);
  for (i = 0; i < gt_array_size(s->LTRs); i++)
    gt_genome_node_delete(*(GtGenomeNode**) gt_array_get(s->LTRs, i));
  gt_array_delete(s->LTRs);
  gt_transcript_exons_delete(s->mRNA_exons_forward);
  gt_transcript_exons_delete(s->mRNA_exons_reverse);
  gt_transcript_exons_delete(s->CDS_exons_forward);
  gt_transcript_exons_delete(s->CDS_exons_reverse);
  gt_transcript_counts_delete(s->mRNA_counts_forward);
  gt_transcript_counts_delete(s->mRNA_counts_reverse);
  gt_transcript_counts_delete(s->CDS_counts_forward);
  gt_transcript_counts_delete(s->CDS_counts_reverse);
  gt_bittab_delete(s->real_mRNA_nucleotides_forward);
  gt_bittab_delete(s->pred_mRNA_nucleotides_forward);
  gt_bittab_delete(s->real_mRNA_nucleotides_reverse);
  gt_bittab_delete(s->pred_mRNA_nucleotides_reverse);
  gt_bittab_delete(s->real_CDS_nucleotides_forward);
  gt_bittab_delete(s->pred_CDS_nucleotides_forward);
  gt_bittab_delete(s->real_CDS_nucleotides_reverse);
  gt_bittab_delete(s->pred_CDS_nucleotides_reverse);
  gt_bittab_delete(s->true_mRNA_genes_forward);
  gt_bittab_delete(s->true_mRNA_genes_reverse);
  gt_bittab_delete(s->true_CDS_genes_forward);
  gt_bittab_delete(s->true_CDS_genes_reverse);
  gt_bittab_delete(s->true_mRNA_mRNAs_forward);
  gt_bittab_delete(s->true_mRNA_mRNAs_reverse);
  gt_bittab_delete(s->true_CDS_mRNAs_forward);
  gt_bittab_delete(s->true_CDS_mRNAs_reverse);
  gt_bittab_delete(s->true_LTRs);
  gt_bittab_delete(s->overlapped_genes_forward);
  gt_bittab_delete(s->overlapped_genes_reverse);
  gt_bittab_delete(s->overlapped_mRNAs_forward);
  gt_bittab_delete(s->overlapped_mRNAs_reverse);
  gt_bittab_delete(s->overlapped_LTRs);
  gt_transcript_bittabs_delete(s->mRNA_exon_bittabs_forward);
  gt_transcript_bittabs_delete(s->mRNA_exon_bittabs_reverse);
  gt_transcript_bittabs_delete(s->CDS_exon_bittabs_forward);
  gt_transcript_bittabs_delete(s->CDS_exon_bittabs_reverse);
  gt_transcript_used_exons_delete(s->used_mRNA_exons_forward);
  gt_transcript_used_exons_delete(s->used_mRNA_exons_reverse);
  gt_transcript_used_exons_delete(s->used_CDS_exons_forward);
  gt_transcript_used_exons_delete(s->used_CDS_exons_reverse);
  gt_free(s);
}

GtStreamEvaluator* gt_stream_evaluator_new(GtNodeStream *reference,
                                           GtNodeStream *prediction,
                                           bool nuceval, bool evalLTR,
                                           unsigned long LTRdelta)
{
  GtStreamEvaluator *evaluator = gt_calloc(1, sizeof (GtStreamEvaluator));
  evaluator->reference = gt_node_stream_ref(reference);
  evaluator->prediction = gt_node_stream_ref(prediction);
  evaluator->nuceval = nuceval;
  evaluator->evalLTR = evalLTR;
  evaluator->LTRdelta = LTRdelta;
  evaluator->slots = gt_hashmap_new(GT_HASH_STRING, gt_free_func,
                                    (GtFree) slot_delete);
  evaluator->mRNA_gene_evaluator = gt_evaluator_new();
  evaluator->CDS_gene_evaluator = gt_evaluator_new();
  evaluator->mRNA_mRNA_evaluator = gt_evaluator_new();
  evaluator->CDS_mRNA_evaluator = gt_evaluator_new();
  evaluator->LTR_evaluator = gt_evaluator_new();
  evaluator->mRNA_exon_evaluators = gt_transcript_evaluators_new();
  evaluator->mRNA_exon_evaluators_collapsed = gt_transcript_evaluators_new();
  evaluator->CDS_exon_evaluators = gt_transcript_evaluators_new();
  evaluator->CDS_exon_evaluators_collapsed = gt_transcript_evaluators_new();
  return evaluator;
}

static int set_actuals_and_sort_them(GT_UNUSED void *key, void *value,
                                     void *data, GT_UNUSED GtError *err)
{
  GtStreamEvaluator *se = (GtStreamEvaluator*) data;
  Slot *s = (Slot*) value;

  gt_error_check(err);
  gt_assert(key && value && data);

  /* set actual genes */
  gt_evaluator_add_actual(se->mRNA_gene_evaluator,
                          gt_array_size(s->genes_forward));
  gt_evaluator_add_actual(se->mRNA_gene_evaluator,
                          gt_array_size(s->genes_reverse));
  gt_evaluator_add_actual(se->CDS_gene_evaluator,
                          gt_array_size(s->genes_forward));
  gt_evaluator_add_actual(se->CDS_gene_evaluator,
                          gt_array_size(s->genes_reverse));

  /* set actual mRNAs */
  gt_evaluator_add_actual(se->mRNA_mRNA_evaluator,
                          gt_array_size(s->mRNAs_forward));
  gt_evaluator_add_actual(se->mRNA_mRNA_evaluator,
                          gt_array_size(s->mRNAs_reverse));
  gt_evaluator_add_actual(se->CDS_mRNA_evaluator,
                          gt_array_size(s->mRNAs_forward));
  gt_evaluator_add_actual(se->CDS_mRNA_evaluator,
                          gt_array_size(s->mRNAs_reverse));

  /* set actual LTRs */
  gt_evaluator_add_actual(se->LTR_evaluator, gt_array_size(s->LTRs));

  /* set actual exons (before uniq!) */
  gt_transcript_evaluators_add_actuals(se->mRNA_exon_evaluators,
                                       s->mRNA_exons_forward);
  gt_transcript_evaluators_add_actuals(se->mRNA_exon_evaluators,
                                       s->mRNA_exons_reverse);
  gt_transcript_evaluators_add_actuals(se->CDS_exon_evaluators,
                                       s->CDS_exons_forward);
  gt_transcript_evaluators_add_actuals(se->CDS_exon_evaluators,
                                       s->CDS_exons_reverse);

  /* sort genes */
  gt_genome_nodes_sort(s->genes_forward);
  gt_genome_nodes_sort(s->genes_reverse);

  /* sort mRNAs */
  gt_genome_nodes_sort(s->mRNAs_forward);
  gt_genome_nodes_sort(s->mRNAs_reverse);

  /* sort LTRs */
  gt_genome_nodes_sort(s->LTRs);

  /* sort exons */
  gt_transcript_exons_sort(s->mRNA_exons_forward);
  gt_transcript_exons_sort(s->mRNA_exons_reverse);
  gt_transcript_exons_sort(s->CDS_exons_forward);
  gt_transcript_exons_sort(s->CDS_exons_reverse);

  /* determine true exons */
  s->mRNA_counts_forward =
    gt_transcript_exons_uniq_in_place_count(s->mRNA_exons_forward);
  s->mRNA_counts_reverse =
    gt_transcript_exons_uniq_in_place_count(s->mRNA_exons_reverse);
  s->CDS_counts_forward =
    gt_transcript_exons_uniq_in_place_count(s->CDS_exons_forward);
  s->CDS_counts_reverse =
    gt_transcript_exons_uniq_in_place_count(s->CDS_exons_reverse);

  /* set actual exons for the collapsed case (after uniq!) */
  gt_transcript_evaluators_add_actuals(se->mRNA_exon_evaluators_collapsed,
                                       s->mRNA_exons_forward);
  gt_transcript_evaluators_add_actuals(se->mRNA_exon_evaluators_collapsed,
                                       s->mRNA_exons_reverse);
  gt_transcript_evaluators_add_actuals(se->CDS_exon_evaluators_collapsed,
                                       s->CDS_exons_forward);
  gt_transcript_evaluators_add_actuals(se->CDS_exon_evaluators_collapsed,
                                       s->CDS_exons_reverse);

  /* make sure that the genes are sorted */
  gt_assert(gt_genome_nodes_are_sorted(s->genes_forward));
  gt_assert(gt_genome_nodes_are_sorted(s->genes_reverse));

  /* make sure that the mRNAs are sorted */
  gt_assert(gt_genome_nodes_are_sorted(s->mRNAs_forward));
  gt_assert(gt_genome_nodes_are_sorted(s->mRNAs_reverse));

  /* make sure that the LTRs are sorted */
  gt_assert(gt_genome_nodes_are_sorted(s->LTRs));

  /* make sure that the exons are sorted */
  gt_assert(gt_transcript_exons_are_sorted(s->mRNA_exons_forward));
  gt_assert(gt_transcript_exons_are_sorted(s->mRNA_exons_reverse));
  gt_assert(gt_transcript_exons_are_sorted(s->CDS_exons_forward));
  gt_assert(gt_transcript_exons_are_sorted(s->CDS_exons_reverse));

  /* init true bittabs */
  s->true_mRNA_genes_forward = gt_array_size(s->genes_forward)
                               ? gt_bittab_new(gt_array_size(s->genes_forward))
                               : NULL;
  s->true_mRNA_genes_reverse = gt_array_size(s->genes_reverse)
                               ? gt_bittab_new(gt_array_size(s->genes_reverse))
                               : NULL;
  s->true_CDS_genes_forward = gt_array_size(s->genes_forward)
                              ? gt_bittab_new(gt_array_size(s->genes_forward))
                              : NULL;
  s->true_CDS_genes_reverse = gt_array_size(s->genes_reverse)
                              ? gt_bittab_new(gt_array_size(s->genes_reverse))
                              : NULL;
  s->true_mRNA_mRNAs_forward = gt_array_size(s->mRNAs_forward)
                               ? gt_bittab_new(gt_array_size(s->mRNAs_forward))
                               : NULL;
  s->true_mRNA_mRNAs_reverse = gt_array_size(s->mRNAs_reverse)
                               ? gt_bittab_new(gt_array_size(s->mRNAs_reverse))
                               : NULL;
  s->true_CDS_mRNAs_forward = gt_array_size(s->mRNAs_forward)
                              ? gt_bittab_new(gt_array_size(s->mRNAs_forward))
                              : NULL;
  s->true_CDS_mRNAs_reverse = gt_array_size(s->mRNAs_reverse)
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
    gt_transcript_exons_create_bittabs(s->mRNA_exons_forward);
  s->mRNA_exon_bittabs_reverse =
    gt_transcript_exons_create_bittabs(s->mRNA_exons_reverse);
  s->CDS_exon_bittabs_forward =
    gt_transcript_exons_create_bittabs(s->CDS_exons_forward);
  s->CDS_exon_bittabs_reverse =
    gt_transcript_exons_create_bittabs(s->CDS_exons_reverse);

  return 0;
}

static void add_real_exon(GtTranscriptExons *te, GtRange range,
                          GtFeatureNode *fn)
{
  gt_assert(te);
  gt_array_add(gt_transcript_exons_get_all(te), range);
  switch (gt_feature_node_get_transcriptfeaturetype(fn)) {
    case TRANSCRIPT_FEATURE_TYPE_UNDETERMINED:
      gt_warning("type of feature (single, initial, internal, or terminal) "
                 "given on line %u in file \"%s\" could not be determined, "
                 "because the feature has no Parent attribute. Treating it as "
                 "single.",
                 gt_genome_node_get_line_number((GtGenomeNode*) fn),
                 gt_genome_node_get_filename((GtGenomeNode*) fn));
      /*@fallthrough@*/
    case TRANSCRIPT_FEATURE_TYPE_SINGLE:
      gt_array_add(gt_transcript_exons_get_single(te), range);
      break;
    case TRANSCRIPT_FEATURE_TYPE_INITIAL:
      gt_array_add(gt_transcript_exons_get_initial(te), range);
      break;
    case TRANSCRIPT_FEATURE_TYPE_INTERNAL:
      gt_array_add(gt_transcript_exons_get_internal(te), range);
      break;
    case TRANSCRIPT_FEATURE_TYPE_TERMINAL:
      gt_array_add(gt_transcript_exons_get_terminal(te), range);
      break;
  }
}

static void add_nucleotide_exon(GtBittab *nucleotides, GtRange range,
                                GtRange real_range, unsigned long *FP)
{
  unsigned long i;
  gt_assert(nucleotides);
  for (i = range.start; i <= range.end; i++) {
    if (gt_range_within(&real_range, i)) {
      gt_assert(i >= real_range.start);
      gt_bittab_set_bit(nucleotides, i - real_range.start);
    }
    else {
      gt_assert(FP);
      (*FP)++;
    }
  }
}

static int process_real_feature(GtFeatureNode *fn, void *data,
                                GT_UNUSED GtError *err)
{
  ProcessRealFeatureInfo *info = (ProcessRealFeatureInfo*) data;
  GtGenomeNode *gn_ref;
  GtRange range;

  gt_error_check(err);
  gt_assert(fn && data);

  if (gt_feature_node_has_type(fn, gt_ft_gene)) {
    switch (gt_feature_node_get_strand(fn)) {
      case GT_STRAND_FORWARD:
        gn_ref = gt_genome_node_ref((GtGenomeNode*) fn);
        gt_array_add(info->slot->genes_forward, gn_ref);
        break;
      case GT_STRAND_REVERSE:
        gn_ref = gt_genome_node_ref((GtGenomeNode*) fn);
        gt_array_add(info->slot->genes_reverse, gn_ref);
        break;
      default:
        if (info->verbose) {
          fprintf(stderr, "skipping real gene with unknown orientation "
                  "(line %u)\n",
                  gt_genome_node_get_line_number((GtGenomeNode*) fn));
        }
    }
  }
  else if (gt_feature_node_has_type(fn, gt_ft_mRNA)) {
    switch (gt_feature_node_get_strand(fn)) {
      case GT_STRAND_FORWARD:
        gn_ref = gt_genome_node_ref((GtGenomeNode*) fn);
        gt_array_add(info->slot->mRNAs_forward, gn_ref);
        break;
      case GT_STRAND_REVERSE:
        gn_ref = gt_genome_node_ref((GtGenomeNode*) fn);
        gt_array_add(info->slot->mRNAs_reverse, gn_ref);
        break;
      default:
        if (info->verbose) {
          fprintf(stderr, "skipping real mRNA with unknown orientation "
                  "(line %u)\n",
                  gt_genome_node_get_line_number((GtGenomeNode*) fn));
        }
    }
  }
  else if (gt_feature_node_has_type(fn, gt_ft_LTR_retrotransposon)) {
    gn_ref = gt_genome_node_ref((GtGenomeNode*) fn);
    gt_array_add(info->slot->LTRs, gn_ref);
  }
  else if (gt_feature_node_has_type(fn, gt_ft_CDS)) {
    range = gt_genome_node_get_range((GtGenomeNode*) fn);
    switch (gt_feature_node_get_strand(fn)) {
      case GT_STRAND_FORWARD:
        add_real_exon(info->slot->CDS_exons_forward, range, fn);
        /* nucleotide level */
        if (info->nuceval) {
          add_nucleotide_exon(info->slot->real_CDS_nucleotides_forward, range,
                              info->slot->real_range, NULL);
        }
        break;
      case GT_STRAND_REVERSE:
        add_real_exon(info->slot->CDS_exons_reverse, range, fn);
        /* nucleotide level */
        if (info->nuceval) {
          add_nucleotide_exon(info->slot->real_CDS_nucleotides_reverse, range,
                              info->slot->real_range, NULL);
        }
        break;
      default:
        if (info->verbose) {
          fprintf(stderr, "skipping real CDS exon with unknown orientation "
                  "(line %u)\n",
                  gt_genome_node_get_line_number((GtGenomeNode*) fn));
        }
    }
  }
  else if (gt_feature_node_has_type(fn, gt_ft_exon)) {
    range = gt_genome_node_get_range((GtGenomeNode*) fn);
    switch (gt_feature_node_get_strand(fn)) {
      case GT_STRAND_FORWARD:
        add_real_exon(info->slot->mRNA_exons_forward, range, fn);
        /* nucleotide level */
        if (info->nuceval) {
          add_nucleotide_exon(info->slot->real_mRNA_nucleotides_forward,
                              range, info->slot->real_range, NULL);
        }
        break;
      case GT_STRAND_REVERSE:
        add_real_exon(info->slot->mRNA_exons_reverse, range, fn);
        /* nucleotide level */
        if (info->nuceval) {
          add_nucleotide_exon(info->slot->real_mRNA_nucleotides_reverse,
                              range, info->slot->real_range, NULL);
        }
        break;
      default:
        if (info->verbose) {
          fprintf(stderr, "skipping real mRNA exon with unknown orientation "
                  "(line %u)\n",
                  gt_genome_node_get_line_number((GtGenomeNode*) fn));
        }
    }
  }
  return 0;
}

typedef struct {
  GtArray *exons;
  const char *feature_type;
} StoreExonFeatureInfo;

static int store_exon(GtFeatureNode *fn, void *data, GT_UNUSED GtError *err)
{
  StoreExonFeatureInfo *info = (StoreExonFeatureInfo*) data;
  GtRange range;
  gt_error_check(err);
  gt_assert(fn && info);
  if (gt_feature_node_has_type(fn, info->feature_type)) {
    range = gt_genome_node_get_range((GtGenomeNode*) fn);
    gt_array_add(info->exons, range);
  }
  return 0;
}

static bool mRNAs_are_equal(GtGenomeNode *gn_1, GtGenomeNode *gn_2,
                            const char *feature_type)
{
  GtArray *exons_1, *exons_2;
  StoreExonFeatureInfo info;
  bool equal;
  GT_UNUSED int had_err;

  gt_assert(gn_1 && gn_2 && feature_type);

  /* init */
  exons_1 = gt_array_new(sizeof (GtRange));
  exons_2 = gt_array_new(sizeof (GtRange));

  /* get exon ranges */
  info.exons = exons_1;
  info.feature_type = feature_type;
  had_err = gt_feature_node_traverse_children(gt_feature_node_cast(gn_1), &info,
                                              store_exon, false, NULL);
  gt_assert(!had_err); /* cannot happen, store_exon() is sane */
  info.exons = exons_2;
  had_err = gt_feature_node_traverse_children(gt_feature_node_cast(gn_2), &info,
                                              store_exon, false, NULL);
  gt_assert(!had_err); /* cannot happen, store_exon() is sane */

  /* sort exon ranges */
  gt_ranges_sort(exons_1);
  gt_ranges_sort(exons_2);

  /* compare exon ranges */
  equal = gt_ranges_are_equal(exons_1, exons_2);

  /* free */
  gt_array_delete(exons_1);
  gt_array_delete(exons_2);

  return equal;
}

typedef struct {
  GtArray *exons,
          *mRNAs;
  const char *feature_type;
} StoreGeneFeatureInfo;

static int store_gene_feature(GtFeatureNode *fn, void *data,
                              GT_UNUSED GtError *err)
{
  StoreGeneFeatureInfo *info = (StoreGeneFeatureInfo*) data;
  GtRange range;
  gt_error_check(err);
  gt_assert(fn && info);
  if (gt_feature_node_has_type(fn, gt_ft_mRNA))
    gt_array_add(info->mRNAs, fn);
  else if (gt_feature_node_has_type(fn, info->feature_type)) {
    range = gt_genome_node_get_range((GtGenomeNode*) fn);
    gt_array_add(info->exons, range);
  }
  return 0;
}

static bool genes_are_equal(GtGenomeNode *gn_1, GtGenomeNode *gn_2,
                            const char *feature_type)
{
  GtArray *exons_1, *exons_2, *mRNAs_1, *mRNAs_2;
  StoreGeneFeatureInfo info;
  unsigned long i;
  bool equal;
  GT_UNUSED int had_err;

  gt_assert(gn_1 && gn_2 && feature_type);

  /* init */
  exons_1 = gt_array_new(sizeof (GtRange));
  exons_2 = gt_array_new(sizeof (GtRange));
  mRNAs_1 = gt_array_new(sizeof (GtGenomeNode*));
  mRNAs_2 = gt_array_new(sizeof (GtGenomeNode*));

  /* get (direct) gene features */
  info.exons = exons_1;
  info.mRNAs = mRNAs_1;
  info.feature_type = feature_type;
  had_err = gt_feature_node_traverse_direct_children(gt_feature_node_cast(gn_1),
                                                     &info, store_gene_feature,
                                                     NULL);
  gt_assert(!had_err); /* cannot happen, store_gene_feature() is sane */
  info.exons = exons_2;
  info.mRNAs = mRNAs_2;
  had_err = gt_feature_node_traverse_direct_children(gt_feature_node_cast(gn_2),
                                                     &info, store_gene_feature,
                                                     NULL);
  gt_assert(!had_err); /* cannot happen, store_gene_feature() is sane */

  /* sort exon ranges */
  gt_ranges_sort(exons_1);
  gt_ranges_sort(exons_2);

  /* compare exon ranges */
  equal = gt_ranges_are_equal(exons_1, exons_2);

  /* compare mRNAs, if necessary */
  if (equal && gt_array_size(mRNAs_1) == gt_array_size(mRNAs_2)) {
    /* sort mRNAs */
    gt_genome_nodes_sort(mRNAs_1);
    gt_genome_nodes_sort(mRNAs_2);
    for (i = 0; i < gt_array_size(mRNAs_1); i++) {
      gt_assert(equal);
      equal = mRNAs_are_equal(*(GtGenomeNode**) gt_array_get(mRNAs_1, i),
                              *(GtGenomeNode**) gt_array_get(mRNAs_2, i),
                              feature_type);
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

static void store_predicted_exon(GtTranscriptEvaluators *te, GtFeatureNode *fn)
{
  gt_assert(te && fn);
  gt_evaluator_add_predicted(gt_transcript_evaluators_get_all(te), 1);
  switch (gt_feature_node_get_transcriptfeaturetype(fn)) {
    case TRANSCRIPT_FEATURE_TYPE_UNDETERMINED:
      gt_warning("type of feature (single, initial, internal, or terminal) "
                 "given on line %u in file \"%s\" could not be determined, "
                 "because the feature has no Parent attribute. Treating it as "
                 "single.",
                 gt_genome_node_get_line_number((GtGenomeNode*) fn),
                 gt_genome_node_get_filename((GtGenomeNode*) fn));
      /*@fallthrough@*/
    case TRANSCRIPT_FEATURE_TYPE_SINGLE:
      gt_evaluator_add_predicted(gt_transcript_evaluators_get_single(te), 1);
    break;
    case TRANSCRIPT_FEATURE_TYPE_INITIAL:
      gt_evaluator_add_predicted(gt_transcript_evaluators_get_initial(te), 1);
    break;
    case TRANSCRIPT_FEATURE_TYPE_INTERNAL:
      gt_evaluator_add_predicted(gt_transcript_evaluators_get_internal(te), 1);
    break;
    case TRANSCRIPT_FEATURE_TYPE_TERMINAL:
      gt_evaluator_add_predicted(gt_transcript_evaluators_get_terminal(te), 1);
    break;
  }
}

/* adds exon only if necessary */
static void add_predicted_collapsed(GtDlist *used_exons,
                                    GtRange *predicted_range,
                                    GtEvaluator *exon_evaluator_collapsed)
{
  GtRange *used_range;
  if (!gt_dlist_find(used_exons, predicted_range)) {
    used_range = gt_malloc(sizeof (GtRange));
    used_range->start = predicted_range->start;
    used_range->end = predicted_range->end;
    gt_dlist_add(used_exons, used_range);
    gt_evaluator_add_predicted(exon_evaluator_collapsed, 1);
  }
}

static void store_predicted_exon_collapsed(GtTranscriptUsedExons *used_exons,
                                           GtRange *predicted_range,
                                           GtTranscriptEvaluators *te,
                                           GtFeatureNode *fn)
{
  add_predicted_collapsed(gt_transcript_used_exons_get_all(used_exons),
                          predicted_range,
                          gt_transcript_evaluators_get_all(te));
  switch (gt_feature_node_get_transcriptfeaturetype(fn)) {
    case TRANSCRIPT_FEATURE_TYPE_UNDETERMINED:
      /* we do not show a warning here, because store_predicted_exon() has been
         called before and already shown one */
    case TRANSCRIPT_FEATURE_TYPE_SINGLE:
      add_predicted_collapsed(gt_transcript_used_exons_get_single(used_exons),
                              predicted_range,
                              gt_transcript_evaluators_get_single(te));
      break;
    case TRANSCRIPT_FEATURE_TYPE_INITIAL:
      add_predicted_collapsed(gt_transcript_used_exons_get_initial(used_exons),
                              predicted_range,
                              gt_transcript_evaluators_get_initial(te));
      break;
    case TRANSCRIPT_FEATURE_TYPE_INTERNAL:
      add_predicted_collapsed(gt_transcript_used_exons_get_internal(used_exons),
                              predicted_range,
                              gt_transcript_evaluators_get_internal(te));
      break;
    case TRANSCRIPT_FEATURE_TYPE_TERMINAL:
      add_predicted_collapsed(gt_transcript_used_exons_get_terminal(used_exons),
                              predicted_range,
                              gt_transcript_evaluators_get_terminal(te));
      break;
  }
}

static void mark_and_show_false_exon(GtFeatureNode *fn, bool exondiff)
{
  gt_feature_node_mark(fn); /* mark false exons */
  if (exondiff) {
    gt_gff3_output_leading(fn, NULL);
    printf(".\n");
  }
}

static void determine_true_exon(GtFeatureNode *fn, GtStrand predicted_strand,
                                bool exondiff, bool exondiffcollapsed,
                                GtRange *predicted_range,
                                GtArray *exons_forward,
                                GtArray *exons_reverse,
                                GtArray *true_exons_forward,
                                GtArray *true_exons_reverse,
                                GtBittab *true_exons_forward_collapsed,
                                GtBittab *true_exons_reverse_collapsed,
                                GtEvaluator *exon_evaluator,
                                GtEvaluator *exon_evaluator_collapsed)
{
  GtRange *actual_range;
  unsigned long num, *ctr_ptr;

  if ((actual_range = bsearch(predicted_range,
                              predicted_strand == GT_STRAND_FORWARD
                              ? gt_array_get_space(exons_forward)
                              : gt_array_get_space(exons_reverse),
                              predicted_strand == GT_STRAND_FORWARD
                              ? gt_array_size(exons_forward)
                              : gt_array_size(exons_reverse), sizeof (GtRange),
                              (GtCompare) gt_range_compare))) {
    if (predicted_strand == GT_STRAND_FORWARD) {
      num = actual_range - (GtRange*) gt_array_get_space(exons_forward);
      ctr_ptr = gt_array_get(true_exons_forward, num);
      if (*ctr_ptr) {
        (*ctr_ptr)--;
        gt_evaluator_add_true(exon_evaluator);
      }
      else
        mark_and_show_false_exon(fn, exondiff);
      if (true_exons_forward_collapsed &&
          !gt_bittab_bit_is_set(true_exons_forward_collapsed, num)) {
        gt_bittab_set_bit(true_exons_forward_collapsed, num);
        gt_evaluator_add_true(exon_evaluator_collapsed);
      }
    }
    else {
      num = actual_range - (GtRange*) gt_array_get_space(exons_reverse);
      ctr_ptr = gt_array_get(true_exons_reverse, num);
      if (*ctr_ptr) {
        (*ctr_ptr)--;
        gt_evaluator_add_true(exon_evaluator);
      }
      else
        mark_and_show_false_exon(fn, exondiff);
      if (true_exons_reverse_collapsed &&
          !gt_bittab_bit_is_set(true_exons_reverse_collapsed, num)) {
        gt_bittab_set_bit(true_exons_reverse_collapsed, num);
        gt_evaluator_add_true(exon_evaluator_collapsed);
      }
    }
  }
  else {
    mark_and_show_false_exon(fn, exondiff);
    if (exondiffcollapsed) {
      gt_gff3_output_leading(fn, NULL);
      printf(".\n");
    }
  }
}

static void store_true_exon(GtFeatureNode *fn, GtStrand predicted_strand,
                            GtRange *predicted_range, bool exondiff,
                            bool exondiffcollapsed,
                            GtTranscriptExons *exons_forward,
                            GtTranscriptExons *exons_reverse,
                            GtTranscriptCounts *counts_forward,
                            GtTranscriptCounts *counts_reverse,
                            GtTranscriptBittabs *exon_bittabs_forward,
                            GtTranscriptBittabs *exon_bittabs_reverse,
                            GtTranscriptEvaluators *exon_evaluators,
                            GtTranscriptEvaluators *exon_evaluators_collapsed)
{
  gt_assert(fn && predicted_range && exons_forward && exons_reverse);
  determine_true_exon(fn, predicted_strand, exondiff, exondiffcollapsed,
                      predicted_range,
                      gt_transcript_exons_get_all(exons_forward),
                      gt_transcript_exons_get_all(exons_reverse),
                      gt_transcript_counts_get_all(counts_forward),
                      gt_transcript_counts_get_all(counts_reverse),
                      gt_transcript_bittabs_get_all(exon_bittabs_forward),
                      gt_transcript_bittabs_get_all(exon_bittabs_reverse),
                      gt_transcript_evaluators_get_all(exon_evaluators),
                      gt_transcript_evaluators_get_all(
                                                    exon_evaluators_collapsed));
  switch (gt_feature_node_get_transcriptfeaturetype(fn)) {
    case TRANSCRIPT_FEATURE_TYPE_UNDETERMINED:
    case TRANSCRIPT_FEATURE_TYPE_SINGLE:
      determine_true_exon(fn, predicted_strand, exondiff, exondiffcollapsed,
                          predicted_range,
                          gt_transcript_exons_get_single(exons_forward),
                          gt_transcript_exons_get_single(exons_reverse),
                          gt_transcript_counts_get_single(counts_forward),
                          gt_transcript_counts_get_single(counts_reverse),
                          gt_transcript_bittabs_get_single(
                                                      exon_bittabs_forward),
                          gt_transcript_bittabs_get_single(
                                                      exon_bittabs_reverse),
                          gt_transcript_evaluators_get_single(exon_evaluators),
                          gt_transcript_evaluators_get_single(
                            exon_evaluators_collapsed));
      break;
    case TRANSCRIPT_FEATURE_TYPE_INITIAL:
      determine_true_exon(fn, predicted_strand, exondiff, exondiffcollapsed,
                          predicted_range,
                          gt_transcript_exons_get_initial(exons_forward),
                          gt_transcript_exons_get_initial(exons_reverse),
                          gt_transcript_counts_get_initial(counts_forward),
                          gt_transcript_counts_get_initial(counts_reverse),
                          gt_transcript_bittabs_get_initial(
                                                      exon_bittabs_forward),
                          gt_transcript_bittabs_get_initial(
                                                      exon_bittabs_reverse),
                          gt_transcript_evaluators_get_initial(exon_evaluators),
                          gt_transcript_evaluators_get_initial(
                            exon_evaluators_collapsed));
      break;
    case TRANSCRIPT_FEATURE_TYPE_INTERNAL:
      determine_true_exon(fn, predicted_strand, exondiff, exondiffcollapsed,
                          predicted_range,
                          gt_transcript_exons_get_internal(exons_forward),
                          gt_transcript_exons_get_internal(exons_reverse),
                          gt_transcript_counts_get_internal(counts_forward),
                          gt_transcript_counts_get_internal(counts_reverse),
                          gt_transcript_bittabs_get_internal(
                                                      exon_bittabs_forward),
                          gt_transcript_bittabs_get_internal(
                                                      exon_bittabs_reverse),
                          gt_transcript_evaluators_get_internal(
                                                      exon_evaluators),
                          gt_transcript_evaluators_get_internal(
                            exon_evaluators_collapsed));
      break;
    case TRANSCRIPT_FEATURE_TYPE_TERMINAL:
      determine_true_exon(fn, predicted_strand, exondiff, exondiffcollapsed,
                          predicted_range,
                          gt_transcript_exons_get_terminal(exons_forward),
                          gt_transcript_exons_get_terminal(exons_reverse),
                          gt_transcript_counts_get_terminal(counts_forward),
                          gt_transcript_counts_get_terminal(counts_reverse),
                          gt_transcript_bittabs_get_terminal(
                                                        exon_bittabs_forward),
                          gt_transcript_bittabs_get_terminal(
                                                        exon_bittabs_reverse),
                          gt_transcript_evaluators_get_terminal(
                                                        exon_evaluators),
                          gt_transcript_evaluators_get_terminal(
                            exon_evaluators_collapsed));
      break;
  }
}

typedef bool (*FeaturesAreEqualFunc)(GtGenomeNode *gn_1, GtGenomeNode *gn_2,
                                     const char *feature_type);

static void compare_features(GtArray *real_genome_nodes, GtFeatureNode *fn,
                             GtArray *genes_forward, GtArray *genes_reverse,
                             GtBittab *true_genes_forward,
                             GtBittab *true_genes_reverse,
                             GtEvaluator *gene_evaluator,
                             FeaturesAreEqualFunc features_are_equal,
                             const char *feature_type)
{
  GtStrand predicted_strand;
  GtGenomeNode **real_gn;
  unsigned long i, num;
  gt_assert(real_genome_nodes && fn && genes_forward && genes_reverse);
  gt_assert(gene_evaluator && features_are_equal && feature_type);
  predicted_strand = gt_feature_node_get_strand(fn);
  for (i = 0; i < gt_array_size(real_genome_nodes); i++) {
    real_gn = *(GtGenomeNode***) gt_array_get(real_genome_nodes, i);
    if (features_are_equal((GtGenomeNode*) fn, *real_gn, feature_type)) {
      if (predicted_strand == GT_STRAND_FORWARD) {
        num = real_gn - (GtGenomeNode**) gt_array_get_space(genes_forward);
        if (!gt_bittab_bit_is_set(true_genes_forward, num)) {
          gt_bittab_set_bit(true_genes_forward, num);
          gt_evaluator_add_true(gene_evaluator);
          /*@loopbreak@*/
          break;
        }
      }
      else {
        num = real_gn - (GtGenomeNode**) gt_array_get_space(genes_reverse);
        if (!gt_bittab_bit_is_set(true_genes_reverse, num)) {
          gt_bittab_set_bit(true_genes_reverse, num);
          gt_evaluator_add_true(gene_evaluator);
          /*@loopbreak@*/
          break;
        }
      }
    }
  }
}

static int process_predicted_feature(GtFeatureNode *fn, void *data,
                                     GT_UNUSED GtError *err)
{
  ProcessPredictedFeatureInfo *info = (ProcessPredictedFeatureInfo*) data;
  GtRange predicted_range;
  unsigned long i, num;
  GtStrand predicted_strand;
  GtArray *real_genome_nodes;
  GtGenomeNode **real_gn;

  gt_error_check(err);
  gt_assert(fn && data);

  predicted_range = gt_genome_node_get_range((GtGenomeNode*) fn);
  predicted_strand = gt_feature_node_get_strand(fn);
  real_genome_nodes = gt_array_new(sizeof (GtGenomeNode**));

  if (gt_feature_node_has_type(fn, gt_ft_gene)) {
    /* store predicted gene */
    gt_evaluator_add_predicted(info->mRNA_gene_evaluator, 1);
    gt_evaluator_add_predicted(info->CDS_gene_evaluator, 1);
    /* determine true gene */
    switch (predicted_strand) {
      case GT_STRAND_FORWARD:
      case GT_STRAND_REVERSE:
        gt_bsearch_all_mark(real_genome_nodes, &fn,
                            predicted_strand == GT_STRAND_FORWARD
                            ? gt_array_get_space(info->slot->genes_forward)
                            : gt_array_get_space(info->slot->genes_reverse),
                            predicted_strand == GT_STRAND_FORWARD
                            ? gt_array_size(info->slot->genes_forward)
                            : gt_array_size(info->slot->genes_reverse),
                            sizeof (GtGenomeNode*),
                            (GtCompareWithData)
                            gt_genome_node_compare_with_data,
                            NULL,
                            predicted_strand == GT_STRAND_FORWARD
                            ? info->slot->overlapped_genes_forward
                            : info->slot->overlapped_genes_reverse);
        if (gt_array_size(real_genome_nodes)) {
          /* gene(s) with the same range found -> check if they are equal */
          compare_features(real_genome_nodes, fn, info->slot->genes_forward,
                           info->slot->genes_reverse,
                           info->slot->true_mRNA_genes_forward,
                           info->slot->true_mRNA_genes_reverse,
                           info->mRNA_gene_evaluator, genes_are_equal,
                           gt_ft_exon);
          compare_features(real_genome_nodes, fn, info->slot->genes_forward,
                           info->slot->genes_reverse,
                           info->slot->true_CDS_genes_forward,
                           info->slot->true_CDS_genes_reverse,
                           info->CDS_gene_evaluator, genes_are_equal,
                           gt_ft_CDS);
        }
        else {
          /* no gene with the same range found -> check if this is a wrong
             gene */
          if (!gt_feature_node_overlaps_nodes_mark(fn,
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
                  "(line %u)\n",
                  gt_genome_node_get_line_number((GtGenomeNode*) fn));
        }
    }
  }
  else if (gt_feature_node_has_type(fn, gt_ft_mRNA)) {
    /* store predicted mRNA */
    gt_evaluator_add_predicted(info->mRNA_mRNA_evaluator, 1);
    gt_evaluator_add_predicted(info->CDS_mRNA_evaluator, 1);
    /* determine true mRNA */
    switch (predicted_strand) {
      case GT_STRAND_FORWARD:
      case GT_STRAND_REVERSE:
        gt_bsearch_all_mark(real_genome_nodes, &fn,
                            predicted_strand == GT_STRAND_FORWARD
                            ? gt_array_get_space(info->slot->mRNAs_forward)
                            : gt_array_get_space(info->slot->mRNAs_reverse),
                            predicted_strand == GT_STRAND_FORWARD
                            ? gt_array_size(info->slot->mRNAs_forward)
                            : gt_array_size(info->slot->mRNAs_reverse),
                            sizeof (GtGenomeNode*),
                            (GtCompareWithData)
                            gt_genome_node_compare_with_data,
                            NULL,
                            predicted_strand == GT_STRAND_FORWARD
                            ? info->slot->overlapped_mRNAs_forward
                            : info->slot->overlapped_mRNAs_reverse);
        if (gt_array_size(real_genome_nodes)) {
          /* mRNA(s) with the same range found -> check if they are equal */
          compare_features(real_genome_nodes, fn, info->slot->mRNAs_forward,
                           info->slot->mRNAs_reverse,
                           info->slot->true_mRNA_mRNAs_forward,
                           info->slot->true_mRNA_mRNAs_reverse,
                           info->mRNA_mRNA_evaluator, mRNAs_are_equal,
                           gt_ft_exon);
          compare_features(real_genome_nodes, fn, info->slot->mRNAs_forward,
                           info->slot->mRNAs_reverse,
                           info->slot->true_CDS_mRNAs_forward,
                           info->slot->true_CDS_mRNAs_reverse,
                           info->CDS_mRNA_evaluator, mRNAs_are_equal,
                           gt_ft_CDS);
        }
        else {
          /* no mRNA with the same range found -> check if this is a wrong
             mRNA */
          if (!gt_feature_node_overlaps_nodes_mark(fn,
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
                  "(line %u)\n",
                  gt_genome_node_get_line_number((GtGenomeNode*) fn));
        }
    }
  }
  else if (gt_feature_node_has_type(fn, gt_ft_LTR_retrotransposon)) {
    /* store predicted LTR */
    gt_evaluator_add_predicted(info->LTR_evaluator, 1);
    /* determine true LTR */
    gt_bsearch_all_mark(real_genome_nodes, &fn,
                        gt_array_get_space(info->slot->LTRs),
                        gt_array_size(info->slot->LTRs),
                        sizeof (GtGenomeNode*),
                        (GtCompareWithData) gt_genome_node_compare_delta,
                        &info->LTRdelta,
                        info->slot->overlapped_LTRs);

    if (gt_array_size(real_genome_nodes)) {
      for (i = 0; i < gt_array_size(real_genome_nodes); i++) {
        real_gn = *(GtGenomeNode***) gt_array_get(real_genome_nodes, i);
        num = real_gn - (GtGenomeNode**) gt_array_get_space(info->slot->LTRs);
        if (!gt_bittab_bit_is_set(info->slot->true_LTRs, num)) {
          gt_bittab_set_bit(info->slot->true_LTRs, num);
          gt_evaluator_add_true(info->LTR_evaluator);
          /*@loopbreak@*/
          break;
        }
      }
    }
    else {
      /* no LTR with the same range found -> check if this is a wrong LTR */
      if (!gt_feature_node_overlaps_nodes_mark(fn, info->slot->LTRs,
                                               info->slot->overlapped_LTRs)) {
        (*info->wrong_LTRs)++;
      }
    }
  }
  else if (gt_feature_node_has_type(fn, gt_ft_exon)) {
    /* store predicted exon (mRNA level)*/
    store_predicted_exon(info->mRNA_exon_evaluators, fn);

    /* store predicted exon (mRNA level, collapsed) */
    store_predicted_exon_collapsed(predicted_strand == GT_STRAND_FORWARD
                                   ? info->slot->used_mRNA_exons_forward
                                   : info->slot->used_mRNA_exons_reverse,
                                   &predicted_range,
                                   info->mRNA_exon_evaluators_collapsed, fn);

    /* determine true exon (mRNA level)*/
    switch (predicted_strand) {
      case GT_STRAND_FORWARD:
      case GT_STRAND_REVERSE:
        store_true_exon(fn, predicted_strand, &predicted_range,
                        info->exondiff, info->exondiffcollapsed,
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
                  "(line %u)\n",
                  gt_genome_node_get_line_number((GtGenomeNode*) fn));
        }
    }
  }
  else if (gt_feature_node_has_type(fn, gt_ft_CDS)) {
    /* store predicted exon (CDS level)*/
    store_predicted_exon(info->CDS_exon_evaluators, fn);

    /* store predicted exon (CDS level, collapsed) */
    store_predicted_exon_collapsed(predicted_strand == GT_STRAND_FORWARD
                                   ? info->slot->used_CDS_exons_forward
                                   : info->slot->used_CDS_exons_reverse,
                                   &predicted_range,
                                   info->CDS_exon_evaluators_collapsed, fn);

    /* determine true exon (CDS level) */
    switch (predicted_strand) {
      case GT_STRAND_FORWARD:
      case GT_STRAND_REVERSE:
        store_true_exon(fn, predicted_strand, &predicted_range,
                        info->exondiff, info->exondiffcollapsed,
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
                  "(line %u)\n",
                  gt_genome_node_get_line_number((GtGenomeNode*) fn));
        }
      }
  }
  gt_array_delete(real_genome_nodes);
  return 0;
}

static int determine_missing_features(GT_UNUSED void *key, void *value,
                                      void *data, GT_UNUSED GtError *err)
{
  GtStreamEvaluator *se = (GtStreamEvaluator*) data;
  Slot *slot = (Slot*) value;
  gt_error_check(err);
  gt_assert(key && value && data);
  if (slot->overlapped_genes_forward) {
    se->missing_genes +=
      gt_bittab_size(slot->overlapped_genes_forward) -
      gt_bittab_count_set_bits(slot->overlapped_genes_forward);
  }
  if (slot->overlapped_genes_reverse) {
    se->missing_genes +=
      gt_bittab_size(slot->overlapped_genes_reverse) -
      gt_bittab_count_set_bits(slot->overlapped_genes_reverse);
  }
  if (slot->overlapped_mRNAs_forward) {
    se->missing_mRNAs +=
      gt_bittab_size(slot->overlapped_mRNAs_forward) -
      gt_bittab_count_set_bits(slot->overlapped_mRNAs_forward);
  }
  if (slot->overlapped_mRNAs_reverse) {
    se->missing_mRNAs +=
      gt_bittab_size(slot->overlapped_mRNAs_reverse) -
      gt_bittab_count_set_bits(slot->overlapped_mRNAs_reverse);
  }
  if (slot->overlapped_LTRs) {
    se->missing_LTRs  += gt_bittab_size(slot->overlapped_LTRs) -
                         gt_bittab_count_set_bits(slot->overlapped_LTRs);
  }
  return 0;
}

static void add_nucleotide_values(NucEval *nucleotides, GtBittab *real,
                                  GtBittab *pred, GtBittab *tmp,
                                  const char *level)
{
  gt_assert(nucleotides && real && pred && tmp);
  if (gt_log_enabled()) {
    gt_log_log("%s", level);
    gt_log_log("reference:");
    gt_bittab_show(real, gt_log_fp());
    gt_log_log("prediction:");
    gt_bittab_show(pred, gt_log_fp());
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

static int compute_nucleotides_values(GT_UNUSED void *key, void *value,
                                      void *data, GT_UNUSED GtError *err)
{
  GtStreamEvaluator *se = (GtStreamEvaluator*) data;
  Slot *slot = (Slot*) value;
  GtBittab *tmp;
  gt_error_check(err);
  gt_assert(key && value && data);
  /* add ``out of range'' FPs */
  se->mRNA_nucleotides.FP += slot->FP_mRNA_nucleotides_forward;
  se->mRNA_nucleotides.FP += slot->FP_mRNA_nucleotides_reverse;
  se->CDS_nucleotides.FP  += slot->FP_CDS_nucleotides_forward;
  se->CDS_nucleotides.FP  += slot->FP_CDS_nucleotides_reverse;
  /* add other values */
  tmp = gt_bittab_new(gt_range_length(&slot->real_range));
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

int gt_stream_evaluator_evaluate(GtStreamEvaluator *se, bool verbose,
                                 bool exondiff, bool exondiffcollapsed,
                                 GtNodeVisitor *nv, GtError *err)
{
  GtGenomeNode *gn;
  GtFeatureNode *fn;
  Slot *slot;
  ProcessRealFeatureInfo real_info;
  ProcessPredictedFeatureInfo predicted_info;
  int had_err;

  gt_error_check(err);
  gt_assert(se);

  /* init */
  real_info.nuceval = se->nuceval;
  real_info.verbose = verbose;
  predicted_info.nuceval = se->nuceval;
  predicted_info.verbose = verbose;
  predicted_info.exondiff = exondiff;
  predicted_info.exondiffcollapsed = exondiffcollapsed;
  predicted_info.LTRdelta = se->LTRdelta;
  predicted_info.mRNA_gene_evaluator = se->mRNA_gene_evaluator;
  predicted_info.CDS_gene_evaluator = se->CDS_gene_evaluator;
  predicted_info.mRNA_mRNA_evaluator = se->mRNA_mRNA_evaluator;
  predicted_info.CDS_mRNA_evaluator = se->CDS_mRNA_evaluator;
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

  /* process the reference stream completely */
  while (!(had_err = gt_node_stream_next(se->reference, &gn, err)) && gn) {
    if (gt_region_node_try_cast(gn)) {
      /* each sequence region gets its own ``slot'' */
      if (!(slot = gt_hashmap_get(se->slots,
                                  gt_str_get(gt_genome_node_get_seqid(gn))))) {
        slot = slot_new(se->nuceval, gt_genome_node_get_range(gn));
        gt_hashmap_add(se->slots,
                       gt_cstr_dup(gt_str_get(gt_genome_node_get_seqid(gn))),
                       slot);
      }
      gt_assert(slot);
    }
    /* we consider only genome features */
    if ((fn = gt_feature_node_try_cast(gn))) {
      /* each sequence must have its own ``slot'' at this point */
      slot = gt_hashmap_get(se->slots,
                            gt_str_get(gt_genome_node_get_seqid(gn)));
      gt_assert(slot);
      /* store the exons */
      real_info.slot = slot;
      gt_feature_node_determine_transcripttypes(fn);
      had_err = gt_feature_node_traverse_children(fn, &real_info,
                                                  process_real_feature, false,
                                                  NULL);
      gt_assert(!had_err); /* cannot happen, process_real_feature() is sane */
    }
    if (nv)
      gt_genome_node_accept(gn, nv, err);
    gt_genome_node_delete(gn);
  }

  /* set the actuals and sort them */
  if (!had_err) {
    had_err = gt_hashmap_foreach(se->slots, set_actuals_and_sort_them, se,
                                 NULL);
    gt_assert(!had_err); /* set_actuals_and_sort_them() is sane */
  }

  /* process the prediction stream */
  if (!had_err) {
    while (!(had_err = gt_node_stream_next(se->prediction, &gn, err)) &&
           gn) {
      /* we consider only genome features */
      if ((fn = gt_feature_node_try_cast(gn))) {
        /* get (real) slot */
        slot = gt_hashmap_get(se->slots,
                              gt_str_get(gt_genome_node_get_seqid(gn)));
        if (slot) {
          predicted_info.slot = slot;
          gt_feature_node_determine_transcripttypes(fn);
          had_err = gt_feature_node_traverse_children(fn, &predicted_info,
                                                      process_predicted_feature,
                                                      false, NULL);
          gt_assert(!had_err); /* cannot happen, process_predicted_feature() is
                               sane */
        }
        else {
          /* we got no (real) slot */
          gt_warning("sequence id \"%s\" (with predictions) not given in "
                     "reference", gt_str_get(gt_genome_node_get_seqid(gn)));
        }
      }
      if (nv)
        had_err = gt_genome_node_accept(gn, nv, err);
      gt_genome_node_delete(gn);
    }
  }

  /* determine the missing mRNAs */
  if (!had_err) {
    had_err = gt_hashmap_foreach(se->slots, determine_missing_features, se,
                                 NULL);
    gt_assert(!had_err); /* determine_missing_features() is sane */
  }

  /* compute the nucleotides values */
  if (!had_err && se->nuceval) {
    had_err = gt_hashmap_foreach(se->slots, compute_nucleotides_values, se,
                                 NULL);
    gt_assert(!had_err); /* compute_nucleotides_values() is sane */
  }

  return had_err;
}

static void show_transcript_values(GtTranscriptEvaluators *te,
                                   const char *level,
                                   const char *additional_info, GtFile *outfp)
{
  gt_assert(te);

  gt_file_xprintf(outfp, "exon sensitivity (%s level, all%s): ", level,
                  additional_info);
  gt_evaluator_show_sensitivity(gt_transcript_evaluators_get_all(te), outfp);
  gt_file_xfputc('\n', outfp);

  gt_file_xprintf(outfp, "exon specificity (%s level, all%s): ", level,
                  additional_info);
  gt_evaluator_show_specificity(gt_transcript_evaluators_get_all(te), outfp);
  gt_file_xfputc('\n', outfp);

  gt_file_xprintf(outfp, "exon sensitivity (%s level, single%s): ", level,
                  additional_info);
  gt_evaluator_show_sensitivity(gt_transcript_evaluators_get_single(te), outfp);
  gt_file_xfputc('\n', outfp);

  gt_file_xprintf(outfp, "exon specificity (%s level, single%s): ", level,
                  additional_info );
  gt_evaluator_show_specificity(gt_transcript_evaluators_get_single(te), outfp);
  gt_file_xfputc('\n', outfp);

  gt_file_xprintf(outfp, "exon sensitivity (%s level, initial%s): ", level,
                  additional_info);
  gt_evaluator_show_sensitivity(gt_transcript_evaluators_get_initial(te),
                                outfp);
  gt_file_xfputc('\n', outfp);

  gt_file_xprintf(outfp, "exon specificity (%s level, initial%s): ", level,
                  additional_info);
  gt_evaluator_show_specificity(gt_transcript_evaluators_get_initial(te),
                                outfp);
  gt_file_xfputc('\n', outfp);

  gt_file_xprintf(outfp, "exon sensitivity (%s level, internal%s): ", level,
                  additional_info);
  gt_evaluator_show_sensitivity(gt_transcript_evaluators_get_internal(te),
                                outfp);
  gt_file_xfputc('\n', outfp);

  gt_file_xprintf(outfp, "exon specificity (%s level, internal%s): ", level,
                  additional_info);
  gt_evaluator_show_specificity(gt_transcript_evaluators_get_internal(te),
                                outfp);
  gt_file_xfputc('\n', outfp);

  gt_file_xprintf(outfp, "exon sensitivity (%s level, terminal%s): ", level,
                  additional_info);
  gt_evaluator_show_sensitivity(gt_transcript_evaluators_get_terminal(te),
                                outfp);
  gt_file_xfputc('\n', outfp);

  gt_file_xprintf(outfp, "exon specificity (%s level, terminal%s): ", level,
                  additional_info);
  gt_evaluator_show_specificity(gt_transcript_evaluators_get_terminal(te),
                                outfp);
  gt_file_xfputc('\n', outfp);
}

static void show_nucleotide_values(NucEval *nucleotides, const char *level,
                                   GtFile *outfp)
{
  double sensitivity = 1.0, specificity = 1.0;
  gt_assert(nucleotides && level);
  if (nucleotides->TP || nucleotides->FN) {
    sensitivity = (double) nucleotides->TP /
                  (nucleotides->TP + nucleotides->FN);
  }
  if (nucleotides->TP || nucleotides->FP) {
    specificity = (double) nucleotides->TP /
                  (nucleotides->TP + nucleotides->FP);
  }
  gt_file_xprintf(outfp, "nucleotide sensitivity (%s level): %6.2f%% "
                  "(TP=%lu/(TP=%lu + FN=%lu))\n", level, sensitivity * 100.0,
                  nucleotides->TP, nucleotides->TP, nucleotides->FN);
  gt_file_xprintf(outfp, "nucleotide specificity (%s level): %6.2f%% "
                  "(TP=%lu/(TP=%lu + FP=%lu))\n", level, specificity * 100.0,
                  nucleotides->TP, nucleotides->TP, nucleotides->FP);
}

void gt_stream_evaluator_show(GtStreamEvaluator *se, GtFile *outfp)
{
  gt_assert(se);

  if (!se->evalLTR) {
    /* gene level */
    gt_file_xprintf(outfp, "gene sensitivity (mRNA level): ");
    gt_evaluator_show_sensitivity(se->mRNA_gene_evaluator, outfp);
    gt_file_xprintf(outfp, " (missing genes: %lu)\n", se->missing_genes);

    gt_file_xprintf(outfp, "gene specificity (mRNA level): ");
    gt_evaluator_show_specificity(se->mRNA_gene_evaluator, outfp);
    gt_file_xprintf(outfp, " (wrong genes: %lu)\n", se->wrong_genes);

    gt_file_xprintf(outfp, "gene sensitivity (CDS level): ");
    gt_evaluator_show_sensitivity(se->CDS_gene_evaluator, outfp);
    gt_file_xprintf(outfp, " (missing genes: %lu)\n", se->missing_genes);

    gt_file_xprintf(outfp, "gene specificity (CDS level): ");
    gt_evaluator_show_specificity(se->CDS_gene_evaluator, outfp);
    gt_file_xprintf(outfp, " (wrong genes: %lu)\n", se->wrong_genes);

    /* mRNA level */
    gt_file_xprintf(outfp, "mRNA sensitivity (mRNA level): ");
    gt_evaluator_show_sensitivity(se->mRNA_mRNA_evaluator, outfp);
    gt_file_xprintf(outfp, " (missing mRNAs: %lu)\n", se->missing_mRNAs);

    gt_file_xprintf(outfp, "mRNA specificity (mRNA level): ");
    gt_evaluator_show_specificity(se->mRNA_mRNA_evaluator, outfp);
    gt_file_xprintf(outfp, " (wrong mRNAs: %lu)\n", se->wrong_mRNAs);

    gt_file_xprintf(outfp, "mRNA sensitivity (CDS level): ");
    gt_evaluator_show_sensitivity(se->CDS_mRNA_evaluator, outfp);
    gt_file_xprintf(outfp, " (missing mRNAs: %lu)\n", se->missing_mRNAs);

    gt_file_xprintf(outfp, "mRNA specificity (CDS level): ");
    gt_evaluator_show_specificity(se->CDS_mRNA_evaluator, outfp);
    gt_file_xprintf(outfp, " (wrong mRNAs: %lu)\n", se->wrong_mRNAs);

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
    gt_file_xprintf(outfp, "LTR_retrotransposon sensitivity: ");
    gt_evaluator_show_sensitivity(se->LTR_evaluator, outfp);
    gt_file_xprintf(outfp, " (missing LTRs: %lu)\n", se->missing_LTRs);

    gt_file_xprintf(outfp, "LTR_retrotransposon specificity: ");
    gt_evaluator_show_specificity(se->LTR_evaluator, outfp);
    gt_file_xprintf(outfp, " (wrong LTRs: %lu)\n", se->wrong_LTRs);
  }
}

void gt_stream_evaluator_delete(GtStreamEvaluator *se)
{
  if (!se) return;
  gt_node_stream_delete(se->reference);
  gt_node_stream_delete(se->prediction);
  gt_hashmap_delete(se->slots);
  gt_evaluator_delete(se->mRNA_gene_evaluator);
  gt_evaluator_delete(se->CDS_gene_evaluator);
  gt_evaluator_delete(se->mRNA_mRNA_evaluator);
  gt_evaluator_delete(se->CDS_mRNA_evaluator);
  gt_evaluator_delete(se->LTR_evaluator);
  gt_transcript_evaluators_delete(se->mRNA_exon_evaluators);
  gt_transcript_evaluators_delete(se->mRNA_exon_evaluators_collapsed);
  gt_transcript_evaluators_delete(se->CDS_exon_evaluators);
  gt_transcript_evaluators_delete(se->CDS_exon_evaluators_collapsed);
  gt_free(se);
}
