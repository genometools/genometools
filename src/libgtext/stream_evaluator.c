/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include <libgtext/bsearch.h>
#include <libgtext/evaluator.h>
#include <libgtext/gff3_output.h>
#include <libgtext/stream_evaluator.h>
#include <libgtext/transcript_evaluators.h>
#include <libgtext/transcript_exons.h>
#include <libgtext/transcript_used_exons.h>

struct StreamEvaluator {
  GenomeStream *reality,
               *prediction;
  Hashtable *real_features; /* sequence name -> feature type hash */
  Evaluator *gene_evaluator,
            *mRNA_evaluator;
  TranscriptEvaluators *mRNA_exon_evaluators,
                       *mRNA_exon_evaluators_collapsed,
                       *CDS_exon_evaluators,
                       *CDS_exon_evaluators_collapsed;
  unsigned long missing_genes,
                wrong_genes,
                missing_mRNAs,
                wrong_mRNAs;
};

typedef struct {
  Array *genes_forward,
        *genes_reverse,
        *mRNAs_forward,
        *mRNAs_reverse;
  TranscriptExons *mRNA_exons_forward,
                  *mRNA_exons_reverse,
                  *CDS_exons_forward,
                  *CDS_exons_reverse;
  TranscriptCounts *mRNA_counts_forward,
                   *mRNA_counts_reverse,
                   *CDS_counts_forward,
                   *CDS_counts_reverse;
  Bittab *true_genes_forward,
         *true_genes_reverse,
         *true_mRNAs_forward,
         *true_mRNAs_reverse,
         *overlapped_genes_forward,
         *overlapped_genes_reverse,
         *overlapped_mRNAs_forward,
         *overlapped_mRNAs_reverse;
  TranscriptBittabs *mRNA_exon_bittabs_forward,
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
  bool verbose;
} Process_real_feature_data;

typedef struct {
  Slot *slot;
  bool verbose,
       exondiff;
  Evaluator *gene_evaluator,
            *mRNA_evaluator;
  TranscriptEvaluators *mRNA_exon_evaluators,
                       *mRNA_exon_evaluators_collapsed,
                       *CDS_exon_evaluators,
                       *CDS_exon_evaluators_collapsed;
  unsigned long *wrong_genes,
                *wrong_mRNAs;
} Process_predicted_feature_info;

static Slot* slot_new(Env *env)
{
  Slot *s = env_ma_calloc(env, 1, sizeof (Slot));
  s->genes_forward = array_new(sizeof (GenomeNode*), env);
  s->genes_reverse = array_new(sizeof (GenomeNode*), env);
  s->mRNAs_forward = array_new(sizeof (GenomeNode*), env);
  s->mRNAs_reverse = array_new(sizeof (GenomeNode*), env);
  s->mRNA_exons_forward = transcript_exons_new(env);
  s->mRNA_exons_reverse = transcript_exons_new(env);
  s->CDS_exons_forward = transcript_exons_new(env);
  s->CDS_exons_reverse = transcript_exons_new(env);
  s->used_mRNA_exons_forward = transcript_used_exons_new(env);
  s->used_mRNA_exons_reverse = transcript_used_exons_new(env);
  s->used_CDS_exons_forward = transcript_used_exons_new(env);
  s->used_CDS_exons_reverse = transcript_used_exons_new(env);
  return s;
}

static void slot_delete(Slot *s, Env *env)
{
  unsigned long i;
  assert(s);
  for (i = 0; i < array_size(s->genes_forward); i++)
    genome_node_rec_delete(*(GenomeNode**) array_get(s->genes_forward, i), env);
  array_delete(s->genes_forward, env);
  for (i = 0; i < array_size(s->genes_reverse); i++)
    genome_node_rec_delete(*(GenomeNode**) array_get(s->genes_reverse, i), env);
  array_delete(s->genes_reverse, env);
  for (i = 0; i < array_size(s->mRNAs_forward); i++)
    genome_node_rec_delete(*(GenomeNode**) array_get(s->mRNAs_forward, i), env);
  array_delete(s->mRNAs_forward, env);
  for (i = 0; i < array_size(s->mRNAs_reverse); i++)
    genome_node_rec_delete(*(GenomeNode**) array_get(s->mRNAs_reverse, i), env);
  array_delete(s->mRNAs_reverse, env);
  transcript_exons_delete(s->mRNA_exons_forward, env);
  transcript_exons_delete(s->mRNA_exons_reverse, env);
  transcript_exons_delete(s->CDS_exons_forward, env);
  transcript_exons_delete(s->CDS_exons_reverse, env);
  transcript_counts_delete(s->mRNA_counts_forward, env);
  transcript_counts_delete(s->mRNA_counts_reverse, env);
  transcript_counts_delete(s->CDS_counts_forward, env);
  transcript_counts_delete(s->CDS_counts_reverse, env);
  bittab_delete(s->true_genes_forward, env);
  bittab_delete(s->true_genes_reverse, env);
  bittab_delete(s->true_mRNAs_forward, env);
  bittab_delete(s->true_mRNAs_reverse, env);
  bittab_delete(s->overlapped_genes_forward, env);
  bittab_delete(s->overlapped_genes_reverse, env);
  bittab_delete(s->overlapped_mRNAs_forward, env);
  bittab_delete(s->overlapped_mRNAs_reverse, env);
  transcript_bittabs_delete(s->mRNA_exon_bittabs_forward, env);
  transcript_bittabs_delete(s->mRNA_exon_bittabs_reverse, env);
  transcript_bittabs_delete(s->CDS_exon_bittabs_forward, env);
  transcript_bittabs_delete(s->CDS_exon_bittabs_reverse, env);
  transcript_used_exons_delete(s->used_mRNA_exons_forward, env);
  transcript_used_exons_delete(s->used_mRNA_exons_reverse, env);
  transcript_used_exons_delete(s->used_CDS_exons_forward, env);
  transcript_used_exons_delete(s->used_CDS_exons_reverse, env);
  env_ma_free(s, env);
}

StreamEvaluator* stream_evaluator_new(GenomeStream *reality,
                                      GenomeStream *prediction, Env *env)
{
  StreamEvaluator *evaluator = env_ma_malloc(env, sizeof (StreamEvaluator));
  evaluator->reality = reality;
  evaluator->prediction = prediction;
  evaluator->real_features = hashtable_new(HASH_STRING, NULL,
                                           (FreeFunc) slot_delete, env);
  evaluator->gene_evaluator = evaluator_new(env);
  evaluator->mRNA_evaluator = evaluator_new(env);
  evaluator->mRNA_exon_evaluators = transcript_evaluators_new(env);
  evaluator->mRNA_exon_evaluators_collapsed = transcript_evaluators_new(env);
  evaluator->CDS_exon_evaluators = transcript_evaluators_new(env);
  evaluator->CDS_exon_evaluators_collapsed = transcript_evaluators_new(env);
  evaluator->missing_genes = 0;
  evaluator->wrong_genes = 0;
  evaluator->missing_mRNAs = 0;
  evaluator->wrong_mRNAs = 0;
  return evaluator;
}

static int set_actuals_and_sort_them(void *key, void *value, void *data,
                                     Env *env)
{
  StreamEvaluator *se = (StreamEvaluator*) data;
  Slot *s = (Slot*) value;

  env_error_check(env);
  assert(key && value && data);

  /* set actual genes */
  evaluator_add_actual(se->gene_evaluator, array_size(s->genes_forward));
  evaluator_add_actual(se->gene_evaluator, array_size(s->genes_reverse));

  /* set actual mRNAs */
  evaluator_add_actual(se->mRNA_evaluator, array_size(s->mRNAs_forward));
  evaluator_add_actual(se->mRNA_evaluator, array_size(s->mRNAs_reverse));

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

  /* sort exons */
  transcript_exons_sort(s->mRNA_exons_forward);
  transcript_exons_sort(s->mRNA_exons_reverse);
  transcript_exons_sort(s->CDS_exons_forward);
  transcript_exons_sort(s->CDS_exons_reverse);

  /* determine true exons */
  s->mRNA_counts_forward =
    transcript_exons_uniq_in_place_count(s->mRNA_exons_forward, env);
  s->mRNA_counts_reverse =
    transcript_exons_uniq_in_place_count(s->mRNA_exons_reverse, env);
  s->CDS_counts_forward =
    transcript_exons_uniq_in_place_count(s->CDS_exons_forward, env);
  s->CDS_counts_reverse =
    transcript_exons_uniq_in_place_count(s->CDS_exons_reverse, env);

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

  /* make sure that the exons are sorted */
  assert(transcript_exons_are_sorted(s->mRNA_exons_forward));
  assert(transcript_exons_are_sorted(s->mRNA_exons_reverse));
  assert(transcript_exons_are_sorted(s->CDS_exons_forward));
  assert(transcript_exons_are_sorted(s->CDS_exons_reverse));

  /* init true bittabs */
  s->true_genes_forward = array_size(s->genes_forward)
                          ? bittab_new(array_size(s->genes_forward), env)
                          : NULL;
  s->true_genes_reverse = array_size(s->genes_reverse)
                          ? bittab_new(array_size(s->genes_reverse), env)
                          : NULL;
  s->true_mRNAs_forward = array_size(s->mRNAs_forward)
                          ? bittab_new(array_size(s->mRNAs_forward), env)
                          : NULL;
  s->true_mRNAs_reverse = array_size(s->mRNAs_reverse)
                          ? bittab_new(array_size(s->mRNAs_reverse), env)
                          : NULL;

  /* init overlap bittabs */
  s->overlapped_genes_forward = array_size(s->genes_forward)
                                ? bittab_new(array_size(s->genes_forward), env)
                                : NULL;
  s->overlapped_genes_reverse = array_size(s->genes_reverse)
                                ? bittab_new(array_size(s->genes_reverse), env)
                                : NULL;
  s->overlapped_mRNAs_forward = array_size(s->mRNAs_forward)
                                ? bittab_new(array_size(s->mRNAs_forward), env)
                                : NULL;
  s->overlapped_mRNAs_reverse = array_size(s->mRNAs_reverse)
                                ? bittab_new(array_size(s->mRNAs_reverse), env)
                                : NULL;

  /* init bittabs (for collapsed exons) */
  s->mRNA_exon_bittabs_forward =
    transcript_exons_create_bittabs(s->mRNA_exons_forward, env);
  s->mRNA_exon_bittabs_reverse =
    transcript_exons_create_bittabs(s->mRNA_exons_reverse, env);
  s->CDS_exon_bittabs_forward =
    transcript_exons_create_bittabs(s->CDS_exons_forward, env);
  s->CDS_exon_bittabs_reverse =
    transcript_exons_create_bittabs(s->CDS_exons_reverse, env);

  return 0;
}

static void add_exon(TranscriptExons *te, Range range, GenomeFeature *gf,
                     Env *env)
{
  assert(te);
  array_add(transcript_exons_get_all(te), range, env);
  switch (genome_feature_get_transcriptfeaturetype(gf)) {
    case TRANSCRIPT_FEATURE_TYPE_SINGLE:
      array_add(transcript_exons_get_single(te), range, env);
      break;
    case TRANSCRIPT_FEATURE_TYPE_INITIAL:
      array_add(transcript_exons_get_initial(te), range, env);
      break;
    case TRANSCRIPT_FEATURE_TYPE_INTERNAL:
      array_add(transcript_exons_get_internal(te), range, env);
      break;
    case TRANSCRIPT_FEATURE_TYPE_TERMINAL:
      array_add(transcript_exons_get_terminal(te), range, env);
      break;
    case TRANSCRIPT_FEATURE_TYPE_UNDETERMINED: assert(0);
  }
}

static int process_real_feature(GenomeNode *gn, void *data, Env *env)
{
  Process_real_feature_data *process_real_feature_data =
    (Process_real_feature_data*) data;
  GenomeNode *gn_ref;
  GenomeFeature *gf;
  Range range;

  env_error_check(env);
  assert(gn && data);
  gf = (GenomeFeature*) gn;

  switch (genome_feature_get_type(gf)) {
    case gft_gene:
      switch (genome_feature_get_strand(gf)) {
        case STRAND_FORWARD:
          gn_ref = genome_node_rec_ref(gn, env);
          array_add(process_real_feature_data->slot->genes_forward, gn_ref,
                    env);
          break;
        case STRAND_REVERSE:
          gn_ref = genome_node_rec_ref(gn, env);
          array_add(process_real_feature_data->slot->genes_reverse, gn_ref,
                    env);
          break;
        default:
          if (process_real_feature_data->verbose) {
            fprintf(stderr, "skipping real gene with unknown orientation "
                    "(line %lu)\n", genome_node_get_line_number(gn));
          }
      }
      break;
    case gft_mRNA:
      switch (genome_feature_get_strand(gf)) {
        case STRAND_FORWARD:
          gn_ref = genome_node_rec_ref(gn, env);
          array_add(process_real_feature_data->slot->mRNAs_forward, gn_ref,
                    env);
          break;
        case STRAND_REVERSE:
          gn_ref = genome_node_rec_ref(gn, env);
          array_add(process_real_feature_data->slot->mRNAs_reverse, gn_ref,
                    env);
          break;
        default:
          if (process_real_feature_data->verbose) {
            fprintf(stderr, "skipping real mRNA with unknown orientation "
                    "(line %lu)\n", genome_node_get_line_number(gn));
          }
      }
      break;
    case gft_CDS:
      range = genome_node_get_range(gn);
      switch (genome_feature_get_strand(gf)) {
        case STRAND_FORWARD:
          add_exon(process_real_feature_data->slot->CDS_exons_forward, range,
                   gf, env);
          break;
        case STRAND_REVERSE:
          add_exon(process_real_feature_data->slot->CDS_exons_reverse, range,
                   gf, env);
          break;
        default:
          if (process_real_feature_data->verbose) {
            fprintf(stderr, "skipping real CDS exon with unknown orientation "
                    "(line %lu)\n", genome_node_get_line_number(gn));
          }
      }
      break;
    case gft_exon:
      range = genome_node_get_range(gn);
      switch (genome_feature_get_strand(gf)) {
        case STRAND_FORWARD:
          add_exon(process_real_feature_data->slot->mRNA_exons_forward, range,
                   gf, env);
          break;
        case STRAND_REVERSE:
          add_exon(process_real_feature_data->slot->mRNA_exons_reverse, range,
                   gf, env);
          break;
        default:
          if (process_real_feature_data->verbose) {
            fprintf(stderr, "skipping real mRNA exon with unknown orientation "
                    "(line %lu)\n", genome_node_get_line_number(gn));
          }
      }
      break;
    default:
      assert(1); /* shut up compiler */
  }
  return 0;
}

static int store_exon(GenomeNode *gn, void *data, Env *env)
{
  Array *exons = (Array*) data;
  Range range;
  GenomeFeature *gf;
  env_error_check(env);
  gf = genome_node_cast(genome_feature_class(), gn);
  assert(gf && exons);
  if (genome_feature_get_type(gf) == gft_exon) {
    range = genome_node_get_range(gn);
    array_add(exons, range, env);
  }
  return 0;
}

static bool mRNAs_are_equal(GenomeNode *gn_1, GenomeNode *gn_2, Env *env)
{
  Array *exons_1, *exons_2;
  bool equal;
  int has_err;

  assert(gn_1 && gn_2);

  /* init */
  exons_1 = array_new(sizeof (Range), env);
  exons_2 = array_new(sizeof (Range), env);

  /* get exon ranges */
  has_err = genome_node_traverse_children(gn_1, exons_1, store_exon, false,
                                          env);
  assert(!has_err); /* cannot happen, store_exon() is sane */
  has_err = genome_node_traverse_children(gn_2, exons_2, store_exon, false,
                                          env);
  assert(!has_err); /* cannot happen, store_exon() is sane */

  /* sort exon ranges */
  ranges_sort(exons_1);
  ranges_sort(exons_2);

  /* compare exon ranges */
  equal = ranges_are_equal(exons_1, exons_2);

  /* free */
  array_delete(exons_1, env);
  array_delete(exons_2, env);

  return equal;
}

typedef struct {
  Array *exons,
        *mRNAs;
} Store_gene_feature_info;

static int store_gene_feature(GenomeNode *gn, void *data, Env *env)
{
  GenomeFeature *gf;
  Store_gene_feature_info *info = (Store_gene_feature_info*) data;
  Range range;
  env_error_check(env);
  gf = genome_node_cast(genome_feature_class(), gn);
  assert(gf && info);
  switch (genome_feature_get_type(gf)) {
    case gft_mRNA:
      array_add(info->mRNAs, gf, env);
      break;
    case gft_exon:
      range = genome_node_get_range(gn);
      array_add(info->exons, range, env);
      break;
    default:
      assert(1);
  }
  return 0;
}

static bool genes_are_equal(GenomeNode *gn_1, GenomeNode *gn_2, Env *env)
{
  Array *exons_1, *exons_2, *mRNAs_1, *mRNAs_2;
  Store_gene_feature_info info;
  unsigned long i;
  bool equal;
  int has_err;

  /* init */
  exons_1 = array_new(sizeof (Range), env);
  exons_2 = array_new(sizeof (Range), env);
  mRNAs_1 = array_new(sizeof (GenomeNode*), env);
  mRNAs_2 = array_new(sizeof (GenomeNode*), env);

  /* get (direct) gene features */
  info.exons = exons_1;
  info.mRNAs = mRNAs_1;
  has_err = genome_node_traverse_direct_children(gn_1, &info,
                                                 store_gene_feature, env);
  assert(!has_err); /* cannot happen, store_gene_feature() is sane */
  info.exons = exons_2;
  info.mRNAs = mRNAs_2;
  has_err = genome_node_traverse_direct_children(gn_2, &info,
                                                 store_gene_feature, env);
  assert(!has_err); /* cannot happen, store_gene_feature() is sane */

  /* sort exon ranges */
  ranges_sort(exons_1);
  ranges_sort(exons_2);

  /* compare exon ranges */
  equal = ranges_are_equal(exons_1, exons_2);

  /* compare mRNAs, if necessary */
  if (equal && array_size(mRNAs_1) == array_size(mRNAs_2)) {
    /* sort mRNAs */
    genome_nodes_sort(mRNAs_1);
    genome_nodes_sort(mRNAs_2);
    for (i = 0; i < array_size(mRNAs_1); i++) {
      assert(equal);
      equal = mRNAs_are_equal(*(GenomeNode**) array_get(mRNAs_1, i),
                              *(GenomeNode**) array_get(mRNAs_2, i), env);
      if (!equal)
        break;
    }
  }

  /* free */
  array_delete(exons_1, env);
  array_delete(exons_2, env);
  array_delete(mRNAs_1, env);
  array_delete(mRNAs_2, env);

  return equal;
}

static void store_predicted_exon(TranscriptEvaluators *te, GenomeFeature *gf)
{
  assert(te && gf);
  evaluator_add_predicted(transcript_evaluators_get_all(te), 1);
  switch (genome_feature_get_transcriptfeaturetype(gf)) {
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
    case TRANSCRIPT_FEATURE_TYPE_UNDETERMINED: assert(0);
  }
}

/* adds exon only if necessary */
static void add_predicted_collapsed(Dlist *used_exons, Range *predicted_range,
                                    Evaluator *exon_evaluator_collapsed,
                                    Env *env)
{
  Range *used_range;
  if (!dlist_find(used_exons, predicted_range)) {
    used_range = env_ma_malloc(env, sizeof (Range));
    used_range->start = predicted_range->start;
    used_range->end = predicted_range->end;
    dlist_add(used_exons, used_range, env);
    evaluator_add_predicted(exon_evaluator_collapsed, 1);
  }
}

static void store_predicted_exon_collapsed(TranscriptUsedExons *used_exons,
                                           Range *predicted_range,
                                           TranscriptEvaluators *te,
                                           GenomeFeature *gf, Env *env)
{
  add_predicted_collapsed(transcript_used_exons_get_all(used_exons),
                          predicted_range, transcript_evaluators_get_all(te),
                          env);
  switch (genome_feature_get_transcriptfeaturetype(gf)) {
    case TRANSCRIPT_FEATURE_TYPE_SINGLE:
      add_predicted_collapsed(transcript_used_exons_get_single(used_exons),
                              predicted_range,
                              transcript_evaluators_get_single(te), env);
      break;
    case TRANSCRIPT_FEATURE_TYPE_INITIAL:
      add_predicted_collapsed(transcript_used_exons_get_initial(used_exons),
                              predicted_range,
                              transcript_evaluators_get_initial(te), env);
          break;
    case TRANSCRIPT_FEATURE_TYPE_INTERNAL:
      add_predicted_collapsed(transcript_used_exons_get_internal(used_exons),
                              predicted_range,
                              transcript_evaluators_get_internal(te), env);
          break;
    case TRANSCRIPT_FEATURE_TYPE_TERMINAL:
      add_predicted_collapsed(transcript_used_exons_get_terminal(used_exons),
                              predicted_range,
                              transcript_evaluators_get_terminal(te), env);
          break;
    case TRANSCRIPT_FEATURE_TYPE_UNDETERMINED: assert(0);
  }
}

static void determine_true_exon(GenomeNode *gn, Strand predicted_strand,
                                bool exondiff, Range *predicted_range,
                                Array *exons_forward,
                                Array *exons_reverse,
                                Array *true_exons_forward,
                                Array *true_exons_reverse,
                                Bittab *true_exons_forward_collapsed,
                                Bittab *true_exons_reverse_collapsed,
                                Evaluator *exon_evaluator,
                                Evaluator *exon_evaluator_collapsed)
{
  Range *actual_range;
  unsigned long num, *ctr_ptr;

  if ((actual_range = bsearch(predicted_range,
                              predicted_strand == STRAND_FORWARD
                              ? array_get_space(exons_forward)
                              : array_get_space(exons_reverse),
                              predicted_strand == STRAND_FORWARD
                              ? array_size(exons_forward)
                              : array_size(exons_reverse), sizeof (Range),
                              (Compare) range_compare_ptr))) {
    if (predicted_strand == STRAND_FORWARD) {
      num = actual_range - (Range*) array_get_space(exons_forward);
      ctr_ptr = array_get(true_exons_forward, num);
      if (*ctr_ptr) {
        (*ctr_ptr)--;
        evaluator_add_true(exon_evaluator);
      }
      if (true_exons_forward_collapsed &&
          !bittab_bit_is_set(true_exons_forward_collapsed, num)) {
        bittab_set_bit(true_exons_forward_collapsed, num);
        evaluator_add_true(exon_evaluator_collapsed);
      }
    }
    else {
      num = actual_range - (Range*) array_get_space(exons_reverse);
      ctr_ptr = array_get(true_exons_reverse, num);
      if (*ctr_ptr) {
        (*ctr_ptr)--;
        evaluator_add_true(exon_evaluator);
      }
      if (true_exons_reverse_collapsed &&
          !bittab_bit_is_set(true_exons_reverse_collapsed, num)) {
        bittab_set_bit(true_exons_reverse_collapsed, num);
        evaluator_add_true(exon_evaluator_collapsed);
      }
    }
  }
  else if (exondiff) {
    printf("> ");
    gff3_output_leading((GenomeFeature*) gn, NULL);
    printf(".\n");
  }
}

static void store_true_exon(GenomeNode *gn, Strand predicted_strand,
                            Range *predicted_range, bool exondiff,
                            TranscriptExons *exons_forward,
                            TranscriptExons *exons_reverse,
                            TranscriptCounts *counts_forward,
                            TranscriptCounts *counts_reverse,
                            TranscriptBittabs *exon_bittabs_forward,
                            TranscriptBittabs *exon_bittabs_reverse,
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
  switch (genome_feature_get_transcriptfeaturetype((GenomeFeature*) gn)) {
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
    case TRANSCRIPT_FEATURE_TYPE_UNDETERMINED: assert(0);
  }
}

static int process_predicted_feature(GenomeNode *gn, void *data, Env *env)
{
  Process_predicted_feature_info *info = (Process_predicted_feature_info*) data;
  Range predicted_range;
  unsigned long i, num;
  Strand predicted_strand;
  Array *real_genome_nodes;
  GenomeNode **real_gn;

  env_error_check(env);
  assert(gn && data);

  predicted_range = genome_node_get_range(gn);
  predicted_strand = genome_feature_get_strand((GenomeFeature*) gn);
  real_genome_nodes = array_new(sizeof (GenomeNode**), env);

  switch (genome_feature_get_type((GenomeFeature*) gn)) {
    case gft_gene:
      /* store predicted gene */
      evaluator_add_predicted(info->gene_evaluator, 1);
      /* determine true gene */
      switch (predicted_strand) {
        case STRAND_FORWARD:
        case STRAND_REVERSE:
          bsearch_all_mark(real_genome_nodes, &gn,
                           predicted_strand == STRAND_FORWARD
                           ? array_get_space(info->slot->genes_forward)
                           : array_get_space(info->slot->genes_reverse),
                           predicted_strand == STRAND_FORWARD
                           ? array_size(info->slot->genes_forward)
                           : array_size(info->slot->genes_reverse),
                           sizeof (GenomeNode*), (Compare) genome_node_compare,
                           predicted_strand == STRAND_FORWARD
                           ? info->slot->overlapped_genes_forward
                           : info->slot->overlapped_genes_reverse, env);
          if (array_size(real_genome_nodes)) {
            /* gene(s) with the same range found -> check if they are equal */
            for (i = 0; i < array_size(real_genome_nodes); i++) {
              real_gn = *(GenomeNode***) array_get(real_genome_nodes, i);
              if (genes_are_equal(gn, *real_gn, env)) {
                if (predicted_strand == STRAND_FORWARD) {
                  num = real_gn - (GenomeNode**)
                        array_get_space(info->slot->genes_forward);
                  if (!bittab_bit_is_set(info->slot->true_genes_forward, num)) {
                    bittab_set_bit(info->slot->true_genes_forward, num);
                    evaluator_add_true(info->gene_evaluator);
                    /*@loopbreak@*/
                    break;
                  }
                }
                else {
                  num = real_gn - (GenomeNode**)
                        array_get_space(info->slot->genes_reverse);
                  if (!bittab_bit_is_set(info->slot->true_genes_reverse, num)) {
                    bittab_set_bit(info->slot->true_genes_reverse, num);
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
            if (!genome_node_overlaps_nodes_mark(gn,
                                      predicted_strand == STRAND_FORWARD
                                      ? info->slot->genes_forward
                                      : info->slot->genes_reverse,
                                      predicted_strand == STRAND_FORWARD
                                      ? info->slot->overlapped_genes_forward
                                      : info->slot->overlapped_genes_reverse)) {
              (*info->wrong_genes)++;
            }
          }
          break;
        default:
          if (info->verbose) {
            fprintf(stderr, "skipping predicted gene with unknown orientation "
                    "(line %lu)\n", genome_node_get_line_number(gn));
          }
      }
      break;
    case gft_mRNA:
      /* store predicted mRNA */
      evaluator_add_predicted(info->mRNA_evaluator, 1);
      /* determine true mRNA */
      switch (predicted_strand) {
        case STRAND_FORWARD:
        case STRAND_REVERSE:
          bsearch_all_mark(real_genome_nodes, &gn,
                           predicted_strand == STRAND_FORWARD
                           ? array_get_space(info->slot->mRNAs_forward)
                           : array_get_space(info->slot->mRNAs_reverse),
                           predicted_strand == STRAND_FORWARD
                           ? array_size(info->slot->mRNAs_forward)
                           : array_size(info->slot->mRNAs_reverse),
                           sizeof (GenomeNode*), (Compare) genome_node_compare,
                           predicted_strand == STRAND_FORWARD
                           ? info->slot->overlapped_mRNAs_forward
                           : info->slot->overlapped_mRNAs_reverse, env);
          if (array_size(real_genome_nodes)) {
            /* mRNA(s) with the same range found -> check if they are equal */
            for (i = 0; i < array_size(real_genome_nodes); i++) {
              real_gn = *(GenomeNode***) array_get(real_genome_nodes, i);
              if (mRNAs_are_equal(gn, *real_gn, env)) {
                if (predicted_strand == STRAND_FORWARD) {
                  num = real_gn - (GenomeNode**)
                        array_get_space(info->slot->mRNAs_forward);
                  if (!bittab_bit_is_set(info->slot->true_mRNAs_forward, num)) {
                    bittab_set_bit(info->slot->true_mRNAs_forward, num);
                    evaluator_add_true(info->mRNA_evaluator);
                    /*@loopbreak@*/
                    break;
                  }
                }
                else {
                  num = real_gn - (GenomeNode**)
                        array_get_space(info->slot->mRNAs_reverse);
                  if (!bittab_bit_is_set(info->slot->true_mRNAs_reverse, num)) {
                    bittab_set_bit(info->slot->true_mRNAs_reverse, num);
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
            if (!genome_node_overlaps_nodes_mark(gn,
                                      predicted_strand == STRAND_FORWARD
                                      ? info->slot->mRNAs_forward
                                      : info->slot->mRNAs_reverse,
                                      predicted_strand == STRAND_FORWARD
                                      ? info->slot->overlapped_mRNAs_forward
                                      : info->slot->overlapped_mRNAs_reverse)) {
              (*info->wrong_mRNAs)++;
            }
          }
          break;
        default:
          if (info->verbose) {
            fprintf(stderr, "skipping predicted mRNA with unknown orientation "
                    "(line %lu)\n", genome_node_get_line_number(gn));
          }
      }
    break;
    case gft_exon:
      /* store predicted exon (mRNA level)*/
      store_predicted_exon(info->mRNA_exon_evaluators, (GenomeFeature*) gn);

      /* store predicted exon (mRNA level, collapsed) */
      store_predicted_exon_collapsed(predicted_strand == STRAND_FORWARD
                                     ? info->slot->used_mRNA_exons_forward
                                     : info->slot->used_mRNA_exons_reverse,
                                     &predicted_range,
                                     info->mRNA_exon_evaluators_collapsed,
                                     (GenomeFeature*) gn, env);

      /* determine true exon (mRNA level)*/
      switch (predicted_strand) {
        case STRAND_FORWARD:
        case STRAND_REVERSE:
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
          break;
        default:
          if (info->verbose) {
            fprintf(stderr, "skipping predicted exon with unknown orientation "
                    "(line %lu)\n", genome_node_get_line_number(gn));
          }
      }
      break;
    case gft_CDS:
      /* store predicted exon (CDS level)*/
      store_predicted_exon(info->CDS_exon_evaluators, (GenomeFeature*) gn);

      /* store predicted exon (CDS level, collapsed) */
      store_predicted_exon_collapsed(predicted_strand == STRAND_FORWARD
                                     ? info->slot->used_CDS_exons_forward
                                     : info->slot->used_CDS_exons_reverse,
                                     &predicted_range,
                                     info->CDS_exon_evaluators_collapsed,
                                     (GenomeFeature*) gn, env);

      /* determine true exon (CDS level) */
      switch (predicted_strand) {
        case STRAND_FORWARD:
        case STRAND_REVERSE:
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
          break;
        default:
          if (info->verbose) {
            fprintf(stderr, "skipping predicted exon with unknown orientation "
                    "(line %lu)\n", genome_node_get_line_number(gn));
          }
        }
      break;
    default:
      assert(1); /* shut up compiler */
  }
  array_delete(real_genome_nodes, env);
  return 0;
}

int determine_missing_features(void *key, void *value, void *data, Env *env)
{
  StreamEvaluator *se = (StreamEvaluator*) data;
  Slot *slot = (Slot*) value;
  env_error_check(env);
  assert(key && value && data);
  if (slot->overlapped_genes_forward) {
    se->missing_genes += bittab_size(slot->overlapped_genes_forward) -
                         bittab_count_set_bits(slot->overlapped_genes_forward);
  }
  if (slot->overlapped_genes_reverse) {
    se->missing_genes += bittab_size(slot->overlapped_genes_reverse) -
                         bittab_count_set_bits(slot->overlapped_genes_reverse);
  }
  if (slot->overlapped_mRNAs_forward) {
    se->missing_mRNAs += bittab_size(slot->overlapped_mRNAs_forward) -
                         bittab_count_set_bits(slot->overlapped_mRNAs_forward);
  }
  if (slot->overlapped_mRNAs_reverse) {
    se->missing_mRNAs += bittab_size(slot->overlapped_mRNAs_reverse) -
                         bittab_count_set_bits(slot->overlapped_mRNAs_reverse);
  }
  return 0;
}

int stream_evaluator_evaluate(StreamEvaluator *se, bool verbose, bool exondiff,
                              Env *env)
{
  GenomeNode *gn;
  SequenceRegion *sr;
  GenomeFeature *gf;
  Slot *slot;
  Process_real_feature_data process_real_feature_data;
  Process_predicted_feature_info info;
  int has_err;

  assert(se);
  env_error_check(env);

  /* init */
  process_real_feature_data.verbose = verbose;
  info.verbose = verbose;
  info.exondiff = exondiff;
  info.gene_evaluator = se->gene_evaluator;
  info.mRNA_evaluator = se->mRNA_evaluator;
  info.mRNA_exon_evaluators = se->mRNA_exon_evaluators;
  info.mRNA_exon_evaluators_collapsed = se->mRNA_exon_evaluators_collapsed;
  info.CDS_exon_evaluators = se->CDS_exon_evaluators;
  info.CDS_exon_evaluators_collapsed = se->CDS_exon_evaluators_collapsed;
  info.wrong_genes = &se->wrong_genes;
  info.wrong_mRNAs = &se->wrong_mRNAs;

  /* process the reality stream completely */
  while (!(has_err = genome_stream_next_tree(se->reality, &gn, env)) && gn) {
    sr = genome_node_cast(sequence_region_class(), gn);
    if (sr) {
      /* each sequence region gets its own ``slot'' */
      if (!(slot = hashtable_get(se->real_features,
                                 str_get(genome_node_get_seqid(gn))))) {

        slot = slot_new(env);
        hashtable_add(se->real_features, str_get(genome_node_get_seqid(gn)),
                      slot, env);
      }
      assert(slot);
    }
    gf = genome_node_cast(genome_feature_class(), gn);
    /* we consider only genome features */
    if (gf) {
      /* each sequence must have its own ``slot'' at this point */
      slot = hashtable_get(se->real_features,
                           str_get(genome_node_get_seqid(gn)));
      assert(slot);
      /* store the exons */
      process_real_feature_data.slot = slot;
      genome_feature_determine_transcripttypes(gf, env);
      has_err = genome_node_traverse_children(gn, &process_real_feature_data,
                                              process_real_feature, false,
                                              env);
      assert(!has_err); /* cannot happen, process_real_feature() is sane */
    }
    genome_node_rec_delete(gn, env);
  }

  /* set the actuals and sort them */
  if (!has_err) {
    has_err = hashtable_foreach(se->real_features, set_actuals_and_sort_them,
                                se, env);
  }

  /* process the prediction stream */
  if (!has_err) {
    while (!(has_err = genome_stream_next_tree(se->prediction, &gn, env)) &&
           gn) {
      gf = genome_node_cast(genome_feature_class(), gn);
      /* we consider only genome features */
      if (gf) {
        /* get (real) slot */
        slot = hashtable_get(se->real_features,
                             str_get(genome_node_get_seqid(gn)));
        if (slot) {
          info.slot = slot;
          genome_feature_determine_transcripttypes(gf, env);
          has_err = genome_node_traverse_children(gn, &info,
                                                  process_predicted_feature,
                                                  false, env);
          assert(!has_err); /* cannot happen, process_predicted_feature() is
                               sane */
        }
        else {
          /* we got no (real) slot */
          warning("sequence id \"%s\" (with predictions) not given in "
                  "``reality''", str_get(genome_node_get_seqid(gn)));
        }
      }
      genome_node_rec_delete(gn, env);
    }
  }

  /* determine the missing mRNAs */
  if (!has_err) {
    has_err = hashtable_foreach(se->real_features, determine_missing_features,
                                se, env);
  }

  return has_err;
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

void stream_evaluator_show(StreamEvaluator *se, FILE *outfp)
{
  assert(se);

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
}

void stream_evaluator_delete(StreamEvaluator *se, Env *env)
{
  if (!se) return;
  genome_stream_delete(se->reality, env);
  genome_stream_delete(se->prediction, env);
  hashtable_delete(se->real_features, env);
  evaluator_delete(se->gene_evaluator, env);
  evaluator_delete(se->mRNA_evaluator, env);
  transcript_evaluators_delete(se->mRNA_exon_evaluators, env);
  transcript_evaluators_delete(se->mRNA_exon_evaluators_collapsed, env);
  transcript_evaluators_delete(se->CDS_exon_evaluators, env);
  transcript_evaluators_delete(se->CDS_exon_evaluators_collapsed, env);
  env_ma_free(se, env);
}
