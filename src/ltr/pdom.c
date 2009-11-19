/*
  Copyright (c) 2008-2009 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2008-2009 Center for Bioinformatics, University of Hamburg

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

#ifdef HAVE_HMMER

#include <string.h>
#include <ctype.h>
#include <float.h>
#include <pthread.h>
#include "core/codon.h"
#include "core/log.h"
#include "core/ma.h"
#include "core/mathsupport.h"
#include "core/translator.h"
#include "core/undef.h"
#include "core/unused_api.h"
#include "extended/globalchaining.h"
#include "extended/reverse.h"
#include "ltr/pdom.h"

/* HMMER related includes */
#include "structs.h"
#include "globals.h"
#include "squid.h"
#include "funcs.h"

/* number of tophit_s structs to preallocate in HMMER back-end */
#define MAX_TOPHITS 50

struct GtPdomFinder {
  GtStrArray *hmm_files;
  GtArray *models;
  GtLTRElement *elem;
  double glob_eval_cutoff;
  unsigned int nof_threads,
               chain_max_gap_length;
};

struct GtPdomModel {
  struct plan7_s *model;
  struct threshold_s thresh;
  unsigned long reference_count;
};

struct GtPdomSingleHit {
  GtPhase phase;
  GtRange range;
  GtStr *aa_seq_model,
        *mline,
        *aa_seq_matched;
  double eval;
  GtPdomModelHit *mhit;
  unsigned long reference_count;
};

struct GtPdomModelHit {
  struct tophit_s *hits_fwd, *hits_rev;
  GtStrand strand;
  GtArray *best_chain;
  unsigned long reference_count;
  GtLTRElement *elem;
};

struct GtPdomResults {
  GtHashmap *domains;
  double combined_e_value_fwd,
         combined_e_value_rev;
  bool empty;
};

typedef struct GtPdomDomainTraverseInfo {
  GtPdomIteratorFunc func;
  void *data;
} GtPdomDomainTraverseInfo;

typedef struct GtPdomSharedMem {
  pthread_t *thread;
  int nof_threads;
  GtPdomFinder *gpf;
  unsigned long next_hmm;
  GtArray *hmms;
  GtStr **fwd, **rev;
  GtPdomResults *results;
  pthread_mutex_t in_lock, out_lock;
} GtPdomSharedMem;

/* --------------- GtPdomSingleHit ---------------------------------- */

void gt_pdom_single_hit_format_alignment(const GtPdomSingleHit *sh,
                                         unsigned long width,
                                         GtStr *dest)
{
  unsigned long pos,
                match_start,
                match_end,
                i,
                len, modellen, matchlen;
  char *buffer, cpbuf[BUFSIZ];
  const char *model,
             *mline,
             *matched;
  gt_assert(sh && width > 0 && dest && width < BUFSIZ);
  buffer = gt_calloc(width+1, sizeof (char));
  model = gt_str_get(sh->aa_seq_model);
  mline = gt_str_get(sh->mline);
  matched = gt_str_get(sh->aa_seq_matched);
  gt_str_reset(dest);

  modellen = gt_str_length(sh->aa_seq_model);
  matchlen = gt_str_length(sh->aa_seq_matched);
  len = MAX(modellen, matchlen)-3; /* -3 to remove the '<-*' */

  match_end = sh->range.start - 1;
  if (width > len)
    width = len;
  for (pos = 3; pos + width < len; pos += width) /* 3 to remove the '*->' */
  {
    match_start = match_end + 1;
    for (i = pos; i < pos + width && matched[i] != '\0'; i++)
      if (!isgap(matched[i])) match_end++;
    strncpy(buffer, model+pos, width);
    snprintf(cpbuf, BUFSIZ, "%6s %s\n", " ", buffer);
    gt_str_append_cstr(dest, cpbuf);
    strncpy(buffer, mline+pos, width);
    snprintf(cpbuf, BUFSIZ, "%6s %s\n", " ", buffer);
    gt_str_append_cstr(dest, cpbuf);
    strncpy(buffer, matched+pos, width);
    snprintf(cpbuf, BUFSIZ, "%6lu %s %6lu\n\n",
                            match_start,
                            buffer,
                            match_end);
    gt_str_append_cstr(dest, cpbuf);
  }
  if (pos < len)
  {
    unsigned long rest_len;
    rest_len = len - pos;
    match_start = match_end + 1;
    for (i = pos; i < len && matched[i] != '\0'; i++)
      if (!isgap(matched[i])) match_end++;
    strncpy(buffer, model+pos, rest_len);
    buffer[rest_len]='\0';
    snprintf(cpbuf, BUFSIZ, "%6s %s\n", " ", buffer);
    gt_str_append_cstr(dest, cpbuf);
    strncpy(buffer, mline+pos, rest_len);
    snprintf(cpbuf, BUFSIZ, "%6s %s\n", " ", buffer);
    gt_str_append_cstr(dest, cpbuf);
    strncpy(buffer, matched+pos, rest_len);
    snprintf(cpbuf, BUFSIZ, "%6lu %s %6lu\n\n",
                            match_start,
                            buffer,
                            match_end);
    gt_str_append_cstr(dest, cpbuf);
  }
  gt_free(buffer);
}

void gt_pdom_single_hit_get_aaseq(const GtPdomSingleHit *sh,
                                  GtStr *dest)
{
  const char *matched;
  unsigned long i;
  gt_assert(sh && dest);
  gt_str_reset(dest);

  matched = gt_str_get(sh->aa_seq_matched);
  for (i = 0; i < gt_str_length(sh->aa_seq_matched); i++)
  {
    if (!isgap(matched[i]))
      gt_str_append_char(dest, toupper(matched[i]));
  }
}

static GtPdomSingleHit* gt_pdom_single_hit_new(struct hit_s *singlehit,
                                               GtPdomModelHit *mhit)
{
  GtPdomSingleHit *hit;
  gt_assert(singlehit);
  hit = gt_calloc(1, sizeof (GtPdomSingleHit));
  hit->phase = gt_phase_get(singlehit->name[0]);
  hit->range.start = singlehit->sqfrom-1;
  hit->range.end = singlehit->sqto-1;
  hit->mhit = mhit;
  hit->eval = singlehit->pvalue;
  if (singlehit->ali && singlehit->ali->aseq)
    hit->aa_seq_matched = gt_str_new_cstr(singlehit->ali->aseq);
  if (singlehit->ali && singlehit->ali->model)
    hit->aa_seq_model = gt_str_new_cstr(singlehit->ali->model);
  if (singlehit->ali && singlehit->ali->mline)
    hit->mline = gt_str_new_cstr(singlehit->ali->mline);
  gt_assert(hit->range.start <= hit->range.end);

  return hit;
}

GtPdomSingleHit* gt_pdom_single_hit_ref(GtPdomSingleHit *sh)
{
  gt_assert(sh);
  sh->reference_count++;
  return sh;
}

GtPhase gt_pdom_single_hit_get_phase(const GtPdomSingleHit *singlehit)
{
  gt_assert(singlehit);
  return singlehit->phase;
}

GtRange gt_pdom_single_hit_get_range(const GtPdomSingleHit *singlehit)
{
  GtRange retrng;
  gt_assert(singlehit);
  switch (singlehit->mhit->strand)
  {
    case GT_STRAND_FORWARD:
    default:
      retrng.start = singlehit->mhit->elem->leftLTR_5
                        + (singlehit->range.start)*GT_CODON_LENGTH
                        + (unsigned long) singlehit->phase;
      retrng.end   =  retrng.start
                        + gt_range_length(&singlehit->range)*GT_CODON_LENGTH;
      break;
    case GT_STRAND_REVERSE:
      retrng.start = singlehit->mhit->elem->rightLTR_3
                        - (singlehit->range.end+1)*GT_CODON_LENGTH
                        - (unsigned long) singlehit->phase;
      retrng.end   =  retrng.start
                        + gt_range_length(&singlehit->range)*GT_CODON_LENGTH;
      break;
  }
  return retrng;
}

double gt_pdom_single_hit_get_evalue(const GtPdomSingleHit *singlehit)
{
  gt_assert(singlehit);
  return singlehit->eval;
}

void gt_pdom_single_hit_delete(void *value)
{
  GtPdomSingleHit *sh;
  if (!value) return;
  sh = (GtPdomSingleHit*) value;
  if (sh->reference_count) {
    sh->reference_count--;
    return;
  }
  gt_str_delete(sh->aa_seq_matched);
  gt_str_delete(sh->aa_seq_model);
  gt_str_delete(sh->mline);
  gt_free(sh);
}

/* --------------- GtPdomModelHit ---------------------------------- */

static GtPdomModelHit* gt_pdom_model_hit_new(GtLTRElement *elem)
{
  GtPdomModelHit *hit;
  hit = gt_calloc(1, sizeof (GtPdomModelHit));
  hit->hits_fwd   = AllocTophits(MAX_TOPHITS);
  hit->hits_rev   = AllocTophits(MAX_TOPHITS);
  hit->strand     = GT_STRAND_UNKNOWN;
  hit->elem       = elem;
  hit->best_chain = gt_array_new(sizeof (GtPdomSingleHit*));
  return hit;
}

GtPdomModelHit *gt_pdom_model_hit_ref(GtPdomModelHit *mh)
{
  gt_assert(mh);
  mh->reference_count++;
  return mh;
}

void gt_pdom_model_hit_add_single_hit(GtPdomModelHit *mh,
                                      GtPdomSingleHit *sh)
{
  gt_assert(mh && mh->best_chain && sh);
  gt_array_add(mh->best_chain, sh);
}

GtStrand gt_pdom_model_hit_get_best_strand(const GtPdomModelHit *h)
{
  gt_assert(h);
  return h->strand;
}

unsigned long gt_pdom_model_hit_best_chain_length(const GtPdomModelHit *h)
{
  gt_assert(h);
  return gt_array_size(h->best_chain);
}

GtPdomSingleHit* gt_pdom_model_hit_best_single_hit(const GtPdomModelHit *h,
                                                   unsigned long nth)
{
  gt_assert(h);
  return *(GtPdomSingleHit**) gt_array_get(h->best_chain, nth);
}

void gt_pdom_model_hit_delete(void *value)
{
  GtPdomModelHit *hit;
  if (!value) return;
  hit = (GtPdomModelHit*) value;
  if (hit->reference_count) {
    hit->reference_count--;
    return;
  }
  if (hit->hits_fwd)
    FreeTophits(hit->hits_fwd);
  if (hit->hits_rev)
    FreeTophits(hit->hits_rev);
  if (hit->best_chain)
  {
    unsigned long i;
    for (i=0;i<gt_array_size(hit->best_chain);i++)
    {
      gt_pdom_single_hit_delete(*(GtPdomSingleHit**)
                                 gt_array_get(hit->best_chain, i));
    }
    gt_array_delete(hit->best_chain);
  }
  gt_free(hit);
}

/* --------------- GtPdomModel ------------------------------------ */

static void print_cutoffs(GtPdomModel *m)
{
  if (!m) return;
  gt_log_log("cutoffs for model %s are", m->model->name);
  gt_log_log("globT = %f    domT = %f", m->thresh.globT, m->thresh.domT);
  gt_log_log("globE = %f    domE = %f", m->thresh.globE, m->thresh.domE);
}

int gt_pdom_model_set_trusted_cutoff(GtPdomModel *m, GtError *err)
{
  int had_err = 0;
  gt_assert(m);
  gt_error_check(err);
  gt_log_log("setting trusted cutoffs for model %s", m->model->name);
  m->thresh.autocut = CUT_TC;
  if (!SetAutocuts(&m->thresh, m->model)) {
    gt_error_set(err, "HMM %s does not contain a TC cutoff!", m->model->name);
    had_err = -1;
  }
  print_cutoffs(m);
  return had_err;
}

int gt_pdom_model_set_gathering_cutoff(GtPdomModel *m, GtError *err)
{
  int had_err = 0;
  gt_assert(m);
  gt_error_check(err);
  gt_log_log("setting gathering cutoffs for model %s", m->model->name);
  m->thresh.autocut = CUT_GA;
  if (!SetAutocuts(&m->thresh, m->model)) {
    gt_error_set(err, "HMM %s does not contain a GA cutoff!", m->model->name);
    had_err = -1;
  }
  print_cutoffs(m);
  return had_err;
}

void gt_pdom_model_reset_cutoffs(GtPdomModel *m)
{
  gt_assert(m);
  m->thresh.globT = m->thresh.domT = -FLT_MAX;
  m->thresh.domE = m->thresh.globE = FLT_MAX;
  m->thresh.autocut = CUT_NONE;
}

static GtPdomModel* gt_pdom_model_new(struct plan7_s *model)
{
  GtPdomModel *newmodel;
  gt_assert(model);
  newmodel = gt_calloc(1, sizeof (GtPdomModel));
  newmodel->model = model;
  /* initialize with no cutoffs */
  gt_pdom_model_reset_cutoffs(newmodel);
  return newmodel;
}

int gt_pdom_model_set_evalue_cutoff(GtPdomModel *m, double eval_cutoff,
                                    GT_UNUSED GtError *err)
{
  int had_err = 0;
  gt_assert(m);
  gt_error_check(err);
  gt_log_log("setting evalue cutoff for model %s to %f", m->model->name,
                                                         eval_cutoff);
  m->thresh.globE   = eval_cutoff;
  print_cutoffs(m);
  return had_err;
}

GtPdomModel *gt_pdom_model_ref(GtPdomModel *m)
{
  gt_assert(m);
  m->reference_count++;
  return m;
}

const char* gt_pdom_model_get_name(const GtPdomModel *model)
{
  gt_assert(model);
  return model->model->name;
}

const char* gt_pdom_model_get_acc(const GtPdomModel *model)
{
  gt_assert(model);
  return model->model->acc;
}

struct threshold_s* gt_pdom_model_get_thresholds(const GtPdomModel *model)
{
  gt_assert(model);
  return (struct threshold_s*) &model->thresh;
}

void gt_pdom_model_delete(void *value)
{
  GtPdomModel *model;
  if (!value) return;
  model = (GtPdomModel*) value;
  if (model->reference_count) {
    model->reference_count--;
    return;
  }
  FreePlan7(model->model);
  gt_free(model);
}

/* --------------- GtPdomResults ---------------------------------- */

static GtPdomResults* gt_pdom_results_new(void)
{
  GtPdomResults *res;
  res = gt_calloc(1, sizeof (GtPdomResults));
  res->domains = gt_hashmap_new(HASH_DIRECT, gt_pdom_model_delete,
                                gt_pdom_model_hit_delete);
  res->empty = TRUE;
  res->combined_e_value_fwd = res->combined_e_value_rev = 0.0;
  return res;
}

bool gt_pdom_results_empty(GtPdomResults *res)
{
  gt_assert(res);
  return res->empty;
}

double gt_pdom_results_get_combined_evalue_fwd(GtPdomResults *result)
{
  gt_assert(result);
  return result->combined_e_value_fwd;
}

double gt_pdom_results_get_combined_evalue_rev(GtPdomResults *result)
{
  gt_assert(result);
  return result->combined_e_value_rev;
}

void gt_pdom_results_delete(GtPdomResults *res)
{
  if (!res) return;
  gt_hashmap_delete(res->domains);
  gt_free(res);
}

/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */

static int visit_domain(void *key, void *value, void *data,
                        GT_UNUSED GtError *err)
{
  GtPdomModel *model = (GtPdomModel *) key;
  GtPdomModelHit *hit = (GtPdomModelHit*) value;
  GtPdomDomainTraverseInfo *info = (GtPdomDomainTraverseInfo*) data;
  return info->func(model, hit, info->data, err);
}

int gt_pdom_results_foreach_domain_hit(GtPdomResults *results,
                                       GtPdomIteratorFunc func,
                                       void *data,
                                       GtError *err)
{
  GtPdomDomainTraverseInfo info;
  info.func = func;
  info.data = data;
  return gt_hashmap_foreach(results->domains,
                            visit_domain,
                            &info,
                            err);
}

void gt_hmmer_search(struct plan7_s *hmm,
                     char *seq,
                     unsigned long seqlen,
                     char *seqdesc,
                     struct threshold_s *thresh,
                     struct tophit_s *ghit,
                     struct tophit_s *dhit,
                     pthread_mutex_t *lock)
{
  struct dpmatrix_s *mx;
  struct p7trace_s *tr;
  unsigned char *dsq;
  float sc;
  int retval;
  double pvalue, evalue;

  mx = CreatePlan7Matrix(1, hmm->M, 25, 0);

  dsq = DigitizeSequence(seq, seqlen);

  if (P7ViterbiSpaceOK(seqlen, hmm->M, mx))
    sc = P7Viterbi(dsq, seqlen, hmm, mx, &tr);
  else
    sc = P7SmallViterbi(dsq, seqlen, hmm, mx, &tr);

  sc -= TraceScoreCorrection(hmm, tr, dsq);

  if ((retval = pthread_mutex_lock(lock)) != 0)
  {
    fprintf(stderr, "pthread_mutex_lock failure: %s\n", strerror(retval));
    exit(EXIT_FAILURE);
  }
  pvalue = PValue(hmm, sc);
  evalue = thresh->Z ? (double) thresh->Z * pvalue : (double) pvalue;

  if (sc >= thresh->globT && evalue <= thresh->globE)
  {
    sc = PostprocessSignificantHit(ghit, dhit,
                                   tr, hmm, dsq, seqlen,
                                   seqdesc, NULL, NULL,
                                   false, sc, true, thresh, FALSE);
  }
  if ((retval = pthread_mutex_unlock(lock)) != 0)
  {
    fprintf(stderr, "pthread_mutex_unlock failure: %s\n", strerror(retval));
    exit(EXIT_FAILURE);
  }
  P7FreeTrace(tr);
  free(dsq);
  FreePlan7Matrix(mx);
}

static void chainproc(GtChain *c, GtFragment *f,
                      GT_UNUSED unsigned long nof_frags,
                      GT_UNUSED unsigned long gap_length, void *data)
{
  unsigned long i;
  GtPdomModelHit *hit;
  hit = (GtPdomModelHit*) data;

  gt_log_log("resulting chain has %ld GtFragments, score %ld",
             gt_chain_size(c),
             gt_chain_get_score(c));
  for (i=0;i<gt_chain_size(c);i++)
  {
    GtPdomSingleHit *sh;
    GtFragment frag;
    frag = f[gt_chain_get_fragnum(c, i)];
    gt_log_log("(%lu %lu) (%lu %lu), %p", frag.startpos1, frag.endpos1,
                                          frag.startpos2, frag.endpos2,
                                          frag.data);
    sh = gt_pdom_single_hit_new(frag.data, hit);
    gt_pdom_model_hit_add_single_hit(hit, sh);
  }
  gt_log_log("\n");
}

static int gt_fragcmp(const void *frag1, const void *frag2)
{
  GtFragment *f1 = (GtFragment*) frag1;
  GtFragment *f2 = (GtFragment*) frag2;
  if (f1->startpos2 == f2->startpos2)
    return 0;
  else return (f1->startpos2 < f2->startpos2 ? -1 : 1);
}

void* gt_pdom_per_domain_worker_thread(void *data)
{
  GtPdomSharedMem *shared;
  GtPdomModel *hmm;
  int retval;
  struct tophit_s *ghit = NULL, *hits = NULL;
  bool best_fwd = TRUE;
  unsigned long i;
  GtFragment *frags;
  struct threshold_s* thresh;
  GtPdomModelHit *hit;

  shared = (GtPdomSharedMem*) data;

  for (;;)
  {
    /* Lock input */
    if ((retval = pthread_mutex_lock(&shared->in_lock)) != 0)
    {
      fprintf(stderr, "Failed to lock: %s\n", strerror(retval));
      exit(EXIT_FAILURE);
    }
    /* Have all HMMs been distributed? If so, we are done here. */
    if (shared->next_hmm == gt_array_size(shared->gpf->models))
    {
      if ((retval = pthread_mutex_unlock(&shared->in_lock)) != 0)
      {
        fprintf(stderr, "Failed to unlock: %s\n", strerror(retval));
        exit(EXIT_FAILURE);
      }
      pthread_exit(NULL);
    }

    /* Get work from HMM list */
    hmm = *(GtPdomModel**) gt_array_get(shared->gpf->models,
                                        shared->next_hmm);
    shared->next_hmm++;

    /* work claimed, release input */
    if ((retval = pthread_mutex_unlock(&shared->in_lock)) != 0)
    {
      fprintf(stderr, "pthread_mutex_unlock failure: %s\n", strerror(retval));
      exit(EXIT_FAILURE);
    }

    hit = gt_pdom_model_hit_new(shared->gpf->elem);
    ghit = AllocTophits(MAX_TOPHITS);
    thresh = gt_pdom_model_get_thresholds(hmm);

    gt_hmmer_search(hmm->model,
                    gt_str_get(shared->fwd[0]),
                    gt_str_length(shared->fwd[0]),
                    "0+",
                    thresh,
                    ghit,
                    hit->hits_fwd,
                    &shared->out_lock);
    gt_hmmer_search(hmm->model,
                    gt_str_get(shared->fwd[1]),
                    gt_str_length(shared->fwd[1]),
                    "1+",
                    thresh,
                    ghit,
                    hit->hits_fwd,
                    &shared->out_lock);
    gt_hmmer_search(hmm->model,
                    gt_str_get(shared->fwd[2]),
                    gt_str_length(shared->fwd[2]),
                    "2+",
                    thresh,
                    ghit,
                    hit->hits_fwd,
                    &shared->out_lock);
    gt_hmmer_search(hmm->model,
                    gt_str_get(shared->rev[0]),
                    gt_str_length(shared->rev[0]),
                    "0-",
                    thresh,
                    ghit,
                    hit->hits_rev,
                    &shared->out_lock);
    gt_hmmer_search(hmm->model,
                    gt_str_get(shared->rev[1]),
                    gt_str_length(shared->rev[1]),
                    "1-",
                    thresh,
                    ghit,
                    hit->hits_rev,
                    &shared->out_lock);
    gt_hmmer_search(hmm->model,
                    gt_str_get(shared->rev[2]),
                    gt_str_length(shared->rev[2]),
                    "2-",
                    thresh,
                    ghit,
                    hit->hits_rev,
                    &shared->out_lock);

    FullSortTophits(hit->hits_fwd);
    FullSortTophits(hit->hits_rev);

    /* check if there were any hits */
    if (hit->hits_fwd->num > 0 || hit->hits_rev->num > 0)
    {
      shared->results->empty = FALSE;
      if (hit->hits_fwd->num > 0)
      {
        if (hit->hits_rev->num > 0)
        {
          if (gt_double_compare(hit->hits_fwd->hit[0]->score,
                                hit->hits_rev->hit[0]->score) < 0)
            best_fwd = FALSE;
        }
        else best_fwd = TRUE;
      }
      else best_fwd = FALSE;

      /* determine best-scoring strand */
      hits = (best_fwd ? hit->hits_fwd : hit->hits_rev);
      gt_assert(hits);

      /* no need to chain if there is only one hit */
      if (hits->num > 1)
      {
        /* create GtFragment set for chaining */
        frags = (GtFragment*) gt_calloc(hits->num, sizeof (GtFragment));
        for (i=0;i<hits->num;i++)
        {
          frags[i].startpos1 = hits->hit[i]->hmmfrom;
          frags[i].endpos1   = hits->hit[i]->hmmto;
          frags[i].startpos2 = hits->hit[i]->sqfrom;
          frags[i].endpos2   = hits->hit[i]->sqto;
          /* let weight(f) be targetlength(f) multiplied by hmmerscore(f) */
          frags[i].weight    = (hits->hit[i]->sqto - hits->hit[i]->sqfrom + 1)
                                 * hits->hit[i]->score;
          frags[i].data      = hits->hit[i];
        }

        /* sort GtFragments by position */
        qsort(frags, hits->num, sizeof (GtFragment), gt_fragcmp);
        for (i=0;i<hits->num;i++)
        {
          gt_log_log("(%lu %lu) (%lu %lu) %p", frags[i].startpos1,
                                               frags[i].endpos1,
                                               frags[i].startpos2,
                                               frags[i].endpos2,
                                               frags[i].data);
        }
        gt_log_log("chaining %d frags", hits->num);
        /* do chaining */
        gt_globalchaining_max(frags, hits->num,
                              shared->gpf->chain_max_gap_length,
                              chainproc, hit);
        gt_free(frags);
      }
      else
      {
        GtPdomSingleHit *sh = gt_pdom_single_hit_new(hits->hit[0],
                                                     hit);
        gt_pdom_model_hit_add_single_hit(hit, sh);
      }

      /* Lock results, we want to write to the result hashtable */
      if ((retval = pthread_mutex_lock(&(shared->out_lock))) != 0)
      {
        fprintf(stderr, "Failed to lock: %s\n", strerror(retval));
        exit(EXIT_FAILURE);
      }

      /* register results */
      gt_hashmap_add(shared->results->domains,
                     gt_pdom_model_ref(hmm),
                     hit);
      if (best_fwd)
      {
        shared->results->combined_e_value_fwd
          += log(hit->hits_fwd->hit[0]->pvalue);
        hit->strand = GT_STRAND_FORWARD;
      }
      else
      {
        shared->results->combined_e_value_rev
          += log(hit->hits_rev->hit[0]->pvalue);
        hit->strand = GT_STRAND_REVERSE;
       }

      /* unlock results */
      if ((retval = pthread_mutex_unlock(&(shared->out_lock))) != 0)
      {
        fprintf(stderr, "Failed to unlock: %s\n", strerror(retval));
        exit(EXIT_FAILURE);
      }
    }
    else
      gt_pdom_model_hit_delete(hit);

    FreeTophits(ghit);
  }
}

static GtPdomSharedMem* gt_pdom_run_threads(GtArray *hmms, int nof_threads,
                                            GtStr **fwd, GtStr **rev,
                                            GtPdomResults *results,
                                            GtPdomFinder *gpf)
{
  int retval, i;
  GtPdomSharedMem *shared;
  pthread_attr_t attr;

  gt_assert(hmms && nof_threads > 0 && fwd && rev && results && gpf);

  shared = gt_calloc(1, sizeof (GtPdomSharedMem));

  shared->nof_threads = nof_threads;
  shared->thread = gt_calloc(nof_threads, sizeof (pthread_t));
  shared->hmms = hmms;
  shared->fwd = fwd;
  shared->rev = rev;
  shared->next_hmm = 0;
  shared->results = results;
  shared->gpf = gpf;

  if ((retval = pthread_mutex_init(&shared->in_lock, NULL)) != 0)
  {
    fprintf(stderr, "Could not initialize lock! %s\n", strerror(retval));
    exit(EXIT_FAILURE);
  }
  if ((retval = pthread_mutex_init(&shared->out_lock, NULL)) != 0)
  {
    fprintf(stderr, "Could not initialize lock! %s\n", strerror(retval));
    exit(EXIT_FAILURE);
  }

  pthread_attr_init(&attr);

  for (i = 0; i < nof_threads; i++)
  {
    if ((retval = pthread_create(&(shared->thread[i]), &attr,
           gt_pdom_per_domain_worker_thread, (void *) shared)) != 0)
    {
      fprintf(stderr,"Failed to create thread %d, return code %d\n", i, retval);
      exit(EXIT_FAILURE);
    }
  }
  pthread_attr_destroy(&attr);
  return shared;
}

static void gt_pdom_free_shared(GtPdomSharedMem *shared)
{
  gt_free(shared->thread);
  gt_free(shared);
}

GtPdomResults* gt_pdom_finder_find(GtPdomFinder *gpf, const char *seq,
                                   const char *rev_seq, GtLTRElement *element,
                                   GtError *err)
{
  GtStr *fwd[3], *rev[3];
  char translated;
  unsigned long seqlen;
  GtPdomSharedMem *shared = NULL;
  GtPdomResults *results = NULL;
  GtTranslator *tr;
  int i,
      had_err = 0;
  unsigned int frame;
  gt_assert(seq && rev_seq && strlen(seq) == strlen(rev_seq) && element);
  gt_error_check(err);

  seqlen = gt_ltrelement_length(element);
  gpf->elem = element;
  results = gt_pdom_results_new();
  tr = gt_translator_new();

  for (i=0;i<3;i++)
  {
    fwd[i] = gt_str_new();
    rev[i] = gt_str_new();
  }

  /* create translations */
  had_err = gt_translator_start(tr, seq, seqlen, &translated,
                                &frame, err);
  while (!had_err && translated)
  {
    gt_str_append_char(fwd[frame], translated);
    had_err = gt_translator_next(tr, &translated, &frame, NULL);
  }
  if (!had_err)
  {
    had_err = gt_translator_start(tr, rev_seq, seqlen, &translated,
                                  &frame, err);
    while (!had_err && translated)
    {
      gt_str_append_char(rev[frame], translated);
      had_err = gt_translator_next(tr, &translated, &frame, NULL);
    }
  }

  if (!had_err)
  {
    for (i=0;i<3;i++)
    {
      gt_assert(gt_str_length(fwd[i]) == (seqlen - i) / 3);
      gt_assert(gt_str_length(rev[i]) == (seqlen - i) / 3);
    }

    /* start worker threads */
    shared = gt_pdom_run_threads(gpf->models, gpf->nof_threads,
                                 fwd, rev, results, gpf);

    /* continue when all threads are done */
    for (i = 0; i < shared->nof_threads; i++)
    {
      if (pthread_join(shared->thread[i], NULL) != 0)
      {
        fprintf(stderr, "Could not join threads!");
        exit(EXIT_FAILURE);
      }
    }

    /* cleanup */
    gt_pdom_free_shared(shared);
  }

  SqdClean();
  for (i=0;i<3;i++)
  {
    gt_str_delete(fwd[i]);
    gt_str_delete(rev[i]);
  }

  gt_translator_delete(tr);
  return results;
}

/* --------------- GtPdomFinder ------------------------------------- */

static int gt_pdom_finder_load_files(GtPdomFinder *gpf, GtStrArray *files,
                                     GtPdomCutoff cutoff, GtError *err)
{
  int had_err = 0;
  GtPdomModel *gpm = NULL;
  unsigned long i;
  gt_assert(gpf && files && err);
  for (i=0;i<gt_str_array_size(files);i++)
  {
    struct plan7_s *hmm;
    HMMFILE *hmmfp;
    char *hmmfile = (char*) gt_str_array_get(files, i);
    gt_log_log("trying to load HMM file '%s'...", hmmfile);
    if ((hmmfp = HMMFileOpen(hmmfile, "HMMERDB")) == NULL)
    {
      gt_error_set(err, "Failed to open HMM file '%s'", hmmfile);
      had_err = -1;
      if (hmmfp) HMMFileClose(hmmfp);
      break;
    }
    if (!had_err && !HMMFileRead(hmmfp, &hmm))
    {
      gt_error_set(err, "Failed to read any HMMs from file '%s'", hmmfile);
      had_err = -1;
      if (hmmfp) HMMFileClose(hmmfp);
      break;
    }
    if (!had_err && hmm == NULL)
    {
      had_err = -1;
      gt_error_set(err, "HMM file '%s' corrupt or in incorrect format?",
                   hmmfile);
      if (hmmfp) HMMFileClose(hmmfp);
      break;
    }
    if (!had_err)
    {
      gpm = gt_pdom_model_new(hmm);
      if (cutoff == GT_PHMM_CUTOFF_TC) {
        had_err = gt_pdom_model_set_trusted_cutoff(gpm, err);
      }
      if (!had_err && (cutoff == GT_PHMM_CUTOFF_GA)) {
        had_err = gt_pdom_model_set_gathering_cutoff(gpm, err);
      }
      if (!had_err) {
        had_err = gt_pdom_model_set_evalue_cutoff(gpm, gpf->glob_eval_cutoff,
                                                  err);
      }
      if (had_err) {
        if (hmmfp) HMMFileClose(hmmfp);
        break;
      }
    }
    if (!had_err) {
      P7Logoddsify(hmm, true);
      gt_array_add(gpf->models, gpm);
    }
    if (hmmfp) HMMFileClose(hmmfp);
  }
  return had_err;
}

static void gt_pdom_clear_hmms(GtArray *hmms)
{
  unsigned long i;
  if (!hmms) return;
  for (i=0;i<gt_array_size(hmms);i++)
  {
    gt_pdom_model_delete(*(GtPdomModel**) gt_array_get(hmms,i));
  }
  gt_array_delete(hmms);
}

GtPdomFinder* gt_pdom_finder_new(GtStrArray *hmmfiles, double eval_cutoff,
                                 unsigned int nof_threads,
                                 unsigned int chain_max_gap_length,
                                 GtPdomCutoff cutoff,
                                 GtError *e)
{
  int had_err = 0;
  GtPdomFinder *gpf;
  gt_assert(hmmfiles && e && nof_threads);
  gpf = gt_calloc(1, sizeof (GtPdomFinder));
  gpf->nof_threads = nof_threads;
  gpf->hmm_files = gt_str_array_ref(hmmfiles);
  gpf->chain_max_gap_length = chain_max_gap_length;
  gpf->models = gt_array_new(sizeof (struct GtPdomModel*));
  gpf->glob_eval_cutoff = eval_cutoff;
  had_err = gt_pdom_finder_load_files(gpf, hmmfiles, cutoff, e);
  if (had_err)
  {
    gt_pdom_clear_hmms(gpf->models);
    gt_str_array_delete(gpf->hmm_files);
    gt_free(gpf);
    return NULL;
  } else return gpf;
}

void gt_pdom_finder_delete(GtPdomFinder *gpf)
{
  if (!gpf) return;
  gt_pdom_clear_hmms(gpf->models);
  gt_str_array_delete(gpf->hmm_files);
  gt_free(gpf);
}

#endif
