/*
  Copyright (c) 2008-2010 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2008-2010 Center for Bioinformatics, University of Hamburg

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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <float.h>
#include <unistd.h>

#include "core/codon_api.h"
#include "core/codon_iterator_simple_api.h"
#include "core/log.h"
#include "core/ma.h"
#include "core/mathsupport.h"
#include "core/minmax.h"
#include "core/thread_api.h"
#include "core/translator.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "extended/globalchaining.h"
#include "extended/reverse_api.h"
#include "ltr/pdom.h"

/* HMMER3 related includes */
#include "p7_config.h"
#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_sq.h"
#include "esl_ssi.h"
#include "esl_sqio.h"
#include "esl_stopwatch.h"
#include "p7_config.h"
#ifdef HMMER_THREADS
#include "esl_threads.h"
#include "esl_workqueue.h"
#endif
#include "hmmer.h"

struct GtPdomFinder {
  GtStrArray *hmm_files;
  GtArray *models;
  GtLTRElement *elem;
  ESL_GETOPTS *getopts;
  double glob_eval_cutoff;
  bool reportall;
  unsigned int chain_max_gap_length;
};

struct GtPdomModel {
  P7_HMM *model;
  ESL_ALPHABET *abc;
  unsigned long reference_count;
};

struct GtPdomSingleHit {
  GtPhase phase;
  GtRange range;
  GtStr *aa_seq_model,
        *mline,
        *aa_seq_matched;
  double eval;
  bool chained;
  GtPdomModelHit *mhit;
  unsigned long reference_count;
  GtArray *chains;
};

struct GtPdomModelHit {
  P7_TOPHITS *hits_fwd, *hits_rev;
  GtStrand strand;
  GtArray *members;
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
  GtPdomFinder *gpf;
  unsigned long next_hmm;
  GtArray *hmms;
  GtStr **fwd, **rev;
  GtPdomResults *results;
  GtMutex *in_lock, *out_lock;
} GtPdomSharedMem;

/* --------------- GtPdomSingleHit ---------------------------------- */

#define gt_pdom_isgap(c) ((c) == ' ' || (c) == '.' || (c) == '_' \
                            || (c) == '-' || (c) == '~')

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
  len = MAX(modellen, matchlen);

  match_end = sh->range.start - 1;
  if (width > len)
    width = len;
  for (pos = 0; pos + width < len; pos += width)
  {
    match_start = match_end + 1;
    for (i = pos; i < pos + width && matched[i] != '\0'; i++)
      if (!gt_pdom_isgap(matched[i])) match_end++;
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
      if (!gt_pdom_isgap(matched[i])) match_end++;
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
    if (!gt_pdom_isgap(matched[i])) {
      /* replace stop codons by 'X'es */
      if (matched[i] == '*') {
        gt_str_append_char(dest, 'X');
      } else {
        gt_str_append_char(dest, toupper(matched[i]));
      }
    }
  }
}

static GtPdomSingleHit* gt_pdom_single_hit_new(P7_DOMAIN *singlehit,
                                               GtPdomModelHit *mhit)
{
  GtPdomSingleHit *hit;
  gt_assert(singlehit && singlehit->ad);
  hit = gt_calloc(1, sizeof (GtPdomSingleHit));
  hit->phase = gt_phase_get(singlehit->ad->sqname[0]);
  hit->range.start = singlehit->ad->sqfrom-1;
  hit->range.end = singlehit->ad->sqto-1;
  gt_log_log("create single hit %lu-%lu", hit->range.start, hit->range.end);
  gt_assert(gt_range_length(&(hit->range)) > 0);
  hit->mhit = mhit;
  hit->eval = singlehit->pvalue;
  if (singlehit->ad && singlehit->ad->aseq)
    hit->aa_seq_matched = gt_str_new_cstr(singlehit->ad->aseq);
  if (singlehit->ad && singlehit->ad->model)
    hit->aa_seq_model = gt_str_new_cstr(singlehit->ad->model);
  if (singlehit->ad && singlehit->ad->mline)
    hit->mline = gt_str_new_cstr(singlehit->ad->mline);
  gt_assert(hit->range.start <= hit->range.end);
  hit->chains = gt_array_new(sizeof (unsigned long));
  hit->chained = false;
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
                        + (singlehit->range.start) * GT_CODON_LENGTH
                        + (unsigned long) singlehit->phase;
      retrng.end   =  retrng.start
                        + (gt_range_length(&singlehit->range) - 1)
                           * GT_CODON_LENGTH;
      break;
    case GT_STRAND_REVERSE:
      retrng.start = singlehit->mhit->elem->rightLTR_3
                        - (singlehit->range.end+1) * GT_CODON_LENGTH
                        - (unsigned long) singlehit->phase;
      retrng.end   =  retrng.start
                        + (gt_range_length(&singlehit->range) - 1)
                           * GT_CODON_LENGTH;
      break;
  }
  gt_log_log("range: (%lu-%lu, phase %d) -> %lu-%lu (len %lu) strand: %c\n",
             singlehit->range.start, singlehit->range.end,
             singlehit->phase, retrng.start, retrng.end,
             gt_range_length(&retrng),
             GT_STRAND_CHARS[singlehit->mhit->strand]);
  return retrng;
}

double gt_pdom_single_hit_get_evalue(const GtPdomSingleHit *singlehit)
{
  gt_assert(singlehit);
  return singlehit->eval;
}

void gt_pdom_single_hit_set_chained(GtPdomSingleHit *singlehit,
                                    unsigned long cno)
{
  gt_assert(singlehit);
  gt_array_add(singlehit->chains, cno);
  gt_log_log("adding single hit %lu-%lu to chain %lu", singlehit->range.start,
                                                       singlehit->range.end,
                                                       cno);
  singlehit->chained = true;
}

bool gt_pdom_single_hit_is_chained(GtPdomSingleHit *singlehit)
{
  gt_assert(singlehit);
  return singlehit->chained;
}

GtArray* gt_pdom_single_hit_get_chains(GtPdomSingleHit *singlehit)
{
  gt_assert(singlehit);
  return singlehit->chains;
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
  gt_array_delete(sh->chains);
  gt_free(sh);
}

/* --------------- GtPdomModelHit ---------------------------------- */

static GtPdomModelHit* gt_pdom_model_hit_new(GtLTRElement *elem)
{
  GtPdomModelHit *hit;
  hit = gt_calloc(1, sizeof (GtPdomModelHit));
  hit->hits_fwd   = p7_tophits_Create();
  hit->hits_rev   = p7_tophits_Create();
  hit->strand     = GT_STRAND_UNKNOWN;
  hit->elem       = elem;
  hit->members = gt_array_new(sizeof (GtPdomSingleHit*));
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
  gt_assert(mh && mh->members && sh);
  gt_array_add(mh->members, sh);
}

GtStrand gt_pdom_model_hit_get_best_strand(const GtPdomModelHit *h)
{
  gt_assert(h);
  return h->strand;
}

unsigned long gt_pdom_model_hit_num_of_single_hits(const GtPdomModelHit *h)
{
  gt_assert(h);
  return gt_array_size(h->members);
}

GtPdomSingleHit* gt_pdom_model_hit_single_hit(const GtPdomModelHit *h,
                                          unsigned long nth)
{
  gt_assert(h);
  return *(GtPdomSingleHit**) gt_array_get(h->members, nth);
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
    p7_tophits_Destroy(hit->hits_fwd);
  if (hit->hits_rev)
    p7_tophits_Destroy(hit->hits_rev);
  if (hit->members)
  {
    unsigned long i;
    for (i=0;i<gt_array_size(hit->members);i++)
    {
      gt_pdom_single_hit_delete(*(GtPdomSingleHit**)
                                 gt_array_get(hit->members, i));
    }
    gt_array_delete(hit->members);
  }
  gt_free(hit);
}

/* --------------- GtPdomModel ------------------------------------ */

static GtPdomModel* gt_pdom_model_new(P7_HMM *model)
{
  GtPdomModel *newmodel;
  gt_assert(model);
  newmodel = gt_calloc(1, sizeof (GtPdomModel));
  newmodel->model = model;
  return newmodel;
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

void gt_pdom_model_delete(void *value)
{
  GtPdomModel *model;
  if (!value) return;
  model = (GtPdomModel*) value;
  if (model->reference_count) {
    model->reference_count--;
    return;
  }
  if (model->abc) esl_alphabet_Destroy(model->abc);
  p7_hmm_Destroy(model->model);
  gt_free(model);
}

/* --------------- GtPdomResults ---------------------------------- */

static GtPdomResults* gt_pdom_results_new(void)
{
  GtPdomResults *res;
  res = gt_calloc(1, sizeof (GtPdomResults));
  res->domains = gt_hashmap_new(GT_HASH_DIRECT, gt_pdom_model_delete,
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

static inline void gt_hmmer_search(P7_PIPELINE *pli,
                                   char *in_seq,
                                   ESL_ALPHABET *abc,
                                   P7_OPROFILE *om,
                                   P7_BG *bg,
                                   char *seqdesc,
                                   P7_TOPHITS *hits,
                                   GT_UNUSED GtMutex *lock)
{
  ESL_SQ *sq = NULL;
  /* create and encode sequence */
  sq = esl_sq_CreateFrom(seqdesc, in_seq, "\0", "\0", NULL);
  esl_sq_Digitize(abc, sq);
  /* configure for new sequence */
  p7_pli_NewSeq(pli, sq);
  p7_bg_SetLength(bg, sq->n);
  p7_oprofile_ReconfigLength(om, sq->n);
  /* process this sequence */
  p7_Pipeline(pli, om, bg, sq, hits);
  /* cleanup */
  esl_sq_Destroy(sq);  /* XXX: reuse me? */
  p7_pipeline_Reuse(pli);
}

static void chainproc(GtChain *c, GtFragment *f,
                      GT_UNUSED unsigned long nof_frags,
                      GT_UNUSED unsigned long gap_length, void *data)
{
  unsigned long i,
                *chainno = (unsigned long*) data;

  gt_log_log("resulting chain has %ld GtFragments, score %ld",
             gt_chain_size(c),
             gt_chain_get_score(c));
  for (i = 0; i < gt_chain_size(c); i++)
  {
    GtFragment frag;
    frag = f[gt_chain_get_fragnum(c, i)];
    gt_log_log("(%lu %lu) (%lu %lu)", frag.startpos1, frag.endpos1,
                                          frag.startpos2, frag.endpos2);
    gt_pdom_single_hit_set_chained((GtPdomSingleHit*) frag.data, *chainno);
  }
  (*chainno)++;
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

static void* gt_pdom_per_domain_worker_thread(void *data)
{
  GtPdomSharedMem *shared;
  GtPdomModel *hmm;
  P7_PIPELINE *pli;
  P7_TOPHITS *hits = NULL;
  bool best_fwd = TRUE;
  GtFragment *frags;
  GtPdomModelHit *hit;

  shared = (GtPdomSharedMem*) data;

  P7_BG *bg = NULL;
  P7_PROFILE *gm = NULL;
  P7_OPROFILE *om = NULL;

  for (;;)
  {
    unsigned long h = 0,
                  d = 0,
                  nof_domains = 0;
    P7_DOMAIN *lastdom = NULL;

    gt_mutex_lock(shared->in_lock);
    /* Have all HMMs been distributed? If so, we are done here. */
    if (shared->next_hmm == gt_array_size(shared->gpf->models))
    {
      gt_mutex_unlock(shared->in_lock);
      return NULL;
    }

    /* Get work from HMM list */
    hmm = *(GtPdomModel**) gt_array_get(shared->gpf->models,
                                        shared->next_hmm);
    shared->next_hmm++;

    /* work claimed, release input */
    gt_mutex_unlock(shared->in_lock);

    /* prepare HMM for searching */
    pli = p7_pipeline_Create(shared->gpf->getopts, hmm->model->M, 400,
                             p7_SEARCH_SEQS);

    bg = p7_bg_Create(hmm->abc);
    gm = p7_profile_Create(hmm->model->M, hmm->abc);
    om = p7_oprofile_Create(hmm->model->M, hmm->abc);
    p7_ProfileConfig(hmm->model, bg, gm, 400, p7_LOCAL);
    p7_oprofile_Convert(gm, om);
    p7_pli_NewModel(pli, om, bg);

    hit = gt_pdom_model_hit_new(shared->gpf->elem);

    gt_hmmer_search(pli,
                    gt_str_get(shared->fwd[0]),
                    hmm->abc,
                    om,
                    bg,
                    "0+",
                    hit->hits_fwd,
                    shared->out_lock);
    gt_hmmer_search(pli,
                    gt_str_get(shared->fwd[1]),
                    hmm->abc,
                    om,
                    bg,
                    "1+",
                    hit->hits_fwd,
                    shared->out_lock);
    gt_hmmer_search(pli,
                    gt_str_get(shared->fwd[2]),
                    hmm->abc,
                    om,
                    bg,
                    "2+",
                    hit->hits_fwd,
                    shared->out_lock);
    gt_hmmer_search(pli,
                    gt_str_get(shared->rev[0]),
                    hmm->abc,
                    om,
                    bg,
                    "0-",
                    hit->hits_rev,
                    shared->out_lock);
    gt_hmmer_search(pli,
                    gt_str_get(shared->rev[1]),
                    hmm->abc,
                    om,
                    bg,
                    "1-",
                    hit->hits_rev,
                    shared->out_lock);
    gt_hmmer_search(pli,
                    gt_str_get(shared->rev[2]),
                    hmm->abc,
                    om,
                    bg,
                    "2-",
                    hit->hits_rev,
                    shared->out_lock);

    p7_tophits_Sort(hit->hits_fwd);
    p7_tophits_Sort(hit->hits_rev);
    p7_tophits_Threshold(hit->hits_fwd, pli);
    p7_tophits_Threshold(hit->hits_rev, pli);

    p7_pipeline_Destroy(pli);
    p7_oprofile_Destroy(om);
    p7_profile_Destroy(gm);
    p7_bg_Destroy(bg);

    /* check if there were any hits */
    if (hit->hits_fwd->nreported > 0 || hit->hits_rev->nreported > 0)
    {
      gt_log_log("N: %llu/%llu, nreported: %llu/%llu, nincluded: %llu/%llu",
                  (unsigned long long) hit->hits_fwd->N,
                  (unsigned long long) hit->hits_rev->N,
                  (unsigned long long) hit->hits_fwd->nreported,
                  (unsigned long long) hit->hits_rev->nreported,
                  (unsigned long long) hit->hits_fwd->nincluded,
                  (unsigned long long) hit->hits_rev->nincluded);
      shared->results->empty = FALSE;
      if (hit->hits_fwd->nreported > 0)
      {
        if (hit->hits_rev->nreported > 0)
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

      /* determine number of domain hits, track single domain hits */
      for (h = 0; h < hits->N; h++) {
        if (hits->hit[h]->flags & p7_IS_REPORTED) {
          for (d = 0; d < hits->hit[h]->ndom; d++) {
            if (hits->hit[h]->dcl[d].is_reported) {
              nof_domains++;
              lastdom = &(hits->hit[h]->dcl[d]);
              gt_log_log("domain %s, from: %lu, to %lu, on %s (%lu/%d)",
                         hits->hit[h]->dcl[d].ad->hmmname,
                         hits->hit[h]->dcl[d].ad->sqfrom,
                         hits->hit[h]->dcl[d].ad->sqto,
                         hits->hit[h]->dcl[d].ad->sqname,
                         d+1,
                         hits->hit[h]->ndom);
            }
          }
        }
      }

      /* no need to chain if there is only one domain hit */
      if (nof_domains > 1)
      {
        unsigned long i = 0;
        gt_log_log("nof_domains: %lu", nof_domains);
        /* create GtFragment set for chaining */
        frags = (GtFragment*) gt_calloc(nof_domains, sizeof (GtFragment));
        for (h = 0; h < hits->N; h++) {
          if (hits->hit[h]->flags & p7_IS_REPORTED) {
            for (d = 0; d < hits->hit[h]->ndom; d++) {
              if (hits->hit[h]->dcl[d].is_reported) {
                P7_DOMAIN *dom = &(hits->hit[h]->dcl[d]);
                GtPdomSingleHit *sh = gt_pdom_single_hit_new(dom, hit);
                gt_pdom_model_hit_add_single_hit(hit, sh);
                frags[i].startpos1 = dom->ad->hmmfrom;
                frags[i].endpos1   = dom->ad->hmmto;
                frags[i].startpos2 = dom->ad->sqfrom;
                frags[i].endpos2   = dom->ad->sqto;
                /* let weight(f) be targetlength(f) multiplied by bitscore(f) */
                frags[i].weight    = (dom->ad->sqto - dom->ad->sqfrom + 1)
                                         * dom->bitscore;
                frags[i].data      = sh;
                i++;
              }
            }
          }
        }

        /* sort GtFragments by position */
        qsort(frags, nof_domains, sizeof (GtFragment), gt_fragcmp);
        for (i=0;i<nof_domains;i++)
        {
          gt_log_log("(%lu %lu) (%lu %lu)", frags[i].startpos1,
                                               frags[i].endpos1,
                                               frags[i].startpos2,
                                               frags[i].endpos2);
        }
        gt_log_log("chaining %lu frags", nof_domains);

        unsigned long chainno = 0;
        /* do chaining */
        gt_globalchaining_max(frags, nof_domains,
                              shared->gpf->chain_max_gap_length,
                              chainproc, &chainno);
        gt_free(frags);
      }
      else
      {
        /* single domain hit, no chaining necessary */
        if (lastdom != NULL) {
          GtPdomSingleHit *sh;
          sh = gt_pdom_single_hit_new(lastdom, hit);
          gt_pdom_model_hit_add_single_hit(hit, sh);
          gt_pdom_single_hit_set_chained(sh, 0UL);
        }
      }

      /* Lock results, we want to write to the result hashtable */
      gt_mutex_lock(shared->out_lock);

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
      gt_mutex_unlock(shared->out_lock);
    }
    else
      gt_pdom_model_hit_delete(hit);

  }
}

static void gt_pdom_run_threads(GtArray *hmms,
                                GtStr **fwd, GtStr **rev,
                                GtPdomResults *results,
                                GtPdomFinder *gpf)
{
  GtPdomSharedMem *shared;
  gt_assert(hmms && fwd && rev && results && gpf);

  shared = gt_calloc(1, sizeof (GtPdomSharedMem));

  shared->hmms = hmms;
  shared->fwd = fwd;
  shared->rev = rev;
  shared->next_hmm = 0;
  shared->results = results;
  shared->gpf = gpf;
  shared->in_lock = gt_mutex_new();
  shared->out_lock = gt_mutex_new();

  gt_multithread(gt_pdom_per_domain_worker_thread, shared, NULL);

  gt_mutex_delete(shared->in_lock);
  gt_mutex_delete(shared->out_lock);
  gt_free(shared);
}

GtPdomResults* gt_pdom_finder_find(GtPdomFinder *gpf, const char *seq,
                                   const char *rev_seq, GtLTRElement *element,
                                   GtError *err)
{
  GtStr *fwd[3], *rev[3];
  char translated;
  unsigned long seqlen;
  GtCodonIterator *ci;
  GtPdomResults *results = NULL;
  GtTranslatorStatus status;
  GtTranslator *tr;
  int i,
      had_err = 0;
  unsigned int frame;
  gt_assert(seq && rev_seq && strlen(seq) == strlen(rev_seq) && element);
  gt_error_check(err);

  seqlen = gt_ltrelement_length(element);
  gpf->elem = element;
  results = gt_pdom_results_new();

  for (i=0;i<3;i++)
  {
    fwd[i] = gt_str_new();
    rev[i] = gt_str_new();
  }

  ci = gt_codon_iterator_simple_new(seq, seqlen, NULL);
  gt_assert(ci);
  tr = gt_translator_new(ci);
  /* create translations */
  status = gt_translator_next(tr, &translated, &frame, err);
  while (status == GT_TRANSLATOR_OK && translated)
  {
    gt_str_append_char(fwd[frame], translated);
    status = gt_translator_next(tr, &translated, &frame, NULL);
  }
  if (status == GT_TRANSLATOR_ERROR) had_err = -1;

  if (!had_err)
  {
    gt_codon_iterator_delete(ci);
    ci = gt_codon_iterator_simple_new(rev_seq, seqlen, NULL);
    gt_translator_set_codon_iterator(tr, ci);
    status = gt_translator_next(tr, &translated, &frame, err);
    while (status == GT_TRANSLATOR_OK && translated)
    {
      gt_str_append_char(rev[frame], translated);
      status = gt_translator_next(tr, &translated, &frame, NULL);
    }
    if (status == GT_TRANSLATOR_ERROR) had_err = -1;
  }

  if (!had_err)
  {
    for (i=0;i<3;i++)
    {
      gt_assert(gt_str_length(fwd[i]) == (seqlen - i) / 3);
      gt_assert(gt_str_length(rev[i]) == (seqlen - i) / 3);
    }

    /* start worker threads */
    gt_pdom_run_threads(gpf->models, fwd, rev, results, gpf);
  }

  for (i=0;i<3;i++)
  {
    gt_str_delete(fwd[i]);
    gt_str_delete(rev[i]);
  }

  gt_codon_iterator_delete(ci);
  gt_translator_delete(tr);
  return results;
}

/* --------------- GtPdomFinder ------------------------------------- */

static int gt_pdom_finder_load_files(GtPdomFinder *gpf, GtStrArray *files,
                                     GtError *err)
{
  int had_err = 0;
  GtPdomModel *gpm = NULL;
  unsigned long i;
  gt_assert(gpf && files && err);
  for (i=0;i<gt_str_array_size(files);i++)
  {
    P7_HMM *hmm       = NULL;
    P7_HMMFILE *hmmfp = NULL;
    ESL_ALPHABET *abc = NULL;
    char *hmmfile = (char*) gt_str_array_get(files, i);
    gt_log_log("trying to load HMM file '%s'...", hmmfile);
    if (p7_hmmfile_Open(hmmfile, NULL, &hmmfp) != eslOK)
    {
      gt_error_set(err, "Failed to open HMM file '%s'", hmmfile);
      had_err = -1;
      break;
    }
    if (!had_err && p7_hmmfile_Read(hmmfp, &abc, &hmm) != eslOK)
    {
      gt_error_set(err, "Failed to read any HMMs from file '%s'", hmmfile);
      had_err = -1;
      if (hmmfp) p7_hmmfile_Close(hmmfp);
      break;
    }
    if (!had_err && hmm == NULL)
    {
      had_err = -1;
      gt_error_set(err, "HMM file '%s' corrupt or in incorrect format?",
                   hmmfile);
      if (hmmfp) p7_hmmfile_Close(hmmfp);
      break;
    }
    if (!had_err)
    {
      gpm = gt_pdom_model_new(hmm);
      gpm->abc = abc;
      if (had_err) {
        if (hmmfp) p7_hmmfile_Close(hmmfp);
        break;
      }
    }
    if (!had_err) {
      gt_log_log("HMM file '%s' loaded.", hmmfile);
      gt_array_add(gpf->models, gpm);
    }
    if (hmmfp) p7_hmmfile_Close(hmmfp);
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

static ESL_OPTIONS gt_pdom_hmmer3_options[] = {
  /* name           type         default   env  range   toggles   reqs
   * incomp help docgroup*/
  { "-h", eslARG_NONE, FALSE, NULL, NULL, NULL,  NULL,
    NULL, "show brief help on version and usage", 0 },
  { "-E", eslARG_REAL, "10.0", NULL, "x>0", NULL,
    NULL,  "--cut_ga,--cut_nc,--cut_tc",
    "E-value cutoff for reporting significant sequence hits", 0 },
 { "-T", eslARG_REAL, FALSE, NULL, "x>0", NULL, NULL,
   "--cut_ga,--cut_nc,--cut_tc",
   "bit score cutoff for reporting significant sequence hits", 0 },
 { "-Z", eslARG_REAL, FALSE, NULL, "x>0", NULL, NULL,
   NULL, "set # of comparisons done, for E-value calculation", 0 },
 { "--domE", eslARG_REAL,"1000.0", NULL, "x>0", NULL, NULL,
   "--cut_ga,--cut_nc,--cut_tc",
   "E-value cutoff for reporting individual domains", 0 },
 { "--domT", eslARG_REAL, FALSE, NULL, "x>0", NULL, NULL,
   "--cut_ga,--cut_nc,--cut_tc",
   "bit score cutoff for reporting individual domains", 0 },
 { "--domZ", eslARG_REAL, FALSE, NULL, "x>0", NULL, NULL, NULL,
   "set # of significant seqs, for domain E-value calculation", 0 },
 { "--cut_ga", eslARG_NONE, FALSE, NULL, NULL, NULL, NULL,
   "--seqE,--seqT,--domE,--domT",
   "use GA gathering threshold bit score cutoffs in <hmmfile>", 0 },
 { "--cut_nc", eslARG_NONE, FALSE, NULL, NULL, NULL, NULL,
   "--seqE,--seqT,--domE,--domT",
   "use NC noise threshold bit score cutoffs in <hmmfile>", 0 },
 { "--cut_tc", eslARG_NONE, FALSE, NULL, NULL, NULL, NULL,
   "--seqE,--seqT,--domE,--domT",
   "use TC trusted threshold bit score cutoffs in <hmmfile>", 0 },
 { "--max", eslARG_NONE, FALSE, NULL, NULL, NULL, NULL,
   "--F1,--F2,--F3", "Turn all heuristic filters off (less speed, more power)",
   0 },
 { "--F1", eslARG_REAL, "0.02", NULL, NULL, NULL, NULL,
   "--max", "Stage 1 (MSV) threshold: promote hits w/ P <= F1", 0 },
 { "--F2", eslARG_REAL, "1e-3", NULL, NULL, NULL, NULL,
   "--max", "Stage 2 (Vit) threshold: promote hits w/ P <= F2", 0 },
 { "--F3", eslARG_REAL, "1e-5", NULL, NULL, NULL, NULL,
   "--max", "Stage 3 (Fwd) threshold: promote hits w/ P <= F3", 0 },
 { "--nobias", eslARG_NONE, NULL, NULL, NULL, NULL, NULL,
   "--max", "turn off composition bias filter", 0 },
 { "--nonull2", eslARG_NONE, NULL, NULL, NULL, NULL, NULL, NULL,
   "turn off biased composition score corrections", 0 },
 { "--incE", eslARG_REAL,  "0.01", NULL, "x>0",     NULL,  NULL,  NULL,
   "consider sequences <= this E-value threshold as significant", 5 },
 { "--incdomE",    eslARG_REAL,  "0.01", NULL, "x>0",     NULL,  NULL,  NULL,
   "consider domains <= this E-value threshold as significant", 5 },
 { "--incT",       eslARG_REAL,   FALSE, NULL, NULL,      NULL,    NULL,  NULL,
   "consider sequences >= this score threshold as significant", 5 },
 { "--incdomT",    eslARG_REAL,   FALSE, NULL,  NULL,     NULL,  NULL,  NULL,
   "consider domains >= this score threshold as significant", 5 },
 { "--seed", eslARG_INT, "42", NULL, "n>=0", NULL, NULL, NULL,
 "set RNG seed to <n> (if 0: one-time arbitrary seed)", 0 },
 { "--acc", eslARG_NONE, FALSE, NULL, NULL, NULL, NULL, NULL,
   "output target accessions instead of names if possible", 0 },
 { "--noali", eslARG_NONE, FALSE, NULL, NULL, NULL,  NULL,  NULL,
   "don't output alignments, so output is smaller", 2 },
 {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

GtPdomFinder* gt_pdom_finder_new(GtStrArray *hmmfiles, double eval_cutoff,
                                 unsigned int chain_max_gap_length,
                                 GtPdomCutoff cutoff, GtError *err)
{
  int had_err = 0;
  GtPdomFinder *gpf;
  GtStr *cmd;   /* for HMMER3 parameterization, we need to construct a
                   virtual command line for Easel to parse, see below */
  gt_assert(hmmfiles && err);
  gpf = gt_calloc(1, sizeof (GtPdomFinder));
  gpf->getopts = esl_getopts_Create(gt_pdom_hmmer3_options);
  gpf->hmm_files = gt_str_array_ref(hmmfiles);
  gpf->chain_max_gap_length = chain_max_gap_length;
  gpf->models = gt_array_new(sizeof (struct GtPdomModel*));
  gpf->glob_eval_cutoff = eval_cutoff;
  had_err = gt_pdom_finder_load_files(gpf, hmmfiles, err);

  if (!had_err) {
    cmd = gt_str_new_cstr("--noali ");
    switch (cutoff) {
      case GT_PHMM_CUTOFF_GA:
        gt_str_append_cstr(cmd, "--cut_ga ");
        break;
      case GT_PHMM_CUTOFF_TC:
        gt_str_append_cstr(cmd, "--cut_tc ");
        break;
      case GT_PHMM_CUTOFF_NONE:
        gt_str_append_cstr(cmd, "-E ");
        gt_str_append_double(cmd, eval_cutoff, 100);
        gt_str_append_cstr(cmd, " ");
        break;
    }
    had_err = esl_opt_ProcessSpoof(gpf->getopts, gt_str_get(cmd));
    gt_str_delete(cmd);
  }

  p7_FLogsumInit();
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
  if (gpf->getopts) esl_getopts_Destroy(gpf->getopts);
  gt_str_array_delete(gpf->hmm_files);
  gt_free(gpf);
}

#endif
