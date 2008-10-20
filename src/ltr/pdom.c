/*
  Copyright (c) 2008 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2008 Center for Bioinformatics, University of Hamburg

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
#include "core/log.h"
#include "core/ma.h"
#include "core/mathsupport.h"
#include "core/translate.h"
#include "core/undef.h"
#include "core/unused_api.h"
#include "extended/globalchaining.h"
#include "extended/reverse.h"
#include "ltr/pdom.h"

/* HMMER related includes */
#include "globals.h"
#include "squid.h"
#include "funcs.h"

/* number of tophit_s structs to preallocate in HMMER back-end */
#define MAX_TOPHITS 50

struct GtPdomHit {
  struct tophit_s *hits_fwd, *hits_rev;
  GtStrand strand;
  GtArray *best_chain;
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
  unsigned long next_hmm;
  GtArray *hmms;
  char *fwd_fr1, *fwd_fr2, *fwd_fr3,
       *rev_fr1, *rev_fr2, *rev_fr3;
  GtPdomResults *results;
  GtPdomOptions *opts;
  pthread_mutex_t in_lock, out_lock;
} GtPdomSharedMem;

GtPdomHit* gt_pdom_hit_new(void)
{
  GtPdomHit *hit;
  hit = gt_calloc(1, sizeof (GtPdomHit));
  hit->hits_fwd   = AllocTophits(MAX_TOPHITS);
  hit->hits_rev   = AllocTophits(MAX_TOPHITS);
  hit->strand     = GT_STRAND_UNKNOWN;
  hit->best_chain = gt_array_new(sizeof (struct hit_s*));
  return hit;
}

GtStrand gt_pdom_hit_get_strand(const GtPdomHit *h)
{
  gt_assert(h);
  return h->strand;
}

GtArray* gt_pdom_hit_get_best_chain(const GtPdomHit *h)
{
  gt_assert(h);
  return h->best_chain;
}

void gt_pdom_hit_delete(void *value)
{
  GtPdomHit *hit;
  if (!value) return;
  hit = (GtPdomHit*) value;
  if (hit->hits_fwd)
    FreeTophits(hit->hits_fwd);
  if (hit->hits_rev)
    FreeTophits(hit->hits_rev);
  if (hit->best_chain)
    gt_array_delete(hit->best_chain);
  gt_free(hit);
}

static GtPdomResults* gt_pdom_results_new(void)
{
  GtPdomResults *res;
  res = gt_calloc(1, sizeof (GtPdomResults));
  res->domains = gt_hashmap_new(HASH_DIRECT,
                                NULL,
                                gt_pdom_hit_delete);
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

static int visit_domain(void *key, void *value, void *data,
                        GT_UNUSED GtError *err)
{
  struct plan7_s *model = (struct plan7_s *) key;
  GtPdomHit *hit = (GtPdomHit*) value;
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

int gt_pdom_load_hmm_files(GtPdomOptions *opts, GtError *err)
{
  unsigned long i;
  int had_err = 0;

  gt_assert(opts && err);

  for (i=0;i<gt_str_array_size(opts->hmm_files);i++)
  {
    struct plan7_s *hmm;
    HMMFILE *hmmfp;
    char *hmmfile = (char*) gt_str_array_get(opts->hmm_files, i);
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
      P7Logoddsify(hmm, true);
      gt_array_add(opts->plan7_ts, hmm);
    }
    if (!had_err && !SetAutocuts(&opts->thresh, hmm))
    {
      gt_error_set(err, "HMM %s did not contain the GA, TC, or NC "
                        "cutoffs you needed",
                   hmm->name);
      had_err = -1;
      if (hmmfp) HMMFileClose(hmmfp);
      break;
    }
    if (hmmfp) HMMFileClose(hmmfp);
  }
  gt_log_log("Loaded %lu HMM model(s)",
             gt_array_size(opts->plan7_ts));
  return had_err;
}

static void chainproc(GtChain *c, Fragment *f, void *data)
{
  unsigned long i;
  GtPdomHit *hit;
  hit = (GtPdomHit*) data;

  gt_log_log("resulting chain has %ld fragments, score %ld",
             gt_chain_size(c),
             gt_chain_get_score(c));
  for (i=0;i<gt_chain_size(c);i++)
  {
    Fragment frag = f[gt_chain_get_fragnum(c, i)];
    gt_log_log("(%lu %lu) (%lu %lu), %p", frag.startpos1, frag.endpos1,
                                          frag.startpos2, frag.endpos2,
                                          frag.data);
    gt_array_add(hit->best_chain, frag.data);
  }
  gt_log_log("\n");
}

static int gt_fragcmp(const void *frag1, const void *frag2)
{
  Fragment *f1 = (Fragment*) frag1;
  Fragment *f2 = (Fragment*) frag2;
  if (f1->startpos2 == f2->startpos2)
    return 0;
  else return (f1->startpos2 < f2->startpos2 ? -1 : 1);
}

void* gt_pdom_per_domain_worker_thread(void *data)
{
  GtPdomSharedMem *shared;
  struct plan7_s *hmm;
  int retval;
  struct tophit_s *ghit = NULL, *hits = NULL;
  bool best_fwd = TRUE;
  unsigned long i;
  Fragment *frags;
  GtPdomHit *hit;

  shared = (GtPdomSharedMem *) data;

  for (;;)
  {
    /* Lock input */
    if ((retval = pthread_mutex_lock(&shared->in_lock)) != 0)
    {
      fprintf(stderr, "Failed to lock: %s\n", strerror(retval));
      exit(EXIT_FAILURE);
    }
    /* Have all HMMs been distributed? If so, we are done here. */
    if (shared->next_hmm == gt_array_size(shared->opts->plan7_ts))
    {
      if ((retval = pthread_mutex_unlock(&shared->in_lock)) != 0)
      {
        fprintf(stderr, "Failed to unlock: %s\n", strerror(retval));
        exit(EXIT_FAILURE);
      }
      pthread_exit(NULL);
    }

    /* Get work from HMM list */
    hmm = *(struct plan7_s**) gt_array_get(shared->opts->plan7_ts,
                                           shared->next_hmm);
    shared->next_hmm++;

    /* work claimed, release input */
    if ((retval = pthread_mutex_unlock(&shared->in_lock)) != 0)
    {
      fprintf(stderr, "pthread_mutex_unlock failure: %s\n", strerror(retval));
      exit(EXIT_FAILURE);
    }

    hit = gt_malloc(sizeof (GtPdomHit));
    ghit = AllocTophits(MAX_TOPHITS);
    hit->hits_fwd = AllocTophits(MAX_TOPHITS);
    hit->hits_rev = AllocTophits(MAX_TOPHITS);
    hit->best_chain = gt_array_new(sizeof (struct hit_s*));

    gt_hmmer_search(hmm,shared->fwd_fr1,strlen(shared->fwd_fr1),"0+",
                 &shared->opts->thresh, ghit, hit->hits_fwd, &shared->out_lock);
    gt_hmmer_search(hmm,shared->fwd_fr2,strlen(shared->fwd_fr2),"1+",
                 &shared->opts->thresh, ghit, hit->hits_fwd, &shared->out_lock);
    gt_hmmer_search(hmm,shared->fwd_fr3,strlen(shared->fwd_fr3),"2+",
                 &shared->opts->thresh, ghit, hit->hits_fwd, &shared->out_lock);
    gt_hmmer_search(hmm,shared->rev_fr1,strlen(shared->rev_fr1),"0-",
                 &shared->opts->thresh, ghit, hit->hits_rev, &shared->out_lock);
    gt_hmmer_search(hmm,shared->rev_fr2,strlen(shared->rev_fr2),"1-",
                 &shared->opts->thresh, ghit, hit->hits_rev, &shared->out_lock);
    gt_hmmer_search(hmm,shared->rev_fr3,strlen(shared->rev_fr3),"2-",
                 &shared->opts->thresh, ghit, hit->hits_rev, &shared->out_lock);

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
        /* create fragment set for chaining */
        frags = (Fragment*) gt_calloc(hits->num, sizeof (Fragment));
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

        /* sort fragments by position */
        qsort(frags, hits->num, sizeof (Fragment), gt_fragcmp);
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
                              shared->opts->chain_max_gap_length,
                              chainproc, hit);
        gt_free(frags);
      }
      else
      {
        gt_array_add(hit->best_chain, hits->hit[0]);
      }

      /* Lock results, we want to write to the result hashtable */
      if ((retval = pthread_mutex_lock(&(shared->out_lock))) != 0)
      {
        fprintf(stderr, "Failed to lock: %s\n", strerror(retval));
        exit(EXIT_FAILURE);
      }

      /* register results */
      gt_hashmap_add(shared->results->domains, hmm, hit);
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
      gt_pdom_hit_delete(hit);

    FreeTophits(ghit);
  }
}

static GtPdomSharedMem* gt_pdom_run_threads(GtArray *hmms, int nof_threads,
                                          char *fwd_fr1,char *fwd_fr2,
                                            char *fwd_fr3,
                                          char *rev_fr1,char *rev_fr2,
                                            char *rev_fr3,
                                          GtPdomResults *results,
                                          GtPdomOptions *opts)
{
  int retval, i;
  GtPdomSharedMem *shared;
  pthread_attr_t attr;

  gt_assert(hmms && nof_threads > 0 && *fwd_fr1 && *fwd_fr2 && *fwd_fr3
          && *rev_fr1 && rev_fr2 && rev_fr3 && results && opts);

  shared = gt_calloc(1, sizeof (GtPdomSharedMem));

  shared->nof_threads = nof_threads;
  shared->thread = gt_calloc(nof_threads, sizeof (pthread_t));
  shared->hmms = hmms;
  shared->fwd_fr1 = fwd_fr1;
  shared->fwd_fr2 = fwd_fr2;
  shared->fwd_fr3 = fwd_fr3;
  shared->rev_fr1 = rev_fr1;
  shared->rev_fr2 = rev_fr2;
  shared->rev_fr3 = rev_fr3;
  shared->next_hmm = 0;
  shared->results = results;
  shared->opts = opts;

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

GtPdomResults* gt_pdom_find(const char *seq, const char *rev_seq,
                             GtLTRElement *element, GtPdomOptions *opts)
{
  char *fwd_fr1, *fwd_fr2, *fwd_fr3,
       *rev_fr1, *rev_fr2, *rev_fr3;
  unsigned long seqlen = gt_ltrelement_length(element);
  GtPdomSharedMem *shared = NULL;
  GtPdomResults *results = NULL;
  int i;

  gt_assert(seq && rev_seq && element && opts);

  results = gt_pdom_results_new();

  /* create translations */
  gt_translate_all_frames(&fwd_fr1, &fwd_fr2, &fwd_fr3,     seq, seqlen);
  gt_translate_all_frames(&rev_fr1, &rev_fr2, &rev_fr3, rev_seq, seqlen);

  /* start worker threads */
  shared = gt_pdom_run_threads(opts->plan7_ts, opts->nof_threads,
                               fwd_fr1, fwd_fr2, fwd_fr3,
                               rev_fr1, rev_fr2, rev_fr3,
                               results, opts);

  /* continue when all threads are done */
  for (i = 0; i < shared->nof_threads; i++)
  {
    if (pthread_join(shared->thread[i],NULL) != 0)
    {
      fprintf(stderr, "Could not join threads!");
      exit(EXIT_FAILURE);
    }
  }

  gt_pdom_free_shared(shared);

  SqdClean();
  gt_free(fwd_fr1); gt_free(fwd_fr2); gt_free(fwd_fr3);
  gt_free(rev_fr1); gt_free(rev_fr2); gt_free(rev_fr3);

  return results;
}

void gt_pdom_clear_hmms(GtArray *hmms)
{
  unsigned long i;
  if (!hmms) return;
  for (i=0;i<gt_array_size(hmms);i++)
  {
    FreePlan7(*(struct plan7_s**) gt_array_get(hmms,i));
  }
  gt_array_delete(hmms);
}

#endif
