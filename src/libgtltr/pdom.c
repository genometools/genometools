/*
  Copyright (c) 2008 Sascha Steinbiss <ssteinbiss@zbh.uni-hamburg.de>
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

#include <string.h>
#include <ctype.h>
#include <float.h>
#include <pthread.h>
#include "libgtcore/codon.h"
#include "libgtcore/log.h"
#include "libgtcore/ma.h"
#include "libgtcore/translate.h"
#include "libgtcore/unused.h"
#include "libgtext/globalchaining.h"
#include "libgtext/reverse.h"
#include "libgtltr/pdom.h"

/* HMMER related includes */
#include "globals.h"
#include "squid.h"
#include "funcs.h"

typedef struct pdom_shared {
  pthread_t *thread;
  int nof_threads;
  unsigned long next_hmm;
  Array *hmms;
  char *fwd_fr1, *fwd_fr2, *fwd_fr3,
       *rev_fr1, *rev_fr2, *rev_fr3;
  PdomResults *results;
  PdomOptions *opts;
  pthread_mutex_t in_lock, out_lock;
} pdom_shared_s;

void hmmer_search(struct plan7_s *hmm,
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
  int rtn;
  double pvalue, evalue;

  mx = CreatePlan7Matrix(1, hmm->M, 25, 0);

  dsq = DigitizeSequence(seq, seqlen);

  if (P7ViterbiSpaceOK(seqlen, hmm->M, mx))
    sc = P7Viterbi(dsq, seqlen, hmm, mx, &tr);
  else
    sc = P7SmallViterbi(dsq, seqlen, hmm, mx, &tr);

  sc -= TraceScoreCorrection(hmm, tr, dsq);

  if ((rtn = pthread_mutex_lock(lock)) != 0)
  {
    fprintf(stderr, "pthread_mutex_lock failure: %s\n", strerror(rtn));
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
  if ((rtn = pthread_mutex_unlock(lock)) != 0)
  {
    fprintf(stderr, "pthread_mutex_unlock failure: %s\n", strerror(rtn));
    exit(EXIT_FAILURE);
  }
  P7FreeTrace(tr);
  free(dsq);
  FreePlan7Matrix(mx);
}

void pdom_convert_frame_position(Range *rng, int frame)
{
  rng->start = (rng->start - 1)*CODONLENGTH + frame;
  rng->end   = (rng->end   - 1)*CODONLENGTH + frame;
}

int pdom_load_hmm_files(PdomOptions *opts, Error *err)
{
  unsigned long i;
  int had_err = 0;

  assert(opts && err);

  for (i=0;i<strarray_size(opts->hmm_files);i++)
  {
    HMMFILE *hmmfp;
    struct plan7_s *hmm;
    char *hmmfile = (char*)strarray_get(opts->hmm_files, i);
    if ((hmmfp = HMMFileOpen(hmmfile, "HMMERDB")) == NULL)
    {
      error_set(err, "Failed to open HMM file '%s'", hmmfile);
      had_err = -1;
      break;
    }
    if (!had_err && !HMMFileRead(hmmfp, &hmm))
    {
      error_set(err, "Failed to read any HMMs from file '%s'", hmmfile);
      had_err = -1;
      break;
    }
    if (!had_err && hmm == NULL)
    {
      had_err = -1;
      error_set(err, "HMM file '%s' corrupt or in incorrect format?", hmmfile);
      break;
    }
    if (!had_err)
    {
      P7Logoddsify(hmm, true);
      array_add(opts->plan7_ts, hmm);
    }
    if (!SetAutocuts(&opts->thresh, hmm))
    {
      fprintf(stderr,"HMM %s did not contain the GA, TC, or NC "
                     "cutoffs you needed",
                     hmm->name);
      exit(EXIT_FAILURE);
    }
    if (hmmfp) HMMFileClose(hmmfp);
  }
  log_log("Loaded %lu HMM model(s)",
          array_size(opts->plan7_ts));
  return had_err;
}

static void chainproc(Chain *c, Fragment *f, void *data)
{
  unsigned long i;
  PdomHit *hit;
  hit = (PdomHit *) data;

  log_log("resulting chain has %ld fragments, score %ld", chain_size(c),
                                                          chain_get_score(c));
  for (i=0;i<chain_size(c);i++)
  {
    Fragment frag = f[chain_get_fragnum(c, i)];
    log_log("(%lu %lu) (%lu %lu), %p", frag.startpos1, frag.endpos1,
                                       frag.startpos2, frag.endpos2, frag.data);
    array_add(hit->best_chain, frag.data);
  }
  log_log("\n");
}

static int fragcmp(const void *frag1, const void *frag2)
{
  Fragment *f1 = (Fragment*) frag1;
  Fragment *f2 = (Fragment*) frag2;
  if (f1->startpos2 == f2->startpos2)
    return 0;
  else return (f1->startpos2 < f2->startpos2 ? -1 : 1);
}

void* pdom_per_domain_worker_thread(void *data)
{
  pdom_shared_s *shared;
  struct plan7_s *hmm;
  int rtn;
  struct tophit_s *ghit, *hits;
  bool best_fwd = TRUE;
  unsigned long i;
  Fragment *frags;
  PdomHit *hit;

  shared = (pdom_shared_s *) data;

  for (;;)
  {
    /* Lock input */
    if ((rtn = pthread_mutex_lock(&shared->in_lock)) != 0)
    {
      fprintf(stderr, "Failed to lock: %s\n", strerror(rtn));
      exit(EXIT_FAILURE);
    }
    /* Have all HMMs been distributed? If so, we are done here. */
    if (shared->next_hmm == array_size(shared->opts->plan7_ts))
    {
      if ((rtn = pthread_mutex_unlock(&shared->in_lock)) != 0)
      {
        fprintf(stderr, "Failed to unlock: %s\n", strerror(rtn));
        exit(EXIT_FAILURE);
      }
      pthread_exit(NULL);
    }

    /* Get work from HMM list */
    hmm = *(struct plan7_s**) array_get(shared->opts->plan7_ts,
                                        shared->next_hmm);
    shared->next_hmm++;

    /* work claimed, release input */
    if ((rtn = pthread_mutex_unlock(&shared->in_lock)) != 0)
    {
      fprintf(stderr, "pthread_mutex_unlock failure: %s\n", strerror(rtn));
      exit(EXIT_FAILURE);
    }

    hit = ma_malloc(sizeof (PdomHit));
    ghit = AllocTophits(20);
    hit->hits_fwd = AllocTophits(20);
    hit->hits_rev = AllocTophits(20);
    hit->best_chain = array_new(sizeof(struct hit_s*));

    hmmer_search(hmm,shared->fwd_fr1,strlen(shared->fwd_fr1),"0+",
                 &shared->opts->thresh,ghit, hit->hits_fwd, &shared->out_lock);
    hmmer_search(hmm,shared->fwd_fr2,strlen(shared->fwd_fr2),"1+",
                 &shared->opts->thresh,ghit, hit->hits_fwd, &shared->out_lock);
    hmmer_search(hmm,shared->fwd_fr3,strlen(shared->fwd_fr3),"2+",
                 &shared->opts->thresh,ghit, hit->hits_fwd, &shared->out_lock);
    hmmer_search(hmm,shared->rev_fr1,strlen(shared->rev_fr1),"0-",
                 &shared->opts->thresh,ghit, hit->hits_rev, &shared->out_lock);
    hmmer_search(hmm,shared->rev_fr2,strlen(shared->rev_fr2),"1-",
                 &shared->opts->thresh,ghit, hit->hits_rev, &shared->out_lock);
    hmmer_search(hmm,shared->rev_fr3,strlen(shared->rev_fr3),"2-",
                 &shared->opts->thresh,ghit, hit->hits_rev, &shared->out_lock);

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
          if (!(hit->hits_fwd->hit[0]->score > hit->hits_rev->hit[0]->score))
            best_fwd = FALSE;
        }
      }
      else
        best_fwd = FALSE;

      /* determine best-scoring strand */
      hits = (best_fwd ? hit->hits_fwd : hit->hits_rev);

      /* no need to chain if there is only one hit */
      if (hits->num > 1)
      {
        /* create fragment set for chaining */
        frags = (Fragment*) ma_calloc(hits->num, sizeof (Fragment));
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
        qsort(frags, hits->num, sizeof (Fragment), fragcmp);
        for (i=0;i<hits->num;i++)
        {
          log_log("(%lu %lu) (%lu %lu) %p", frags[i].startpos1,
                                         frags[i].endpos1,
                                         frags[i].startpos2,
                                         frags[i].endpos2,
                                         frags[i].data);
        }
        log_log("chaining %d frags", hits->num);
        /* do chaining */
        globalchaining_max(frags, hits->num,
                           shared->opts->chain_max_gap_length,
                           chainproc, hit);
        ma_free(frags);
      }
      else
        array_add(hit->best_chain, hits->hit[0]);

      /* Lock results, we want to write to the result hashtable */
      if ((rtn = pthread_mutex_lock(&(shared->out_lock))) != 0)
      {
        fprintf(stderr, "Failed to lock: %s\n", strerror(rtn));
        exit(EXIT_FAILURE);
      }

      /* register results */
      hashtable_add(shared->results->domains, hmm, hit);
      if (best_fwd)
      {
        shared->results->combined_e_value_fwd *= hit->hits_fwd->hit[0]->pvalue;
        hit->strand = STRAND_FORWARD;
      }
      else
      {
        shared->results->combined_e_value_rev *= hit->hits_rev->hit[0]->pvalue;
        hit->strand = STRAND_REVERSE;
      }

      /* unlock results */
      if ((rtn = pthread_mutex_unlock(&(shared->out_lock))) != 0)
      {
        fprintf(stderr, "Failed to unlock: %s\n", strerror(rtn));
        exit(EXIT_FAILURE);
      }
    }
    else
      pdom_clear_domain_hit(hit);

    FreeTophits(ghit);
  }
}

static pdom_shared_s* pdom_run_threads(Array *hmms, int nof_threads,
                                  char *fwd_fr1,char *fwd_fr2,char *fwd_fr3,
                                  char *rev_fr1,char *rev_fr2,char *rev_fr3,
                                  PdomResults *results, PdomOptions *opts)
{
  int rtn, i;
  pdom_shared_s *shared;
  pthread_attr_t attr;

  assert(hmms && nof_threads > 0 && *fwd_fr1 && *fwd_fr2 && *fwd_fr3
          && *rev_fr1 && rev_fr2 && rev_fr3 && results && opts);

  shared = ma_malloc(sizeof (pdom_shared_s));

  shared->nof_threads = nof_threads;
  shared->thread = ma_malloc(sizeof (pthread_t) * nof_threads);
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

  if ((rtn = pthread_mutex_init(&shared->in_lock, NULL)) != 0)
  {
    fprintf(stderr, "Could not initialize lock! %s\n", strerror(rtn));
    exit(EXIT_FAILURE);
  }
  if ((rtn = pthread_mutex_init(&shared->out_lock, NULL)) != 0)
  {
    fprintf(stderr, "Could not initialize lock! %s\n", strerror(rtn));
    exit(EXIT_FAILURE);
  }

  pthread_attr_init(&attr);

  for (i = 0; i < nof_threads; i++)
  {
    if ((rtn = pthread_create(&(shared->thread[i]), &attr,
           pdom_per_domain_worker_thread, (void *) shared)) != 0)
    {
      fprintf(stderr,"Failed to create thread %d, return code %d\n", i, rtn);
      exit(EXIT_FAILURE);
    }
  }
  pthread_attr_destroy(&attr);
  return shared;
}

static void pdom_free_shared(pdom_shared_s *shared)
{
  ma_free(shared->thread);
  ma_free(shared);
}

void pdom_find(const char *seq, const char *rev_seq, LTRElement *element,
               PdomResults *results, PdomOptions *opts)
{
  char *fwd_fr1, *fwd_fr2, *fwd_fr3,
       *rev_fr1, *rev_fr2, *rev_fr3;
  unsigned long seqlen = ltrelement_length(element);
  pdom_shared_s *shared;
  int i;

  assert(seq && rev_seq && element && results && opts);

  results->empty = TRUE;
  results->combined_e_value_fwd = results->combined_e_value_rev = 1.0;

  /* create translations */
  translate_all_frames(&fwd_fr1,&fwd_fr2,&fwd_fr3,    seq,seqlen);
  translate_all_frames(&rev_fr1,&rev_fr2,&rev_fr3,rev_seq,seqlen);

  /* start worker threads */
  shared = pdom_run_threads(opts->plan7_ts, opts->nof_threads,
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

  pdom_free_shared(shared);

  SqdClean();
  ma_free(fwd_fr1);ma_free(fwd_fr2);ma_free(fwd_fr3);
  ma_free(rev_fr1);ma_free(rev_fr2);ma_free(rev_fr3);
}

void pdom_clear_domain_hit(void *value)
{
  PdomHit *hit;
  if (!value) return;
  hit = (PdomHit*) value;
  if (hit->hits_fwd)
    FreeTophits(hit->hits_fwd);
  if (hit->hits_rev)
    FreeTophits(hit->hits_rev);
  if (hit->best_chain)
    array_delete(hit->best_chain);
  ma_free(hit);
}

void pdom_clear_hmms(Array *hmms)
{
  unsigned long i;
  if (!hmms) return;
  for (i=0;i<array_size(hmms);i++)
  {
    FreePlan7(*(struct plan7_s**) array_get(hmms,i));
  }
  array_delete(hmms);
}
