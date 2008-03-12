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
#include "libgtcore/ma.h"
#include "libgtcore/translate.h"
#include "libgtcore/unused.h"
#include "libgtext/reverse.h"
#include "libgtltr/pdom.h"

/* HMMER related includes */
#include "globals.h"
#include "squid.h"
#include "funcs.h"

void hmmer_search(struct plan7_s *hmm,
                  char *seq,
                  unsigned long seqlen,
                  char *seqdesc,
                  struct threshold_s *thresh,
                  struct tophit_s *ghit,
                  struct tophit_s *dhit)
{
  struct dpmatrix_s *mx;
  struct p7trace_s *tr;
  unsigned char   *dsq;
  float  sc;
  double pvalue, evalue;

  mx = CreatePlan7Matrix(1, hmm->M, 25, 0);

  dsq = DigitizeSequence(seq, seqlen);

  if (P7ViterbiSpaceOK(seqlen, hmm->M, mx))
    sc = P7Viterbi(dsq, seqlen, hmm, mx, &tr);
  else
    sc = P7SmallViterbi(dsq, seqlen, hmm, mx, &tr);

  sc -= TraceScoreCorrection(hmm, tr, dsq);

  pvalue = PValue(hmm, sc);
  evalue = thresh->Z ? (double) thresh->Z * pvalue : (double) pvalue;

  if (sc >= thresh->globT && evalue <= thresh->globE)
  {
    sc = PostprocessSignificantHit(ghit, dhit,
                                   tr, hmm, dsq, seqlen,
                                   seqdesc, NULL, NULL,
                                   false, sc, true, thresh, FALSE);
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

int pdom_domain_report_hits(void *key, void *value, UNUSED void *data,
                            UNUSED Error *err)
{
  struct plan7_s *model = (struct plan7_s *) key;
  PdomHit *hit = (PdomHit*) value;
  Range rng;
  int frame = atoi(hit->best_hit->name);
  rng.start = hit->best_hit->sqfrom;
  rng.end = hit->best_hit->sqto;
  pdom_convert_frame_position(&rng, frame);
  printf("    Pdom: \t%s \t%g \t(%c, %lu, %lu)\n", model->name,
                                               hit->best_hit->pvalue,
                                               hit->best_hit->name[1],
                                               rng.start,
                                               rng.end);
  return 0;
}

int load_hmm_files(StrArray *files, Array *models, Error *err)
{
  unsigned long i;
  int had_err = 0;
  for(i=0;i<strarray_size(files);i++)
  {
    HMMFILE *hmmfp;
    struct plan7_s *hmm;
    char *hmmfile = (char*)strarray_get(files, i);
    if ((hmmfp = HMMFileOpen(hmmfile, "HMMERDB")) == NULL)
    {
      error_set(err, "failed to open HMM file '%s'", hmmfile);
      had_err = -1;
      break;
    }
    if (!had_err && !HMMFileRead(hmmfp, &hmm))
    {
      error_set(err, "failed to read any HMMs from file '%s'", hmmfile);
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
      array_add(models, hmm);
    }
    /*if (!SetAutocuts(thresh, hmm))
      Die("HMM %s did not contain the GA, TC, or NC cutoffs you needed",
          hmm->name);*/
    if(hmmfp) HMMFileClose(hmmfp);
  }
  return had_err;
}

void pdom_find(const char *seq, const char *rev_seq, LTRElement *element,
               PdomResults *results, PdomOptions *opts)
{
  struct threshold_s thresh;
  char *fwd_fr1, *fwd_fr2, *fwd_fr3,
       *rev_fr1, *rev_fr2, *rev_fr3;
  unsigned long i,
                seqlen = ltrelement_length(element);

  thresh.globT   = -FLT_MAX;
  thresh.domT    = -FLT_MAX;
  thresh.domE    = FLT_MAX;
  thresh.autocut = CUT_NONE;
  thresh.Z       = 1;
  thresh.globE   = opts->evalue_cutoff;

  results->empty = TRUE;

  /* create translations */
  translate_all_frames(&fwd_fr1,&fwd_fr2,&fwd_fr3,    seq,seqlen);
  translate_all_frames(&rev_fr1,&rev_fr2,&rev_fr3,rev_seq,seqlen);

  /* search for all pHMMs in the current sequence.
   * TODO: we may save time by only looking at the inner region w/o LTRs */
  for(i=0;i<array_size(opts->plan7_ts);i++)
  {
    struct plan7_s *hmm = *(struct plan7_s**) array_get(opts->plan7_ts,i);
    struct tophit_s *ghit;
    PdomHit *hit = ma_malloc(sizeof(PdomHit));
    ghit = AllocTophits(20);
    hit->hits_fwd = AllocTophits(20);
    hit->hits_rev = AllocTophits(20);
    hit->best_hit = NULL;

    hmmer_search(hmm,fwd_fr1,strlen(fwd_fr1),"0+",&thresh,ghit,hit->hits_fwd);
    hmmer_search(hmm,fwd_fr2,strlen(fwd_fr2),"1+",&thresh,ghit,hit->hits_fwd);
    hmmer_search(hmm,fwd_fr3,strlen(fwd_fr3),"2+",&thresh,ghit,hit->hits_fwd);
    hmmer_search(hmm,rev_fr1,strlen(rev_fr1),"0-",&thresh,ghit,hit->hits_rev);
    hmmer_search(hmm,rev_fr2,strlen(rev_fr2),"1-",&thresh,ghit,hit->hits_rev);
    hmmer_search(hmm,rev_fr3,strlen(rev_fr3),"2-",&thresh,ghit,hit->hits_rev);

    FullSortTophits(hit->hits_fwd);
    FullSortTophits(hit->hits_rev);
    if(hit->hits_fwd->num > 0 || hit->hits_rev->num > 0)
    {
      results->empty = FALSE;
      if (hit->hits_fwd->num > 0)
      {
        if (hit->hits_rev->num > 0)
        {
          if (hit->hits_fwd->hit[0]->score > hit->hits_rev->hit[0]->score)
            hit->best_hit = hit->hits_fwd->hit[0];
          else
            hit->best_hit = hit->hits_rev->hit[0];
        }
        else
          hit->best_hit = hit->hits_fwd->hit[0];
      }
      else
        hit->best_hit = hit->hits_rev->hit[0];
      hashtable_add(results->domains, hmm, hit);
    }
    else pdom_clear_domain_hit(hit);
    FreeTophits(ghit);
  }
  SqdClean();
  ma_free(fwd_fr1);ma_free(fwd_fr2);ma_free(fwd_fr3);
  ma_free(rev_fr1);ma_free(rev_fr2);ma_free(rev_fr3);
}

void pdom_clear_domain_hit(void *value)
{
  if (!value) return;
  PdomHit *hit = (PdomHit*) value;
  FreeTophits(hit->hits_fwd);
  FreeTophits(hit->hits_rev);
  ma_free(hit);
}

void pdom_clear_hmms(Array *hmms)
{
  unsigned long i;
  for(i=0;i<array_size(hmms);i++)
  {
    FreePlan7(*(struct plan7_s**) array_get(hmms,i));
  }
  array_delete(hmms);
}
