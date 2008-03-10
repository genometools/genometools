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

#include <assert.h>
#include <math.h>
#include "libgtcore/ma.h"
#include "libgtcore/minmax.h"
#include "libgtcore/xansi.h"
#include "libgtcore/array.h"
#include "libgtltr/ppt.h"

HMM* ppt_hmm_new(const Alpha *alpha)
{
  HMM *hmm;

  hmm = hmm_new(PPT_NOF_STATES, alpha_size(alpha));

  /* set emission probabilities */
  hmm_set_emission_probability(hmm, PPT_OUT,  alpha_encode(alpha, 'G'), 0.25);
  hmm_set_emission_probability(hmm, PPT_OUT,  alpha_encode(alpha, 'A'), 0.25);
  hmm_set_emission_probability(hmm, PPT_OUT,  alpha_encode(alpha, 'C'), 0.25);
  hmm_set_emission_probability(hmm, PPT_OUT,  alpha_encode(alpha, 'T'), 0.25);
  hmm_set_emission_probability(hmm, PPT_IN,  alpha_encode(alpha, 'G'),  0.485);
  hmm_set_emission_probability(hmm, PPT_IN,  alpha_encode(alpha, 'A'),  0.485);
  hmm_set_emission_probability(hmm, PPT_IN,  alpha_encode(alpha, 'C'),  0.015);
  hmm_set_emission_probability(hmm, PPT_IN,  alpha_encode(alpha, 'T'),  0.015);
  hmm_set_emission_probability(hmm, PPT_UBOX, alpha_encode(alpha, 'G'), 0.03);
  hmm_set_emission_probability(hmm, PPT_UBOX, alpha_encode(alpha, 'A'), 0.03);
  hmm_set_emission_probability(hmm, PPT_UBOX, alpha_encode(alpha, 'C'), 0.03);
  hmm_set_emission_probability(hmm, PPT_UBOX, alpha_encode(alpha, 'T'), 0.91);
  hmm_set_emission_probability(hmm, PPT_N,    alpha_encode(alpha, 'G'), 0.00);
  hmm_set_emission_probability(hmm, PPT_N,    alpha_encode(alpha, 'A'), 0.00);
  hmm_set_emission_probability(hmm, PPT_N,    alpha_encode(alpha, 'C'), 0.00);
  hmm_set_emission_probability(hmm, PPT_N,    alpha_encode(alpha, 'T'), 0.00);
  hmm_set_emission_probability(hmm, PPT_N,    alpha_encode(alpha, 'N'), 1.00);

  /* set transition probabilities */
  hmm_set_transition_probability(hmm, PPT_OUT, PPT_IN,   0.05);
  hmm_set_transition_probability(hmm, PPT_OUT, PPT_N,    0.05);
  hmm_set_transition_probability(hmm, PPT_OUT, PPT_UBOX, 0.05);
  hmm_set_transition_probability(hmm, PPT_UBOX, PPT_OUT, 0.05);
  hmm_set_transition_probability(hmm, PPT_UBOX, PPT_N,   0.05);
  hmm_set_transition_probability(hmm, PPT_UBOX, PPT_IN,  0.05);
  hmm_set_transition_probability(hmm, PPT_IN, PPT_UBOX,  0.05);
  hmm_set_transition_probability(hmm, PPT_IN, PPT_OUT,   0.05);
  hmm_set_transition_probability(hmm, PPT_IN, PPT_N,     0.05);
  hmm_set_transition_probability(hmm, PPT_N, PPT_UBOX,   0.05);
  hmm_set_transition_probability(hmm, PPT_N, PPT_OUT,    0.05);
  hmm_set_transition_probability(hmm, PPT_N, PPT_IN,     0.05);
  hmm_set_missing_transition_probabilities(hmm);
  assert(hmm_is_valid(hmm));

  return hmm;
}

double ppt_score(unsigned long posdiff, unsigned int width)
{
  double score, optscore;
  optscore = pow(width,4);
  score = (-pow(posdiff,4))+optscore;
  return score/optscore;
}

static unsigned long score_hits(Array *results, unsigned long seqlen,
                                unsigned long ltrlen, unsigned int radius)
{
  unsigned long i = 0,
                highest_index = UNDEF_ULONG;
  double highest_score = 0.0;

  /* score and report PPTs */
  for (i=0;i<array_size(results);i++)
  {
    PPT_Hit *hit = *(PPT_Hit**) array_get(results,i);
    if (hit->state == PPT_IN)
    {
      hit->score = ppt_score(abs((seqlen-ltrlen-1)
                                 -(seqlen-ltrlen-1)-radius+hit->end),
                            radius);
      if (i>0)
      {
        PPT_Hit *uhit = *(PPT_Hit**) array_get(results,i-1);
        if (uhit->state == PPT_UBOX)
          /* this PPT has a U-box, handle accordingly */
          hit->ubox = uhit;
      }
      if (hit->score > highest_score)
      {
        highest_index = i;
        highest_score = hit->score;
      }
    }
  }
  return highest_index;
}

static void group_hits(unsigned int *decoded, Array *results, PPTOptions *o,
                       unsigned long radius, Strand strand)
{
  PPT_Hit *cur_hit;
  unsigned long i = 0;

  /* group hits into stretches */
  cur_hit = ma_malloc(sizeof (PPT_Hit));
  cur_hit->start = 0UL;
  cur_hit->score = 0.0;
  cur_hit->strand = strand;
  cur_hit->ubox = NULL;
  for (i=0;i<2*radius-1;i++)
  {
    cur_hit->state = decoded[i];
    cur_hit->end=i;
    if (decoded[i+1] != decoded[i] || i+2==2*radius)
    {
      if ((cur_hit->state == PPT_IN
           && cur_hit->end-cur_hit->start+1 >= o->ppt_minlen)
          || (cur_hit->state == PPT_UBOX
           && cur_hit->end-cur_hit->start+1 >= o->ubox_minlen))
        array_add(results, cur_hit);
      else {
        cur_hit->state = PPT_OUT;
        array_add(results, cur_hit);
      }
      if (i+2!=2*radius)
      {
        cur_hit = ma_malloc(sizeof (PPT_Hit));
        cur_hit->start = i+1;
        cur_hit->score = 0.0;
        cur_hit->strand = strand;
        cur_hit->ubox = NULL;
      }
    }
  }
  cur_hit->end++;
}

void ppt_find(const char *seq,
              const char *rev_seq,
              LTRElement *element,
              PPTResults *results,
              PPTOptions *o)
{
  assert(seq && element && results && o);
  unsigned int *encoded_seq=NULL,
               *decoded=NULL;
  const Alpha *alpha = alpha_new_dna();
  HMM *hmm = ppt_hmm_new(alpha);
  Array *results_fwd = array_new(sizeof(PPT_Hit*)),
        *results_rev = array_new(sizeof(PPT_Hit*));
  unsigned long i = 0,
                radius = 0,
                seqlen = ltrelement_length(element),
                highest_fwd,
                highest_rev,
                ltrlen;
  PPT_Hit *tmp;

  results->best_hit = NULL;

  /* do PPT finding on forward strand
   * -------------------------------- */
  ltrlen = ltrelement_rightltrlen(element);
  /* make sure that we do not cross the LTR boundary */
  radius = MIN(o->radius, ltrlen);
  /* encode sequence */
  encoded_seq = ma_malloc(sizeof (unsigned int) * seqlen);
  for (i=0;i<seqlen;i++)
  {
    encoded_seq[i] = alpha_encode(alpha, seq[i]);
  }
  /* use Viterbi algorithm to decode emissions within radius */
  decoded = ma_malloc(sizeof (unsigned int) * 2*radius+2);
  hmm_decode(hmm, decoded, encoded_seq+seqlen-ltrlen-radius, 2*radius);
  group_hits(decoded, results_fwd, o, radius, STRAND_FORWARD);
  highest_fwd = score_hits(results_fwd, seqlen, ltrlen, radius);
  results->hits_fwd = results_fwd;
  if (highest_fwd != UNDEF_ULONG)
  {
    tmp = *(PPT_Hit**) array_get(results_fwd, highest_fwd);
    results->best_hit = tmp;
  }
  /* do PPT finding on reverse strand
   * -------------------------------- */
  ltrlen = ltrelement_leftltrlen(element);
  /* make sure that we do not cross the LTR boundary */
  radius = MIN(o->radius, ltrlen);
  /* encode sequence */
  for (i=0;i<seqlen;i++)
  {
    encoded_seq[i] = alpha_encode(alpha, rev_seq[i]);
  }
  /* use Viterbi algorithm to decode emissions within radius */
  hmm_decode(hmm, decoded, encoded_seq+seqlen-ltrlen-radius, 2*radius);
  group_hits(decoded, results_rev, o, radius, STRAND_REVERSE);
  highest_rev = score_hits(results_rev, seqlen, ltrlen, radius);
  results->hits_rev = results_rev;
  if (highest_rev != UNDEF_ULONG)
  {
    tmp = *(PPT_Hit**) array_get(results_rev, highest_rev);
    if (!results->best_hit || tmp->score > results->best_hit->score)
      results->best_hit = tmp;
  }

  ma_free(encoded_seq);
  ma_free(decoded);
  alpha_delete((Alpha*) alpha);
  hmm_delete(hmm);
}

void ppt_clear_results(PPTResults *results)
{
  unsigned long i;
  for (i=0;i<array_size(results->hits_fwd);i++)
  {
    ma_free(*(PPT_Hit**) array_get(results->hits_fwd,i));
  }
  array_delete(results->hits_fwd);
  for (i=0;i<array_size(results->hits_rev);i++)
  {
    ma_free(*(PPT_Hit**) array_get(results->hits_rev,i));
  }
  array_delete(results->hits_rev);
}
