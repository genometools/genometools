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
  hmm_set_emission_probability(hmm, PPT_OUT, alpha_encode(alpha, 'g'), 0.25);
  hmm_set_emission_probability(hmm, PPT_OUT, alpha_encode(alpha, 'a'), 0.25);
  hmm_set_emission_probability(hmm, PPT_OUT, alpha_encode(alpha, 'c'), 0.25);
  hmm_set_emission_probability(hmm, PPT_OUT, alpha_encode(alpha, 't'), 0.25);
  hmm_set_emission_probability(hmm, PPT_IN, alpha_encode(alpha, 'g'), 0.485);
  hmm_set_emission_probability(hmm, PPT_IN, alpha_encode(alpha, 'a'), 0.485);
  hmm_set_emission_probability(hmm, PPT_IN, alpha_encode(alpha, 'c'), 0.015);
  hmm_set_emission_probability(hmm, PPT_IN, alpha_encode(alpha, 't'), 0.015);
  hmm_set_emission_probability(hmm, PPT_UBOX, alpha_encode(alpha, 'g'), 0.02);
  hmm_set_emission_probability(hmm, PPT_UBOX, alpha_encode(alpha, 'a'), 0.02);
  hmm_set_emission_probability(hmm, PPT_UBOX, alpha_encode(alpha, 'c'), 0.02);
  hmm_set_emission_probability(hmm, PPT_UBOX, alpha_encode(alpha, 't'), 0.94);

  /* set transition probabilities */
  hmm_set_transition_probability(hmm, PPT_IN, PPT_OUT, 0.05);
  hmm_set_transition_probability(hmm, PPT_OUT, PPT_IN, 0.05);
  hmm_set_transition_probability(hmm, PPT_OUT, PPT_UBOX, 0.05);
  hmm_set_transition_probability(hmm, PPT_UBOX, PPT_OUT, 0.05);
  hmm_set_transition_probability(hmm, PPT_UBOX, PPT_IN, 0.05);
  hmm_set_transition_probability(hmm, PPT_IN, PPT_UBOX, 0.05);
  hmm_set_missing_transition_probabilities(hmm);
  assert(hmm_is_valid(hmm));

  return hmm;
}

/* simple quadratic score, please refine me! */
double ppt_score(unsigned long posdiff, unsigned int width)
{
  double score, optscore;
  optscore = 0.5*pow(width,4);
  score = (0.5*-pow(posdiff,4))+optscore;
  return score/optscore;
}

unsigned long ppt_find(const char *seq,
                       unsigned long seqlen,
                       unsigned long ltrlen,
                       Array *results,
                       LTRharvestoptions *lo,
                       double score_func(unsigned long, unsigned int),
                       Strand strand)
{
  assert(seq && seqlen>0 && ltrlen>0 && results && score_func && lo
          && strand != STRAND_UNKNOWN && strand != STRAND_BOTH);
  unsigned int *encoded_seq=NULL,
               *decoded=NULL;
  const Alpha *alpha = alpha_new_dna();
  HMM *hmm = ppt_hmm_new(alpha);
  PPT_Hit *cur_hit;
  unsigned long i = 0,
                radius = 0,
                highest_index = UNDEF_ULONG;
  double highest_score = 0.0;

  /* encode sequence */
  encoded_seq = ma_malloc(sizeof (unsigned int) * seqlen);
  for (i=0;i<seqlen;i++)
  {
    encoded_seq[i] = alpha_encode(alpha, seq[i]);
  }

  /* make sure that we do not cross the LTR boundary */
  radius = MIN(lo->ppt_radius, ltrlen);

  /* use Viterbi algorithm to decode emissions within radius */
  decoded = ma_malloc(sizeof (unsigned int) * 2*radius+2);
  hmm_decode(hmm, decoded, encoded_seq+seqlen-ltrlen-radius, 2*radius);

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
           && cur_hit->end-cur_hit->start+1 >= lo->ppt_minlen)
          || (cur_hit->state == PPT_UBOX
           && cur_hit->end-cur_hit->start+1 >= lo->ubox_minlen))
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

  /* score and report PPTs */
  for (i=0;i<array_size(results);i++)
  {
    unsigned long ltrlen;
    PPT_Hit *hit = *(PPT_Hit**) array_get(results,i);
    if (hit->state == PPT_IN)
    {
      hit->score = score_func(abs((seqlen-ltrlen-1)
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
        highest_index = i;
    }
  }

  ma_free(encoded_seq);
  ma_free(decoded);
  alpha_delete((Alpha*) alpha);
  hmm_delete(hmm);
  return highest_index;
}
