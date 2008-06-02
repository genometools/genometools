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
#include "libgtcore/ensure.h"
#include "libgtcore/ma.h"
#include "libgtcore/mathsupport.h"
#include "libgtcore/minmax.h"
#include "libgtcore/xansi.h"
#include "libgtcore/array.h"
#include "libgtext/reverse.h"
#include "libgtltr/ppt.h"

/* This enumeration defines the states in the PPT detection HMM. */
typedef enum {
  PPT_IN,
  PPT_OUT,
  PPT_UBOX,
  PPT_N,
  PPT_NOF_STATES
} PPT_States;

struct PPTHit {
  Range rng;
  double score;
  PPT_States state;
  PPTHit *ubox;
  Strand strand;
  PPTResults *res;
};

struct PPTResults {
  Array *hits;
  PPTHit *best_hit;
  LTRElement *elem;
  PPTOptions *opts;
};

static PPTResults* ppt_results_new()
{
  PPTResults *res = ma_calloc(1, sizeof (PPTResults));
  assert(res);
  return res;
}

static PPTHit* ppt_hit_new(Strand strand, PPTResults *r)
{
  PPTHit *h = ma_calloc(1, sizeof (PPTHit));
  assert(h);
  h->strand = strand;
  h->res = r;
  h->score = 0.0;
  return h;
}

Range ppt_hit_get_coords(PPTHit *h)
{
  Range r;
  assert(h);
  r.start = h->rng.start;
  r.end = h->rng.end;
  ltrelement_offset2pos(h->res->elem,
                        &r,
                        h->res->opts->radius,
                        OFFSET_BEGIN_RIGHT_LTR,
                        h->strand);
  return r;
}

PPTHit* ppt_hit_get_ubox(PPTHit *h)
{
  assert(h);
  return (h->ubox ? h->ubox : NULL);
}

Strand ppt_hit_get_strand(PPTHit *h)
{
  assert(h);
  return h->strand;
}

unsigned long ppt_results_get_number_of_hits(PPTResults *r)
{
  assert(r);
  return array_size(r->hits);
}

PPTHit* ppt_results_get_ranked_hit(PPTResults *r, unsigned long i)
{
  assert(r);
  return *(PPTHit**) array_get(r->hits, i);
}

static HMM* ppt_hmm_new(const Alpha *alpha)
{
  HMM *hmm;

  assert(alpha);

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

static double ppt_score(unsigned long radius, unsigned long end)
{
  double ret;
  unsigned long r2;
  r2 = radius * radius;
  ret = ((double) r2 - pow(abs((double) radius - (double) end),2))/(double) r2;
  return ret;
}

static int ppt_hit_cmp(const void *h1, const void *h2)
{
  PPTHit* ph1 = *(PPTHit**) h1;
  PPTHit* ph2 = *(PPTHit**) h2;
  assert(h1 && h2);
  return double_compare(ph2->score, ph1->score);
}

static bool ubox_ok(PPTHit *hit, Range uboxlen)
{
  assert(hit);
  return (hit->state == PPT_UBOX
           && hit->rng.end-hit->rng.start+1 >= uboxlen.start
           && hit->rng.end-hit->rng.start+1 <= uboxlen.end);
}

static bool ppt_ok(PPTHit *hit, Range pptlen)
{
  assert(hit);
  return (hit->state == PPT_IN
           && hit->rng.end-hit->rng.start+1 >= pptlen.start
           && hit->rng.end-hit->rng.start+1 <= pptlen.end);
}

static void group_hits(unsigned int *decoded, PPTResults *results,
                       unsigned long radius, Strand strand)
{
  PPTHit *cur_hit = NULL,
         *tmp = NULL;
  unsigned long i = 0,
                ltrlen = 0;

  assert(decoded && results && strand != STRAND_UNKNOWN);

  /* group hits into stretches */
  cur_hit = ppt_hit_new(strand, results);
  for (i=0;i<2*radius-1;i++)
  {
    cur_hit->state = decoded[i];
    cur_hit->rng.end=i;
    if (decoded[i+1] != decoded[i] || i+2==2*radius)
    {
      switch (cur_hit->state)
      {
        case PPT_UBOX:
          if (ubox_ok(cur_hit, results->opts->ubox_len))
            tmp = cur_hit;
          else
          {
            ma_free(tmp);
            tmp = NULL;
            ma_free(cur_hit);
            cur_hit = NULL;
          }
          break;
        case PPT_IN:
          if (ppt_ok(cur_hit, results->opts->ppt_len))
          {
            ltrlen = (cur_hit->strand == STRAND_FORWARD ?
                             ltrelement_rightltrlen(results->elem) :
                             ltrelement_leftltrlen(results->elem));
            cur_hit->score = ppt_score(radius, cur_hit->rng.end);
            array_add(results->hits, cur_hit);
            if (tmp)
            {
              /* this PPT has a U-box, handle accordingly */
              cur_hit->ubox = tmp;
              tmp = NULL;
            }
          }
          else
          {
            ma_free(tmp);
            tmp = NULL;
            ma_free(cur_hit);
            cur_hit = NULL;
          }
          break;
        default:
          ma_free(tmp);
          tmp = NULL;
          ma_free(cur_hit);
          cur_hit = NULL;
          break;
      }
      if (i+2!=2*radius)
      {
        cur_hit = ppt_hit_new(strand,results);
        cur_hit->rng.start = i+1;
      }
    }
  }
  if (cur_hit)
    cur_hit->rng.end++;
  ma_free(tmp);
  tmp = NULL;
}

PPTResults* ppt_find(const char *seq,
                     const char *rev_seq,
                     LTRElement *element,
                     PPTOptions *o)
{
  unsigned int *encoded_seq=NULL,
               *decoded=NULL;
  const Alpha *alpha;
  HMM *hmm;
  PPTResults *results = NULL;
  unsigned long i = 0,
                radius = 0,
                seqlen = ltrelement_length(element),
                ltrlen = 0;

  assert(seq && rev_seq && element && o);

  results = ppt_results_new();
  results->elem = element;
  results->opts = o;
  results->hits = array_new(sizeof (PPTHit*));

  alpha = alpha_new_dna();
  hmm = ppt_hmm_new(alpha);

  results->best_hit = NULL;

  /* do PPT finding on forward strand
   * -------------------------------- */
  ltrlen = ltrelement_rightltrlen(element);
  /* make sure that we do not cross the LTR boundary */
  radius = MIN(o->radius, ltrlen-1);
  /* encode sequence */
  encoded_seq = ma_malloc(sizeof (unsigned int) * seqlen);
  for (i=0;i<seqlen;i++)
  {
    encoded_seq[i] = alpha_encode(alpha, seq[i]);
  }
  /* use Viterbi algorithm to decode emissions within radius */
  decoded = ma_malloc(sizeof (unsigned int) * 2*radius+1);
  hmm_decode(hmm, decoded, encoded_seq+seqlen-ltrlen-radius+1, 2*radius);
  group_hits(decoded, results, radius, STRAND_FORWARD);
  /* radius length may change in the next strand, so reallocate */
  ma_free(decoded);

  /* do PPT finding on reverse strand
   * -------------------------------- */

  ltrlen = ltrelement_leftltrlen(element);
  /* make sure that we do not cross the LTR boundary */
  radius = MIN(o->radius, ltrlen-1);
  /* encode sequence */
  for (i=0;i<seqlen;i++)
  {
    encoded_seq[i] = alpha_encode(alpha, rev_seq[i]);
  }
  /* use Viterbi algorithm to decode emissions within radius */
  decoded = ma_malloc(sizeof (unsigned int) * 2*radius+1);
  hmm_decode(hmm, decoded, encoded_seq+seqlen-ltrlen-radius, 2*radius);
  group_hits(decoded, results, radius, STRAND_REVERSE);

  /* rank hits by descending score */
  array_sort(results->hits, ppt_hit_cmp);

  ma_free(encoded_seq);
  ma_free(decoded);
  alpha_delete((Alpha*) alpha);
  hmm_delete(hmm);

  return results;
}

void ppt_results_delete(PPTResults *results)
{
  unsigned long i;

  if (!results) return;

  if (results->hits)
  {
    for (i=0;i<array_size(results->hits);i++)
    {
      PPTHit * hit = *(PPTHit**) array_get(results->hits,i);
      if (hit->ubox)
        ma_free(hit->ubox);
      ma_free(hit);
    }
    array_delete(results->hits);
  }
  ma_free(results);
}

int ppt_unit_test(Error *err)
{
  int had_err = 0;
  PPTOptions o;
  PPTResults *rs;
  PPTHit *h;
  LTRElement element;
  Range rng;
  const char *seq = "gatcagtcgactcgatcgactcgatcgactcgagcacggcgacgat"
                    "gctggtcggctaactggggggggaggatcgacttcgactcgacgatcgactcga";
  char *rev_seq = ma_malloc(strlen(seq)*sizeof (char*));

  memcpy(rev_seq, seq, sizeof (char) * strlen(seq));
  reverse_complement(rev_seq,strlen(seq),err);
  element.leftLTR_3 = 21;
  element.leftLTR_5 = 0;
  element.rightLTR_3 = 99;
  element.rightLTR_5 = 72;
  memset(&o, 0, sizeof (PPTOptions));
  o.ppt_len.start = 5;
  o.ppt_len.end = 15;
  o.ubox_len.start = 2;
  o.ubox_len.end = 15;
  o.radius = 12;
  rs = ppt_find(seq, rev_seq, &element, &o);

  ensure(had_err, ppt_results_get_number_of_hits(rs));
  h = ppt_results_get_ranked_hit(rs, 0);
  ensure(had_err, h);
  rng = ppt_hit_get_coords(h);
  ensure(had_err, rng.start == 61);
  ensure(had_err, rng.end == 72);
  ma_free(rev_seq);
  return had_err;
}
