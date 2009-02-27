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

#include <math.h>
#include "core/ensure.h"
#include "core/ma.h"
#include "core/mathsupport.h"
#include "core/minmax.h"
#include "core/xansi.h"
#include "core/array.h"
#include "extended/reverse.h"
#include "ltr/ppt.h"

/* This enumeration defines the states in the PPT detection HMM. */
typedef enum {
  PPT_IN,
  PPT_OUT,
  PPT_UBOX,
  PPT_N,
  PPT_NOF_STATES
} GtPPTStates;

struct GtPPTHit {
  GtRange rng;
  double score;
  GtPPTStates state;
  GtPPTHit *ubox;
  GtStrand strand;
  GtPPTResults *res;
};

struct GtPPTResults {
  GtArray *hits;
  GtLTRElement *elem;
  GtPPTOptions *opts;
};

static GtPPTResults* gt_ppt_results_new(GtLTRElement *elem,
                                        GtPPTOptions *opts)
{
  GtPPTResults *res = gt_calloc(1, sizeof (GtPPTResults));
  res->elem = elem;
  res->opts = opts;
  res->hits = gt_array_new(sizeof (GtPPTHit*));
  return res;
}

static GtPPTHit* gt_ppt_hit_new(GtStrand strand, GtPPTResults *r)
{
  GtPPTHit *h = gt_calloc(1, sizeof (GtPPTHit));
  gt_assert(h);
  h->strand = strand;
  h->res = r;
  h->score = 0.0;
  return h;
}

GtRange gt_ppt_hit_get_coords(const GtPPTHit *h)
{
  GtRange rng;
  gt_assert(h);
  rng.start = h->rng.start;
  rng.end = h->rng.end;
/*  printf("%lu-%lu", rng.start, rng.end); */
  switch (h->strand)
  {
    case GT_STRAND_FORWARD:
    default:
      rng.start = h->res->elem->rightLTR_5 - 1 - h->res->opts->radius + rng.start;
      rng.end = rng.start + (gt_range_length(&h->rng) - 1);
      break;
    case GT_STRAND_REVERSE:
      rng.end = h->res->elem->leftLTR_3 + 1 + h->res->opts->radius - rng.start;
      rng.start = rng.end - (gt_range_length(&h->rng) - 1);
      break;
  }
/*  printf("%lu-%lu", rng.start, rng.end);
  printf("old: %lu, new %lu\n", gt_range_length(&h->rng) , gt_range_length(&rng)); */
  gt_assert(gt_range_length(&rng) == gt_range_length(&h->rng));
  return rng;
}

GtPPTHit* gt_ppt_hit_get_ubox(const GtPPTHit *h)
{
  gt_assert(h);
  return (h->ubox ? h->ubox : NULL);
}

GtStrand gt_ppt_hit_get_strand(const GtPPTHit *h)
{
  gt_assert(h);
  return h->strand;
}

unsigned long gt_ppt_results_get_number_of_hits(GtPPTResults *r)
{
  gt_assert(r);
  return gt_array_size(r->hits);
}

GtPPTHit* gt_ppt_results_get_ranked_hit(GtPPTResults *r, unsigned long i)
{
  gt_assert(r);
  return *(GtPPTHit**) gt_array_get(r->hits, i);
}

GtHMM* gt_ppt_hmm_new(const GtAlpha *alpha, GtPPTOptions *opts)
{
  GtHMM *hmm;
  double non_u_prob = 0.0;

  gt_assert(alpha);

  hmm = gt_hmm_new(PPT_NOF_STATES, gt_alpha_size(alpha));

  /* set emission probabilities */
  gt_hmm_set_emission_probability(hmm, PPT_OUT,
                                  gt_alpha_encode(alpha, 'G'),
                                  opts->bkg_g_prob);
  gt_hmm_set_emission_probability(hmm, PPT_OUT,
                                  gt_alpha_encode(alpha, 'A'),
                                  opts->bkg_a_prob);
  gt_hmm_set_emission_probability(hmm, PPT_OUT,
                                  gt_alpha_encode(alpha, 'C'),
                                  opts->bkg_c_prob);
  gt_hmm_set_emission_probability(hmm, PPT_OUT,
                                  gt_alpha_encode(alpha, 'T'),
                                  opts->bkg_t_prob);
  gt_hmm_set_emission_probability(hmm, PPT_IN,
                                  gt_alpha_encode(alpha, 'G'),
                                  opts->ppt_purine_prob/2);
  gt_hmm_set_emission_probability(hmm, PPT_IN,
                                  gt_alpha_encode(alpha, 'A'),
                                  opts->ppt_purine_prob/2);
  gt_hmm_set_emission_probability(hmm, PPT_IN,
                                  gt_alpha_encode(alpha, 'C'),
                                  opts->ppt_pyrimidine_prob/2);
  gt_hmm_set_emission_probability(hmm, PPT_IN,
                                  gt_alpha_encode(alpha, 'T'),
                                  opts->ppt_pyrimidine_prob/2);
  gt_hmm_set_emission_probability(hmm, PPT_UBOX,
                                  gt_alpha_encode(alpha, 'T'),
                                  opts->ubox_u_prob);
  /* calculate non-U probabilities (still uniform, may be optimised) */
  non_u_prob = (1.0 - (opts->ubox_u_prob)) / ((double) 3);
  gt_hmm_set_emission_probability(hmm, PPT_UBOX,
                                  gt_alpha_encode(alpha, 'G'), non_u_prob);
  gt_hmm_set_emission_probability(hmm, PPT_UBOX,
                                  gt_alpha_encode(alpha, 'A'), non_u_prob);
  gt_hmm_set_emission_probability(hmm, PPT_UBOX,
                                  gt_alpha_encode(alpha, 'C'), non_u_prob);
  gt_hmm_set_emission_probability(hmm, PPT_N,
                                  gt_alpha_encode(alpha, 'G'), 0.00);
  gt_hmm_set_emission_probability(hmm, PPT_N,
                                  gt_alpha_encode(alpha, 'A'), 0.00);
  gt_hmm_set_emission_probability(hmm, PPT_N,
                                  gt_alpha_encode(alpha, 'C'), 0.00);
  gt_hmm_set_emission_probability(hmm, PPT_N,
                                  gt_alpha_encode(alpha, 'T'), 0.00);
  gt_hmm_set_emission_probability(hmm, PPT_N,
                                  gt_alpha_encode(alpha, 'N'), 1.00);

  /* set transition probabilities */
  gt_hmm_set_transition_probability(hmm, PPT_OUT, PPT_IN,   0.05);
  gt_hmm_set_transition_probability(hmm, PPT_OUT, PPT_N,    0.05);
  gt_hmm_set_transition_probability(hmm, PPT_OUT, PPT_UBOX, 0.05);
  gt_hmm_set_transition_probability(hmm, PPT_UBOX, PPT_OUT, 0.05);
  gt_hmm_set_transition_probability(hmm, PPT_UBOX, PPT_N,   0.05);
  gt_hmm_set_transition_probability(hmm, PPT_UBOX, PPT_IN,  0.05);
  gt_hmm_set_transition_probability(hmm, PPT_IN, PPT_UBOX,  0.05);
  gt_hmm_set_transition_probability(hmm, PPT_IN, PPT_OUT,   0.05);
  gt_hmm_set_transition_probability(hmm, PPT_IN, PPT_N,     0.05);
  gt_hmm_set_transition_probability(hmm, PPT_N, PPT_UBOX,   0.05);
  gt_hmm_set_transition_probability(hmm, PPT_N, PPT_OUT,    0.05);
  gt_hmm_set_transition_probability(hmm, PPT_N, PPT_IN,     0.05);
  gt_hmm_set_missing_transition_probabilities(hmm);

  if (!gt_hmm_is_valid(hmm))
  {
    gt_hmm_delete(hmm);
    return NULL;
  } else  return hmm;
}

static double gt_ppt_score(unsigned long radius, unsigned long end)
{
  double ret;
  unsigned long r2;
  r2 = radius * radius;
  ret = ((double) r2 - pow(abs((double) radius - (double) end),2))/(double) r2;
  return ret;
}

static int gt_ppt_hit_cmp(const void *h1, const void *h2)
{
  GtPPTHit* ph1 = *(GtPPTHit**) h1;
  GtPPTHit* ph2 = *(GtPPTHit**) h2;
  gt_assert(h1 && h2);
  return gt_double_compare(ph2->score, ph1->score);
}

static bool gt_ubox_ok(GtPPTHit *hit, GtRange uboxlen)
{
  gt_assert(hit);
  return (hit->state == PPT_UBOX
           && hit->rng.end-hit->rng.start+1 >= uboxlen.start
           && hit->rng.end-hit->rng.start+1 <= uboxlen.end);
}

static bool gt_ppt_ok(GtPPTHit *hit, GtRange pptlen)
{
  gt_assert(hit);
  return (hit->state == PPT_IN
           && hit->rng.end-hit->rng.start+1 >= pptlen.start
           && hit->rng.end-hit->rng.start+1 <= pptlen.end);
}

static void gt_group_hits(unsigned int *decoded, GtPPTResults *results,
                          unsigned long radius, GtStrand strand)
{
  GtPPTHit *cur_hit = NULL,
           *tmp = NULL;
  unsigned long i = 0;

  gt_assert(decoded && results && strand != GT_STRAND_UNKNOWN);

  /* group hits into stretches */
  cur_hit = gt_ppt_hit_new(strand, results);
  for (i=0;i<2*radius-1;i++)
  {
    cur_hit->state = decoded[i];
    cur_hit->rng.end = i;
    if (decoded[i+1] != decoded[i] || i+2==2*radius)
    {
      switch (cur_hit->state)
      {
        case PPT_UBOX:
          if (gt_ubox_ok(cur_hit, results->opts->ubox_len))
            tmp = cur_hit;
          else
          {
            gt_free(tmp);
            tmp = NULL;
            gt_free(cur_hit);
            cur_hit = NULL;
          }
          break;
        case PPT_IN:
          if (gt_ppt_ok(cur_hit, results->opts->ppt_len))
          {
            cur_hit->score = gt_ppt_score(radius, cur_hit->rng.end);
            gt_array_add(results->hits, cur_hit);
            if (tmp)
            {
              /* this PPT has a U-box, handle accordingly */
              cur_hit->ubox = tmp;
              tmp = NULL;
            }
          }
          else
          {
            gt_free(tmp);
            tmp = NULL;
            gt_free(cur_hit);
            cur_hit = NULL;
          }
          break;
        default:
          gt_free(tmp);
          tmp = NULL;
          gt_free(cur_hit);
          cur_hit = NULL;
          break;
      }
      if (i+2!=2*radius)
      {
        cur_hit = gt_ppt_hit_new(strand,results);
        cur_hit->rng.start = i+1;
      }
    }
  }
  if (cur_hit)
    cur_hit->rng.end++;
  gt_free(tmp);
}

GtPPTResults* gt_ppt_find(const char *seq,
                          const char *rev_seq,
                          GtLTRElement *element,
                          GtPPTOptions *o)
{
  unsigned int *encoded_seq=NULL,
               *decoded=NULL;
  const GtAlpha *alpha;
  GtHMM *hmm;
  GtPPTResults *results = NULL;
  unsigned long i = 0,
                radius = 0,
                seqlen = gt_ltrelement_length(element),
                ltrlen = 0;

  gt_assert(seq && rev_seq && element && o);

  results = gt_ppt_results_new(element, o);

  alpha = gt_alpha_new_dna();
  hmm = gt_ppt_hmm_new(alpha, o);

  /* do PPT finding on forward strand
   * -------------------------------- */
  ltrlen = gt_ltrelement_rightltrlen(element);
  /* make sure that we do not cross the LTR boundary */
  radius = MIN(o->radius, ltrlen-1);
  /* encode sequence */
  encoded_seq = gt_malloc(sizeof (unsigned int) * seqlen);
  for (i=0;i<seqlen;i++)
  {
    encoded_seq[i] = gt_alpha_encode(alpha, seq[i]);
  }
  /* use Viterbi algorithm to decode emissions within radius */
  decoded = gt_malloc(sizeof (unsigned int) * (2*radius+1));
  gt_hmm_decode(hmm, decoded,
                encoded_seq + (seqlen-1) - (ltrlen-1) - radius - 1,
                2*radius+1);
  gt_group_hits(decoded, results, radius, GT_STRAND_FORWARD);
  /* radius length may change in the next strand, so reallocate */
  gt_free(decoded);

  /* do PPT finding on reverse strand
   * -------------------------------- */
  ltrlen = gt_ltrelement_leftltrlen(element);
  /* make sure that we do not cross the LTR boundary */
  radius = MIN(o->radius, ltrlen-1);
  /* encode sequence */
  for (i=0;i<seqlen;i++)
  {
    encoded_seq[i] = gt_alpha_encode(alpha, rev_seq[i]);
  }
  /* use Viterbi algorithm to decode emissions within radius */
  decoded = gt_malloc(sizeof (unsigned int) * (2*radius+1));
  gt_hmm_decode(hmm, decoded,
                encoded_seq + (seqlen-1) - (ltrlen-1) - radius - 1,
                2*radius+1);
  gt_group_hits(decoded, results, radius, GT_STRAND_REVERSE);

  /* rank hits by descending score */
  gt_array_sort(results->hits, gt_ppt_hit_cmp);

  gt_free(encoded_seq);
  gt_free(decoded);
  gt_alpha_delete((GtAlpha*) alpha);
  gt_hmm_delete(hmm);

  return results;
}

void gt_ppt_results_delete(GtPPTResults *results)
{
  unsigned long i;

  if (!results) return;

  if (results->hits)
  {
    for (i=0;i<gt_array_size(results->hits);i++)
    {
      GtPPTHit * hit = *(GtPPTHit**) gt_array_get(results->hits,i);
      if (hit->ubox)
        gt_free(hit->ubox);
      gt_free(hit);
    }
    gt_array_delete(results->hits);
  }
  gt_free(results);
}

int gt_ppt_unit_test(GtError *err)
{
  int had_err = 0;
  GtPPTOptions o;
  GtPPTResults *rs;
  GtPPTHit *h;
  GtLTRElement element;
  GtRange rng;
  char *rev_seq,
       *seq,
       tmp[BUFSIZ];
  const char *fullseq =                           "aaaaaaaaaaaaaaaaaaaa"
                    "tatagcactgcatttcgaatatagtttcgaatatagcactgcatttcgaa"
                    "tatagcactgcatttcgaatatagtttcgaatatagcactgcatttcgaa"
                    "aaag"
                        "tcttctttct" /* <- PPT reverse */
                                  "aaaaatatagtttcgaatatagcactgcatttcgaa"
                    "tatagcactgcatttcgaatatagtttcgaatatagcactgcatttcgaa"
                    "tatagcactgcatttcgaatatagtttcgaatatagcactgcatttcgaa"
                    "tatagcactgcatttcgaatatagtttcgaatatagcactgcatttcgaa"
                    "tatagcactgcatttcgaatatagtttcgaatatagcactgcatttcgaa"
                    "tatagcactgcatttcgaatatagtttcgaatatagcactgcatttcgaa"
                    "tatagcactgcatttcgaatatagtttcgaatatagcactgcatttcgaa"
                    "tatagcactgcatttcgaatatagtttcgaatatttt"
                                   /* PPT forward -> */  "gggatagggggag"
                    "tatagcactgcatttcgaatatagtttcgaatatagcactgcatttcgaa"
                    "tatagcactgcatttcgaatatagtttcgaatatagcactgcatttcgaa"
                    "aaaaaaaaaaaaaaaaaaaa";

  seq     = gt_malloc(600 * sizeof (char));
  rev_seq = gt_malloc(600 * sizeof (char));
  memcpy(seq,     fullseq + 20, 600);
  memcpy(rev_seq, fullseq + 20, 600);
  gt_reverse_complement(rev_seq, 600, err);

  element.leftLTR_5 = 20;
  element.leftLTR_3 = 119;
  element.rightLTR_5 = 520;
  element.rightLTR_3 = 619;

  /* run PPT finding */
  memset(&o, 0, sizeof (GtPPTOptions));
  o.ppt_len.start = 5;
  o.ppt_len.end = 15;
  o.ubox_len.start = 2;
  o.ubox_len.end = 15;
  o.radius = 30;
  o.ppt_pyrimidine_prob = PPT_PYRIMIDINE_PROB;
  o.ppt_purine_prob = PPT_PURINE_PROB;
  o.bkg_a_prob = BKG_A_PROB;
  o.bkg_g_prob = BKG_G_PROB;
  o.bkg_t_prob = BKG_T_PROB;
  o.bkg_c_prob = BKG_C_PROB;
  o.ubox_u_prob = UBOX_U_PROB;
  rs = gt_ppt_find(seq, rev_seq, &element, &o);

  ensure(had_err, gt_ppt_results_get_number_of_hits(rs) == 2);
  h = gt_ppt_results_get_ranked_hit(rs, 0);
  ensure(had_err, h);
  rng = gt_ppt_hit_get_coords(h);
  ensure(had_err, rng.start == 507);
  ensure(had_err, rng.end == 519);
  memset(tmp, 0, BUFSIZ);
  memcpy(tmp, fullseq + (rng.start * sizeof (char)),
         (rng.end - rng.start + 1) * sizeof (char));
  ensure(had_err, strcmp(tmp, "gggatagggggag" ) == 0);

  h = gt_ppt_results_get_ranked_hit(rs, 1);
  ensure(had_err, h);
  rng = gt_ppt_hit_get_coords(h);
  ensure(had_err, rng.start == 124);
  ensure(had_err, rng.end == 133);
  memset(tmp, 0, BUFSIZ);
  memcpy(tmp, fullseq + (rng.start * sizeof (char)),
         (rng.end - rng.start + 1) * sizeof (char));
  ensure(had_err, strcmp(tmp, "tcttctttct" ) == 0);

  gt_free(rev_seq);
  gt_free(seq);

  gt_ppt_results_delete(rs);
  return had_err;
}
