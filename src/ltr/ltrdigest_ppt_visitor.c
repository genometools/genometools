/*
  Copyright (c) 2008-2013 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2008-2013 Center for Bioinformatics, University of Hamburg

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
#include <signal.h>
#include <string.h>
#include "core/array_api.h"
#include "core/log.h"
#include "core/ma.h"
#include "core/mathsupport.h"
#include "core/minmax.h"
#include "core/range_api.h"
#include "core/str_api.h"
#include "core/strand_api.h"
#include "core/symbol_api.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "core/warning_api.h"
#include "extended/node_visitor_api.h"
#include "extended/extract_feature_sequence.h"
#include "extended/feature_node.h"
#include "extended/feature_node_iterator_api.h"
#include "extended/feature_type.h"
#include "extended/hmm.h"
#include "extended/reverse_api.h"
#include "ltr/ltrdigest_def.h"
#include "ltr/ltrdigest_ppt_visitor.h"

struct GtLTRdigestPPTVisitor {
  const GtNodeVisitor parent_instance;
  GtRegionMapping *rmap;
  GtStr *tag;
  GtHMM *hmm;
  GtAlphabet *alpha;
  GtRange ppt_len, ubox_len;
  double ppt_pyrimidine_prob,
         ppt_purine_prob,
         bkg_a_prob,
         bkg_g_prob,
         bkg_t_prob,
         bkg_c_prob,
         ubox_u_prob;
  unsigned int radius,
               max_ubox_dist;
};

typedef struct GtPPTHit GtPPTHit;
typedef struct GtPPTResults GtPPTResults;

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
  GtRange leftltrrng,
          rightltrrng;
};

const GtNodeVisitorClass* gt_ltrdigest_ppt_visitor_class(void);

#define gt_ltrdigest_ppt_visitor_cast(GV)\
        gt_node_visitor_cast(gt_ltrdigest_ppt_visitor_class(), GV)

static GtPPTResults* gt_ppt_results_new(GtRange leftltrrng,
                                        GtRange rightltrrng)
{
  GtPPTResults *res = gt_calloc((size_t) 1, sizeof (GtPPTResults));
  res->leftltrrng = leftltrrng;
  res->rightltrrng = rightltrrng;
  res->hits = gt_array_new(sizeof (GtPPTHit*));
  return res;
}

static GtPPTHit* gt_ppt_hit_new(GtStrand strand, GtPPTResults *r)
{
  GtPPTHit *h = gt_calloc((size_t) 1, sizeof (GtPPTHit));
  gt_assert(h);
  h->strand = strand;
  h->res = r;
  h->score = 0.0;
  return h;
}

  GtRange gt_ppt_hit_get_coords(const GtPPTHit *h,
                                     GtLTRdigestPPTVisitor *lv)
{
  GtRange rng;
  gt_assert(h);
  rng.start = h->rng.start;
  rng.end = h->rng.end;
  switch (h->strand)
  {
    case GT_STRAND_FORWARD:
    default:
      rng.start = h->res->rightltrrng.start - 2 - lv->radius
                    + rng.start;
      rng.end = rng.start + (gt_range_length(&h->rng) - 1);
      break;
    case GT_STRAND_REVERSE:
      rng.end = h->res->leftltrrng.end + lv->radius - rng.start;
      rng.start = rng.end - (gt_range_length(&h->rng) - 1);
      break;
  }
  gt_assert(gt_range_length(&rng) == gt_range_length(&h->rng));
  return rng;
}

static GtPPTHit* gt_ppt_hit_get_ubox(const GtPPTHit *h)
{
  gt_assert(h);
  return h->ubox;
}

static GtStrand gt_ppt_hit_get_strand(const GtPPTHit *h)
{
  gt_assert(h);
  return h->strand;
}

static GtUword gt_ppt_results_get_number_of_hits(GtPPTResults *r)
{
  gt_assert(r);
  return gt_array_size(r->hits);
}

static GtPPTHit* gt_ppt_results_get_ranked_hit(GtPPTResults *r, GtUword i)
{
  gt_assert(r);
  return *(GtPPTHit**) gt_array_get(r->hits, i);
}

static GtHMM* gt_ppt_hmm_new(const GtAlphabet *alpha, GtLTRdigestPPTVisitor *lv)
{
  GtHMM *hmm;
  double non_u_prob = 0.0;

  gt_assert(alpha);

  hmm = gt_hmm_new((unsigned int) PPT_NOF_STATES, gt_alphabet_size(alpha));

  /* set emission probabilities */
  gt_hmm_set_emission_probability(hmm, (unsigned int) PPT_OUT,
                                  (unsigned int) gt_alphabet_encode(alpha, 'G'),
                                  lv->bkg_g_prob);
  gt_hmm_set_emission_probability(hmm, (unsigned int) PPT_OUT,
                                  (unsigned int) gt_alphabet_encode(alpha, 'A'),
                                  lv->bkg_a_prob);
  gt_hmm_set_emission_probability(hmm, (unsigned int) PPT_OUT,
                                  (unsigned int) gt_alphabet_encode(alpha, 'C'),
                                  lv->bkg_c_prob);
  gt_hmm_set_emission_probability(hmm, (unsigned int) PPT_OUT,
                                  (unsigned int) gt_alphabet_encode(alpha, 'T'),
                                  lv->bkg_t_prob);
  gt_hmm_set_emission_probability(hmm, (unsigned int) PPT_IN,
                                  (unsigned int) gt_alphabet_encode(alpha, 'G'),
                                  lv->ppt_purine_prob/2);
  gt_hmm_set_emission_probability(hmm, (unsigned int) PPT_IN,
                                  (unsigned int) gt_alphabet_encode(alpha, 'A'),
                                  lv->ppt_purine_prob/2);
  gt_hmm_set_emission_probability(hmm, (unsigned int) PPT_IN,
                                  (unsigned int) gt_alphabet_encode(alpha, 'C'),
                                  lv->ppt_pyrimidine_prob/2);
  gt_hmm_set_emission_probability(hmm, (unsigned int) PPT_IN,
                                  (unsigned int) gt_alphabet_encode(alpha, 'T'),
                                  lv->ppt_pyrimidine_prob/2);
  gt_hmm_set_emission_probability(hmm, (unsigned int) PPT_UBOX,
                                  (unsigned int) gt_alphabet_encode(alpha, 'T'),
                                  lv->ubox_u_prob);
  /* calculate non-U probabilities (still uniform, may be optimised) */
  non_u_prob = (1.0 - (lv->ubox_u_prob)) / ((double) 3);
  gt_hmm_set_emission_probability(hmm, (unsigned int) PPT_UBOX,
                                  (unsigned int) gt_alphabet_encode(alpha, 'G'),
                                  non_u_prob);
  gt_hmm_set_emission_probability(hmm, (unsigned int) PPT_UBOX,
                                  (unsigned int) gt_alphabet_encode(alpha, 'A'),
                                  non_u_prob);
  gt_hmm_set_emission_probability(hmm, (unsigned int) PPT_UBOX,
                                  (unsigned int) gt_alphabet_encode(alpha, 'C'),
                                  non_u_prob);
  gt_hmm_set_emission_probability(hmm, (unsigned int) PPT_N,
                                  (unsigned int) gt_alphabet_encode(alpha, 'G'),
                                  0.00);
  gt_hmm_set_emission_probability(hmm, (unsigned int) PPT_N,
                                  (unsigned int) gt_alphabet_encode(alpha, 'A'),
                                  0.00);
  gt_hmm_set_emission_probability(hmm, (unsigned int) PPT_N,
                                  (unsigned int) gt_alphabet_encode(alpha, 'C'),
                                  0.00);
  gt_hmm_set_emission_probability(hmm, (unsigned int) PPT_N,
                                  (unsigned int) gt_alphabet_encode(alpha, 'T'),
                                  0.00);
  gt_hmm_set_emission_probability(hmm, (unsigned int) PPT_N,
                                  (unsigned int) gt_alphabet_encode(alpha, 'N'),
                                  1.00);

  /* set transition probabilities */
  gt_hmm_set_transition_probability(hmm, (unsigned int) PPT_OUT,
                                    (unsigned int) PPT_IN, 0.05);
  gt_hmm_set_transition_probability(hmm, (unsigned int) PPT_OUT,
                                    (unsigned int) PPT_N, 0.05);
  gt_hmm_set_transition_probability(hmm, (unsigned int) PPT_OUT,
                                    (unsigned int) PPT_UBOX, 0.05);
  gt_hmm_set_transition_probability(hmm, (unsigned int) PPT_UBOX,
                                    (unsigned int) PPT_OUT, 0.05);
  gt_hmm_set_transition_probability(hmm, (unsigned int) PPT_UBOX,
                                    (unsigned int) PPT_N, 0.05);
  gt_hmm_set_transition_probability(hmm, (unsigned int) PPT_UBOX,
                                    (unsigned int) PPT_IN, 0.05);
  gt_hmm_set_transition_probability(hmm, (unsigned int) PPT_IN,
                                    (unsigned int) PPT_UBOX, 0.05);
  gt_hmm_set_transition_probability(hmm, (unsigned int) PPT_IN,
                                    (unsigned int) PPT_OUT, 0.05);
  gt_hmm_set_transition_probability(hmm, (unsigned int) PPT_IN,
                                    (unsigned int) PPT_N, 0.05);
  gt_hmm_set_transition_probability(hmm, (unsigned int) PPT_N,
                                    (unsigned int) PPT_UBOX, 0.05);
  gt_hmm_set_transition_probability(hmm, (unsigned int) PPT_N,
                                    (unsigned int) PPT_OUT, 0.05);
  gt_hmm_set_transition_probability(hmm, (unsigned int) PPT_N,
                                    (unsigned int) PPT_IN, 0.05);
  gt_hmm_set_missing_transition_probabilities(hmm);

  if (!gt_hmm_is_valid(hmm))
  {
    gt_hmm_delete(hmm);
    return NULL;
  } else
    return hmm;
}

static inline double gt_ppt_score(GtUword radius, GtUword end)
{
  double ret;
  GtUword r2;
  r2 = radius * radius;
  ret = ((double) r2
           - pow(fabs((double) radius - (double) end), 2.0))/(double) r2;
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

static void gt_group_hits(GtLTRdigestPPTVisitor *lv,
                          unsigned int *decoded, GtPPTResults *results,
                          GtUword radius, GT_UNUSED const char *seq,
                          GtStrand strand)
{
  GtPPTHit *cur_hit = NULL,
           *potential_ubox = NULL;
  GtUword i = 0;

  gt_assert(decoded && results && strand != GT_STRAND_UNKNOWN);

  /* group hits into stretches */
  cur_hit = gt_ppt_hit_new(strand, results);
  for (i = 0; i < 2 * radius - 1; i++)
  {
    gt_assert(cur_hit != NULL);
    cur_hit->state = (GtPPTStates) decoded[i];
    cur_hit->rng.end = i;
    if (decoded[i+1] != decoded[i] || i + 2 == 2 * radius)
    {
      switch (cur_hit->state)
      {
        case PPT_UBOX:
          if (gt_ubox_ok(cur_hit, lv->ubox_len)) {
            if (potential_ubox != NULL) {
              gt_free(potential_ubox);
              potential_ubox = NULL;
            }
            potential_ubox = cur_hit;
            cur_hit = NULL;
          }
          else
          {
            gt_free(cur_hit);
            cur_hit = NULL;
          }
          break;
        case PPT_IN:
          if (gt_ppt_ok(cur_hit, lv->ppt_len))
          {
            cur_hit->score = gt_ppt_score(radius, cur_hit->rng.end);
            gt_array_add(results->hits, cur_hit);
            if (potential_ubox != NULL)
            {
              if (cur_hit->rng.start - potential_ubox->rng.end
                    <= (GtUword) lv->max_ubox_dist)
              {
                /* this PPT has a U-box, handle accordingly */
                cur_hit->ubox = potential_ubox;
              } else {
                gt_free(potential_ubox);
              }
              potential_ubox = NULL;
            }
            cur_hit = NULL;
          }
          else
          {
            if (potential_ubox != NULL) {
              gt_free(potential_ubox);
              potential_ubox = NULL;
            }
            gt_free(cur_hit);
            cur_hit = NULL;
          }
          break;
        default:
          if (potential_ubox != NULL) {
                gt_free(potential_ubox);
                potential_ubox = NULL;
          }
          gt_free(cur_hit);
          cur_hit = NULL;
          break;
      }
      if (i+2!=2*radius)
      {
        /* we assume that the current hit has been processed */
        gt_assert(cur_hit == NULL);
        cur_hit = gt_ppt_hit_new(strand, results);
        cur_hit->rng.start = i+1;
      }
    }
  }
  if (cur_hit != NULL)
    cur_hit->rng.end++;
  gt_free(potential_ubox);
}

static GtPPTResults* gt_ppt_find(GtLTRdigestPPTVisitor *v,
                                 const char *seq, const char *rev_seq,
                                 GtUword seqlen,
                                 GtRange rightltrrng,
                                 GtRange leftltrrng)
{
  unsigned int *encoded_seq = NULL,
               *decoded = NULL;
  GtPPTResults *results = NULL;
  GtUword i = 0,
                radius = 0,
                ltrlen = 0;

  gt_assert(seq && rev_seq && v);

  results = gt_ppt_results_new(leftltrrng, rightltrrng);

    /* do PPT finding on forward strand
     -------------------------------- */
  ltrlen = gt_range_length(&rightltrrng);
  /* make sure that we do not cross the LTR boundary */
  radius = MIN((GtUword) v->radius, ltrlen-1);
  /* encode sequence */
  encoded_seq = gt_malloc(sizeof (unsigned int) * seqlen);
  for (i = 0; i < seqlen; i++) {
    encoded_seq[i] = (unsigned int) gt_alphabet_encode(v->alpha, seq[i]);
  }
  /* use Viterbi algorithm to decode emissions within radius */
  decoded = gt_malloc(sizeof (unsigned int) * (2 * radius + 1));
  gt_hmm_decode(v->hmm, decoded,
                encoded_seq + (seqlen-1) - (ltrlen-1) - radius - 1,
                (unsigned int) (2 * radius + 1));
  gt_group_hits(v, decoded, results, radius,
                seq + (seqlen - 1) - (ltrlen - 1) - radius - 1,
                GT_STRAND_FORWARD);
  /* radius length may change in the next strand, so reallocate */
  gt_free(decoded);

  /* do PPT finding on reverse strand
     -------------------------------- */
  ltrlen = gt_range_length(&leftltrrng);
  /* make sure that we do not cross the LTR boundary */
  radius = MIN((GtUword) v->radius, ltrlen - 1);
  /* encode sequence */
  for (i = 0; i < seqlen; i++) {
    encoded_seq[i] = (unsigned int) gt_alphabet_encode(v->alpha, rev_seq[i]);
  }
  /* use Viterbi algorithm to decode emissions within radius */
  decoded = gt_malloc(sizeof (unsigned int) * (2 * radius + 1));
  gt_hmm_decode(v->hmm, decoded,
                encoded_seq + (seqlen - 1) - (ltrlen - 1) - radius - 1,
                (unsigned int) (2 * radius + 1));
  gt_group_hits(v, decoded, results, radius,
                rev_seq + (seqlen - 1) - (ltrlen - 1) - radius - 1,
                GT_STRAND_REVERSE);

  /* rank hits by descending score */
  gt_array_sort(results->hits, gt_ppt_hit_cmp);

  gt_free(encoded_seq);
  gt_free(decoded);

  return results;
}

void gt_ppt_results_delete(GtPPTResults *results)
{
  GtUword i;
  if (results == NULL) return;

  if (results->hits != NULL)
  {
    for (i=0;i<gt_array_size(results->hits);i++)
    {
      GtPPTHit * hit = *(GtPPTHit**) gt_array_get(results->hits,i);
      if (hit->ubox != NULL)
        gt_free(hit->ubox);
      gt_free(hit);
    }
    gt_array_delete(results->hits);
  }
  gt_free(results);
}

static void ppt_attach_results_to_gff3(GtLTRdigestPPTVisitor *lv,
                                       GtPPTResults *results,
                                       GtFeatureNode *mainnode,
                                       GtStrand *canonical_strand)
{
  GtRange ppt_range;
  GtUword i = 0;
  GtGenomeNode *gf;
  GtPPTHit *hit = gt_ppt_results_get_ranked_hit(results, i++),
           *ubox;
  if (*canonical_strand == GT_STRAND_UNKNOWN)
    *canonical_strand = hit->strand;
  else {
    /* find best-scoring PPT on the given canonical strand */
    while (hit->strand != *canonical_strand
             && i < gt_ppt_results_get_number_of_hits(results)) {
      gt_log_log("dropping PPT because of nonconsistent strand: %s\n",
                 gt_feature_node_get_attribute(mainnode, "ID"));
      hit = gt_ppt_results_get_ranked_hit(results, i++);
    }
    /* if there is none, do not report a PPT */
    if (hit->strand != *canonical_strand)
      return;
  }
  ppt_range = gt_ppt_hit_get_coords(hit, lv);
  ppt_range.start++; ppt_range.end++;  /* GFF3 is 1-based */
  gf = gt_feature_node_new(gt_genome_node_get_seqid((GtGenomeNode*) mainnode),
                           gt_ft_RR_tract,
                           ppt_range.start,
                           ppt_range.end,
                           hit->strand);
  gt_feature_node_set_source((GtFeatureNode*) gf, lv->tag);
  gt_feature_node_set_strand(mainnode, hit->strand);
  gt_feature_node_add_child(mainnode, (GtFeatureNode*) gf);
  if ((ubox = gt_ppt_hit_get_ubox(hit)) != NULL) {
    GtRange ubox_range = gt_ppt_hit_get_coords(ubox, lv);
    ubox_range.start++; ubox_range.end++;
    gf = gt_feature_node_new(gt_genome_node_get_seqid((GtGenomeNode*) mainnode),
                             gt_ft_U_box,
                             ubox_range.start,
                             ubox_range.end,
                             gt_ppt_hit_get_strand(ubox));
    gt_feature_node_set_source((GtFeatureNode*) gf, lv->tag);
    gt_feature_node_set_strand(mainnode, gt_ppt_hit_get_strand(ubox));
    gt_feature_node_add_child(mainnode, (GtFeatureNode*) gf);
  }
}

static int gt_ltrdigest_ppt_visitor_feature_node(GtNodeVisitor *nv,
                                                 GtFeatureNode *fn,
                                                 GtError *err)
{
  GT_UNUSED GtLTRdigestPPTVisitor *lv;
  GtFeatureNodeIterator *fni;
  GtRange leftltrrng = { 0, 0 }, rightltrrng = { 0, 0 };
  bool seen_left = false;
  GtFeatureNode *curnode = NULL,
                *ltr_retrotrans = NULL;
  int had_err = 0;
  lv = gt_ltrdigest_ppt_visitor_cast(nv);
  gt_assert(lv);
  gt_error_check(err);

  /* traverse annotation subgraph and find LTR element */
  fni = gt_feature_node_iterator_new(fn);
  while (!had_err && (curnode = gt_feature_node_iterator_next(fni))) {
    if (strcmp(gt_feature_node_get_type(curnode),
               gt_ft_LTR_retrotransposon) == 0) {
      ltr_retrotrans = curnode;
    }
    if (strcmp(gt_feature_node_get_type(curnode),
               gt_ft_long_terminal_repeat) == 0) {
      if (seen_left)   /* XXX */
        rightltrrng = gt_genome_node_get_range((GtGenomeNode*) curnode);
      else {
        leftltrrng = gt_genome_node_get_range((GtGenomeNode*) curnode);
        seen_left = true;
      }
    }
  }
  gt_feature_node_iterator_delete(fni);

  if (!had_err && ltr_retrotrans != NULL) {
    char *rev_seq;
    GtRange rng;
    GtStrand canonical_strand = gt_feature_node_get_strand(ltr_retrotrans);
    GtPPTResults *res;
    GtUword seqlen;
    GtStr *seq = NULL;
    rng = gt_genome_node_get_range((GtGenomeNode*) ltr_retrotrans);

    if (gt_range_length(&rng) < (GtUword) 10UL) {
      gt_warning("LTR_retrotransposon (%s, line %u) is too short for "
                 "PPT detection (" GT_WU " nt), skipped this step",
            gt_genome_node_get_filename((GtGenomeNode*) ltr_retrotrans),
            gt_genome_node_get_line_number((GtGenomeNode*) ltr_retrotrans),
            gt_range_length(&rng));
      return had_err;
    }

    seq = gt_str_new();
    seqlen = gt_range_length(&rng);

    had_err = gt_extract_feature_sequence(seq, (GtGenomeNode*) ltr_retrotrans,
                                          gt_symbol(gt_ft_LTR_retrotransposon),
                                          false, NULL, NULL, lv->rmap, err);

    if (!had_err) {
      rev_seq = gt_malloc((size_t) seqlen * sizeof (char));
      strncpy(rev_seq, gt_str_get(seq), (size_t) seqlen * sizeof (char));
      (void) gt_reverse_complement(rev_seq, seqlen, NULL);

      res = gt_ppt_find(lv, gt_str_get(seq), rev_seq, seqlen, rightltrrng,
                        leftltrrng);
      if (gt_ppt_results_get_number_of_hits(res) > 0) {
        ppt_attach_results_to_gff3(lv, res, ltr_retrotrans, &canonical_strand);
      }
      gt_ppt_results_delete(res);
      gt_free(rev_seq);
    }
    gt_str_delete(seq);
  }

  return had_err;
}

void gt_ltrdigest_ppt_visitor_free(GtNodeVisitor *nv)
{
  GT_UNUSED GtLTRdigestPPTVisitor *lv;
  if (!nv) return;
  lv = gt_ltrdigest_ppt_visitor_cast(nv);
  gt_str_delete(lv->tag);
  gt_alphabet_delete(lv->alpha);
  gt_hmm_delete(lv->hmm);
}

const GtNodeVisitorClass* gt_ltrdigest_ppt_visitor_class(void)
{
  static const GtNodeVisitorClass *nvc = NULL;
  if (!nvc) {
    nvc = gt_node_visitor_class_new(sizeof (GtLTRdigestPPTVisitor),
                                   gt_ltrdigest_ppt_visitor_free,
                                   NULL,
                                   gt_ltrdigest_ppt_visitor_feature_node,
                                   NULL,
                                   NULL,
                                   NULL);
  }
  return nvc;
}

GtNodeVisitor* gt_ltrdigest_ppt_visitor_new(GtRegionMapping *rmap,
                                            GtRange ppt_len,
                                            GtRange ubox_len,
                                            double ppt_pyrimidine_prob,
                                            double ppt_purine_prob,
                                            double bkg_a_prob,
                                            double bkg_g_prob,
                                            double bkg_t_prob,
                                            double bkg_c_prob,
                                            double ubox_u_prob,
                                            unsigned int radius,
                                            unsigned int max_ubox_dist,
                                            GT_UNUSED GtError *err)
{
  GtNodeVisitor *nv = NULL;
  GtLTRdigestPPTVisitor *lv;
  GT_UNUSED int had_err = 0, i;
  gt_assert(rmap);
  nv = gt_node_visitor_create(gt_ltrdigest_ppt_visitor_class());
  lv = gt_ltrdigest_ppt_visitor_cast(nv);
  lv->ppt_len = ppt_len;
  lv->ubox_len = ubox_len;
  lv->ppt_pyrimidine_prob = ppt_pyrimidine_prob;
  lv->ppt_purine_prob = ppt_purine_prob;
  lv->bkg_a_prob = bkg_a_prob;
  lv->bkg_g_prob = bkg_g_prob;
  lv->bkg_t_prob = bkg_t_prob;
  lv->bkg_c_prob = bkg_c_prob;
  lv->ubox_u_prob = ubox_u_prob;
  lv->radius = radius;
  lv->rmap = rmap;
  lv->max_ubox_dist = max_ubox_dist;
  lv->tag = gt_str_new_cstr(GT_LTRDIGEST_TAG);
  lv->alpha = gt_alphabet_new_dna();
  lv->hmm = gt_ppt_hmm_new(lv->alpha, lv);
  if (!lv->hmm) {
    gt_node_visitor_delete(nv);
    gt_error_set(err, "PPT HMM parameters are not valid!");
    had_err = -1;
    nv = NULL;
  }
  return nv;
}
