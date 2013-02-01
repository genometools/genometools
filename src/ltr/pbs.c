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
#include "core/ma.h"
#include "core/mathsupport.h"
#include "core/minmax.h"
#include "core/xansi_api.h"
#include "core/array.h"
#include "core/bioseq.h"
#include "core/ensure.h"
#include "core/fa.h"
#include "core/fileutils_api.h"
#include "extended/reverse_api.h"
#include "extended/swalign.h"
#include "ltr/pbs.h"

struct GtPBSHit {
  unsigned long start,
                end,
                edist,
                offset,
                tstart,
                alilen;
  GtStrand strand;
  double score;
  const char *trna;
  GtPBSResults *res;
};

struct GtPBSResults {
  GtArray *hits;
  GtLTRElement *elem;
  GtPBSOptions *opts;
};

static GtPBSHit* gt_pbs_hit_new(unsigned long alilen, GtStrand strand,
                                const char *tRNA, unsigned long tstart,
                                unsigned long start, unsigned long end,
                                unsigned long offset, unsigned long edist,
                                double score, GtPBSResults *r)
{
  GtPBSHit *hit = gt_malloc(sizeof (GtPBSHit));
  hit->alilen  = alilen;
  hit->strand  = strand;
  hit->trna    = tRNA;
  hit->tstart  = tstart;
  hit->start   = start;
  hit->end     = end;
  hit->offset  = offset;
  hit->edist   = edist;
  hit->score   = score;
  hit->res     = r;
  return hit;
}

static GtPBSResults* gt_pbs_results_new(GtLTRElement *elem,
                                        GtPBSOptions *opts)
{
  GtPBSResults *res = gt_calloc((size_t) 1, sizeof (GtPBSResults));
  res->elem = elem;
  res->opts = opts;
  res->hits = gt_array_new(sizeof (GtPBSHit*));
  return res;
}

GtRange gt_pbs_hit_get_coords(const GtPBSHit *h)
{
  GtRange rng;
  gt_assert(h && h->end >= h->start);
  rng.start = h->start;
  rng.end = h->end;
  switch (h->strand)
  {
    case GT_STRAND_FORWARD:
    default:
      rng.start = h->res->elem->leftLTR_3 + 1 - h->res->opts->radius
                    + rng.start;
      rng.end = rng.start + (h->end - h->start);
      break;
    case GT_STRAND_REVERSE:
      rng.end = h->res->elem->rightLTR_5 - 1 + h->res->opts->radius - rng.start;
      rng.start = rng.end - (h->end - h->start);
      break;
  }
  gt_assert(gt_range_length(&rng) == (h->end - h->start + 1));
  return rng;
}

unsigned long gt_pbs_hit_get_offset(const GtPBSHit *h)
{
  gt_assert(h);
  return h->offset;
}

double gt_pbs_hit_get_score(const GtPBSHit *h)
{
  gt_assert(h);
  return h->score;
}

unsigned long gt_pbs_hit_get_edist(const GtPBSHit *h)
{
  gt_assert(h);
  return h->edist;
}

const char* gt_pbs_hit_get_trna(const GtPBSHit *h)
{
  gt_assert(h);
  return h->trna;
}

unsigned long gt_pbs_hit_get_tstart(const GtPBSHit *h)
{
  gt_assert(h);
  return h->tstart;
}

GtStrand gt_pbs_hit_get_strand(const GtPBSHit *h)
{
  gt_assert(h);
  return h->strand;
}

unsigned long gt_pbs_hit_get_alignment_length(const GtPBSHit *h)
{
  gt_assert(h);
  return h->alilen;
}

unsigned long gt_pbs_results_get_number_of_hits(const GtPBSResults *r)
{
  gt_assert(r);
  return gt_array_size(r->hits);
}

GtPBSHit* gt_pbs_results_get_ranked_hit(const GtPBSResults *r, unsigned long i)
{
  gt_assert(r);
  return *(GtPBSHit**) gt_array_get(r->hits, i);
}

static GtScoreFunction* gt_dna_scorefunc_new(GtAlphabet *a, int match,
                                             int mismatch, int insertion,
                                             int deletion)
{
  GtScoreMatrix *sm = gt_score_matrix_new(a);
  GtScoreFunction *sf = gt_score_function_new(sm, insertion, deletion);
  unsigned int m,n;

  for (m=0;m<gt_alphabet_size(a);m++)
  {
    for (n=0;n<gt_alphabet_size(a);n++)
    {
      gt_score_matrix_set_score(sm, m, n, (n==m ? match : mismatch));
    }
  }
  /* make N-N a mismatch! */
  gt_score_matrix_set_score(sm, (unsigned int) gt_alphabet_encode(a, 'n'),
                            (unsigned int) gt_alphabet_encode(a, 'n'),
                            mismatch);
  return sf;
}

static double gt_pbs_score_func(unsigned long edist, unsigned long offset,
                                unsigned long alilen, unsigned long trnalen,
                                unsigned long trna_offset)
{
  double penalties;
  if (edist == 0 || offset == 0)
    penalties = 1.0;
  else
    penalties = (double) edist * (double) offset;
  return ((double) alilen *
             (double) ((trnalen - trna_offset)/(double) trnalen))
          /penalties;
}

static inline unsigned long gt_ulongabs(unsigned long nr1, unsigned long nr2)
{
  if (nr1 == nr2) return 0UL;
  if (nr1 > nr2) {
    return nr1 - nr2;
  } else {
    return nr2 - nr1;
  }
}

static void gt_pbs_add_hit(GtArray *hitlist, GtAlignment *ali, GtPBSOptions *o,
                           unsigned long trna_seqlen, const char *desc,
                           GtStrand strand, GtPBSResults *r)
{
  unsigned long dist;
  GtPBSHit *hit;
  unsigned long offset, alilen;
  GtRange urange, vrange;
  gt_assert(hitlist && desc);

  if (!ali) return;
  gt_assert(ali);

  dist = gt_alignment_eval(ali);
  urange = gt_alignment_get_urange(ali);
  vrange = gt_alignment_get_vrange(ali);
  offset = gt_ulongabs((unsigned long) o->radius, urange.start);
  alilen = gt_ulongabs(urange.end, urange.start)+1;

  if (dist <= (unsigned long) o->max_edist
        && offset <= o->offsetlen.end
        && offset >= o->offsetlen.start
        && alilen <= o->alilen.end
        && alilen >= o->alilen.start
        && vrange.start <= o->trnaoffsetlen.end
        && vrange.start >= o->trnaoffsetlen.start)
  {
    hit = gt_pbs_hit_new(alilen,
                         strand,
                         desc,
                         vrange.start,
                         urange.start,
                         urange.end,
                         offset,
                         dist,
                         gt_pbs_score_func(dist,
                                           offset,
                                           urange.end-urange.start+1,
                                           trna_seqlen,
                                           vrange.start),
                         r);
    gt_array_add(hitlist, hit);
  }
}

static int gt_pbs_hit_compare(const void *h1, const void *h2)
{
  GtPBSHit *hp1 = *(GtPBSHit**) h1;
  GtPBSHit *hp2 = *(GtPBSHit**) h2;

  return (gt_double_compare(hp2->score, hp1->score));
}

GtPBSResults* gt_pbs_find(const char *seq,
                          const char *rev_seq,
                          GtLTRElement *element,
                          GtPBSOptions *o,
                          GtError *err)
{
  GtSeq *seq_forward, *seq_rev;
  GtPBSResults *results;
  unsigned long j;
  GtAlignment *ali;
  GtAlphabet *a = gt_alphabet_new_dna();
  GtScoreFunction *sf = gt_dna_scorefunc_new(a,
                                             o->ali_score_match,
                                             o->ali_score_mismatch,
                                             o->ali_score_insertion,
                                             o->ali_score_deletion);

  gt_assert(seq && rev_seq && sf && a && element);

  results = gt_pbs_results_new(element, o);

  seq_forward = gt_seq_new(seq + (gt_ltrelement_leftltrlen(element))
                               - (o->radius),
                           (unsigned long) (2*o->radius + 1),
                           a);

  seq_rev     = gt_seq_new(rev_seq + (gt_ltrelement_rightltrlen(element))
                                   - (o->radius),
                           (unsigned long) (2*o->radius + 1),
                           a);

  for (j=0;j<gt_bioseq_number_of_sequences(o->trna_lib);j++)
  {
    GtSeq *trna_seq, *trna_from3;
    char *trna_from3_full;
    unsigned long trna_seqlen;

    trna_seq = gt_bioseq_get_seq(o->trna_lib, j);
    trna_seqlen = gt_seq_length(trna_seq);

    trna_from3_full = gt_calloc((size_t) trna_seqlen, sizeof (char));
    memcpy(trna_from3_full, gt_seq_get_orig(trna_seq),
           sizeof (char)*trna_seqlen);
    (void) gt_reverse_complement(trna_from3_full, trna_seqlen, err);
    trna_from3 = gt_seq_new_own(trna_from3_full, trna_seqlen, a);

    ali = gt_swalign(seq_forward, trna_from3, sf);
    gt_pbs_add_hit(results->hits, ali, o, trna_seqlen,
                   gt_seq_get_description(trna_seq), GT_STRAND_FORWARD,
                   results);
    gt_alignment_delete(ali);

    ali = gt_swalign(seq_rev, trna_from3, sf);
    gt_pbs_add_hit(results->hits, ali, o, trna_seqlen,
                   gt_seq_get_description(trna_seq), GT_STRAND_REVERSE,
                   results);
    gt_alignment_delete(ali);

    gt_seq_delete(trna_seq);
    gt_seq_delete(trna_from3);
  }
  gt_seq_delete(seq_forward);
  gt_seq_delete(seq_rev);
  gt_score_function_delete(sf);
  gt_alphabet_delete(a);
  gt_array_sort(results->hits, gt_pbs_hit_compare);
  return results;
}

void gt_pbs_results_delete(GtPBSResults *results)
{
    unsigned long i;
    if (!results) return;
    if (results->hits != NULL)
    {
      for (i=0;i<gt_array_size(results->hits);i++)
      {
        GtPBSHit *hit = *(GtPBSHit**) gt_array_get(results->hits, i);
        gt_free(hit);
      }
      gt_array_delete(results->hits);
    }
    gt_free(results);
}

int gt_pbs_unit_test(GtError *err)
{
  int had_err = 0;
  GtLTRElement element;
  GtPBSOptions o;
  GtStr *tmpfilename;
  FILE *tmpfp;
  GtPBSResults *res;
  GtPBSHit *hit;
  double score1, score2;
  GtRange rng;
  char *rev_seq,
       *seq,
       tmp[BUFSIZ];
  const char *fullseq =                           "aaaaaaaaaaaaaaaaaaaa"
                    "tatagcactgcatttcgaatatagtttcgaatatagcactgcatttcgaa"
                    "tatagcactgcatttcgaatatagtttcgaatatagcactgcatttcgaa"
                    "acatactaggatgctag" /* <- PBS forward */
                                     "aatatagtttcgaatatagcactgcatttcgaa"
                    "tatagcactgcatttcgaatatagtttcgaatatagcactgcatttcgaa"
                    "tatagcactgcatttcgaatatagtttcgaatatagcactgcatttcgaa"
                    "tatagcactgcatttcgaatatagtttcgaatatagcactgcatttcgaa"
                    "tatagcactgcatttcgaatatagtttcgaatatagcactgcatttcgaa"
                    "tatagcactgcatttcgaatatagtttcgaatatagcactgcatttcgaa"
                    "tatagcactgcatttcgaatatagtttcgaatatagcactgcatttcgaa"
                    "tatagcactgcatttcgaatatagtttcgaatatag"
                                   /* PBS reverse -> */ "gatcctaaggctac"
                    "tatagcactgcatttcgaatatagtttcgaatatagcactgcatttcgaa"
                    "tatagcactgcatttcgaatatagtttcgaatatagcactgcatttcgaa"
                    "aaaaaaaaaaaaaaaaaaaa";

  /* notice previous errors */
  gt_error_check(err);

  /* create temporary tRNA library file */
  tmpfilename = gt_str_new();
  tmpfp = gt_xtmpfp(tmpfilename);
  fprintf(tmpfp, ">test1\nccccccccccccccctagcatcctagtatgtccc\n"
                 ">test2\ncccccccccgatcctagggctaccctttc\n");
  gt_fa_xfclose(tmpfp);
  gt_ensure(had_err, gt_file_exists(gt_str_get(tmpfilename)));

  /* setup testing parameters */
  o.radius = 30U;
  o.max_edist = 1U;
  o.alilen.start = 11UL;
  o.alilen.end = 30UL;
  o.offsetlen.start = 0UL;
  o.offsetlen.end = 5UL;
  o.trnaoffsetlen.start = 0UL;
  o.trnaoffsetlen.end =  40UL;
  o.ali_score_match = 5;
  o.ali_score_mismatch = -10;
  o.ali_score_insertion = o.ali_score_deletion = -20;
  o.trna_lib = gt_bioseq_new(gt_str_get(tmpfilename), err);
  gt_ensure(had_err, gt_bioseq_number_of_sequences(o.trna_lib) == 2UL);

  element.leftLTR_5 = 20UL;
  element.leftLTR_3 = 119UL;
  element.rightLTR_5 = 520UL;
  element.rightLTR_3 = 619UL;

  /* setup sequences */
  seq     = gt_malloc((size_t) 600 * sizeof (char));
  rev_seq = gt_malloc((size_t) 600 * sizeof (char));
  memcpy(seq,     fullseq + 20, (size_t) 600);
  memcpy(rev_seq, fullseq + 20, (size_t) 600);
  gt_ensure(had_err,
            !gt_reverse_complement(rev_seq, (unsigned long) 600, NULL));

  /* try to find PBS in sequences */
  res = gt_pbs_find(seq, rev_seq, &element, &o, err);
  gt_ensure(had_err, res != NULL);
  gt_ensure(had_err, gt_pbs_results_get_number_of_hits(res) == 2UL);

  /* check first hit on forward strand */
  hit = gt_pbs_results_get_ranked_hit(res, 0UL);
  gt_ensure(had_err, hit != NULL);
  gt_ensure(had_err, gt_pbs_hit_get_alignment_length(hit) == 17UL);
  gt_ensure(had_err, gt_pbs_hit_get_edist(hit) == 0UL);
  gt_ensure(had_err, gt_pbs_hit_get_offset(hit) == 0UL);
  gt_ensure(had_err, gt_pbs_hit_get_tstart(hit) == 3UL);
  gt_ensure(had_err, strcmp(gt_pbs_hit_get_trna(hit), "test1") == 0);
  rng = gt_pbs_hit_get_coords(hit);
  gt_ensure(had_err, rng.start == 120UL);
  gt_ensure(had_err, rng.end == 136UL);
  score1 = gt_pbs_hit_get_score(hit);
  gt_ensure(had_err, gt_pbs_hit_get_strand(hit) == GT_STRAND_FORWARD);
  memset(tmp, 0, BUFSIZ-1);
  memcpy(tmp, fullseq + (rng.start * sizeof (char)),
         (size_t) ((rng.end - rng.start + 1) * sizeof (char)));
  gt_ensure(had_err, strcmp(tmp, "acatactaggatgctag" ) == 0);

  /* check second hit on reverse strand */
  hit = gt_pbs_results_get_ranked_hit(res, 1UL);
  gt_ensure(had_err, hit != NULL);
  gt_ensure(had_err, gt_pbs_hit_get_alignment_length(hit) == 14UL);
  gt_ensure(had_err, gt_pbs_hit_get_edist(hit) == 1UL);
  gt_ensure(had_err, gt_pbs_hit_get_offset(hit) == 0UL);
  gt_ensure(had_err, gt_pbs_hit_get_tstart(hit) == 6UL);
  gt_ensure(had_err, strcmp(gt_pbs_hit_get_trna(hit), "test2") == 0);
  rng = gt_pbs_hit_get_coords(hit);
  gt_ensure(had_err, rng.start == 506UL);
  gt_ensure(had_err, rng.end == 519UL);
  score2 = gt_pbs_hit_get_score(hit);
  gt_ensure(had_err, gt_double_compare(score1, score2) > 0);
  gt_ensure(had_err, gt_pbs_hit_get_strand(hit) == GT_STRAND_REVERSE);
  memset(tmp, 0, BUFSIZ-1);
  memcpy(tmp, fullseq + (rng.start * sizeof (char)),
         (size_t) ((rng.end - rng.start + 1) * sizeof (char)));
  gt_ensure(had_err, strcmp(tmp, "gatcctaaggctac" ) == 0);

  /* clean up */
  gt_xremove(gt_str_get(tmpfilename));
  gt_ensure(had_err, !gt_file_exists(gt_str_get(tmpfilename)));
  gt_str_delete(tmpfilename);
  gt_bioseq_delete(o.trna_lib);
  gt_free(rev_seq);
  gt_free(seq);
  gt_pbs_results_delete(res);

  return had_err;
}
