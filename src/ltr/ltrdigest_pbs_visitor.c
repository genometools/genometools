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

#include <string.h>
#include "core/array_api.h"
#include "core/ensure.h"
#include "core/fa.h"
#include "core/fileutils_api.h"
#include "core/log.h"
#include "core/ma.h"
#include "core/mathsupport.h"
#include "core/minmax.h"
#include "core/range_api.h"
#include "core/score_function.h"
#include "core/str_api.h"
#include "core/strand_api.h"
#include "core/symbol_api.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "extended/alignment.h"
#include "extended/extract_feature_sequence.h"
#include "extended/feature_node.h"
#include "extended/feature_node_iterator_api.h"
#include "extended/feature_type.h"
#include "extended/node_visitor_api.h"
#include "extended/reverse_api.h"
#include "extended/swalign.h"
#include "ltr/ltrdigest_def.h"
#include "ltr/ltrdigest_pbs_visitor.h"

struct GtLTRdigestPBSVisitor{
  const GtNodeVisitor parent_instance;
  GtRegionMapping *rmap;
  GtStr *tag;
  GtFeatureNode *ltr_retrotrans;
  GtUword leftLTR_3,
                rightLTR_5,
                leftltrlen,
                rightltrlen;
  unsigned int radius,
               max_edist;
  GtRange alilen,
          offsetlen,
          trnaoffsetlen;
  int ali_score_match,
      ali_score_mismatch,
      ali_score_insertion,
      ali_score_deletion;
  GtBioseq *trna_lib;
};

typedef struct {
  GtArray *hits;
} GtPBSResults;

typedef struct {
  GtUword start,
                end,
                edist,
                offset,
                tstart,
                alilen;
  GtStrand strand;
  double score;
  const char *trna;
  GtPBSResults *res;
} GtPBSHit;

static GtPBSHit* gt_pbs_hit_new(GtUword alilen, GtStrand strand,
                                const char *tRNA, GtUword tstart,
                                GtUword start, GtUword end,
                                GtUword offset,
                                GtUword edist, double score,
                                GtPBSResults *r)
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

static GtPBSResults* gt_pbs_results_new(void)
{
  GtPBSResults *res = gt_calloc((size_t) 1, sizeof (GtPBSResults));
  res->hits = gt_array_new(sizeof (GtPBSHit*));
  return res;
}

static GtRange gt_pbs_hit_get_coords(GtLTRdigestPBSVisitor *lv,
                                     const GtPBSHit *h)
{
  GtRange rng;
  gt_assert(h && h->end >= h->start);
  rng.start = h->start;
  rng.end = h->end;
  switch (h->strand)
  {
    case GT_STRAND_FORWARD:
    default:
      rng.start = lv->leftLTR_3 + 1 - lv->radius
                    + rng.start;
      rng.end = rng.start + (h->end - h->start);
      break;
    case GT_STRAND_REVERSE:
      rng.end = lv->rightLTR_5 - 1 + lv->radius - rng.start;
      rng.start = rng.end - (h->end - h->start);
      break;
  }
  gt_assert(gt_range_length(&rng) == (h->end - h->start + 1));
  return rng;
}

static GtUword gt_pbs_hit_get_offset(const GtPBSHit *h)
{
  gt_assert(h);
  return h->offset;
}

static double gt_pbs_hit_get_score(const GtPBSHit *h)
{
  gt_assert(h);
  return h->score;
}

static GtUword gt_pbs_hit_get_edist(const GtPBSHit *h)
{
  gt_assert(h);
  return h->edist;
}

static const char* gt_pbs_hit_get_trna(const GtPBSHit *h)
{
  gt_assert(h);
  return h->trna;
}

static GtUword gt_pbs_hit_get_tstart(const GtPBSHit *h)
{
  gt_assert(h);
  return h->tstart;
}

static GtStrand gt_pbs_hit_get_strand(const GtPBSHit *h)
{
  gt_assert(h);
  return h->strand;
}

static GtUword gt_pbs_results_get_number_of_hits(const GtPBSResults *r)
{
  gt_assert(r);
  return gt_array_size(r->hits);
}

static GtPBSHit* gt_pbs_results_get_ranked_hit(const GtPBSResults *r,
                                               GtUword i)
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

static double gt_pbs_score_func(GtUword edist, GtUword offset,
                                GtUword alilen, GtUword trnalen,
                                GtUword trna_offset)
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

static inline GtUword gt_ulongabs(GtUword nr1, GtUword nr2)
{
  if (nr1 == nr2) return 0UL;
  if (nr1 > nr2)
    return nr1 - nr2;
  else
    return nr2 - nr1;
}

static void gt_pbs_add_hit(GtLTRdigestPBSVisitor *lv, GtArray *hitlist,
                           GtAlignment *ali, GtUword trna_seqlen,
                           const char *desc, GtStrand strand, GtPBSResults *r)
{
  GtUword dist;
  GtPBSHit *hit;
  GtUword offset, alilen;
  GtRange urange, vrange;
  gt_assert(hitlist && desc);

  if (!ali) return;
  gt_assert(ali);

  dist = gt_alignment_eval(ali);
  urange = gt_alignment_get_urange(ali);
  vrange = gt_alignment_get_vrange(ali);
  offset = gt_ulongabs((GtUword) lv->radius, urange.start);
  alilen = gt_ulongabs(urange.end, urange.start)+1;

  if (dist <= (GtUword) lv->max_edist
        && offset <= lv->offsetlen.end
        && offset >= lv->offsetlen.start
        && alilen <= lv->alilen.end
        && alilen >= lv->alilen.start
        && vrange.start <= lv->trnaoffsetlen.end
        && vrange.start >= lv->trnaoffsetlen.start)
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
                                           urange.end-urange.start + 1,
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

static GtPBSResults* gt_pbs_find(GtLTRdigestPBSVisitor *lv, const char *seq,
                          const char *rev_seq, GtError *err)
{
  GtSeq *seq_forward, *seq_rev;
  GtPBSResults *results;
  GtUword j;
  GtAlignment *ali;
  GtAlphabet *a = gt_alphabet_new_dna();
  GtScoreFunction *sf;
  gt_assert(lv && seq && rev_seq);

  sf = gt_dna_scorefunc_new(a, lv->ali_score_match, lv->ali_score_mismatch,
                            lv->ali_score_insertion, lv->ali_score_deletion);

  results = gt_pbs_results_new();

  seq_forward = gt_seq_new(seq + (lv->leftltrlen)
                               - (lv->radius),
                           (GtUword) (2 * lv->radius + 1),
                           a);

  seq_rev     = gt_seq_new(rev_seq + (lv->rightltrlen)
                                   - (lv->radius),
                           (GtUword) (2 * lv->radius + 1),
                           a);

  for (j = 0; j < gt_bioseq_number_of_sequences(lv->trna_lib); j++)
  {
    GtSeq *trna_seq, *trna_from3;
    char *trna_from3_full;
    GtUword trna_seqlen;

    trna_seq = gt_bioseq_get_seq(lv->trna_lib, j);
    trna_seqlen = gt_seq_length(trna_seq);

    trna_from3_full = gt_calloc((size_t) trna_seqlen, sizeof (char));
    memcpy(trna_from3_full, gt_seq_get_orig(trna_seq),
           sizeof (char) * trna_seqlen);
    (void) gt_reverse_complement(trna_from3_full, trna_seqlen, err);
    trna_from3 = gt_seq_new_own(trna_from3_full, trna_seqlen, a);

    ali = gt_swalign(seq_forward, trna_from3, sf);
    gt_pbs_add_hit(lv, results->hits, ali, trna_seqlen,
                   gt_seq_get_description(trna_seq), GT_STRAND_FORWARD,
                   results);
    gt_alignment_delete(ali);

    ali = gt_swalign(seq_rev, trna_from3, sf);
    gt_pbs_add_hit(lv, results->hits, ali, trna_seqlen,
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

static void gt_pbs_results_delete(GtPBSResults *results)
{
    GtUword i;
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

static void pbs_attach_results_to_gff3(GtLTRdigestPBSVisitor *lv,
                                       GtPBSResults *results,
                                       GtFeatureNode *mainnode,
                                       GtStrand *canonical_strand)
{
  GtRange pbs_range;
  GtUword i = 0;
  GtGenomeNode *gf;
  char buffer[BUFSIZ];
  GtPBSHit* hit = gt_pbs_results_get_ranked_hit(results, i++);

  gt_log_log("attaching to %p: canonical %c this is %c", mainnode,
           GT_STRAND_CHARS[*canonical_strand], GT_STRAND_CHARS[hit->strand]);
  if (*canonical_strand == GT_STRAND_UNKNOWN)
    *canonical_strand = hit->strand;
  else {
    /* find best-scoring PBS on the given canonical strand */
    while (hit->strand != *canonical_strand
             && i < gt_pbs_results_get_number_of_hits(results)) {
      gt_log_log("dropping PBS because of nonconsistent strand: %s\n",
                 gt_feature_node_get_attribute(mainnode, "ID"));
      hit = gt_pbs_results_get_ranked_hit(results, i++);
    }
    /* if there is none, do not report a PBS */
    if (hit->strand != *canonical_strand)
      return;
  }
  gt_log_log("final strand %c", GT_STRAND_CHARS[hit->strand]);
  pbs_range = gt_pbs_hit_get_coords(lv, hit);
  pbs_range.start++; pbs_range.end++;  /* GFF3 is 1-based */
  gf = gt_feature_node_new(gt_genome_node_get_seqid((GtGenomeNode*) mainnode),
                           gt_ft_primer_binding_site,
                           pbs_range.start,
                           pbs_range.end,
                           hit->strand);
  gt_feature_node_set_source((GtFeatureNode*) gf, lv->tag);
  gt_feature_node_set_score((GtFeatureNode*) gf,
                            (float) gt_pbs_hit_get_score(hit));
  if (gt_pbs_hit_get_trna(hit) != NULL) {
    gt_feature_node_add_attribute((GtFeatureNode*) gf, "trna",
                                   gt_pbs_hit_get_trna(hit));
  }
  gt_feature_node_set_strand(mainnode, gt_pbs_hit_get_strand(hit));
  (void) snprintf(buffer, BUFSIZ-1, ""GT_WU"", gt_pbs_hit_get_tstart(hit));
  gt_feature_node_add_attribute((GtFeatureNode*) gf, "trnaoffset", buffer);
  (void) snprintf(buffer, BUFSIZ-1, ""GT_WU"", gt_pbs_hit_get_offset(hit));
  gt_feature_node_add_attribute((GtFeatureNode*) gf, "pbsoffset", buffer);
  (void) snprintf(buffer, BUFSIZ-1, ""GT_WU"", gt_pbs_hit_get_edist(hit));
  gt_feature_node_add_attribute((GtFeatureNode*) gf, "edist", buffer);

  gt_feature_node_add_child(mainnode, (GtFeatureNode*) gf);
}

const GtNodeVisitorClass* gt_ltrdigest_pbs_visitor_class(void);

#define gt_ltrdigest_pbs_visitor_cast(GV)\
        gt_node_visitor_cast(gt_ltrdigest_pbs_visitor_class(), GV)

static int gt_ltrdigest_pbs_visitor_feature_node(GtNodeVisitor *nv,
                                                 GtFeatureNode *fn,
                                                 GtError *err)
{
  GT_UNUSED GtLTRdigestPBSVisitor *lv;
  GtFeatureNodeIterator *fni;
  GtRange leftltrrng, rightltrrng;
  bool seen_left = false;
  GtFeatureNode *curnode = NULL;
  int had_err = 0;
  lv = gt_ltrdigest_pbs_visitor_cast(nv);
  gt_assert(lv);
  gt_error_check(err);

  /* traverse annotation subgraph and find LTR element */
  fni = gt_feature_node_iterator_new(fn);
  while (!had_err && (curnode = gt_feature_node_iterator_next(fni))) {
    if (strcmp(gt_feature_node_get_type(curnode),
               gt_ft_LTR_retrotransposon) == 0) {
      lv->ltr_retrotrans = curnode;
    }
    if (strcmp(gt_feature_node_get_type(curnode),
               gt_ft_long_terminal_repeat) == 0) {
      if (seen_left) {
        rightltrrng = gt_genome_node_get_range((GtGenomeNode*) curnode);
        lv->rightltrlen = gt_range_length(&rightltrrng);
        lv->rightLTR_5 = rightltrrng.start - 1;
      } else {
        leftltrrng = gt_genome_node_get_range((GtGenomeNode*) curnode);
        lv->leftltrlen = gt_range_length(&leftltrrng);
        lv->leftLTR_3 = leftltrrng.end - 1;
        seen_left = true;
      }
    }
  }
  gt_feature_node_iterator_delete(fni);

  if (!had_err && lv->ltr_retrotrans != NULL) {
    char *rev_seq;
    GtRange rng;
    GtStrand canonical_strand = gt_feature_node_get_strand(lv->ltr_retrotrans);
    GtUword seqlen;
    GtStr *seq = gt_str_new();
    rng = gt_genome_node_get_range((GtGenomeNode*) lv->ltr_retrotrans);
    seqlen = gt_range_length(&rng);

    had_err = gt_extract_feature_sequence(seq,
                                          (GtGenomeNode*) lv->ltr_retrotrans,
                                          gt_symbol(gt_ft_LTR_retrotransposon),
                                          false, NULL, NULL, lv->rmap, err);

    if (!had_err) {
      GtPBSResults *pbs_results = NULL;
      rev_seq = gt_malloc((size_t) (seqlen * sizeof (char)));
      strncpy(rev_seq, gt_str_get(seq), (size_t) seqlen * sizeof (char));
      (void) gt_reverse_complement(rev_seq, seqlen, NULL);

      pbs_results = gt_pbs_find(lv, gt_str_get(seq), (const char*) rev_seq,
                                err);
       if (gt_pbs_results_get_number_of_hits(pbs_results) > 0)
       {
        pbs_attach_results_to_gff3(lv, pbs_results, lv->ltr_retrotrans,
                                   &canonical_strand);
       }
       gt_pbs_results_delete(pbs_results);

      gt_free(rev_seq);
    }
    gt_str_delete(seq);
  }

  return had_err;
}

static void gt_ltrdigest_pbs_visitor_free(GtNodeVisitor *nv)
{
  GT_UNUSED GtLTRdigestPBSVisitor *lv;
  if (!nv) return;
  lv = gt_ltrdigest_pbs_visitor_cast(nv);
  gt_str_delete(lv->tag);
}

const GtNodeVisitorClass* gt_ltrdigest_pbs_visitor_class(void)
{
  static const GtNodeVisitorClass *nvc = NULL;
  if (!nvc) {
    nvc = gt_node_visitor_class_new(sizeof (GtLTRdigestPBSVisitor),
                                   gt_ltrdigest_pbs_visitor_free,
                                   NULL,
                                   gt_ltrdigest_pbs_visitor_feature_node,
                                   NULL,
                                   NULL,
                                   NULL);
  }
  return nvc;
}

GtNodeVisitor* gt_ltrdigest_pbs_visitor_new(GtRegionMapping *rmap,
                                            unsigned int radius,
                                            unsigned int max_edist,
                                            GtRange alilen,
                                            GtRange offsetlen,
                                            GtRange trnaoffsetlen,
                                            int ali_score_match,
                                            int ali_score_mismatch,
                                            int ali_score_insertion,
                                            int ali_score_deletion,
                                            GtBioseq *trna_lib,
                                            GT_UNUSED GtError *err)
{
  GtNodeVisitor *nv = NULL;
  GtLTRdigestPBSVisitor *lv;
  gt_assert(rmap && trna_lib);
  nv = gt_node_visitor_create(gt_ltrdigest_pbs_visitor_class());
  lv = gt_ltrdigest_pbs_visitor_cast(nv);
  lv->tag = gt_str_new_cstr(GT_LTRDIGEST_TAG);
  lv->rmap = rmap;
  lv->radius = radius;
  lv->max_edist = max_edist;
  lv->alilen = alilen;
  lv->offsetlen = offsetlen;
  lv->trnaoffsetlen = trnaoffsetlen;
  lv->ali_score_match = ali_score_match;
  lv->ali_score_mismatch = ali_score_mismatch;
  lv->ali_score_insertion = ali_score_insertion;
  lv->ali_score_deletion = ali_score_deletion;
  lv->trna_lib = trna_lib;
  return nv;
}

int gt_ltrdigest_pbs_visitor_unit_test(GT_UNUSED GtError *err)
{
  int had_err = 0;
  #if 0
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
  gt_ensure(gt_file_exists(gt_str_get(tmpfilename)));

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
  gt_ensure(gt_bioseq_number_of_sequences(o.trna_lib) == 2UL);

  element.leftLTR_5 = 20UL;
  element.leftLTR_3 = 119UL;
  element.rightLTR_5 = 520UL;
  element.rightLTR_3 = 619UL;

  /* setup sequences */
  seq     = gt_malloc((size_t) 600 * sizeof (char));
  rev_seq = gt_malloc((size_t) 600 * sizeof (char));
  memcpy(seq,     fullseq + 20, (size_t) 600);
  memcpy(rev_seq, fullseq + 20, (size_t) 600);
  gt_ensure(
            !gt_reverse_complement(rev_seq, (GtUword) 600, NULL));

  /* try to find PBS in sequences */
  res = gt_pbs_find(seq, rev_seq, &element, &o, err);
  gt_ensure(res != NULL);
  gt_ensure(gt_pbs_results_get_number_of_hits(res) == 2UL);

  /* check first hit on forward strand */
  hit = gt_pbs_results_get_ranked_hit(res, 0UL);
  gt_ensure(hit != NULL);
  gt_ensure(gt_pbs_hit_get_alignment_length(hit) == 17UL);
  gt_ensure(gt_pbs_hit_get_edist(hit) == 0UL);
  gt_ensure(gt_pbs_hit_get_offset(hit) == 0UL);
  gt_ensure(gt_pbs_hit_get_tstart(hit) == 3UL);
  gt_ensure(strcmp(gt_pbs_hit_get_trna(hit), "test1") == 0);
  rng = gt_pbs_hit_get_coords(hit);
  gt_ensure(rng.start == 120UL);
  gt_ensure(rng.end == 136UL);
  score1 = gt_pbs_hit_get_score(hit);
  gt_ensure(gt_pbs_hit_get_strand(hit) == GT_STRAND_FORWARD);
  memset(tmp, 0, BUFSIZ-1);
  memcpy(tmp, fullseq + (rng.start * sizeof (char)),
         (size_t) ((rng.end - rng.start + 1) * sizeof (char)));
  gt_ensure(strcmp(tmp, "acatactaggatgctag" ) == 0);

  /* check second hit on reverse strand */
  hit = gt_pbs_results_get_ranked_hit(res, 1UL);
  gt_ensure(hit != NULL);
  gt_ensure(gt_pbs_hit_get_alignment_length(hit) == 14UL);
  gt_ensure(gt_pbs_hit_get_edist(hit) == 1UL);
  gt_ensure(gt_pbs_hit_get_offset(hit) == 0UL);
  gt_ensure(gt_pbs_hit_get_tstart(hit) == 6UL);
  gt_ensure(strcmp(gt_pbs_hit_get_trna(hit), "test2") == 0);
  rng = gt_pbs_hit_get_coords(hit);
  gt_ensure(rng.start == 506UL);
  gt_ensure(rng.end == 519UL);
  score2 = gt_pbs_hit_get_score(hit);
  gt_ensure(gt_double_compare(score1, score2) > 0);
  gt_ensure(gt_pbs_hit_get_strand(hit) == GT_STRAND_REVERSE);
  memset(tmp, 0, BUFSIZ-1);
  memcpy(tmp, fullseq + (rng.start * sizeof (char)),
         (size_t) ((rng.end - rng.start + 1) * sizeof (char)));
  gt_ensure(strcmp(tmp, "gatcctaaggctac" ) == 0);

  /* clean up */
  gt_xremove(gt_str_get(tmpfilename));
  gt_ensure(!gt_file_exists(gt_str_get(tmpfilename)));
  gt_str_delete(tmpfilename);
  gt_bioseq_delete(o.trna_lib);
  gt_free(rev_seq);
  gt_free(seq);
  gt_pbs_results_delete(res);
#endif

  return had_err;
}
