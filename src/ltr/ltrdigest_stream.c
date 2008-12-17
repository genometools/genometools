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

#include <string.h>
#include "core/codon.h"
#include "core/log.h"
#include "core/ma.h"
#include "core/mathsupport.h"
#include "core/range.h"
#include "core/str.h"
#include "core/unused_api.h"
#include "extended/node_stream_rep.h"
#include "extended/feature_node.h"
#include "extended/feature_node_iterator_api.h"
#include "extended/reverse.h"
#include "ltr/ltrdigest_def.h"
#include "ltr/ltrdigest_stream.h"
#include "ltr/ltr_visitor.h"

struct GtLTRdigestStream {
  const GtNodeStream parent_instance;
  GtNodeStream *in_stream;
  GtBioseq *bioseq;
  GtPBSOptions *pbs_opts;
  GtPPTOptions *ppt_opts;
#ifdef HAVE_HMMER
  GtPdomOptions *pdom_opts;
#endif
  GtLTRVisitor *lv;
  GtStr *ltrdigest_tag;
  int tests_to_run;
  GtLTRElement element;
};

#define gt_ltrdigest_stream_cast(GS)\
        gt_node_stream_cast(gt_ltrdigest_stream_class(), GS)

static inline void convert_frame_position(GtRange *rng, int frame)
{
  rng->start = (rng->start - 1)*GT_CODON_LENGTH + frame;
  rng->end   = (rng->end   - 1)*GT_CODON_LENGTH + frame;
}

#ifdef HAVE_HMMER
static int pdom_hit_attach_gff3(struct plan7_s *model, GtPdomHit *hit,
                                void *data, GT_UNUSED GtError *err)
{
  unsigned long i;
  GtRange rng;
  GtLTRdigestStream  *ls = (GtLTRdigestStream *) data;
  GtStrand strand = gt_pdom_hit_get_strand(hit);
  GtArray *best_chain = gt_pdom_hit_get_best_chain(hit);

  /* do not use the hits on the non-predicted strand
      -- maybe identify nested elements ? */
  if (strand != gt_feature_node_get_strand(ls->element.mainnode))
    return 0;

  for (i=0;i<gt_array_size(best_chain);i++)
  {
    GtGenomeNode *gf;
    struct hit_s *singlehit = *(struct hit_s **) gt_array_get(best_chain, i);
    GtPhase frame = gt_phase_get(singlehit->name[0]);
    rng.start = singlehit->sqfrom; rng.end = singlehit->sqto;
    convert_frame_position(&rng, frame);
    gt_ltrelement_offset2pos(&ls->element, &rng, 0,
                             GT_OFFSET_BEGIN_LEFT_LTR,
                             strand);
    gf = gt_feature_node_new(gt_genome_node_get_seqid((GtGenomeNode*)
                                                      ls->element.mainnode),
                             GT_PDOM_TYPE,
                             rng.start,
                             rng.end,
                             strand);
    gt_feature_node_set_source((GtFeatureNode*) gf, ls->ltrdigest_tag);
    gt_feature_node_set_phase((GtFeatureNode*) gf, frame);
    gt_feature_node_add_attribute((GtFeatureNode*) gf,"pfamname", model->name);
    gt_feature_node_add_attribute((GtFeatureNode*) gf,"pfamid", model->acc);
    gt_feature_node_add_child(ls->element.mainnode, (GtFeatureNode*) gf);
  }
  return 0;
}
#endif

static void pbs_attach_results_to_gff3(GtPBSResults *results,
                                       GtLTRElement *element,
                                       GtStrand *canonical_strand,
                                       GtStr *tag)
{
  GtRange pbs_range;
  GtGenomeNode *gf;
  unsigned long i = 0;
  char buffer[BUFSIZ];
  GtPBSHit* hit = gt_pbs_results_get_ranked_hit(results, i++);
  if (*canonical_strand == GT_STRAND_UNKNOWN)
    *canonical_strand = gt_pbs_hit_get_strand(hit);
  else
  {
    /* do we have to satisfy a strand constraint?
     * then find best-scoring PBS on the given canonical strand */
    while (gt_pbs_hit_get_strand(hit) != *canonical_strand
             && i < gt_pbs_results_get_number_of_hits(results))
    {
      gt_log_log("dropping PBS because of nonconsistent strand: %s\n",
                 gt_feature_node_get_attribute(element->mainnode, "ID"));
      hit = gt_pbs_results_get_ranked_hit(results, i++);
    }
    /* if there is none, do not report a PBS */
    if (gt_pbs_hit_get_strand(hit) != *canonical_strand)
      return;
  }
  pbs_range = gt_pbs_hit_get_coords(hit);
  gf = gt_feature_node_new(gt_genome_node_get_seqid((GtGenomeNode*)
                                                    element->mainnode),
                           GT_PBS_TYPE,
                           pbs_range.start,
                           pbs_range.end,
                           gt_pbs_hit_get_strand(hit));
  gt_feature_node_set_source((GtFeatureNode*) gf, tag);
  gt_feature_node_set_score((GtFeatureNode*) gf, gt_pbs_hit_get_score(hit));
  gt_feature_node_add_attribute((GtFeatureNode*) gf,"trna",
                                 gt_pbs_hit_get_trna(hit));
  gt_feature_node_set_strand(element->mainnode, gt_pbs_hit_get_strand(hit));
  (void) snprintf(buffer, BUFSIZ-1, "%lu", gt_pbs_hit_get_tstart(hit));
  gt_feature_node_add_attribute((GtFeatureNode*) gf,"trnaoffset", buffer);
  (void) snprintf(buffer, BUFSIZ-1, "%lu", gt_pbs_hit_get_offset(hit));
  gt_feature_node_add_attribute((GtFeatureNode*) gf,"pbsoffset", buffer);
  (void) snprintf(buffer, BUFSIZ-1, "%lu", gt_pbs_hit_get_edist(hit));
  gt_feature_node_add_attribute((GtFeatureNode*) gf,"edist", buffer);
  gt_feature_node_add_child(element->mainnode, (GtFeatureNode*) gf);
}

static void ppt_attach_results_to_gff3(GtPPTResults *results,
                                       GtLTRElement *element,
                                       GtStrand *canonical_strand,
                                       GtStr *tag)
{
  GtRange ppt_range;
  unsigned long i = 0;
  GtGenomeNode *gf;
  GtPPTHit* hit = gt_ppt_results_get_ranked_hit(results, i++);
  if (*canonical_strand == GT_STRAND_UNKNOWN)
    *canonical_strand = gt_ppt_hit_get_strand(hit);
  else
  {
    /* find best-scoring PPT on the given canonical strand */
    while (gt_ppt_hit_get_strand(hit) != *canonical_strand
             && i < gt_ppt_results_get_number_of_hits(results))
    {
      gt_log_log("dropping PPT because of nonconsistent strand: %s\n",
                 gt_feature_node_get_attribute(element->mainnode, "ID"));
      hit = gt_ppt_results_get_ranked_hit(results, i++);
    }
    /* if there is none, do not report a PPT */
    if (gt_ppt_hit_get_strand(hit) != *canonical_strand)
      return;
  }
  ppt_range = gt_ppt_hit_get_coords(hit);

  gf = gt_feature_node_new(gt_genome_node_get_seqid((GtGenomeNode*)
                                                    element->mainnode),
                           GT_PPT_TYPE,
                           ppt_range.start,
                           ppt_range.end,
                           gt_ppt_hit_get_strand(hit));
  gt_feature_node_set_source((GtFeatureNode*) gf, tag);
  gt_feature_node_set_strand(element->mainnode, gt_ppt_hit_get_strand(hit));
  gt_feature_node_add_child(element->mainnode, (GtFeatureNode*) gf);
}

static void run_ltrdigest(GtLTRElement *element, GtSeq *seq,
                          GtLTRdigestStream *ls, GtError *err)
{
  char *rev_seq;
  const char *base_seq = gt_seq_get_orig(seq)+element->leftLTR_5;
  unsigned long seqlen = gt_ltrelement_length(element);
  GtStrand canonical_strand = GT_STRAND_UNKNOWN;

  /* create reverse strand sequence */
  rev_seq = gt_calloc(seqlen+1, sizeof (char));
  memcpy(rev_seq, base_seq, sizeof (char) * seqlen);
  rev_seq[seqlen] = '\0';
  (void) gt_reverse_complement(rev_seq, seqlen, err);

#ifdef HAVE_HMMER
  /* Protein domain finding
   * ----------------------*/
  if (ls->tests_to_run & GT_LTRDIGEST_RUN_PDOM)
  {
    GtPdomResults *pdom_results = NULL;
    pdom_results = gt_pdom_find((const char*) base_seq, (const char*) rev_seq,
                                element, ls->pdom_opts);
    if (!gt_pdom_results_empty(pdom_results))
    {
      /* determine most likely strand from protein domain results */
      if (gt_double_compare(
                     gt_pdom_results_get_combined_evalue_fwd(pdom_results),
                     gt_pdom_results_get_combined_evalue_rev(pdom_results)) < 0)
        canonical_strand = GT_STRAND_FORWARD;
      else
        canonical_strand = GT_STRAND_REVERSE;
      gt_feature_node_set_strand(ls->element.mainnode,
                                    canonical_strand);
      /* create nodes for protein match annotations */
      (void) gt_pdom_results_foreach_domain_hit(pdom_results,
                                                pdom_hit_attach_gff3,
                                                ls,
                                                err);
    }
    gt_pdom_results_delete(pdom_results);
  }
#endif

  /* PPT finding
   * -----------*/
  if (ls->tests_to_run & GT_LTRDIGEST_RUN_PPT)
  {
    GtPPTResults *ppt_results = NULL;
    ppt_results = gt_ppt_find((const char*) base_seq, (const char*) rev_seq,
                            element, ls->ppt_opts);
    if (gt_ppt_results_get_number_of_hits(ppt_results) > 0)
    {
      ppt_attach_results_to_gff3(ppt_results, element, &canonical_strand,
                                 ls->ltrdigest_tag);
    }
    gt_ppt_results_delete(ppt_results);
  }

  /* PBS finding
   * ----------- */
  if (ls->tests_to_run & GT_LTRDIGEST_RUN_PBS)
  {
    GtPBSResults *pbs_results = NULL;
    pbs_results = gt_pbs_find((const char*) base_seq, (const char*) rev_seq,
                              element, ls->pbs_opts, err);
     if (gt_pbs_results_get_number_of_hits(pbs_results) > 0)
     {
      pbs_attach_results_to_gff3(pbs_results, element, &canonical_strand,
                                 ls->ltrdigest_tag);
     }
     gt_pbs_results_delete(pbs_results);
  }

  gt_free(rev_seq);
}

static int gt_ltrdigest_stream_next(GtNodeStream *gs, GtGenomeNode **gn,
                                    GtError *e)
{
  GtLTRdigestStream *ls;
  GtFeatureNode *fn;
  int had_err;

  gt_error_check(e);
  ls = gt_ltrdigest_stream_cast(gs);

  /* initialize this element */
  memset(&ls->element, 0, sizeof (GtLTRElement));

  /* get annotations from parser */
  had_err = gt_node_stream_next(ls->in_stream, gn, e);
  if (!had_err && *gn)
  {
    GtFeatureNodeIterator *gni;
    GtFeatureNode *mygn;

   /* only process feature nodes */
   if (!(fn = gt_feature_node_try_cast(*gn)))
     return 0;

    /* fill LTRElement structure from GFF3 subgraph */
    gni = gt_feature_node_iterator_new(fn);
    for (mygn = fn; mygn; mygn = gt_feature_node_iterator_next(gni))
      (void) gt_genome_node_accept((GtGenomeNode*) mygn,
                                   (GtNodeVisitor*) ls->lv,
                                   e);
    gt_feature_node_iterator_delete(gni);
  }

  if (ls->element.mainnode)
  {
    unsigned long seqid;
    const char *sreg;
    GtSeq *seq;
    GtRange elemrng;

    sreg = gt_str_get(gt_genome_node_get_seqid((GtGenomeNode*)
                                               ls->element.mainnode));
    /* XXX: this may work for LTRharvest, but not everywhere!
     * fix mapping (like in other tools)! */
    (void) sscanf(sreg,"seq%lu", &seqid);
    seq = gt_bioseq_get_seq(ls->bioseq, seqid);

    elemrng = gt_genome_node_get_range((GtGenomeNode*) ls->element.mainnode);
    if (elemrng.end <= gt_seq_length(seq))
      /* run LTRdigest core routine */
      run_ltrdigest(&ls->element, seq, ls, e);
    else
    {
      /* do not process elements whose positions exceed sequence boundaries
       (obviously annotation and sequence do not match!) */
      gt_error_set(e, "Element '%s' exceeds sequence boundaries! (%lu > %lu)",
        gt_feature_node_get_attribute(ls->element.mainnode, "ID"),
        elemrng.end, gt_seq_length(seq));
      had_err = -1;
    }
  }
  if (had_err) {
    gt_genome_node_delete(*gn);
    *gn = NULL;
  }
  return had_err;
}

static void gt_ltrdigest_stream_free(GtNodeStream *gs)
{
  GtLTRdigestStream *ls = gt_ltrdigest_stream_cast(gs);
  gt_node_visitor_delete((GtNodeVisitor*) ls->lv);
  gt_str_delete(ls->ltrdigest_tag);
  gt_node_stream_delete(ls->in_stream);
}

const GtNodeStreamClass* gt_ltrdigest_stream_class(void)
{
  static const GtNodeStreamClass *gsc;
  if (!gsc)
    gsc = gt_node_stream_class_new(sizeof (GtLTRdigestStream),
                                   gt_ltrdigest_stream_free,
                                   gt_ltrdigest_stream_next );
  return gsc;
}

GtNodeStream* gt_ltrdigest_stream_new(GtNodeStream *in_stream,
                                      int tests_to_run,
                                      GtBioseq *bioseq,
                                      GtPBSOptions *pbs_opts,
                                      GtPPTOptions *ppt_opts
#ifdef HAVE_HMMER
           /* HMMER only */ ,GtPdomOptions *pdom_opts
#endif
           )
{
  GtNodeStream *gs;
  GtLTRdigestStream *ls;
  gs = gt_node_stream_create(gt_ltrdigest_stream_class(), true);
  ls = gt_ltrdigest_stream_cast(gs);
  ls->in_stream = gt_node_stream_ref(in_stream);
  ls->ppt_opts = ppt_opts;
  ls->pbs_opts = pbs_opts;
#ifdef HAVE_HMMER
  ls->pdom_opts = pdom_opts;
#endif
  ls->tests_to_run = tests_to_run;
  ls->bioseq = bioseq;
  ls->ltrdigest_tag = gt_str_new_cstr(GT_LTRDIGEST_TAG);
  ls->lv = (GtLTRVisitor*) gt_ltr_visitor_new(&ls->element);
  return gs;
}
