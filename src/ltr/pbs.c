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

#include <assert.h>
#include <math.h>
#include "core/ma.h"
#include "core/dlist.h"
#include "core/mathsupport.h"
#include "core/minmax.h"
#include "core/xansi.h"
#include "core/array.h"
#include "core/bioseq.h"
#include "extended/reverse.h"
#include "extended/swalign.h"
#include "ltr/pbs.h"

GtScoreFunction* gt_dna_scorefunc_new(GtAlpha *a, int match, int mismatch,
                                      int insertion, int deletion)
{
  GtScoreMatrix *sm = gt_score_matrix_new(a);
  GtScoreFunction *sf = gt_score_function_new(sm, insertion, deletion);
  unsigned int m,n;

  for (m=0;m<5;m++)
  {
    for (n=0;n<5;n++)
    {
      gt_score_matrix_set_score(sm, m,n,(n==m ? match : mismatch));
    }
  }
  return sf;
}

double gt_pbs_score_func(unsigned long edist, unsigned long offset,
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

void gt_pbs_add_hit(GtDlist *hitlist, GtAlignment *ali, GtPBSOptions *o,
                    unsigned long trna_seqlen, const char *desc,
                    GtStrand strand)
{
  unsigned long dist;
  GtPBS_Hit *hit;
  GtRange urange, vrange;

  dist = gt_alignment_eval(ali);
  urange = gt_alignment_get_urange(ali);
  vrange = gt_alignment_get_vrange(ali);
  if (dist <= o->max_edist
        && abs(o->radius-urange.start) <= o->offsetlen.end
        && abs(o->radius-urange.start) >= o->offsetlen.start
        && abs(urange.end-urange.start+1) <= o->alilen.end
        && abs(urange.end-urange.start+1) >= o->alilen.start
        && vrange.start <= o->trnaoffsetlen.end
        && vrange.start >= o->trnaoffsetlen.start)
  {
    hit = gt_calloc(1, sizeof (GtPBS_Hit));
    hit->alilen  = abs(urange.end-urange.start)+1;
    hit->strand  = strand;
    hit->trna    = desc;
    hit->tstart  = vrange.start;
    hit->start   = urange.start;
    hit->end     = urange.end;
    hit->offset  = abs(o->radius-urange.start);
    hit->edist   = dist;
    hit->score   = gt_pbs_score_func(dist,
                                     abs(o->radius-urange.start),
                                     urange.end-urange.start+1,
                                     trna_seqlen,
                                     vrange.start);
    gt_dlist_add(hitlist, hit);
  }
}

int gt_pbs_hit_compare(const void *h1, const void *h2)
{
  GtPBS_Hit *hp1 = (GtPBS_Hit*) h1;
  GtPBS_Hit *hp2 = (GtPBS_Hit*) h2;

  return (gt_double_compare(hp2->score,hp1->score));
}

void gt_pbs_find(const char *seq,
                 const char *rev_seq,
                 GtLTRElement *element,
                 GtPBSResults *results,
                 GtPBSOptions *o,
                 GtError *err)
{
  GtSeq *seq_forward, *seq_rev;
  unsigned long j;
  GtAlignment *ali;
  GtAlpha *a = (GtAlpha*) gt_alpha_new_dna();
  GtScoreFunction *sf = gt_dna_scorefunc_new(a,
                                             o->ali_score_match,
                                             o->ali_score_mismatch,
                                             o->ali_score_insertion,
                                             o->ali_score_deletion);

  assert(seq && rev_seq && sf && a && element && results);

  results->hits_fwd = gt_dlist_new(gt_pbs_hit_compare);
  results->hits_rev = gt_dlist_new(gt_pbs_hit_compare);
  results->best_hit = NULL;

  seq_forward = gt_seq_new(seq +
                             element->leftLTR_3-element->leftLTR_5-o->radius,
                           2*o->radius,
                           a);

  seq_rev     = gt_seq_new(rev_seq +
                             element->rightLTR_3-element->rightLTR_5-o->radius,
                           2*o->radius,
                           a);

  for (j=0;j<gt_bioseq_number_of_sequences(o->trna_lib);j++)
  {
    GtSeq *trna_seq, *trna_from3;
    char *trna_from3_full;
    unsigned long trna_seqlen;

    trna_seq = gt_bioseq_get_seq(o->trna_lib, j);
    trna_seqlen = gt_seq_length(trna_seq);

    trna_from3_full = gt_calloc(trna_seqlen, sizeof (char));
    memcpy(trna_from3_full, gt_seq_get_orig(trna_seq),
           sizeof (char)*trna_seqlen);
    (void) gt_reverse_complement(trna_from3_full, trna_seqlen, err);
    trna_from3 = gt_seq_new_own(trna_from3_full, trna_seqlen, a);

    ali = gt_swalign(seq_forward, trna_from3, sf);
    gt_pbs_add_hit(results->hits_fwd,ali,o,trna_seqlen,
                   gt_seq_get_description(trna_seq), GT_STRAND_FORWARD);
    gt_alignment_delete(ali);

    ali = gt_swalign(seq_rev, trna_from3, sf);
    gt_pbs_add_hit(results->hits_rev,ali,o,trna_seqlen,
                gt_seq_get_description(trna_seq), GT_STRAND_REVERSE);
    gt_alignment_delete(ali);

    gt_seq_delete(trna_from3);
  }
  gt_seq_delete(seq_forward);
  gt_seq_delete(seq_rev);
  gt_score_function_delete(sf);
  gt_alpha_delete(a);

  if (gt_dlist_size(results->hits_fwd) > 0)
  {
    GtDlistelem *delem;
    delem = gt_dlist_first(results->hits_fwd);
    results->best_hit = (GtPBS_Hit*) gt_dlistelem_get_data(delem);
  }
  if (gt_dlist_size(results->hits_rev) > 0)
  {
    GtDlistelem *delem;
    GtPBS_Hit *tmp;
    delem = gt_dlist_first(results->hits_rev);
    tmp = (GtPBS_Hit*) gt_dlistelem_get_data(delem);
    if (!results->best_hit
          || gt_double_compare(tmp->score, results->best_hit->score) > 0)
      results->best_hit = tmp;
  }
}

void gt_pbs_clear_results(GtPBSResults *results)
{
    GtDlistelem *delem;

    if (!results) return;

    if (results->hits_fwd)
    {
      for (delem = gt_dlist_first(results->hits_fwd);
          delem;
          delem = gt_dlistelem_next(delem))
      {
        GtPBS_Hit *hit = (GtPBS_Hit*) gt_dlistelem_get_data(delem);
        gt_free(hit);
      }
      gt_dlist_delete(results->hits_fwd);
    }
    if (results->hits_rev)
    {
      for (delem = gt_dlist_first(results->hits_rev);
          delem;
          delem = gt_dlistelem_next(delem))
      {
        GtPBS_Hit *hit = (GtPBS_Hit*) gt_dlistelem_get_data(delem);
        gt_free(hit);
      }
      gt_dlist_delete(results->hits_rev);
    }
}
