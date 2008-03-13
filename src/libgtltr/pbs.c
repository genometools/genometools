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
#include "libgtcore/dlist.h"
#include "libgtcore/minmax.h"
#include "libgtcore/xansi.h"
#include "libgtcore/array.h"
#include "libgtcore/bioseq.h"
#include "libgtext/reverse.h"
#include "libgtext/swalign.h"
#include "libgtltr/pbs.h"

ScoreFunction* dna_scorefunc_new(Alpha *a, int match, int mismatch,
                                 int insertion, int deletion)
{
  ScoreMatrix *sm = score_matrix_new(a);
  ScoreFunction *sf = scorefunction_new(sm, insertion, deletion);
  unsigned int m,n;

  for(m=0;m<5;m++)
  {
    for(n=0;n<5;n++)
    {
      score_matrix_set_score(sm, m,n,(n==m ? match : mismatch));
    }
  }
  return sf;
}

double pbs_score_func(unsigned long edist, unsigned long offset,
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

void pbs_add_hit(Dlist *hitlist, Alignment *ali, PBSOptions *o,
                 unsigned long trna_seqlen, const char *desc,
                 Strand strand)
{
  unsigned long dist;
  PBS_Hit *hit;
  Range urange, vrange;

  dist = alignment_eval(ali);
  urange = alignment_get_urange(ali);
  vrange = alignment_get_vrange(ali);
  if (dist <= o->max_edist
        && abs(o->radius-urange.start) <= o->max_offset
        && abs(urange.end-urange.start+1) >= o->ali_min_len
        && vrange.start <= o->max_offset_trna)
  {
    hit = ma_malloc(sizeof (PBS_Hit));
    hit->alilen  = abs(urange.end-urange.start)+1;
    hit->strand  = strand;
    hit->trna    = desc;
    hit->tstart  = vrange.start;
    hit->start   = urange.start;
    hit->end     = urange.end;
    hit->offset  = abs(o->radius-urange.start);
    hit->edist   = dist;
    hit->score   = pbs_score_func(dist,
                                  abs(o->radius-urange.start),
                                  urange.end-urange.start+1,
                                  trna_seqlen,
                                  vrange.start);
    dlist_add(hitlist, hit);
  }
}

int pbs_hit_compare(const void *h1, const void *h2)
{
  PBS_Hit *hp1 = (PBS_Hit*) h1;
  PBS_Hit *hp2 = (PBS_Hit*) h2;

  if(hp1->score == hp2->score)
    return 0;
  else return (hp1->score > hp2->score ? -1 : 1);
}

void pbs_find(const char *seq,
              const char *rev_seq,
              LTRElement *element,
              PBSResults *results,
              PBSOptions *o,
              Error *err)
{
  Seq *seq_forward, *seq_rev;
  unsigned long j;
  Alignment *ali;
  Alpha *a = (Alpha*) alpha_new_dna();
  ScoreFunction *sf = dna_scorefunc_new(a,
                                        o->ali_score_match,
                                        o->ali_score_mismatch,
                                        o->ali_score_insertion,
                                        o->ali_score_deletion);

  assert(seq && rev_seq && sf && a && element && results);

  results->hits_fwd = dlist_new(pbs_hit_compare);
  results->hits_rev = dlist_new(pbs_hit_compare);
  results->best_hit = NULL;

  seq_forward = seq_new(seq +
                          element->leftLTR_3-element->leftLTR_5-o->radius,
                        2*o->radius,
                        a);

  seq_rev     = seq_new(rev_seq +
                          element->rightLTR_3-element->rightLTR_5-o->radius,
                        2*o->radius,
                        a);

  for(j=0;j<bioseq_number_of_sequences(o->trna_lib);j++)
  {
    Seq *trna_seq, *trna_from3;
    char *trna_from3_full;
    unsigned long trna_seqlen;

    trna_seq = bioseq_get_seq(o->trna_lib, j);
    trna_seqlen = seq_length(trna_seq);

    trna_from3_full = ma_malloc(sizeof (char)*trna_seqlen);
    memcpy(trna_from3_full, seq_get_orig(trna_seq), sizeof(char)*trna_seqlen);
    reverse_complement(trna_from3_full, trna_seqlen, err);
    trna_from3 = seq_new_own(trna_from3_full, trna_seqlen, a);

    ali = swalign(seq_forward, trna_from3, sf);
    pbs_add_hit(results->hits_fwd,ali,o,trna_seqlen,
                seq_get_description(trna_seq), STRAND_FORWARD);
    alignment_delete(ali);

    ali = swalign(seq_rev, trna_from3, sf);
    pbs_add_hit(results->hits_rev,ali,o,trna_seqlen,
                seq_get_description(trna_seq), STRAND_REVERSE);
    alignment_delete(ali);

    seq_delete(trna_from3);
  }
  seq_delete(seq_forward);
  seq_delete(seq_rev);
  scorefunction_delete(sf);
  alpha_delete(a);

  if (dlist_size(results->hits_fwd) > 0)
  {
    Dlistelem *delem;
    delem = dlist_first(results->hits_fwd);
    results->best_hit = (PBS_Hit*) dlistelem_get_data(delem);
  }
  if (dlist_size(results->hits_rev) > 0)
  {
    Dlistelem *delem;
    PBS_Hit *tmp;
    delem = dlist_first(results->hits_rev);
    tmp = (PBS_Hit*) dlistelem_get_data(delem);
    if (!results->best_hit || tmp->score > results->best_hit->score)
      results->best_hit = tmp;
  }
}

void pbs_clear_results(PBSResults *results)
{
    Dlistelem *delem;
    for(delem = dlist_first(results->hits_fwd);
        delem;
        delem = dlistelem_next(delem))
    {
      PBS_Hit *hit = (PBS_Hit*) dlistelem_get_data(delem);
      ma_free(hit);
    }
    dlist_delete(results->hits_fwd);
    for(delem = dlist_first(results->hits_rev);
        delem;
        delem = dlistelem_next(delem))
    {
      PBS_Hit *hit = (PBS_Hit*) dlistelem_get_data(delem);
      ma_free(hit);
    }
    dlist_delete(results->hits_rev);
}
