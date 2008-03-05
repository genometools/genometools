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
#include "libgtltr/pbs.h"
#include "libgtcore/bioseq.h"
#include "libgtext/reverse.h"
#include "libgtext/swalign.h"

static ScoreFunction* dna_scorefunc_new(Alpha *a, int match, int mismatch,
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

static double pbs_score_func(unsigned long edist, unsigned long offset,
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

static void pbs_add_hit(Dlist *hitlist, Alignment *ali, LTRharvestoptions *lo,
                        unsigned long trna_seqlen, const char *desc,
                        Strand strand, unsigned long seqoffset)
{
  unsigned long dist;
  PBS_Hit *hit;
  Range urange, vrange;

  dist = alignment_eval(ali);
  urange = alignment_get_urange(ali);
  vrange = alignment_get_vrange(ali);
  if (dist <= lo->pbs_maxedist
        && abs(lo->pbs_radius-urange.start) <= lo->pbs_maxoffset_5_ltr
        && abs(urange.end-urange.start+1) >= lo->pbs_aliminlen
        && vrange.start <= lo->pbs_maxoffset_trna)
  {
    hit = ma_malloc(sizeof (PBS_Hit));
    hit->strand  = strand;
    hit->trna    = desc;
    hit->alilen  = abs(urange.end-urange.start+1);
    hit->tstart  = vrange.start;
    if (strand == STRAND_FORWARD)
      hit->start   = seqoffset + urange.start;
    else if (strand == STRAND_REVERSE)
      hit->start   = seqoffset - urange.end;
    hit->end     = urange.end;
    hit->offset  = abs(lo->pbs_radius-urange.start);
    hit->edist   = dist;
    hit->score   = pbs_score_func(dist,
                                  abs(lo->pbs_radius-urange.start),
                                  urange.end-urange.start+1,
                                  trna_seqlen,
                                  vrange.start);
    dlist_add(hitlist, hit);
  }
}

static int pbs_hit_compare(const void *h1, const void *h2)
{
  PBS_Hit *hp1 = (PBS_Hit*) h1;
  PBS_Hit *hp2 = (PBS_Hit*) h2;

  if(hp1->score == hp2->score)
    return 0;
  else return (hp1->score > hp2->score ? -1 : 1);
}

PBS_Hit* pbs_find(const char *seq,
                  LTRboundaries *line,
                  unsigned long seqlen,
                  Bioseq *trna_lib,
                  LTRharvestoptions *lo,
                  Error *err)
{
  Seq *seq_forward, *seq_rev;
  unsigned long j, seqoffset, seqoffset_rev;
  char *seq_rev_full;
  Dlist *hitlist;
  Alignment *ali;
  Alpha *a = (Alpha*) alpha_new_dna();
  ScoreFunction *sf = dna_scorefunc_new(a,
                                        lo->pbs_ali_score_match,
                                        lo->pbs_ali_score_mismatch,
                                        lo->pbs_ali_score_insertion,
                                        lo->pbs_ali_score_deletion);

  /* We use a Dlist to maintain an ordered list of high-scoring local
     alignments. If this is not used in the future, a simple maximization
     could suffice. */
  hitlist = dlist_new(pbs_hit_compare);

  /* get reverse complement */
  seq_rev_full = ma_malloc(sizeof(char)*seqlen);
  memcpy(seq_rev_full, seq, sizeof(char)*seqlen);
  reverse_complement(seq_rev_full, seqlen, err);

  seqoffset = line->leftLTR_3-lo->pbs_radius;
  seqoffset_rev = line->rightLTR_5+lo->pbs_radius;

  seq_forward = seq_new(seq+
                          line->leftLTR_3-line->leftLTR_5-lo->pbs_radius+1,
                        2*lo->pbs_radius,
                        a);

  seq_rev     = seq_new(seq_rev_full+
                          line->rightLTR_3-line->rightLTR_5-lo->pbs_radius+1,
                        2*lo->pbs_radius,
                        a);

  for(j=0;j<bioseq_number_of_sequences(trna_lib);j++)
  {
    Seq *trna_seq, *trna_from3;
    char *trna_from3_full;
    unsigned long trna_seqlen;

    trna_seq = bioseq_get_seq(trna_lib, j);
    trna_seqlen = seq_length(trna_seq);

    trna_from3_full = ma_malloc(sizeof (char)*trna_seqlen);
    memcpy(trna_from3_full, seq_get_orig(trna_seq), sizeof(char)*trna_seqlen);
    reverse_complement(trna_from3_full, trna_seqlen, err);
    trna_from3 = seq_new_own(trna_from3_full, trna_seqlen, a);

    ali = swalign(seq_forward, trna_from3, sf);
    pbs_add_hit(hitlist,ali,lo,trna_seqlen,
                seq_get_description(trna_seq), STRAND_FORWARD,
                seqoffset);
    alignment_delete(ali);

    ali = swalign(seq_rev, trna_from3, sf);
    pbs_add_hit(hitlist,ali,lo,trna_seqlen,
                seq_get_description(trna_seq), STRAND_REVERSE,
                seqoffset_rev);
    alignment_delete(ali);

    seq_delete(trna_from3);
  }
  seq_delete(seq_forward);
  seq_delete(seq_rev);
  ma_free(seq_rev_full);
  scorefunction_delete(sf);
  alpha_delete(a);

  if(dlist_size(hitlist) > 0)
  {
    Dlistelem *delem;
    PBS_Hit* result;
    delem = dlist_first(hitlist);
    result = (PBS_Hit*) dlistelem_get_data(delem);
    for(delem=dlistelem_next(delem);delem;delem = dlistelem_next(delem))
    {
      PBS_Hit *hit = (PBS_Hit*) dlistelem_get_data(delem);
      ma_free(hit);
    }
    dlist_delete(hitlist);
    return result;
  }
  else
  {
    dlist_delete(hitlist);
    return NULL;
  }
}
