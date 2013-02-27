/*
  Copyright (c) 2003-2009 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2003-2008 Center for Bioinformatics, University of Hamburg

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
#include "core/codon_api.h"
#include "core/safearith.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "gth/default.h"
#include "gth/align_protein.h"
#include "gth/compute_scores.h"
#include "gth/gthsplicesitescr.h"

static void evalnewexonifpossible(bool proteineop, bool *newexon,
                                  bool *newintron, bool *firstexon,
                                  bool introncutout, GthSplicedSeq *spliced_seq,
                                  Exoninfo *exon, Introninfo *intron,
                                  GthSA *sa,
                                  Traversealignmentstate *travstate,
                                  GtAlphabet *gen_alphabet,
                                  GthDPParam *dp_param,
                                  GthDPOptionsEST *dp_options_est,
                                  const unsigned char *gen_seq_tran,
                                  const unsigned char *ref_seq_tran,
                                  unsigned long gen_dp_start)
{
  unsigned long splicedpos;

  if (*newexon) { /* in this case an intron will be saved */
    exon->leftgenomicexonborder = gen_dp_start + travstate->genomicptr;
    exon->leftreferenceexonborder = gt_safe_cast2ulong(travstate->referenceptr);
    *newexon   = false;
    *newintron = true;
    if (*firstexon)
      *firstexon = false;
    else
    {
      /* save acceptorsiteprobability */
      if (introncutout) {
        splicedpos =
          gth_spliced_seq_orig_to_spliced_pos(spliced_seq,
                  gt_safe_cast2ulong(travstate->genomicptr - 1 + gen_dp_start));
        if (splicedpos == GT_UNDEF_ULONG) {
          /* XXX: no spliced position has been found -> this is an artificially
             introduced intron, use 0.0 as acceptor site probabilty */
          intron->acceptorsiteprobability = 0.0;
        }
        else {
          intron->acceptorsiteprobability = (GthFlt)
                              exp((double) dp_param->log_Pacceptor[splicedpos]);
        }
      }
      else {
        intron->acceptorsiteprobability = (GthFlt) exp((double)
                              dp_param->log_Pacceptor[travstate->genomicptr-1]);
      }

      /* for cDNAs/ESTs: calculationg acceptorsitescore: going forward from here
       */
      if (proteineop)
        intron->acceptorsitescore = UNDEFINED_SPLICE_SITE_SCORE;
      else {
        gthcalcsplicesitescore(&intron->acceptorsitescore, travstate,
                               gen_seq_tran, ref_seq_tran, gen_alphabet,
                               dp_options_est, true);
      }

      /* saving the intron */
      gth_sa_add_intron(sa, intron);
    }
  }
}

static void evalnewintronifpossible(bool proteineop, bool *newexon,
                                    bool *newintron, bool lastintron,
                                    bool introncutout, bool gs2out,
                                    GthSplicedSeq *spliced_seq,
                                    Exoninfo *exon, Introninfo *intron,
                                    GthFlt *singleexonweight,
                                    GthFlt *maxsingleexonweight,
                                    GthFlt *overallexonweight,
                                    GthFlt *maxoverallexonweight,
                                    unsigned long
                                    *cumulativelengthofscoredexons,
                                    GthSA *sa,
                                    Traversealignmentstate *travstate,
                                    GtAlphabet *gen_alphabet,
                                    GthDPParam *dp_param,
                                    GthDPOptionsEST *dp_options_est,
                                    const unsigned char *gen_seq_tran,
                                    const unsigned char *ref_seq_tran,
                                    unsigned long gen_dp_start,
                                    unsigned long scoreminexonlen)
{
  unsigned long genomicexonlength, splicedpos;

  if (*newintron) { /* in this case an exon will be saved */
    exon->rightgenomicexonborder = gen_dp_start + travstate->genomicptr - 1;
    exon->rightreferenceexonborder = gt_safe_cast2ulong(travstate
                                                        ->referenceptr - 1);
    *newintron = false;
    *newexon   = true;
    if (*maxsingleexonweight > 0.0) {
      exon->exonscore = (GthDbl) ((*singleexonweight) /
                                            (*maxsingleexonweight));
    }
    else
      exon->exonscore = 0.0;

    /* for calculating the alignmentscore and the cumulative length of scored
       exons */
    genomicexonlength = exon->rightgenomicexonborder -
                        exon->leftgenomicexonborder + 1;

    /* short exons are not used for the alignment score */
    if (genomicexonlength >= scoreminexonlen) {
      *overallexonweight             += *singleexonweight;
      *maxoverallexonweight          += *maxsingleexonweight;
    }
    /* coverage includes short exons (not for gs2out): */
    if (gs2out) {
      if (genomicexonlength >= scoreminexonlen)
        *cumulativelengthofscoredexons += genomicexonlength;
    }
    else
      *cumulativelengthofscoredexons += genomicexonlength;

    /* saving the exon */
    gt_assert(exon->leftgenomicexonborder <= exon->rightgenomicexonborder);
    gth_sa_add_exon(sa, exon);

    /* resetting scores */
    *singleexonweight    = (GthFlt) 0.0;
    *maxsingleexonweight = (GthFlt) 0.0;

    /* if this is not the last intron, save the donor site stuff.
       if this is the last intron, this function has been called to save the
       last exon. Therefore, no saving of donor site stuff is necessary. */
    if (!lastintron) {
      /* save donorsiteprobability */
      if (introncutout) {
        splicedpos =
          gth_spliced_seq_orig_to_spliced_pos(spliced_seq,
                      gt_safe_cast2ulong(travstate->genomicptr + gen_dp_start));
        if (splicedpos == GT_UNDEF_ULONG) {
          /* XXX: no spliced position has been found -> this is an artificially
             introduced intron, use 0.0 as donor site probabilty */
          intron->donorsiteprobability = 0.0;
        }
        else {
          intron->donorsiteprobability = (GthFlt)
                                 exp((double) dp_param->log_Pdonor[splicedpos]);
        }
      }
      else {
        intron->donorsiteprobability = (GthFlt)
                      exp((double) dp_param->log_Pdonor[travstate->genomicptr]);
      }

      /* for the cDNAs/ESTs: calculationg donorsitescore: going back from here
       */
      if (proteineop)
        intron->donorsitescore = UNDEFINED_SPLICE_SITE_SCORE;
      else {
        gthcalcsplicesitescore(&intron->donorsitescore, travstate, gen_seq_tran,
                               ref_seq_tran, gen_alphabet, dp_options_est,
                               false);
      }
    }
  }
}

typedef struct
{
  bool proteineop,
       newexon,
       newintron,
#ifndef NDEBUG
       process_mismatch,
#endif
       firstexon,
       introncutout,
       gs2out;
  GthSplicedSeq *spliced_seq;
  GthFlt singleexonweight,              /* (GS2 = escr) */
                  maxsingleexonweight,           /* (GS2 = oescr) */
                  overallexonweight,             /* (GS2 = gescr) */
                  maxoverallexonweight;          /* (GS2 = ogescr) */
  unsigned long   cumulativelengthofscoredexons; /* (GS2 = mlght) */
  Exoninfo exon;
  Introninfo intron;
  GthSA *sa;
  GthDPParam *dp_param;
  GthDPOptionsEST *dp_options_est;      /* used if proteineop == false */
  const unsigned char *gen_seq_tran,
                      *ref_seq_tran,
                      *ref_seq_orig;
  const GtTransTable *transtable;
  unsigned long gen_dp_start,
                scoreminexonlen,
                ref_dp_length;
  GtAlphabet *gen_alphabet;     /* used if proteineop == true */
  const GtUchar *gen_alphabet_characters;
  GthDPScoresProtein *dp_scores_protein;
} Computebordersandscoresdata;

/*
  The following function can only be used for traversals in forward direction,
  because in the protein case we refer to upstream genomic positions.
*/

static void computescoresprocmismatchordeletion(Traversealignmentstate *state,
                                                void *data,
                                                GT_UNUSED
                                                unsigned long lengthofeop)
{
  Computebordersandscoresdata *d = (Computebordersandscoresdata*) data;
  unsigned int gen_alphabet_mapsize = gt_alphabet_size(d->gen_alphabet);
  GthDPOptionsEST *dp_options_est = d->dp_options_est;
  unsigned char genomicchar1,
                genomicchar2,
                genomicchar3,
#ifndef NDEBUG
                origreferencechar,
#endif
                codon;

  gt_assert(lengthofeop == 1);

  evalnewexonifpossible(d->proteineop, &d->newexon, &d->newintron,
                        &d->firstexon, d->introncutout, d->spliced_seq,
                        &d->exon, &d->intron, d->sa, state, d->gen_alphabet,
                        d->dp_param, dp_options_est, d->gen_seq_tran,
                        d->ref_seq_tran, d->gen_dp_start);

  if (d->proteineop) {
    if (state->processing_intron_with_1_base_left) {
      genomicchar1 = d->gen_seq_tran[state->firstbaseleftptr];
      genomicchar2 = d->gen_seq_tran[state->genomicptr];
      genomicchar3 = d->gen_seq_tran[state->genomicptr + 1];
    }
    else if (state->processing_intron_with_2_bases_left) {
      genomicchar1 = d->gen_seq_tran[state->firstbaseleftptr];
      genomicchar2 = d->gen_seq_tran[state->secondbaseleftptr];
      genomicchar3 = d->gen_seq_tran[state->genomicptr];
    }
    else {
      genomicchar1 = d->gen_seq_tran[state->genomicptr];
      genomicchar2 = d->gen_seq_tran[state->genomicptr + 1];
      genomicchar3 = d->gen_seq_tran[state->genomicptr + 2];
    }
    codon = gthgetcodon(genomicchar1, genomicchar2, genomicchar3,
                        d->gen_alphabet_characters, d->transtable);

#ifndef NDEBUG
    if (d->process_mismatch) {
      if (state->processing_intron_with_2_bases_left)
        origreferencechar = d->ref_seq_orig[state->referenceptr - 1];
      else
        origreferencechar = d->ref_seq_orig[state->referenceptr];
      /* genomic codon does not equal reference character */
      gt_assert(codon != origreferencechar);
    }
#endif

    d->maxsingleexonweight += GTHGETSCORE(d->dp_scores_protein, genomicchar1,
                                          genomicchar2, genomicchar3, codon);
  }
  else {
    genomicchar1 = d->gen_seq_tran[state->genomicptr];
    /* SK: replaced
       ADDOUTPUTWEIGHT(d->maxsingleexonweight, genomicchar1, genomicchar1);
       by the following */
    ADDOUTPUTWEIGHTIDENTITY(d->maxsingleexonweight, genomicchar1);
  }
}

static void computescoresprocmismatch(Traversealignmentstate *state,
                                      void *data, unsigned long lengthofeop)
{
#ifndef NDEBUG
  Computebordersandscoresdata *d = (Computebordersandscoresdata*) data;
  d->process_mismatch = true;
#endif
  computescoresprocmismatchordeletion(state, data, lengthofeop);
}

static void computescoresprocdeletion(Traversealignmentstate *state,
                                      void *data, unsigned long lengthofeop)
{
#ifndef NDEBUG
  Computebordersandscoresdata *d = (Computebordersandscoresdata*) data;
  d->process_mismatch = false;
#endif
  computescoresprocmismatchordeletion(state, data, lengthofeop);
}

static void computescoresprocmismatchordeletionwithgap(Traversealignmentstate
                                                       *state,
                                                       void *data,
                                                       GT_UNUSED
                                                       unsigned long
                                                       lengthofeop)
{
  Computebordersandscoresdata *d = (Computebordersandscoresdata*) data;
  GthDPOptionsEST *dp_options_est = d->dp_options_est;
  GthFlt score;

  gt_assert(lengthofeop == 1);
  /* a protein edit operation is processed */
  gt_assert(d->proteineop);
  /* we are not processing with 1 base left here */
  gt_assert(!state->processing_intron_with_1_base_left);
  /* we are not processing with 2 bases left here */
  gt_assert(!state->processing_intron_with_2_bases_left);

  evalnewexonifpossible(d->proteineop, &d->newexon, &d->newintron,
                        &d->firstexon, d->introncutout , d->spliced_seq,
                        &d->exon, &d->intron, d->sa, state, d->gen_alphabet,
                        d->dp_param, dp_options_est, d->gen_seq_tran,
                        d->ref_seq_tran, d->gen_dp_start);

  /* the score is the same, no matter how many DASH'es are given */
  score = GTHGETSCORE(d->dp_scores_protein, DASH, DASH, DASH, DASH);

  gt_assert(score < 0.0);

  /* we subtract the negative score here to increase the maxsingleexonweight
     this is a somewhat arbirtarily chosen value, since it does not reflect
     the maximum value which is possible with a proper codon */
  d->maxsingleexonweight -= score;
}

static void computebordersandscoresprocinsertion(Traversealignmentstate *state,
                                                 void *data,
                                                 GT_UNUSED
                                                 unsigned long lengthofeop)
{
  Computebordersandscoresdata *d = (Computebordersandscoresdata*) data;
  unsigned int gen_alphabet_mapsize = gt_alphabet_size(d->gen_alphabet);
  GthDPOptionsEST *dp_options_est = d->dp_options_est;
  unsigned char genomicchar, referencechar;
  GthFlt score;

  gt_assert(lengthofeop == 1);
  /* we are not processing with 1 base left here */
  gt_assert(!state->processing_intron_with_1_base_left);
  /* we are not processing with 2 bases left here */
  gt_assert(!state->processing_intron_with_2_bases_left);

  evalnewexonifpossible(d->proteineop, &d->newexon, &d->newintron,
                        &d->firstexon, d->introncutout, d->spliced_seq,
                        &d->exon, &d->intron, d->sa, state, d->gen_alphabet,
                        d->dp_param, dp_options_est, d->gen_seq_tran,
                        d->ref_seq_tran, d->gen_dp_start);

  referencechar = d->ref_seq_tran[state->referenceptr];
  if (d->proteineop) {
    score = GTHGETSCORE(d->dp_scores_protein, DASH, DASH, DASH, referencechar);

    gt_assert(score < 0.0);

    /* XXX: maybe remove this */
    d->singleexonweight += score;
    /* we subtract the negative score here to increase the maxsingleexonweight
       this is a somewhat arbirtarily chosen value, since it does not reflect
       the maximum value which is possible with a proper codon */
    d->maxsingleexonweight -= score;
    /* XXX: maybe better add the maximum score which can be achieved by
       a match with referencechar */
  }
  else {
    genomicchar   = (unsigned char) DASH;
    ADDOUTPUTWEIGHT(d->singleexonweight, genomicchar, referencechar);
    /* SK: replaced
       ADDOUTPUTWEIGHT(d->maxsingleexonweight, genomicchar, genomicchar);
       by the following */
    ADDOUTPUTWEIGHTIDENTITY(d->maxsingleexonweight, genomicchar);
  }
}

/*
  The following function can only be used for traversals in forward direction,
  because in the protein case we refer to upstream genomic positions.
*/

static void computebordersandscoresprocmatch(Traversealignmentstate *state,
                                             void *data,
                                             unsigned long lengthofeop)
{
  Computebordersandscoresdata *d = (Computebordersandscoresdata*) data;
  unsigned int gen_alphabet_mapsize = gt_alphabet_size(d->gen_alphabet);
  GthDPOptionsEST *dp_options_est = d->dp_options_est;
  unsigned char genomicchar1,
                genomicchar2,
                genomicchar3,
                referencechar,
                origreferencechar;
  GthFlt genomicinterimvalue   = 0.0,
                  referenceinterimvalue = 0.0;

  evalnewexonifpossible(d->proteineop, &d->newexon, &d->newintron,
                        &d->firstexon, d->introncutout, d->spliced_seq,
                        &d->exon, &d->intron, d->sa, state, d->gen_alphabet,
                        d->dp_param, d->dp_options_est, d->gen_seq_tran,
                        d->ref_seq_tran, d->gen_dp_start);

  if (d->proteineop) {
    if (state->processing_intron_with_1_base_left) {
      genomicchar1      = d->gen_seq_tran[state->firstbaseleftptr];
      genomicchar2      = d->gen_seq_tran[state->genomicptr];
      genomicchar3      = d->gen_seq_tran[state->genomicptr + 1];
      origreferencechar = d->ref_seq_orig[state->referenceptr];
    }
    else if (state->processing_intron_with_2_bases_left) {
      genomicchar1      = d->gen_seq_tran[state->firstbaseleftptr];
      genomicchar2      = d->gen_seq_tran[state->secondbaseleftptr];
      genomicchar3      = d->gen_seq_tran[state->genomicptr];
      origreferencechar = d->ref_seq_orig[state->referenceptr - 1];
      /*                                                     ^^^!
         we are processing a match after an intron with two bases left here.
         therefore, the reference pointer points already to the amino acid after
         this match, because the previous amino acid belongs to the two
         nucleotides before the intron.
         now it should be clear why we have to subtract 1 here. */
    }
    else {
      genomicchar1      = d->gen_seq_tran[state->genomicptr];
      genomicchar2      = d->gen_seq_tran[state->genomicptr + 1];
      genomicchar3      = d->gen_seq_tran[state->genomicptr + 2];
      origreferencechar = d->ref_seq_orig[state->referenceptr];
    }

    /* genomic codon equals reference character */
    gt_assert(origreferencechar == gthgetcodon(genomicchar1, genomicchar2,
                                               genomicchar3,
                                               d->gen_alphabet_characters,
                                               d->transtable));

    genomicinterimvalue  = GTHGETSCORE(d->dp_scores_protein, genomicchar1,
                                       genomicchar2, genomicchar3,
                                       origreferencechar);
    genomicinterimvalue *= lengthofeop;
    d->singleexonweight    += genomicinterimvalue;
    d->maxsingleexonweight += genomicinterimvalue;
  }
  else {
    genomicchar1      = d->gen_seq_tran[state->genomicptr];
    /* referenceptr in valid range */
    gt_assert(state->referenceptr >= 0 &&
              state->referenceptr <  (long) d->ref_dp_length);
    referencechar = d->ref_seq_tran[state->referenceptr];
    /* genomic char equals reference char */
    gt_assert(genomicchar1 == referencechar);
    ADDOUTPUTWEIGHT(referenceinterimvalue, genomicchar1, referencechar);
    /* SK: replaced
       ADDOUTPUTWEIGHT(genomicinterimvalue, genomicchar1, genomicchar1);
       by the following */
    ADDOUTPUTWEIGHTIDENTITY(genomicinterimvalue, genomicchar1);
    genomicinterimvalue   *= lengthofeop;
    referenceinterimvalue *= lengthofeop;
    d->singleexonweight  += referenceinterimvalue;
    d->maxsingleexonweight += genomicinterimvalue;
  }
}

static void computebordersandscoresprocintron(Traversealignmentstate *state,
                                              void *data,
                                              GT_UNUSED
                                              unsigned long lengthofeop)
{
  Computebordersandscoresdata *d = (Computebordersandscoresdata*) data;
  evalnewintronifpossible(d->proteineop, &d->newexon, &d->newintron,
                          false, d->introncutout, d->gs2out, d->spliced_seq,
                          &d->exon, &d->intron, &d->singleexonweight,
                          &d->maxsingleexonweight, &d->overallexonweight,
                          &d->maxoverallexonweight,
                          &d->cumulativelengthofscoredexons , d->sa, state,
                          d->gen_alphabet, d->dp_param, d->dp_options_est,
                          d->gen_seq_tran, d->ref_seq_tran, d->gen_dp_start,
                          d->scoreminexonlen);
}

#ifndef NDEBUG
static bool containsintronsorinsertions(bool leading,
                                        Editoperation *alignment,
                                        long alignmentlength,
                                        bool proteineop)
{
  Eoptype eoptype;
  long i;
  bool breakforloop = false;

  /* check for introns or insertions */
  for (i = leading ? alignmentlength - 1 : 0;
           leading ? i >= 0 : i < alignmentlength;
           leading ? i-- : i++) {
    eoptype = gt_editoperation_type(alignment[i], proteineop);

    /* if match, mismatch, or deletion -> break
       if insertion or intron -> return true */
    switch (eoptype) {
      case EOP_TYPE_MATCH:
      case EOP_TYPE_MISMATCH:
      case EOP_TYPE_MISMATCH_WITH_1_GAP:
      case EOP_TYPE_MISMATCH_WITH_2_GAPS:
      case EOP_TYPE_DELETION:
      case EOP_TYPE_DELETION_WITH_1_GAP:
      case EOP_TYPE_DELETION_WITH_2_GAPS:
        breakforloop = true;
        break;
      case EOP_TYPE_INSERTION:
      case EOP_TYPE_INTRON:
      case EOP_TYPE_INTRON_WITH_1_BASE_LEFT:
      case EOP_TYPE_INTRON_WITH_2_BASES_LEFT:
        return true;
      default: gt_assert(0);
    }
    if (breakforloop)
      break;
  }

  /* no introns or insertions found -> return false */
  return false;
}
#endif

#ifndef NDEBUG
static bool containsnoleadingorterminalintronsorinsertions(Editoperation
                                                           *alignment,
                                                           long
                                                           alignmentlength,
                                                           bool proteineop)
{
  /* check for leading introns or insertions */
  if (containsintronsorinsertions(true, alignment, alignmentlength, proteineop))
    return false;

  /* check for terminal introns or insertions */
  if (containsintronsorinsertions(false, alignment, alignmentlength,
      proteineop)) {
    return false;
  }

  return true;
}
#endif

void gth_compute_scores(GthSA *sa,
                        bool proteineop,
                        GthDPParam *dp_param,
                        void *dp_options_est,
                        const unsigned char *gen_seq_tran,
                        const unsigned char *ref_seq_tran,
                        const unsigned char *ref_seq_orig,
                        const GtTransTable *transtable,
                        unsigned long gen_dp_start,
                        unsigned long scoreminexonlen,
                        bool introncutout,
                        bool gs2out,
                        GthSplicedSeq *spliced_seq,
                        unsigned long ref_dp_length,
                        GtAlphabet *gen_alphabet,
                        GtAlphabet *ref_alphabet,
                        GthDPScoresProtein *dp_scores_protein)
{
  Traversealignmentfunctions travfunctions;
  Traversealignmentstate travstate;
  Computebordersandscoresdata data;
  GthFlt score, coverageofgenomicsegment, coverageofreferencesegment;

  gt_assert(!gth_sa_num_of_exons(sa));
  gt_assert(!gth_sa_num_of_introns(sa));

  travfunctions.processmismatch  = computescoresprocmismatch;
  travfunctions.processdeletion  = computescoresprocdeletion;
  travfunctions.processinsertion = computebordersandscoresprocinsertion;
  travfunctions.processmatch     = computebordersandscoresprocmatch;
  travfunctions.processintron    = computebordersandscoresprocintron;
  travfunctions.breakcondition   = NULL;

  /* additional functions for protein edit operations */
  travfunctions.processintron_with_1_base_left  =
    computebordersandscoresprocintron;
  travfunctions.processintron_with_2_bases_left =
    computebordersandscoresprocintron;
  travfunctions.processmismatch_with_1_gap      =
    computescoresprocmismatchordeletionwithgap;
  travfunctions.processmismatch_with_2_gaps     =
    computescoresprocmismatchordeletionwithgap;
  travfunctions.processdeletion_with_1_gap      =
    computescoresprocmismatchordeletionwithgap;
  travfunctions.processdeletion_with_2_gaps     =
    computescoresprocmismatchordeletionwithgap;

  travstate.proteineop = proteineop;
  travstate.processing_intron_with_1_base_left  = false;
  travstate.processing_intron_with_2_bases_left = false;
  travstate.alignment = gth_sa_get_editoperations(sa);
  travstate.alignmentlength =  gth_sa_get_editoperations_length(sa);
  travstate.eopptr       = travstate.alignment + travstate.alignmentlength - 1;
  travstate.genomicptr   = gth_sa_genomiccutoff_start(sa);
  travstate.referenceptr = gth_sa_referencecutoff_start(sa);

  if (travstate.alignmentlength <= 0) {
    /* in this case the alignmentscore is set to 0, which leads to discarding
       this alignment later */
    gth_sa_set_score(sa, 0.0);
    return;
  }

  /* editoperations contain no zero base exons */
  gt_assert(gth_sa_contains_no_zero_base_exons(sa));
  /* editoperations contain no leading or terminal introns or insertions */
  gt_assert(containsnoleadingorterminalintronsorinsertions(travstate.alignment,
                                                        travstate
                                                        .alignmentlength,
                                                        proteineop));
  /* sum of edit operations equals referencelength */
  gt_assert(gt_eops_equal_referencelength(travstate.alignment,
                                       travstate.alignmentlength,
                                       ref_dp_length
                                       - gth_sa_referencecutoff_start(sa)
                                       - gth_sa_referencecutoff_end(sa),
                                       proteineop));

  data.proteineop                     = proteineop;
  data.newexon                        = true;
  data.newintron                      = true;
  data.firstexon                      = true;
  data.introncutout                   = introncutout;
  data.gs2out                         = gs2out;
  data.spliced_seq                    = spliced_seq;
  data.singleexonweight               = (GthFlt) 0.0;
  data.maxsingleexonweight            = (GthFlt) 0.0;
  data.overallexonweight              = (GthFlt) 0.0;
  data.maxoverallexonweight           = (GthFlt) 0.0;
  data.cumulativelengthofscoredexons  = 0;

  data.exon.leftgenomicexonborder     = GT_UNDEF_ULONG;
  data.exon.rightgenomicexonborder    = GT_UNDEF_ULONG;
  data.exon.leftreferenceexonborder   = GT_UNDEF_ULONG;
  data.exon.rightreferenceexonborder  = GT_UNDEF_ULONG;
  data.exon.exonscore                 = GTH_UNDEF_GTHDBL;

  data.intron.donorsiteprobability    = GTH_UNDEF_GTHFLT;
  data.intron.acceptorsiteprobability = GTH_UNDEF_GTHFLT;
  data.intron.donorsitescore          = GTH_UNDEF_GTHDBL;
  data.intron.acceptorsitescore       = GTH_UNDEF_GTHDBL;

  data.sa                             = sa;
  data.dp_param                       = dp_param;
  data.dp_options_est                 = dp_options_est;
  data.gen_seq_tran                   = gen_seq_tran;
  data.ref_seq_tran                   = ref_seq_tran;
  data.ref_seq_orig                   = ref_seq_orig;
  data.transtable                     = transtable;
  data.gen_dp_start                   = gen_dp_start;
  data.scoreminexonlen                = scoreminexonlen;
  data.ref_dp_length                  = ref_dp_length;
  data.gen_alphabet                   = gen_alphabet;
  data.gen_alphabet_characters        = gen_alphabet
                                        ? gt_alphabet_characters(gen_alphabet)
                                        : NULL;
  data.dp_scores_protein              = dp_scores_protein;

  gthtraversealignment(true, &travstate, proteineop, &data, &travfunctions);

  /* this is for saving the last exon */
  evalnewintronifpossible(proteineop, &data.newexon, &data.newintron, true,
                          data.introncutout, data.gs2out, data.spliced_seq,
                          &data.exon, &data.intron, &data.singleexonweight,
                          &data.maxsingleexonweight, &data.overallexonweight,
                          &data.maxoverallexonweight,
                          &data.cumulativelengthofscoredexons, sa, &travstate,
                          gen_alphabet, data.dp_param, data.dp_options_est,
                          data.gen_seq_tran, data.ref_seq_tran,
                          data.gen_dp_start, data.scoreminexonlen);

  /* saving the scores for the whole alignment */
  if (data.maxoverallexonweight > 0.0) {
    score = data.overallexonweight / data.maxoverallexonweight;
    /* XXX: the way the alignmentscore is computed, it is possible to get a
       score > 1.0. Since we don't want this, we cap it */
    if (score > 1.0)
      score = 1.0;
  }
  else
    score = 0.0;
  gth_sa_set_score(sa, score);
  gth_sa_set_cumlen_scored_exons(sa, data.cumulativelengthofscoredexons);

  /* fraction of the gen_dp_length which is scored/weighted */
  coverageofgenomicsegment   = (GthFlt) data.cumulativelengthofscoredexons /
                               (GthFlt) gth_sa_gen_dp_length(sa);
  /* coverage of genomic segment is valid value */
  gt_assert(coverageofgenomicsegment >= 0.0 && coverageofgenomicsegment <= 1.0);

  /* fraction of the referencelength which is scored/weighted */
  coverageofreferencesegment = (GthFlt) data.cumulativelengthofscoredexons /
                               (GthFlt) ((proteineop ? GT_CODON_LENGTH : 1) *
                                         gth_sa_ref_total_length(sa));

  if (coverageofgenomicsegment > coverageofreferencesegment) {
    gth_sa_set_coverage(sa, coverageofgenomicsegment);
    gth_sa_set_highest_cov(sa, true);
  }
  else {
    gth_sa_set_coverage(sa, coverageofreferencesegment);
    gth_sa_set_highest_cov(sa, false);
  }

  /* test the assumption that the coverage is never larger then the default */
  gt_assert(gth_sa_coverage(sa) <= GTH_DEFAULT_MAX_COVERAGE);

  /* compute poly(A) tail position */
  gth_sa_calc_polyAtailpos(sa, ref_seq_tran, ref_alphabet);

  /* determined exons are forward and consecutive */
  gt_assert(gth_sa_exons_are_forward_and_consecutive(sa));
}
