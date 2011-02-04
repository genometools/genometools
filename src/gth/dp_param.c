/*
  Copyright (c) 2003-2011 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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
#include "core/chardef.h"
#include "core/ma_api.h"
#include "core/range.h"
#include "core/safearith.h"
#include "core/unused_api.h"
#include "gth/bssm_param_rep.h"
#include "gth/dp_param.h"
#include "gth/splice_site_model_rep.h"

/* XXX: why A? */
#define SUBSTITUTEWILDCARDWITHA(C)\
        if (C ==  WILDCARD)\
        {\
          C = gen_alphabet_symbolmap['A'];\
        }

#define CHECK_DP_PARAMETER_ALLOCATION(PTR)\
        if ((PTR) == NULL)\
        {\
          gth_dp_param_delete(dp_param);\
          return NULL;\
        }

static void evalsplicesiteprobformodel(GthFlt *prob, bool donorsite,
                                       const unsigned char *gen_seq_tran,
                                       const GtRange *gen_seq_bounds,
                                       unsigned long genpos,
                                       const GtUchar *gen_alphabet_symbolmap,
                                       GthBSSMModel *bssmmodel)
{
  unsigned long pc, /* previous char */
                cc, /* current char */
                d, i, j;
  long startpos, endpos;
  GthDbl pval = 0.5, Tv[3] = { 0.0 }, Fv[4] = { 0.0 };
  gt_assert(bssmmodel);

  /* set start and endpos */
  if (donorsite) {
    startpos = genpos - bssmmodel->window_size_left,
    endpos   = gt_safe_cast2long(genpos + bssmmodel->window_size_right + 1);
  }
  else { /* acceptorsite */
    startpos =  genpos - bssmmodel->window_size_left - 1,
    endpos   =  genpos + bssmmodel->window_size_right;
  }

  if ((startpos >= gt_safe_cast2long(gen_seq_bounds->start)) &&
      (gt_safe_cast2ulong(endpos) <= gen_seq_bounds->end)) {
    pc = gen_seq_tran[startpos];
    SUBSTITUTEWILDCARDWITHA(pc);
    if (bssmmodel->hypothesis_num == HYPOTHESIS2) {
      gt_assert(bssmmodel->hypotables.hypo2table);
      Tv[0] = (double) bssmmodel->hypotables.hypo2table[0][0][pc][0];
      Fv[0] = (double) bssmmodel->hypotables.hypo2table[1][0][pc][0];
    }
    else {
      gt_assert(bssmmodel->hypotables.hypo7table);
      Tv[0] = (double) bssmmodel->hypotables.hypo7table[0][0][pc][0];
      Tv[1] = (double) bssmmodel->hypotables.hypo7table[1][0][pc][0];
      Tv[2] = (double) bssmmodel->hypotables.hypo7table[2][0][pc][0];
      Fv[0] = (double) bssmmodel->hypotables.hypo7table[3][0][pc][0];
      Fv[1] = (double) bssmmodel->hypotables.hypo7table[4][0][pc][0];
      Fv[2] = (double) bssmmodel->hypotables.hypo7table[5][0][pc][0];
      Fv[3] = (double) bssmmodel->hypotables.hypo7table[6][0][pc][0];
    }
    d = 50 - bssmmodel->window_size_left;
    for (i = startpos + 1; i <= gt_safe_cast2ulong(endpos); i++) {
      j = d + (i - startpos);
      cc = gen_seq_tran[i];
      SUBSTITUTEWILDCARDWITHA(cc);
      if (bssmmodel->hypothesis_num == HYPOTHESIS2) {
        Tv[0] = Tv[0] * bssmmodel->hypotables.hypo2table[0][j][pc][cc];
        Fv[0] = Fv[0] * bssmmodel->hypotables.hypo2table[1][j][pc][cc];
      }
      else {
        Tv[0] = Tv[0] * bssmmodel->hypotables.hypo7table[0][j][pc][cc];
        Tv[1] = Tv[1] * bssmmodel->hypotables.hypo7table[1][j][pc][cc];
        Tv[2] = Tv[2] * bssmmodel->hypotables.hypo7table[2][j][pc][cc];
        Fv[0] = Fv[0] * bssmmodel->hypotables.hypo7table[3][j][pc][cc];
        Fv[1] = Fv[1] * bssmmodel->hypotables.hypo7table[4][j][pc][cc];
        Fv[2] = Fv[2] * bssmmodel->hypotables.hypo7table[5][j][pc][cc];
        Fv[3] = Fv[3] * bssmmodel->hypotables.hypo7table[6][j][pc][cc];
      }
      pc = cc;
    }
    if (bssmmodel->hypothesis_num == HYPOTHESIS2)
      pval = Tv[0] / (Tv[0] + Fv[0]);
    else {
      pval = Tv[0] + Tv[1] + Tv[2];
      pval = pval / (pval + Fv[0] + Fv[1] + Fv[2] + Fv[3]);
    }
    /* XXX: prevents problem with fission_yeast.bssm file */
    if (isnan(pval))
      pval = 0.0;
    /* pval is a valid probability */
    gt_assert(pval >= 0.0 && pval <= 1.0);
    pval = 2.0 * (pval - 0.5); /* XXX: ? */
  }
  else
    pval = 0.0;

  /* return probability */
  *prob = (GthFlt) pval;
}

static void evaldonorprob(GthFlt *prob, const unsigned char *gen_seq_tran,
                          const GtRange *gen_seq_bounds, unsigned long genpos,
                          const GtUchar *gen_alphabet_symbolmap,
                          GthBSSMParam *bssm_param)
{
  *prob = (GthFlt) 0.0;
  gt_assert(bssm_param);

  if (genpos < gen_seq_bounds->end) {
    if (bssm_param->gt_donor_model_set &&
        gen_seq_tran[genpos]     == gen_alphabet_symbolmap['G'] &&
        gen_seq_tran[genpos + 1] == gen_alphabet_symbolmap['T']) {
      evalsplicesiteprobformodel(prob, true, gen_seq_tran, gen_seq_bounds,
                                 genpos, gen_alphabet_symbolmap,
                                 &bssm_param->gt_donor_model);
    }
    else if (bssm_param->gc_donor_model_set &&
             gen_seq_tran[genpos]     == gen_alphabet_symbolmap['G'] &&
             gen_seq_tran[genpos + 1] == gen_alphabet_symbolmap['C']) {
      evalsplicesiteprobformodel(prob, true, gen_seq_tran, gen_seq_bounds,
                                 genpos, gen_alphabet_symbolmap,
                                 &bssm_param->gc_donor_model);
    }
  }
}

static void evalacceptorprob(GthFlt *prob, const unsigned char *gen_seq_tran,
                             const GtRange *gen_seq_bounds,
                             unsigned long genpos,
                             const GtUchar *gen_alphabet_symbolmap,
                             GthBSSMParam *bssm_param)
{
  *prob = (GthFlt) 0.0;
  gt_assert(bssm_param);

  if (genpos > 0 &&
      gen_seq_tran[genpos - 1] == gen_alphabet_symbolmap['A'] &&
      gen_seq_tran[genpos]     == gen_alphabet_symbolmap['G']) {
    if (bssm_param->ag_acceptor_model_set) {
      evalsplicesiteprobformodel(prob, false, gen_seq_tran, gen_seq_bounds,
                                 genpos, gen_alphabet_symbolmap,
                                 &bssm_param->ag_acceptor_model);
    }
  }
}

/*
  The following function evaluates the U12-type intron model.
  That is, if a consensus sequence /[AG]TATCCTT  (where / denotes the exon end
  and $[AG]$ indicates A or G) is found in the sequence part of
  <gen_seq_tran> denoted by <left> and <right>, the probabilities are
  updated.
*/

static void evaluateU12intronmodel(GthFlt *log_Pdonor, GthFlt *log_1minusPdonor,
                                   const unsigned char *gen_seq_tran,
                                   unsigned long left,
                                   unsigned long right,
                                   unsigned long probindex,
                                   GT_UNUSED unsigned long probindex_afterwards,
                                   GthFlt log_U12typedonorprob,
                                   GthFlt log1minus_U12typedonorprob,
                                   GthFlt log_U12typedonorprobonemismatch,
                                   GthFlt log1minus_U12typedonorprobonemismatch,
                                   unsigned char Achar,
                                   unsigned char Cchar,
                                   unsigned char Gchar,
                                   unsigned char Tchar)
{
  unsigned long mismatchcounter;
  const unsigned char *sptr, *pptr;

#define CONSENSUS_PATTERN_SIZE  8

#define CHECK_FOR_U12_CONSENSUS_CHAR(C)\
        if (*pptr != (C))\
        {\
          mismatchcounter++;\
          if (mismatchcounter > 1)\
          {\
            continue;\
          }\
        }\
        pptr--

  for (sptr = gen_seq_tran + left + CONSENSUS_PATTERN_SIZE - 1;
       sptr <= gen_seq_tran + right; sptr++, probindex++)
  {
    pptr            = sptr;
    mismatchcounter = 0;

    /* check the last 6 characters of consensus (we are going backwards)
       one mismatch is allowed in this 6 characters */
    CHECK_FOR_U12_CONSENSUS_CHAR(Tchar);
    CHECK_FOR_U12_CONSENSUS_CHAR(Tchar);
    CHECK_FOR_U12_CONSENSUS_CHAR(Cchar);
    CHECK_FOR_U12_CONSENSUS_CHAR(Cchar);
    CHECK_FOR_U12_CONSENSUS_CHAR(Tchar);
    CHECK_FOR_U12_CONSENSUS_CHAR(Achar);

    /* check second character of consensus (has to match exactly) */
    if (*pptr != Tchar)
      continue;
    pptr--;

    /* check first character of consensus (has to match exactly) */
    if (*pptr == Achar || *pptr == Gchar) {
      if (mismatchcounter == 0) {
        if (log_U12typedonorprob > log_Pdonor[probindex]) {
          log_Pdonor[probindex]       = log_U12typedonorprob;
          log_1minusPdonor[probindex] = log1minus_U12typedonorprob;
        }
      }
      else {
        if (log_U12typedonorprobonemismatch> log_Pdonor[probindex]) {
          log_Pdonor[probindex]       = log_U12typedonorprobonemismatch;
          log_1minusPdonor[probindex] = log1minus_U12typedonorprobonemismatch;
        }
      }
    }
  }

  /* probindex is correct */
  gt_assert(probindex + CONSENSUS_PATTERN_SIZE - 1 == probindex_afterwards);
}

static void calculateprobabilities(GtArray *ranges, unsigned long totallength,
                                   GthFlt *log_Pdonor,
                                   GthFlt *log_1minusPdonor,
                                   GthFlt *log_Pacceptor,
                                   GthFlt *log_1minusPacceptor,
                                   const unsigned char *gen_seq_tran,
                                   const GtRange *gen_seq_bounds,
                                   const GtUchar *gen_alphabet_symbolmap,
                                   GthSpliceSiteModel *ssm,
                                   bool gt_donor_model_set,
                                   bool gc_donor_model_set,
                                   bool ag_acceptor_model_set)
{
  const unsigned char Achar = gen_alphabet_symbolmap['A'],
                      Cchar = gen_alphabet_symbolmap['C'],
                      Gchar = gen_alphabet_symbolmap['G'],
                      Tchar = gen_alphabet_symbolmap['T'];
  unsigned long rangeindex,
                startpos,
                endpos,
                genomicindex,
                probindex = 0,
                probindexbackupforU12evaluation,
                numofranges = gt_array_size(ranges);
  unsigned char cc,
                ccminus1,
                ccplus1;
  bool lastgenomicbase = false,
       usegenericmodel = true;

#define SET_OTHER_SPLICE_SITE_PROB(P)\
        if (usegenericmodel)\
        {\
          log_Pdonor[(P)]          = ssm->log_genericothersplicesitep;\
          log_Pacceptor[(P)]       = ssm->log_genericothersplicesitep;\
          log_1minusPdonor[(P)]    = ssm->log1minus_genericothersplicesitep;\
          log_1minusPacceptor[(P)] = ssm->log1minus_genericothersplicesitep;\
        }\
        else\
        {\
          log_Pdonor[(P)]          = ssm->log_nongenericothersplicesitep;\
          log_Pacceptor[(P)]       = ssm->log_nongenericothersplicesitep;\
          log_1minusPdonor[(P)]    = ssm->log1minus_nongenericothersplicesitep;\
          log_1minusPacceptor[(P)] = ssm->log1minus_nongenericothersplicesitep;\
        }

  if (gt_donor_model_set || gc_donor_model_set || ag_acceptor_model_set)
    usegenericmodel = false;

  for (rangeindex = 0; rangeindex < numofranges; rangeindex++) {
    probindexbackupforU12evaluation = probindex;

    if (rangeindex == 0 &&
        ((GtRange*) gt_array_get_first(ranges))->start ==
        gen_seq_bounds->start) {
      /* in this case we consider the first genomic base */
      gt_assert(probindex == 0);
      SET_OTHER_SPLICE_SITE_PROB(0);
      probindex++;
      startpos = ((GtRange*) gt_array_get(ranges, rangeindex))->start + 1;
    }
    else
      startpos = ((GtRange*) gt_array_get(ranges, rangeindex))->start;

    if (rangeindex == numofranges - 1 &&
        ((GtRange*) gt_array_get(ranges, numofranges - 1))->end ==
        gen_seq_bounds->end) {
      /* in this case we consider the last genomic base */
      SET_OTHER_SPLICE_SITE_PROB(totallength - 1);
      lastgenomicbase = true; /* to increase probindex later */
      endpos = ((GtRange*) gt_array_get(ranges, rangeindex))->end - 1;
    }
    else
      endpos = ((GtRange*) gt_array_get(ranges, rangeindex))->end;

    for (genomicindex = startpos; genomicindex <= endpos; genomicindex++) {
      cc       = gen_seq_tran[genomicindex]; /* cc = "current character" */
      ccminus1 = gen_seq_tran[genomicindex - 1];
      ccplus1  = gen_seq_tran[genomicindex + 1];

      if (cc == Gchar) {
        if (ccplus1 == Tchar) {
          /* donor = GT */
          if (!gt_donor_model_set) {
            log_Pdonor[probindex]       = ssm->log_genericGTdonorprob;
            log_1minusPdonor[probindex] = ssm->log1minus_genericGTdonorprob;
          }
          else {
            log_Pdonor[probindex]       = ssm->log_nongenericGTdonorprob;
            log_1minusPdonor[probindex] = ssm->log1minus_nongenericGTdonorprob;
          }
        }
        else if (ccplus1 == Cchar) {
          /* donor = GC */
          if (!gc_donor_model_set) {
            log_Pdonor[probindex]       = ssm->log_genericGCdonorprob;
            log_1minusPdonor[probindex] = ssm->log1minus_genericGCdonorprob;
          }
          else {
            log_Pdonor[probindex]       = ssm->log_nongenericGCdonorprob;
            log_1minusPdonor[probindex] = ssm->log1minus_nongenericGCdonorprob;
          }
        }
        else {
          /* donor = GA | GG */
          if (usegenericmodel) {
            log_Pdonor[probindex]       = ssm->log_genericothersplicesitep;
            log_1minusPdonor[probindex] =
              ssm->log1minus_genericothersplicesitep;
          }
          else {
            log_Pdonor[probindex]       = ssm->log_nongenericothersplicesitep;
            log_1minusPdonor[probindex] =
              ssm->log1minus_nongenericothersplicesitep;
          }
        }
        if (ccminus1 == Achar) {
          if (!ag_acceptor_model_set) {
            /* acceptor = AG */
            log_Pacceptor[probindex]       = ssm->log_genericAGacceptorprob;
            log_1minusPacceptor[probindex] =
              ssm->log1minus_genericAGacceptorprob;
          }
          else {
            log_Pacceptor[probindex]       = ssm->log_nongenericAGacceptorprob;
            log_1minusPacceptor[probindex] =
              ssm->log1minus_nongenericAGacceptorprob;
          }
        }
        else {
          /* acceptor = CG | GG | TG */
          if (usegenericmodel) {
            log_Pacceptor[probindex]       = ssm->log_genericothersplicesitep;
            log_1minusPacceptor[probindex] =
              ssm->log1minus_genericothersplicesitep;
          }
          else {
            log_Pacceptor[probindex]       =
              ssm->log_nongenericothersplicesitep;
            log_1minusPacceptor[probindex] =
              ssm->log1minus_nongenericothersplicesitep;
          }
        }
      }
      else {
        SET_OTHER_SPLICE_SITE_PROB(probindex);

        if (cc == Achar && ccplus1 == Tchar) {
          /* donor = AT */
          if (ssm->useU12intronmodel) {
            log_Pdonor[probindex]          = ssm->log_nongenericATdonorprob;
            log_1minusPdonor[probindex]    =
              ssm->log1minus_nongenericATdonorprob;
          }
          else {
            log_Pdonor[probindex]          = ssm->log_genericATdonorprob;
            log_1minusPdonor[probindex]    =
              ssm->log1minus_genericATdonorprob;
          }
        }
        if (ccminus1 == Achar && cc == Cchar) {
          /* acceptor = AC */
          if (usegenericmodel) {
            log_Pacceptor[probindex]       = ssm->log_genericACacceptorprob;
            log_1minusPacceptor[probindex] =
              ssm->log1minus_genericACacceptorprob;
          }
          else {
            log_Pacceptor[probindex]       = ssm->log_nongenericACacceptorprob;
            log_1minusPacceptor[probindex] =
              ssm->log1minus_nongenericACacceptorprob;
          }
        }
      }
      probindex++;
    }

    if (lastgenomicbase)
      probindex++;

    /* evalutate U12-type intron model */
    if (ssm->useU12intronmodel) {
      evaluateU12intronmodel(log_Pdonor, log_1minusPdonor, gen_seq_tran,
                             ((GtRange*) gt_array_get(ranges, rangeindex))
                             ->start,
                             ((GtRange*) gt_array_get(ranges, rangeindex))
                             ->end, probindexbackupforU12evaluation,
                             probindex, ssm->log_U12typedonorprob,
                             ssm->log1minus_U12typedonorprob,
                             ssm->log_U12typedonorprobonemismatch,
                             ssm->log1minus_U12typedonorprobonemismatch, Achar,
                             Cchar, Gchar, Tchar);
    }
  }

  gt_assert(probindex == totallength);
}

static void filllogvaluesforonestrand(GtArray *ranges,
                                      const unsigned char *gen_seq_tran,
                                      const GtRange *gen_seq_bounds,
                                      GthFlt *log_Pdonor,
                                      GthFlt *log_1minusPdonor,
                                      GthFlt *log_Pacceptor,
                                      GthFlt *log_1minusPacceptor,
                                      const GtUchar *gen_alphabet_symbolmap,
                                      GthSpliceSiteModel *splice_site_model)
{
  unsigned long rangeindex,
       startpos,
       endpos,
       genomicindex,
       probindex   = 0,
       totallength = gt_ranges_total_length(ranges);
  GthFlt donorprob, acceptorprob;
  double log_donorprob, log_acceptorprob;

  if (!splice_site_model->bssm_param) {
    /* the ``generic'' species */
    calculateprobabilities(ranges, totallength, log_Pdonor, log_1minusPdonor,
                           log_Pacceptor, log_1minusPacceptor, gen_seq_tran,
                           gen_seq_bounds, gen_alphabet_symbolmap,
                           splice_site_model, false, false, false);
  }
  else {
    /* all other species */
    calculateprobabilities(ranges, totallength, log_Pdonor, log_1minusPdonor,
                           log_Pacceptor, log_1minusPacceptor, gen_seq_tran,
                           gen_seq_bounds, gen_alphabet_symbolmap,
                           splice_site_model,
                           splice_site_model->bssm_param->gt_donor_model_set,
                           splice_site_model->bssm_param->gc_donor_model_set,
                           splice_site_model->bssm_param
                           ->ag_acceptor_model_set);

    for (rangeindex = 0; rangeindex < gt_array_size(ranges); rangeindex++) {
      startpos = ((GtRange*) gt_array_get(ranges, rangeindex))->start;
      endpos   = ((GtRange*) gt_array_get(ranges, rangeindex))->end;

      for (genomicindex = startpos; genomicindex <= endpos; genomicindex++) {
        evaldonorprob(&donorprob, gen_seq_tran, gen_seq_bounds, genomicindex,
                       gen_alphabet_symbolmap, splice_site_model->bssm_param);
        if (donorprob > 0.0) {
          log_donorprob = log((double) donorprob);
          if (log_donorprob > (double) log_Pdonor[probindex]) {
            log_Pdonor[probindex]       = (GthFlt) log_donorprob;
            log_1minusPdonor[probindex] = (GthFlt)
                                          log(1.0 - (double) donorprob);
          }
        }
        evalacceptorprob(&acceptorprob, gen_seq_tran, gen_seq_bounds,
                         genomicindex, gen_alphabet_symbolmap,
                         splice_site_model->bssm_param);
        if (acceptorprob > 0.0) {
          log_acceptorprob = log((double) acceptorprob);
          if (log_acceptorprob > (double) log_Pacceptor[probindex]) {
            log_Pacceptor[probindex]       = (GthFlt) log_acceptorprob;
            log_1minusPacceptor[probindex] = (GthFlt)
                                             log(1.0 - (double) acceptorprob);
          }
        }
        probindex++;
      }
    }

    gt_assert(probindex == totallength);
  }
}

GthDPParam* dp_param_alloc(GtArray *ranges)
{
  unsigned long totallength;
  GthDPParam *dp_param;

  gt_assert(ranges);

  dp_param = gt_calloc(1, sizeof *dp_param);
  totallength = gt_ranges_total_length(ranges);

  /* try to allocate space for the DP parameter */
  dp_param->log_Pdonor = malloc(sizeof (GthFlt) * totallength);
  CHECK_DP_PARAMETER_ALLOCATION(dp_param->log_Pdonor);

  dp_param->log_1minusPdonor = malloc(sizeof (GthFlt) * totallength);
  CHECK_DP_PARAMETER_ALLOCATION(dp_param->log_1minusPdonor);

  dp_param->log_Pacceptor = malloc(sizeof (GthFlt) * totallength);
  CHECK_DP_PARAMETER_ALLOCATION(dp_param->log_Pacceptor);

  dp_param->log_1minusPacceptor = malloc(sizeof (GthFlt) * totallength);
  CHECK_DP_PARAMETER_ALLOCATION(dp_param->log_1minusPacceptor);

  return dp_param;
}

GthDPParam* gth_dp_param_new(GtArray *ranges,
                             const unsigned char *gen_seq_tran,
                             const GtRange *gen_seq_bounds,
                             GthSpliceSiteModel *splice_site_model,
                             GtAlphabet *gen_alphabet)
{
 GthDPParam *dp_param;
 gt_assert(ranges && gen_seq_tran && splice_site_model && gen_alphabet);
 if ((dp_param = dp_param_alloc(ranges))) {
   filllogvaluesforonestrand(ranges, gen_seq_tran, gen_seq_bounds,
                             dp_param->log_Pdonor,
                             dp_param->log_1minusPdonor,
                             dp_param->log_Pacceptor,
                             dp_param->log_1minusPacceptor,
                             gt_alphabet_symbolmap(gen_alphabet),
                             splice_site_model);
 }
 return dp_param;
}

GthDPParam* gth_dp_param_new_with_range(unsigned long start,
                                        unsigned long end,
                                        const unsigned char *gen_seq_tran,
                                        const GtRange *gen_seq_bounds,
                                        GthSpliceSiteModel *splice_site_model,
                                        GtAlphabet *gen_alphabet)
{
  GthDPParam *dp_param;
  GtRange range;
  GtArray *ranges;
  gt_assert(start <= end);
  gt_assert(gen_seq_tran && splice_site_model && gen_alphabet);
  ranges = gt_array_new(sizeof (GtRange));
  range.start = start;
  range.end = end;
  gt_array_add(ranges, range);
  dp_param = gth_dp_param_new(ranges, gen_seq_tran, gen_seq_bounds,
                              splice_site_model, gen_alphabet);
  gt_array_delete(ranges);
  return dp_param;
}

void gth_dp_param_delete(GthDPParam *dp_param)
{
  if (!dp_param) return;
  free(dp_param->log_Pdonor);
  free(dp_param->log_1minusPdonor);
  free(dp_param->log_Pacceptor);
  free(dp_param->log_1minusPacceptor);
  gt_free(dp_param);
}
