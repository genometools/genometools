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
#include "core/divmodmul.h"
#include "core/safearith.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "gth/align_dna_imp.h"
#include "gth/array2dim_plain.h"
#include "gth/compute_scores.h"
#include "gth/gthenum.h"
#include "gth/gtherror.h"
#include "gth/path_matrix.h"
#include "gth/path_walker.h"

/* IMPORTANT: Definition has to be consistent with DnaRetrace in
   align_dna_imp.h. */
static const char *dna_retracenames[]= {
  SHOWENUMTYPE(DNA_I_N),
  SHOWENUMTYPE(DNA_E_N),
  SHOWENUMTYPE(DNA_E_NM),
  SHOWENUMTYPE(DNA_I_NM),
  SHOWENUMTYPE(DNA_E_M),
  SHOWENUMTYPE(DNA_I_M),
  SHOWENUMTYPE(DNA_NUMOFRETRACE)
};

static const char* dna_showretracenames(DnaRetrace retrace)
{
  gt_assert(retrace <= DNA_NUMOFRETRACE);
  return dna_retracenames[retrace];
}

/* the following function allocates space for the DP tables for cDNAs/ESTs */
static int dp_matrix_init(GthDPMatrix *dpm,
                          unsigned long gen_dp_length,
                          unsigned long ref_dp_length,
                          unsigned long autoicmaxmatrixsize,
                          bool introncutout,
                          GthJumpTable *jump_table,
                          GthStat *stat)
{
  unsigned long t, n, m, matrixsize, sizeofpathtype =  sizeof (GthPath);

  /* XXX: adjust this check for QUARTER_MATRIX case */
  if (DNA_NUMOFSTATES * sizeofpathtype * (gen_dp_length + 1) >=
      (~0)/(ref_dp_length + 1)) {
    /* in this case the matrix would be larger than the addressable memory
       of this machine -> return ERROR_MATRIX_ALLOCATION_FAILED */
    return GTH_ERROR_MATRIX_ALLOCATION_FAILED;
  }

  matrixsize = gt_safe_mult_ulong((GT_DIV2(gen_dp_length + 1) +
                                   GT_MOD2(gen_dp_length + 1)),
                                   ref_dp_length + 1);

  if (!introncutout && autoicmaxmatrixsize > 0) {
    /* in this case the automatic intron cutout technique is enabled
       check if allocated matrix would be larger as specified maximal
       matrix size. If so, return matrix allocation error */
    if (sizeofpathtype * matrixsize * DNA_NUMOFSTATES >
        autoicmaxmatrixsize << 20) {
      return GTH_ERROR_MATRIX_ALLOCATION_FAILED;
    }
  }

  /* allocate space for dpm->path */
  if (jump_table) {
    gth_array2dim_plain_calloc(dpm->path,
                               GT_DIV2(gen_dp_length + 1) +
                               GT_MOD2(gen_dp_length + 1), ref_dp_length + 1);
  }
  else {
    gth_array2dim_plain_malloc(dpm->path,
                               GT_DIV2(gen_dp_length + 1) +
                               GT_MOD2(gen_dp_length + 1), ref_dp_length + 1);
  }
  dpm->path_jt = NULL;
  if (!dpm->path)
    return GTH_ERROR_MATRIX_ALLOCATION_FAILED;

  /* allocate space for dpm->score */
  for (t = DNA_E_STATE; t < DNA_NUMOFSTATES; t++) {
    for (n = 0; n < DNA_NUMOFSCORETABLES; n++) {
      dpm->score[t][n] = gt_malloc(sizeof (GthFlt) * (ref_dp_length + 1));
    }
  }

  /* allocating space for intronstart and exonstart */
  for (n = 0; n < DNA_NUMOFSCORETABLES; n++) {
    dpm->intronstart[n] = gt_calloc(ref_dp_length + 1,
                                    sizeof *dpm->intronstart);
    dpm->exonstart[n] = gt_calloc(ref_dp_length + 1, sizeof *dpm->exonstart);
  }

  /* initialize the DP matrices */
  dpm->path[0][0]  = DNA_E_NM;
  dpm->path[0][0] |= I_STATE_E_N;
  for (m = 1; m <= ref_dp_length; m++) {
    dpm->path[0][m]  = DNA_E_M;
    dpm->path[0][m] |= I_STATE_I_N;
  }

  for (n = 0; n < DNA_NUMOFSCORETABLES; n++) {
    dpm->score[DNA_E_STATE][n][0] = 0.0;
    dpm->score[DNA_I_STATE][n][0] = 0.0;

    for (m = 1; m <= ref_dp_length; m++) {
      dpm->score[DNA_E_STATE][n][m] = (GthFlt) 0.0;
      /* disallow intron status for 5' non-matching cDNA letters: */
      dpm->score[DNA_I_STATE][n][m] = (GthFlt) GTH_MINUSINFINITY;
    }
  }
  dpm->gen_dp_length = gen_dp_length;
  dpm->ref_dp_length = ref_dp_length;

  /* statistics */
  gth_stat_increment_numofbacktracematrixallocations(stat);
  gth_stat_increase_totalsizeofbacktracematricesinMB(stat,
                                           (sizeofpathtype * matrixsize) >> 20);

  return 0;
}

#if 0
static unsigned long dp_matrix_get_reference_length(DPMatrix *dpm)
{
  return dpm->ref_dp_length;
}
#endif

/* the following function evaluates state E_1m for indices 1 and <m> and
   stores a backtrace reference */
static void E_1m(GthDPMatrix *dpm, unsigned char genomicchar,
                 const unsigned char *ref_seq_tran, unsigned long m,
                 GtAlphabet *gen_alphabet, GthDbl log_probies,
                 GthDPOptionsEST *dp_options_est,
                 GthDPOptionsCore *dp_options_core)
{
  GthFlt value, maxvalue;
  GthPath retrace;
  GthDbl rval = 0.0, outputweight = 0.0;
  unsigned char referencechar = ref_seq_tran[m-1];
  unsigned int gen_alphabet_mapsize = gt_alphabet_size(gen_alphabet);

  /* 0. */
  rval = log_probies;
  ADDOUTPUTWEIGHT(rval,genomicchar,referencechar);
  if ((m < dp_options_est->wdecreasedoutput ||
       m > dpm->ref_dp_length - dp_options_est->wdecreasedoutput) &&
       genomicchar == referencechar) {
    ADDOUTPUTWEIGHT(outputweight, genomicchar, referencechar);
    rval -= (outputweight / 2.0);
  }
  maxvalue = (GthFlt) (dpm->score[DNA_E_STATE][0][m-1] + rval);
  retrace  = DNA_E_NM;

  /* 1. */
  outputweight = 0.0;
  rval = log_probies;
  ADDOUTPUTWEIGHT(rval,genomicchar,referencechar);
  if ((m < dp_options_est->wdecreasedoutput ||
       m > dpm->ref_dp_length - dp_options_est->wdecreasedoutput) &&
       genomicchar == referencechar) {
    ADDOUTPUTWEIGHT(outputweight, genomicchar, referencechar);
    rval -= (outputweight / 2.0);
  }
  value = (GthFlt) (dpm->score[DNA_I_STATE][0][m-1] + rval);
  /* intron from intronstart to n-1 => n-1 - intronstart + 1 */
  if (1 - dpm->intronstart[0][m - 1] < dp_options_core->dpminintronlength)
    value -= dp_options_core->shortintronpenalty;
  UPDATEMAX(DNA_I_NM);

  /* 2. */
  rval = log_probies;
  if (m < dpm->ref_dp_length) {
    ADDOUTPUTWEIGHT(rval, genomicchar, (unsigned char) DASH);
  }
  value = (GthFlt) (dpm->score[DNA_E_STATE][0][m] + rval);
  UPDATEMAX(DNA_E_N);

  /* 3. */
  rval = log_probies;
  if (m < dpm->ref_dp_length) {
    ADDOUTPUTWEIGHT(rval, genomicchar, (unsigned char) DASH);
  }
  value = (GthFlt) (dpm->score[DNA_I_STATE][0][m] + rval);
  /* intron from intronstart to n-1 => n-1 - intronstart + 1 */
  if (1 - dpm->intronstart[0][m] < dp_options_core->dpminintronlength)
    value -= dp_options_core->shortintronpenalty;
  UPDATEMAX(DNA_I_N);

  /* 4. */
  rval = log_probies;
  ADDOUTPUTWEIGHT(rval, (unsigned char) DASH, referencechar);
  value = (GthFlt) (dpm->score[DNA_E_STATE][1][m-1] + rval);
  UPDATEMAX(DNA_E_M);

  /* 5. */
  rval = dpm->score[DNA_I_STATE][1][m-1] + log_probies;
  ADDOUTPUTWEIGHT(rval, (unsigned char) DASH, referencechar);
  value = (GthFlt) (dpm->score[DNA_I_STATE][1][m-1] + rval);
  /* intron from intronstart to n => n - intronstart + 1 */
  if (1 - dpm->intronstart[1][m - 1] + 1 < dp_options_core->dpminintronlength)
    value -= dp_options_core->shortintronpenalty;
  UPDATEMAX(DNA_I_M);

  /* save maximum values */
  dpm->score[DNA_E_STATE][1][m] = maxvalue;
  dpm->path[0][m] |= (retrace << 4);

  switch (retrace) {
    case DNA_I_NM:
    case DNA_I_N:
    case DNA_I_M:
      dpm->exonstart[1][m] = 1;
      break;
    case DNA_E_NM:
      dpm->exonstart[1][m] = dpm->exonstart[0][m - 1];
      break;
    case DNA_E_N:
      dpm->exonstart[1][m] = dpm->exonstart[0][m];
      break;
    case DNA_E_M:
      dpm->exonstart[1][m] = dpm->exonstart[1][m - 1];
      break;
    default: gt_assert(0);
  }
}

/* the following function evaluates state I_1m for indices 1 and <m> and
   stores a backtrace reference */
static void I_1m(GthDPMatrix *dpm, unsigned long m, GthDbl log_1minusprobies)
{
  GthFlt value, maxvalue;
  GthPath retrace;

  /* 0. */
  maxvalue = (GthFlt) (dpm->score[DNA_E_STATE][0][m] + log_1minusprobies);
  retrace  = I_STATE_E_N;

  /* 1. */
  value = (GthFlt) (dpm->score[DNA_I_STATE][0][m] + log_1minusprobies);
  UPDATEMAX(I_STATE_I_N);

  /* save maximum values */
  dpm->score[DNA_I_STATE][1][m] = maxvalue;
  dpm->path[0][m] |= (retrace << 4);

  switch (retrace) {
    case I_STATE_E_N:
      /* begin of a new intron */
      dpm->intronstart[1][m] = 1;
      break;
    case I_STATE_I_N:
      /* continue existing intron */
      dpm->intronstart[1][m] = dpm->intronstart[0][m];
      break;
    default: gt_assert(0);
  }
}

/* the following function evaluate the dynamic programming tables */
static void dna_complete_path_matrix(GthDPMatrix *dpm,
                                     const unsigned char *gen_seq_tran,
                                     const unsigned char *ref_seq_tran,
                                     unsigned long genomic_offset,
                                     GtAlphabet *gen_alphabet,
                                     GthDPParam *dp_param,
                                     GthDPOptionsEST *dp_options_est,
                                     GthDPOptionsCore *dp_options_core)
{
  GthFlt value, maxvalue;
  GthPath retrace;
  unsigned long n, m, modn, modnminus1;
  GthDbl rval, outputweight, **outputweights,
         log_probies,          /* initial exon state probability */
         log_1minusprobies;    /* initial intron state probability */
  GthFlt log_probdelgen,       /* deletion in genomic sequence */
         log_1minusprobdelgen;
  unsigned char genomicchar, referencechar;
  unsigned int gen_alphabet_mapsize = gt_alphabet_size(gen_alphabet);

  gt_assert(dpm->gen_dp_length > 1);

  log_probies = (GthDbl) log((double) dp_options_est->probies);
  log_1minusprobies = (GthDbl) log(1.0 - dp_options_est->probies);
  log_probdelgen = (GthFlt) log((double) dp_options_est->probdelgen);
  log_1minusprobdelgen = (GthFlt) log(1.0 - dp_options_est->probdelgen);

  /* precompute outputweights
     XXX: move this to somewhere else, maybe make it smaller */
  gt_array2dim_calloc(outputweights, UCHAR_MAX+1, UCHAR_MAX+1);
  for (n = 0; n <= UCHAR_MAX; n++) {
    for (m = 0; m <= UCHAR_MAX; m++) {
      ADDOUTPUTWEIGHT(outputweights[n][m], n, m);
    }
  }

  if (!genomic_offset) {
    /* handle case for n equals 1 */
    dpm->path[0][0] |= UPPER_E_N;
    dpm->path[0][0] |= UPPER_I_STATE_I_N;

    /* stepping along the cDNA/EST sequence */
    for (m = 1; m <= dpm->ref_dp_length; m++) {
      E_1m(dpm, gen_seq_tran[0], ref_seq_tran, m, gen_alphabet, log_probies,
           dp_options_est, dp_options_core);
      I_1m(dpm, m, log_1minusprobies);
    }
  }

  /* handle all other n's
     stepping along the genomic sequence */
  if (genomic_offset)
    n = genomic_offset + 1;
  else
    n = 2;

  for (; n <= dpm->gen_dp_length; n++) {
    modn = GT_MOD2(n);
    modnminus1 = GT_MOD2(n-1);
    genomicchar = gen_seq_tran[n-1];

    if (modn) {
      dpm->path[GT_DIV2(n)][0] |= UPPER_E_N;
      dpm->path[GT_DIV2(n)][0] |= UPPER_I_STATE_I_N;
    }
    else {
      dpm->path[GT_DIV2(n)][0]  = DNA_E_N;
      dpm->path[GT_DIV2(n)][0] |= I_STATE_I_N;
    }

    /* stepping along the cDNA/EST sequence */
    for (m = 1; m <= dpm->ref_dp_length; m++) {
      referencechar = ref_seq_tran[m-1];

      /* evaluate E_nm */

      /* 0. */
      outputweight = 0.0;
      rval = (GthDbl) (log_1minusprobdelgen + dp_param->log_1minusPdonor[n-1]);
      rval += outputweights[genomicchar][referencechar];
      if ((m < dp_options_est->wdecreasedoutput ||
           m > dpm->ref_dp_length - dp_options_est->wdecreasedoutput) &&
           genomicchar == referencechar) {
        outputweight += outputweights[genomicchar][referencechar];
        rval -= (outputweight / 2.0);
      }
      maxvalue = (GthFlt) (dpm->score[DNA_E_STATE][modnminus1][m-1] + rval);
      retrace  = DNA_E_NM;

      /* 1. */
      outputweight = 0.0;
      rval = (GthDbl) (dp_param->log_Pacceptor[n-2] + log_1minusprobdelgen);
      rval += outputweights[genomicchar][referencechar];
      if ((m < dp_options_est->wdecreasedoutput ||
           m > dpm->ref_dp_length - dp_options_est->wdecreasedoutput) &&
           genomicchar == referencechar) {
        outputweight += outputweights[genomicchar][referencechar];
        rval -= (outputweight / 2.0);
      }
      value = (GthFlt) (dpm->score[DNA_I_STATE][modnminus1][m-1] + rval);
      /* intron from intronstart to n-1 => n-1 - intronstart + 1 */
      if (n - dpm->intronstart[modnminus1][m - 1] <
          dp_options_core->dpminintronlength) {
        value -= dp_options_core->shortintronpenalty;
      }
      UPDATEMAX(DNA_I_NM);

      /* 2. */
      rval = 0.0;
      if (m < dpm->ref_dp_length || n < dp_options_est->wzerotransition)
        rval += (log_1minusprobdelgen + dp_param->log_1minusPdonor[n-1]);
      if (m < dpm->ref_dp_length)
        rval += outputweights[genomicchar][DASH];
      value = (GthFlt) (dpm->score[DNA_E_STATE][modnminus1][m] + rval);
      UPDATEMAX(DNA_E_N);

      /* 3. */
      rval = (GthDbl) (dp_param->log_Pacceptor[n-2] + log_1minusprobdelgen);
      if (m < dpm->ref_dp_length)
        rval += outputweights[genomicchar][DASH];
      value = (GthFlt) (dpm->score[DNA_I_STATE][modnminus1][m] + rval);
      /* intron from intronstart to n-1 => n-1 - intronstart + 1 */
      if (n - dpm->intronstart[modnminus1][m] <
          dp_options_core->dpminintronlength) {
        value -= dp_options_core->shortintronpenalty;
      }
      UPDATEMAX(DNA_I_N);

      /* 4. */
      rval = 0.0;
      if (n < dpm->gen_dp_length || m < dp_options_est->wzerotransition)
        rval = (GthDbl) log_probdelgen;
      if (n < dpm->gen_dp_length)
        rval += outputweights[DASH][referencechar];
      value = (GthFlt) (dpm->score[DNA_E_STATE][modn][m-1] + rval);
      UPDATEMAX(DNA_E_M);

      /* 5. */
      rval = 0.0;
      if (n < dpm->gen_dp_length)
       rval += (dp_param->log_Pacceptor[n-1] + log_probdelgen);
      if (n < dpm->gen_dp_length)
        rval += outputweights[DASH][referencechar];
      value = (GthFlt) (dpm->score[DNA_I_STATE][modn][m-1] + rval);
      /* intron from intronstart to n => n - intronstart + 1 */
      if (n - dpm->intronstart[modn][m - 1] + 1 <
          dp_options_core->dpminintronlength) {
        value -= dp_options_core->shortintronpenalty;
      }
      UPDATEMAX(DNA_I_M);

      /* save maximum values */
      dpm->score[DNA_E_STATE][modn][m] = maxvalue;
      if (modn)
        dpm->path[GT_DIV2(n)][m] |= (retrace << 4);
      else
        dpm->path[GT_DIV2(n)][m]  = retrace;

      switch (retrace) {
        case DNA_I_NM:
        case DNA_I_N:
        case DNA_I_M:
          dpm->exonstart[modn][m] = n;
          break;
        case DNA_E_NM:
          dpm->exonstart[modn][m] = dpm->exonstart[modnminus1][m - 1];
          break;
        case DNA_E_N:
          dpm->exonstart[modn][m] = dpm->exonstart[modnminus1][m];
          break;
        case DNA_E_M:
          dpm->exonstart[modn][m] = dpm->exonstart[modn][m - 1];
          break;
        default: gt_assert(0);
      }

      /* evaluate I_nm */

      /* 0. */
      maxvalue = dpm->score[DNA_E_STATE][modnminus1][m] +
                 (log_1minusprobdelgen + dp_param->log_Pdonor[n-1]);
      if (n - dpm->exonstart[modnminus1][m]
          < dp_options_core->dpminexonlength) {
         maxvalue -= dp_options_core->shortexonpenalty;
      }
      retrace  = I_STATE_E_N;

      /* 1. */
      value = dpm->score[DNA_I_STATE][modnminus1][m];
      if (!dp_options_core->freeintrontrans && m < dpm->ref_dp_length)
        value += dp_param->log_1minusPacceptor[n-2];
      UPDATEMAX(I_STATE_I_N);

      /* save maximum values */
      dpm->score[DNA_I_STATE][modn][m] = maxvalue;
      if (modn)
        dpm->path[GT_DIV2(n)][m] |= (retrace << 4);
      else
        dpm->path[GT_DIV2(n)][m] |= retrace;

      switch (retrace) {
       case I_STATE_E_N:
          /* begin of a new intron */
          dpm->intronstart[modn][m] = n;
          break;
        case I_STATE_I_N:
          /* continue existing intron */
          dpm->intronstart[modn][m] = dpm->intronstart[modnminus1][m];
          break;
        default: gt_assert(0);
      }
    }
  }

  /* free space  */
  gt_array2dim_delete(outputweights);
}

static void dna_include_exon(GthBacktracePath *backtrace_path,
                             unsigned long exonlength)
{
  /* at least one editoperation already saved */
  gt_assert(gth_backtrace_path_length(backtrace_path));

  while (exonlength) {
    /* XXX: The following would be the better solution, but leads to problems
       in the score computation, because the backtrace then ``produces genomic
       bases'' which are not present in the spliced sequence.

       gth_backtrace_path_add_deletion(backtrace_path); */
    gth_backtrace_path_add_intron(backtrace_path);
    exonlength--;
  }
}

static void dna_include_intron(GthBacktracePath *backtrace_path,
                               unsigned long intronlength)
{
  /* at least one editoperation already saved */
  gt_assert(gth_backtrace_path_length(backtrace_path));

  while (intronlength) {
    gth_backtrace_path_add_intron(backtrace_path);
    intronlength--;
  }
}

static int dna_evaltracepath(GthBacktracePath *backtrace_path, GthDPMatrix *dpm,
                             const unsigned char *ref_seq_tran,
                             const unsigned char *gen_seq_tran,
                             DnaStates actualstate, bool introncutout,
                             GthSplicedSeq *spliced_seq, bool comments,
                             bool noicinintroncheck, GthPathMatrix *pm,
                             GtFile *outfp)
{
  unsigned long genptr = dpm->gen_dp_length, last_genptr = 0,
                refptr = dpm->ref_dp_length;
  GthPath pathtype, pathtype_jt = 0;
  bool lower;

  gt_assert(!gth_backtrace_path_length(backtrace_path));

  while ((genptr > 0) || (refptr > 0)) {
    /* here we map the quarter matrix bitvector stuff back on the simple Retrace
       types.  Thereby, no further changes on the backtracing procedure are
       necessary. */
    pathtype = dpm->path[GT_DIV2(genptr)][refptr];
    if (dpm->path_jt)
      pathtype_jt = dpm->path_jt[GT_DIV2(genptr)][refptr];
    lower = (bool) !GT_MOD2(genptr);
    if (lower) {
      if (actualstate == DNA_E_STATE) {
        pathtype  &= LOWER_E_STATE_MASK;
        pathtype_jt  &= LOWER_E_STATE_MASK;
      }
      else {
        pathtype  &= LOWER_I_STATE_MASK;
        pathtype >>= 3;
        pathtype_jt  &= LOWER_I_STATE_MASK;
        pathtype_jt >>= 3;
      }
    }
    else {
      if (actualstate == DNA_E_STATE) {
        pathtype  &= UPPER_E_STATE_MASK;
        pathtype >>= 4;
        pathtype_jt  &= UPPER_E_STATE_MASK;
        pathtype_jt >>= 4;
      }
      else {
        pathtype  &= UPPER_I_STATE_MASK;
        pathtype >>= 7;
        pathtype_jt  &= UPPER_I_STATE_MASK;
        pathtype_jt >>= 7;
      }
    }

    if (dpm->path_jt) {
#if 0
      printf("genptr=%lu, refptr=%lu\n", genptr, refptr);
      printf("pathtype=%d, pathtype_jt=%d\n", pathtype, pathtype_jt);
#endif
      gt_assert(pathtype == pathtype_jt);
    }

    if (introncutout) {
      if (genptr > 0) {
        if (gth_spliced_seq_pos_is_border(spliced_seq, genptr - 1)) {
          /* ensure that introns are only included into already existing introns
          */
          if (pathtype != DNA_I_N ||
              !gth_backtrace_path_last_is_intron(backtrace_path)) {
            if (noicinintroncheck) {
              /* intron cutout in intron check disabled,
                 include exon consisting of deletions */
              if (genptr != last_genptr) {
                last_genptr = genptr;
                dna_include_exon(backtrace_path,
                                 gth_spliced_seq_border_length(spliced_seq,
                                                               genptr - 1));
              }

            }
            else {
              if (comments) {
                gt_file_xprintf(outfp, "%c abort backtracing, intron cutout "
                                   "at p=%s (genpos=%lu (actual strand!))\n",
                                   COMMENTCHAR,
                                   dna_showretracenames((DnaRetrace) pathtype),
                                   spliced_seq->positionmapping[genptr]);
              }
              return GTH_ERROR_CUTOUT_NOT_IN_INTRON;
            }
          }
          else {
            /* include intron */
            dna_include_intron(backtrace_path,
                               gth_spliced_seq_border_length(spliced_seq,
                                                             genptr - 1));
          }
        }
      }
    }

    if (pm) {
      gth_path_matrix_set_max_path(pm, genptr, refptr,
                                   actualstate == DNA_E_STATE);
    }

    switch (actualstate) {
      case DNA_E_STATE:
        /* we are currently in an exon state, the possible pathtypes for exon
           states are handled here */
        switch (pathtype) {
          case DNA_E_NM:
            if (gen_seq_tran[genptr-1] == ref_seq_tran[refptr-1])
              gth_backtrace_path_add_match(backtrace_path, false);
            else
              gth_backtrace_path_add_mismatch(backtrace_path);
            genptr--;
            refptr--;
            break;
          case DNA_I_NM:
            if (gen_seq_tran[genptr-1] == ref_seq_tran[refptr-1])
              gth_backtrace_path_add_match(backtrace_path, false);
            else
              gth_backtrace_path_add_mismatch(backtrace_path);
            genptr--;
            refptr--;
            actualstate = DNA_I_STATE;
            break;
          case DNA_E_N:
            gth_backtrace_path_add_deletion(backtrace_path);
            genptr--;
            break;
          case DNA_I_N:
            gth_backtrace_path_add_deletion(backtrace_path);
            genptr--;
            actualstate = DNA_I_STATE;
            break;
          case DNA_E_M:
            gth_backtrace_path_add_insertion(backtrace_path);
            refptr--;
            break;
          case DNA_I_M:
            gth_backtrace_path_add_insertion(backtrace_path);
            refptr--;
            actualstate = DNA_I_STATE;
            break;
          default: gt_assert(0);
        }
        break;
      case DNA_I_STATE:
        /* we are currently in an intron state, the possible pathtypes for
           intron states are handled here */
        switch (pathtype) {
          case DNA_E_N:
            gth_backtrace_path_add_intron(backtrace_path);
            genptr--;
            actualstate = DNA_E_STATE;
            break;
          case DNA_I_N:
            gth_backtrace_path_add_intron(backtrace_path);
            genptr--;
            break;
          default: gt_assert(0);
        }
        break;
      default: gt_assert(0);
    }
  }

  /* genptr is at start of genseq */
  gt_assert(genptr == 0);
  /* refptr is at start of refseq */
  gt_assert(refptr == 0);

  if (pm)
    gth_path_matrix_set_max_path(pm, genptr, refptr, actualstate);

  return 0;
}

static int dna_find_optimal_path(GthBacktracePath *backtrace_path,
                                 GthDPMatrix *dpm,
                                 const unsigned char *ref_seq_tran,
                                 const unsigned char *gen_seq_tran,
                                 bool introncutout,
                                 GthSplicedSeq *spliced_seq,
                                 bool comments, bool noicinintroncheck,
                                 bool useintron, /* XXX */
                                 GthPathMatrix *pm,
                                 GtFile *outfp)
{
  int rval;
  GthFlt value, maxvalue;
  GthPath retrace;
  DnaStates state;

  if (useintron) {
    retrace = DNA_I_STATE;
  }
  else {
    maxvalue =
      dpm->score[DNA_E_STATE][GT_MOD2(dpm->gen_dp_length)][dpm->ref_dp_length];
    retrace  = DNA_E_STATE;

    for (state = (DnaStates) 1; state < DNA_NUMOFSTATES; state++) {
      value = dpm->score[state][GT_MOD2(dpm->gen_dp_length)]
                        [dpm->ref_dp_length];
      UPDATEMAX(state);
    }
  }

  if ((rval = dna_evaltracepath(backtrace_path, dpm, ref_seq_tran,
                                gen_seq_tran, (DnaStates) retrace,
                                introncutout, spliced_seq, comments,
                                noicinintroncheck, pm, outfp))) {
    return rval;
  }

  gt_assert(gth_backtrace_path_is_valid(backtrace_path));

  return 0;
}

/* the following function frees the DP tables for cDNAs/ESTs */
static void dp_matrix_free(GthDPMatrix *dpm)
{
  unsigned long t, n;

  /* freeing space for dpm->intronstart and dpm->exonstart */
  for (n = 0; n < DNA_NUMOFSCORETABLES; n++) {
    gt_free(dpm->intronstart[n]);
    gt_free(dpm->exonstart[n]);
  }

  /* freeing space for dpm->score */
  for (t = DNA_E_STATE; t < DNA_NUMOFSTATES; t++) {
    for (n = 0; n < DNA_NUMOFSCORETABLES; n++)
      gt_free(dpm->score[t][n]);
  }

  /* freeing space for dpm->path */
  gth_array2dim_plain_delete(dpm->path);
  if (dpm->path_jt)
    gt_array2dim_delete(dpm->path_jt);
}

#if 0
static GthFlt dp_matrix_get_score_e_state(DPMatrix *dpm,
                                                   unsigned long n,
                                                   unsigned long m)
{
  gt_assert(dpm && n < dpm->gen_dp_length && m < dpm->ref_dp_length);
  return dpm->score[DNA_E_STATE][GT_MOD2(n)][m];
}

static GthFlt dp_matrix_get_score_i_state(DPMatrix *dpm,
                                                   unsigned long n,
                                                   unsigned long m)
{
  gt_assert(dpm && n < dpm->gen_dp_length && m < dpm->ref_dp_length);
  return dpm->score[DNA_I_STATE][GT_MOD2(n)][m];
}

static void dp_matrix_set_score_e_state(DPMatrix *dpm, unsigned long n,
                                        unsigned long m, GthFlt score)
{
  gt_assert(dpm && n < dpm->gen_dp_length && m < dpm->ref_dp_length);
  dpm->score[DNA_E_STATE][GT_MOD2(n)][m] = score;
}

static void dp_matrix_set_score_i_state(DPMatrix *dpm, unsigned long n,
                                        unsigned long m, GthFlt score)
{
  gt_assert(dpm && n < dpm->gen_dp_length && m < dpm->ref_dp_length);
  dpm->score[DNA_I_STATE][GT_MOD2(n)][m] = score;
}

static unsigned long dp_matrix_get_intronstart(DPMatrix *dpm, unsigned long n,
                                                              unsigned long m)
{
  gt_assert(dpm && n < dpm->gen_dp_length && m < dpm->ref_dp_length);
  return dpm->intronstart[GT_MOD2(n)][m];
}

static void dp_matrix_set_intronstart(DPMatrix *dpm, unsigned long n,
                                      unsigned long m, unsigned long start)
{
  gt_assert(dpm && n < dpm->gen_dp_length && m < dpm->ref_dp_length);
  dpm->intronstart[GT_MOD2(n)][m] = start;
}

static unsigned long dp_matrix_get_exonstart(DPMatrix *dpm, unsigned long n,
                                                            unsigned long m)
{
  gt_assert(dpm && n < dpm->gen_dp_length && m < dpm->ref_dp_length);
  return dpm->exonstart[GT_MOD2(n)][m];
}

static void dp_matrix_set_exonstart(DPMatrix *dpm, unsigned long n,
                                    unsigned long m, unsigned long start)
{
  gt_assert(dpm && n < dpm->gen_dp_length && m < dpm->ref_dp_length);
  dpm->exonstart[GT_MOD2(n)][m] = start;
}

static DnaRetrace dp_matrix_get_retrace_e_state(DPMatrix *dpm, unsigned long n,
                                                               unsigned long m)
{
  GthPath pathtype;
  gt_assert(dpm && n < dpm->gen_dp_length && m < dpm->ref_dp_length);
  pathtype = dpm->path[GT_DIV2(n)][m];
  if (!GT_MOD2(n)) {
    pathtype  &= LOWER_E_STATE_MASK;
  }
  else {
    pathtype  &= UPPER_E_STATE_MASK;
    pathtype >>= 4;
  }
  return pathtype;
}

static DnaRetrace dp_matrix_get_retrace_i_state(DPMatrix *dpm, unsigned long n,
                                                               unsigned long m)
{
  GthPath pathtype;
  gt_assert(dpm && n < dpm->gen_dp_length && m < dpm->ref_dp_length);
  pathtype = dpm->path[GT_DIV2(n)][m];
  if (!GT_MOD2(n)) {
    pathtype  &= LOWER_I_STATE_MASK;
    pathtype >>= 3;
  }
  else {
    pathtype  &= UPPER_I_STATE_MASK;
    pathtype >>= 7;
  }
  return pathtype;
}

static void dp_matrix_set_retrace_e_state(DPMatrix *dpm, unsigned long n,
                                          unsigned long m, DnaRetrace retrace)
{
  gt_assert(dpm && n < dpm->gen_dp_length && m < dpm->ref_dp_length);
  gt_assert(retrace < DNA_NUMOFRETRACE);
  if (!GT_MOD2(n)) {
    dpm->path[GT_DIV2(n)][m] &= ~LOWER_E_STATE_MASK;
    dpm->path[GT_DIV2(n)][m] |= retrace;
  }
  else {
    dpm->path[GT_DIV2(n)][m] &= ~UPPER_E_STATE_MASK;
    dpm->path[GT_DIV2(n)][m] |= retrace << 4;
  }
}

static void dp_matrix_set_retrace_i_state(DPMatrix *dpm, unsigned long n,
                                          unsigned long m, DnaRetrace retrace)
{
  gt_assert(dpm && n < dpm->gen_dp_length && m < dpm->ref_dp_length);
  gt_assert(retrace == DNA_I_N || retrace == DNA_E_N);
  if (!GT_MOD2(n)) {
    dpm->path[GT_DIV2(n)][m] &= ~LOWER_I_STATE_MASK;
    dpm->path[GT_DIV2(n)][m] |= retrace << 3;
  }
  else {
    dpm->path[GT_DIV2(n)][m] &= ~UPPER_I_STATE_MASK;
    dpm->path[GT_DIV2(n)][m] |= retrace << 7;
  }
}
#endif

#if 0
static void dp_matrix_copy_terminal(DPMatrix *dpm_terminal, DPMatrix *dpm,
                                    unsigned long genomic_offset,
                                    unsigned long genomic_overlap,
                                    unsigned long reference_offset)
{
  GthFlt score;
  unsigned long n, m, start;
  DnaRetrace retrace;
  gt_assert(dpm_terminal && dpm);
  /* copy last score, intronstart, and exonstart row */
  n = genomic_overlap - 1;
  for (m = 0; m < dp_matrix_get_reference_length(dpm_terminal); m++) {
    score = dp_matrix_get_score_e_state(dpm, n + genomic_offset,
                                             m + reference_offset);
    dp_matrix_set_score_e_state(dpm_terminal, n, m, score);
    score = dp_matrix_get_score_i_state(dpm, n + genomic_offset,
                                             m + reference_offset);
    dp_matrix_set_score_i_state(dpm_terminal, n, m, score);
    start = dp_matrix_get_intronstart(dpm, n + genomic_offset,
                                           m + reference_offset);
    dp_matrix_set_intronstart(dpm_terminal, n, m, start);
    start = dp_matrix_get_exonstart(dpm, n + genomic_offset,
                                         m + reference_offset);
    dp_matrix_set_exonstart(dpm_terminal, n, m, start);
  }
  /* copy path matrix */
  for (n = 0; n < genomic_overlap; n++) {
    for (m = 0; m < dp_matrix_get_reference_length(dpm_terminal); m++) {
      retrace = dp_matrix_get_retrace_e_state(dpm, n + genomic_offset,
                                                   m + reference_offset);
      dp_matrix_set_retrace_e_state(dpm_terminal, n, m, retrace);
      retrace = dp_matrix_get_retrace_i_state(dpm, n + genomic_offset,
                                                   m + reference_offset);
      dp_matrix_set_retrace_i_state(dpm_terminal, n, m, retrace);
    }
  }
}
#endif

/* XXX */
#define DETECT_SMALL_EXONS_DP_LENGTH      100000
#define DETECT_SMALL_EXONS_EXTRA_OVERLAP  0

static void detect_small_terminal_exons(GthSA *sa,
                                        unsigned long gen_dp_start,
                                        const unsigned char *gen_seq_tran,
                                        unsigned long gen_dp_length,
                                        const unsigned char *ref_seq_tran,
                                        unsigned long ref_dp_length,
                                        const GtRange *gen_seq_bounds,
                                        bool comments,
                                        GthSpliceSiteModel *splice_site_model,
                                        GtAlphabet *gen_alphabet,
                                        GT_UNUSED GthDPMatrix *dpm,
                                        GthDPOptionsEST *dp_options_est,
                                        GthDPOptionsCore *dp_options_core,
                                        GthStat *stat, GtFile *outfp)
{
  unsigned long gen_dp_start_terminal, gen_dp_length_terminal,
                ref_dp_start_terminal, ref_dp_length_terminal;
  GthDPParam *dp_param_terminal;
  GthDPMatrix dpm_terminal;
  GthBacktracePath *backtrace_path;
  GthDPOptionsCore *dp_options_core_terminal;
  GthDPOptionsEST *dp_options_est_terminal;
  GT_UNUSED int rval; /*XXX*/

  /* determine positions of terminal DP */
  gen_dp_start_terminal = gen_dp_start + gen_dp_length -
                          gth_sa_genomiccutoff_end(sa);
  if (gen_dp_start_terminal >= gen_seq_bounds->end)
    return;
  gen_dp_length_terminal = gth_sa_genomiccutoff_end(sa) +
                             DETECT_SMALL_EXONS_DP_LENGTH;
  ref_dp_start_terminal = 0 + ref_dp_length -
                          gth_sa_referencecutoff_end(sa);
  ref_dp_length_terminal = gth_sa_referencecutoff_end(sa) +
                             DETECT_SMALL_EXONS_EXTRA_OVERLAP;
  if (comments) {
    gt_file_xprintf(outfp, "%c detect_small_terminal_exons():\n",
                       COMMENTCHAR);
    gt_file_xprintf(outfp, "%c gen_dp_start_terminal=%lu\n", COMMENTCHAR,
                       gen_dp_start_terminal);
    gt_file_xprintf(outfp, "%c gen_dp_length_terminal=%lu\n", COMMENTCHAR,
                       gen_dp_length_terminal);
    gt_file_xprintf(outfp, "%c gen_dp_length=%lu\n", COMMENTCHAR,
                       gen_dp_length);
    gt_file_xprintf(outfp, "%c gen_seq_bounds=%lu, %lu\n", COMMENTCHAR,
                       gen_seq_bounds->start, gen_seq_bounds->end);
  }
  if (gen_dp_start_terminal + gen_dp_length_terminal - 1 >
      gen_seq_bounds->end) {
    gen_dp_length_terminal -= gen_dp_start_terminal + gen_dp_length_terminal - 1
                              - gen_seq_bounds->end;
    gt_assert(gen_dp_start_terminal + gen_dp_length_terminal - 1 ==
              gen_seq_bounds->end);
  }
  gt_assert(gen_dp_start_terminal + gen_dp_length_terminal - 1 <=
            gen_seq_bounds->end);
  if (comments) {
    gt_file_xprintf(outfp, "%c gen_dp_start_terminal=%lu\n", COMMENTCHAR,
                       gen_dp_start_terminal);
    gt_file_xprintf(outfp, "%c gen_dp_length_terminal=%lu\n", COMMENTCHAR,
                       gen_dp_length_terminal);
  }

  if (dp_matrix_init(&dpm_terminal, gen_dp_length_terminal,
                     ref_dp_length_terminal, 0, false, NULL, stat)) {
    /* out of memory */
    return;
  }
  dp_param_terminal = gth_dp_param_new_with_range(gen_dp_start_terminal,
                                                  gen_dp_start_terminal +
                                                  gen_dp_length_terminal - 1,
                                                  gen_seq_tran, gen_seq_bounds,
                                                  splice_site_model,
                                                  gen_alphabet);
  if (!dp_param_terminal) {
    /* out of memory */
    dp_matrix_free(&dpm_terminal);
    return;
  }
  dp_options_core_terminal = gth_dp_options_core_clone(dp_options_core);
  dp_options_est_terminal = gth_dp_options_est_clone(dp_options_est);
  /* XXX */
#if 0
  dp_options_est_terminal->wzerotransition = 0; /* XXX */
  dp_options_est_terminal->wdecreasedoutput = 0; /* XXX */
  dp_options_core_terminal->freeintrontrans = true;
#endif
#if 0
  dp_matrix_copy_terminal(&dpm_terminal, dpm,
                          dpm->gen_dp_length - gth_sa_genomiccutoff_end(sa),
                          gth_sa_genomiccutoff_end(sa),
                          ref_dp_length - ref_dp_length_terminal);
#endif
  dna_complete_path_matrix(&dpm_terminal,
                           gen_seq_tran + gen_dp_start_terminal,
                          ref_seq_tran + ref_dp_length - ref_dp_length_terminal,
                           /*gth_sa_genomiccutoff_end(sa), */
                           0,
                           gen_alphabet,
                           dp_param_terminal, dp_options_est_terminal,
                           dp_options_core_terminal);
  backtrace_path = gth_backtrace_path_new(gen_dp_start_terminal,
                                          gen_dp_length_terminal,
                                          ref_dp_start_terminal,
                                          ref_dp_length_terminal);
  gth_backtrace_path_set_alphatype(backtrace_path, DNA_ALPHA);
  rval = dna_find_optimal_path(backtrace_path, &dpm_terminal,
                               ref_seq_tran + ref_dp_length
                               - ref_dp_length_terminal,
                               gen_seq_tran + gen_dp_start_terminal,
                               false, NULL, comments, false, false, NULL,
                               outfp);
  gt_assert(!rval);

  gt_assert(gth_sa_is_valid(sa)); /* XXX */

  gth_sa_cutoff_end(sa);
  gt_assert(gth_sa_is_valid(sa)); /* XXX */

  gth_sa_append(sa, backtrace_path);
  gt_assert(gth_sa_is_valid(sa)); /* XXX */

  gth_backtrace_path_delete(backtrace_path);
  gth_dp_param_delete(dp_param_terminal);
  dp_matrix_free(&dpm_terminal);
  gth_dp_options_core_delete(dp_options_core_terminal);
  gth_dp_options_est_delete(dp_options_est_terminal);
}

static void detect_small_initial_exons(GthSA *sa,
                                       unsigned long *gen_dp_start,
                                       const unsigned char *gen_seq_tran,
                                       unsigned long gen_dp_length,
                                       const unsigned char *ref_seq_tran,
                                       const GtRange *gen_seq_bounds,
                                       bool showeops, bool comments,
                                       GthSpliceSiteModel *splice_site_model,
                                       GtAlphabet *gen_alphabet,
                                       GthDPOptionsEST *dp_options_est,
                                       GthDPOptionsCore *dp_options_core,
                                       GthStat *stat, GtFile *outfp)
{
  unsigned long gen_dp_start_initial, gen_dp_length_initial,
                ref_dp_start_initial, ref_dp_length_initial;
  GthDPParam *dp_param_initial;
  GthDPMatrix dpm_initial;
  GthBacktracePath *backtrace_path;
  GthPathWalker *path_walker = NULL; /* XXX */
  GT_UNUSED int rval; /*XXX*/

  if (comments) {
    gt_file_xprintf(outfp, "%c detect_small_initial_exons():\n",
                       COMMENTCHAR);
    gt_file_xprintf(outfp, "%c gen_dp_start=%lu\n", COMMENTCHAR,
                       *gen_dp_start);
    gt_file_xprintf(outfp, "%c DETECT_SMALL_EXONS_DP_LENGTH=%d\n",
                       COMMENTCHAR, DETECT_SMALL_EXONS_DP_LENGTH);
    gt_file_xprintf(outfp, "%c genomic_cutoff=%lu\n", COMMENTCHAR,
                       gth_sa_genomiccutoff_start(sa));
  }

  /* determine positions of initial DP */
  ref_dp_start_initial = 0; /* XXX */
  ref_dp_length_initial = gth_sa_referencecutoff_start(sa) +
                            DETECT_SMALL_EXONS_EXTRA_OVERLAP;
  if (*gen_dp_start >= DETECT_SMALL_EXONS_DP_LENGTH +
                     gth_sa_genomiccutoff_start(sa)) {
    gen_dp_start_initial = *gen_dp_start - DETECT_SMALL_EXONS_DP_LENGTH;
    gen_dp_length_initial = gth_sa_genomiccutoff_start(sa) +
                              DETECT_SMALL_EXONS_DP_LENGTH;
  }
  else {
    gen_dp_start_initial = 0;
    gen_dp_length_initial = *gen_dp_start + gth_sa_genomiccutoff_start(sa);
  }

  path_walker = gth_path_walker_new(gth_sa_backtrace_path(sa), true);
  gth_path_walker_try_steps(path_walker, 10); /* XXX */
  gen_dp_length_initial += gth_path_walker_gen_distance(path_walker);
  ref_dp_length_initial += gth_path_walker_ref_distance(path_walker);

  if (comments) {
    gt_file_xprintf(outfp, "%c detect_small_initial_exons():\n",
                       COMMENTCHAR);
    gth_path_walker_show(path_walker, outfp);
    gt_file_xprintf(outfp, "%c gen_dp_start=%lu\n", COMMENTCHAR,
                       *gen_dp_start);
    gt_file_xprintf(outfp, "%c gen_dp_start_initial=%lu\n", COMMENTCHAR,
                       gen_dp_start_initial);
    gt_file_xprintf(outfp, "%c gen_dp_length_initial=%lu\n", COMMENTCHAR,
                       gen_dp_length_initial);
    gt_file_xprintf(outfp, "%c gen_dp_length=%lu\n", COMMENTCHAR,
                       gen_dp_length);
  }
  if (gen_dp_start_initial + gen_dp_length_initial - 1 > gen_seq_bounds->end) {
    gen_dp_length_initial -= gen_seq_bounds->end - gen_dp_start_initial +
                             gen_dp_length_initial - 1;
  }
  gt_assert(gen_dp_start_initial + gen_dp_length_initial - 1 <
            gen_seq_bounds->end);

  if (dp_matrix_init(&dpm_initial, gen_dp_length_initial,
                     ref_dp_length_initial, 0, false, NULL, stat)) {
    /* out of memory */
    return;
  }
  dp_param_initial = gth_dp_param_new_with_range(gen_dp_start_initial,
                                                 gen_dp_start_initial +
                                                 gen_dp_length_initial - 1,
                                                 gen_seq_tran, gen_seq_bounds,
                                                 splice_site_model,
                                                 gen_alphabet);
  if (!dp_param_initial) {
    /* out of memory */
    dp_matrix_free(&dpm_initial);
    return;
  }
  dna_complete_path_matrix(&dpm_initial,
                           gen_seq_tran + gen_dp_start_initial, ref_seq_tran, 0,
                           gen_alphabet, dp_param_initial, dp_options_est,
                           dp_options_core);
  backtrace_path = gth_backtrace_path_new(gen_dp_start_initial,
                                          gen_dp_length_initial,
                                          ref_dp_start_initial,
                                          ref_dp_length_initial);
  gth_backtrace_path_set_alphatype(backtrace_path, DNA_ALPHA);
  rval = dna_find_optimal_path(backtrace_path, &dpm_initial, ref_seq_tran,
                               gen_seq_tran + gen_dp_start_initial,
                               false, NULL, comments, false, false, NULL,
                               outfp);
  gt_assert(!rval);
  gt_assert(gth_sa_is_valid(sa)); /* XXX */

  gth_sa_cutoff_start(sa);
  gt_assert(gth_sa_is_valid(sa)); /* XXX */

  gth_sa_cutoff_walked_path(sa, path_walker, showeops, outfp);
  gt_assert(gth_sa_is_valid(sa)); /* XXX */

  if (comments)
    gth_backtrace_path_show(backtrace_path, false, 0, NULL);

  gth_sa_prepend(sa, backtrace_path);
  gt_assert(gth_sa_is_valid(sa)); /* XXX */

  gth_path_walker_delete(path_walker);
  gth_backtrace_path_delete(backtrace_path);
  gth_dp_param_delete(dp_param_initial);
  dp_matrix_free(&dpm_initial);
  *gen_dp_start = gen_dp_start_initial;
}

static void detect_small_exons(GthSA *sa, unsigned long *gen_dp_start,
                               const unsigned char *gen_seq_tran,
                               unsigned long gen_dp_length,
                               const unsigned char *ref_seq_tran,
                               unsigned long ref_dp_length,
                               const GtRange *gen_seq_bounds, bool showeops,
                               bool comments,
                               GthSpliceSiteModel *splice_site_model,
                               GtAlphabet *gen_alphabet, GthDPMatrix *dpm,
                               GthDPOptionsEST *dp_options_est,
                               GthDPOptionsCore *dp_options_core,
                               GthStat *stat, GtFile *outfp)
{
  gt_assert(sa);

  if (gth_sa_referencecutoff_end(sa)) { /* XXX: relax condition */
    /* the cDNA has not been aligned completely -> try to align end again */
    detect_small_terminal_exons(sa, *gen_dp_start, gen_seq_tran, gen_dp_length,
                                ref_seq_tran, ref_dp_length, gen_seq_bounds,
                                comments, splice_site_model, gen_alphabet,
                                dpm, dp_options_est, dp_options_core, stat,
                                outfp);
  }

#if 0
  if (gth_sa_referencecutoff_start(sa)) { /* XXX: relax condition */
#endif
    /* the cDNA has not been aligned completely -> try to align start again */
    detect_small_initial_exons(sa, gen_dp_start, gen_seq_tran, gen_dp_length,
                               ref_seq_tran, gen_seq_bounds, showeops, comments,
                               splice_site_model, gen_alphabet, dp_options_est,
                               dp_options_core, stat, outfp);
#if 0
  }
#endif
}

int gth_align_dna(GthSA *sa,
                  GtArray *gen_ranges,
                  const unsigned char *gen_seq_tran,
                  GT_UNUSED const unsigned char *gen_seq_orig,
                  const unsigned char *ref_seq_tran,
                  const unsigned char *ref_seq_orig,
                  unsigned long ref_dp_length,
                  GtAlphabet *gen_alphabet,
                  GtAlphabet *ref_alphabet,
                  bool introncutout,
                  unsigned long autoicmaxmatrixsize,
                  bool showeops,
                  bool comments,
                  bool gs2out,
                  const GtRange *gen_seq_bounds,
                  GthSpliceSiteModel *splice_site_model,
                  GthDPOptionsCore *dp_options_core,
                  GthDPOptionsEST *dp_options_est,
                  GthDPOptionsPostpro *dp_options_postpro,
                  GthDNACompletePathMatrixJT dna_complete_path_matrix_jt,
                  GthJumpTable *jump_table,
                  unsigned long ref_offset,
                  GthStat *stat,
                  GtFile *outfp)
{
  unsigned long gen_dp_start, gen_dp_end, gen_dp_length;
  GthSplicedSeq *spliced_seq = NULL;
  GthPathMatrix *pm = NULL;
  GthDPParam *dp_param;
  GthDPMatrix dpm;
  int rval;

  gt_assert(gen_ranges);

  /* initialization */
  gen_dp_start  = ((GtRange*) gt_array_get_first(gen_ranges))->start;
  gen_dp_end    = ((GtRange*) gt_array_get_last(gen_ranges))->end;
  gt_assert(gen_dp_start <= gen_dp_end);
  gen_dp_length = gen_dp_end - gen_dp_start + 1;
  dp_param = gth_dp_param_new(gen_ranges, gen_seq_tran, gen_seq_bounds,
                              splice_site_model, gen_alphabet);
  if (!dp_param)
    return GTH_ERROR_DP_PARAMETER_ALLOCATION_FAILED;
  if (introncutout) {
    spliced_seq = gth_spliced_seq_new_with_comments(gen_seq_tran, gen_ranges,
                                                    comments, outfp);
  }
  if ((rval = dp_matrix_init(&dpm,
                             introncutout ? spliced_seq->splicedseqlen
                                          : gen_dp_length,
                             ref_dp_length, autoicmaxmatrixsize, introncutout,
                             jump_table, stat))) {
    gth_dp_param_delete(dp_param);
    gth_spliced_seq_delete(spliced_seq);
    return rval;
  }
  gth_sa_set(sa, DNA_ALPHA, gen_dp_start, gen_dp_length);

  /* calculation */
  if (jump_table) {
    gt_assert(dna_complete_path_matrix_jt);
    dna_complete_path_matrix_jt(&dpm,
                                introncutout ? spliced_seq->splicedseq
                                             : gen_seq_tran + gen_dp_start,
                                ref_seq_tran, 0, gen_alphabet, dp_param,
                                dp_options_est, dp_options_core, jump_table,
                                gen_ranges, ref_dp_length, ref_offset, &pm);
  }
  else {
    dna_complete_path_matrix(&dpm,
                             introncutout ? spliced_seq->splicedseq
                                          : gen_seq_tran + gen_dp_start,
                             ref_seq_tran, 0, gen_alphabet, dp_param,
                             dp_options_est, dp_options_core);
  }

  /* debugging */
  if (!dpm.path_jt &&
      dp_options_core->btmatrixgenrange.start != GT_UNDEF_ULONG) {
    pm = gth_path_matrix_new(dpm.path, dpm.gen_dp_length, dpm.ref_dp_length,
                             &dp_options_core->btmatrixgenrange,
                             &dp_options_core->btmatrixrefrange, NULL);
  }

  /* backtracing */
  if ((rval = dna_find_optimal_path(gth_sa_backtrace_path(sa), &dpm,
                                    ref_seq_tran,
                                    introncutout ? spliced_seq->splicedseq
                                                 : gen_seq_tran + gen_dp_start,
                                    introncutout, spliced_seq, comments,
                                    dp_options_core->noicinintroncheck, false,
                                    pm, outfp))) {
    if (rval == GTH_ERROR_CUTOUT_NOT_IN_INTRON) {
      dp_matrix_free(&dpm);
      gth_dp_param_delete(dp_param);
      gth_spliced_seq_delete(spliced_seq);
    }
    return rval;
  }

  if (pm) {
    gth_path_matrix_show(pm);
    gth_path_matrix_delete(pm);
  }

  /* intron cutout is done after this point */
  if (showeops) {
    gt_file_xprintf(outfp, "showeops (before cutoffs): ");
    gth_backtrace_path_show(gth_sa_backtrace_path(sa), false, 0, outfp);
  }

  /* determine cutoffs if not switched off by command line option */
  gth_sa_determine_cutoffs(sa, dp_options_postpro->leadcutoffsmode,
                               dp_options_postpro->termcutoffsmode,
                               dp_options_postpro->cutoffsminexonlen);

  if (showeops) {
    gt_file_xprintf(outfp, "showeops (after cutoffs): ");
    gth_backtrace_path_show(gth_sa_backtrace_path(sa), false, 0, outfp);
  }

  /* remove zero base exons */
  gth_sa_remove_zero_base_exons(sa, stat);

  /* try to detect small exons */
  if (dp_options_est->detectsmallexons) {
    if (gth_sa_get_editoperations_length(sa)) {
      /* only try to detect small exons if the alignment does not already have
         length 0 */
      detect_small_exons(sa, &gen_dp_start, gen_seq_tran, gen_dp_length,
                         ref_seq_tran, ref_dp_length, gen_seq_bounds, showeops,
                         comments, splice_site_model, gen_alphabet, &dpm,
                         dp_options_est, dp_options_core, stat, outfp);
      gth_sa_determine_cutoffs(sa, dp_options_postpro->leadcutoffsmode,
                                   dp_options_postpro->termcutoffsmode,
                                   dp_options_postpro->cutoffsminexonlen);
      if (showeops) {
        gt_file_xprintf(outfp, "showeops (after small exon detection): ");
        gth_backtrace_path_show(gth_sa_backtrace_path(sa), false, 0, outfp);
      }
    }
  }

  /* compute borders and scores */
  gth_compute_scores(sa, false, dp_param, dp_options_est,
                     gen_seq_tran + gen_dp_start, ref_seq_tran, ref_seq_orig,
                     NULL, gen_dp_start, dp_options_postpro->scoreminexonlen,
                     introncutout, gs2out, spliced_seq, ref_dp_length,
                     gen_alphabet, ref_alphabet, NULL);

  /* free space */
  dp_matrix_free(&dpm);
  gth_dp_param_delete(dp_param);
  gth_spliced_seq_delete(spliced_seq);

  return 0;
}

GthSA* gth_align_dna_simple(GthInput *input,
                            const GtRange *gen_range,
                            unsigned long gen_file_num,
                            unsigned long gen_seq_num,
                            bool gen_strand_forward,
                            unsigned long ref_file_num,
                            unsigned long ref_seq_num,
                            GthSpliceSiteModel *splice_site_model)
{
  GthDPOptionsCore *dp_options_core;
  GthDPOptionsEST *dp_options_est;
  GthDPOptionsPostpro *dp_options_postpro;
  GtRange gen_seq_bounds, ref_seq_bounds;
  GtArray *gen_ranges;
  GthStat *stat;
  GthSA *sa;
  int rval;
  dp_options_core = gth_dp_options_core_new();
  dp_options_est = gth_dp_options_est_new();
  dp_options_postpro = gth_dp_options_postpro_new();
  stat = gth_stat_new();
  gen_seq_bounds = gth_input_get_genomic_range(input, gen_file_num,
                                               gen_seq_num);
  gt_assert(gen_range->start >= gen_seq_bounds.start);
  gt_assert(gen_range->end <= gen_seq_bounds.end);
  gth_input_load_reference_file(input, ref_file_num, true);
  ref_seq_bounds = gth_input_get_reference_range(input, ref_file_num,
                                                 ref_seq_num);
  sa = gth_sa_new_and_set(gen_strand_forward, true, input, gen_file_num,
                          gen_seq_num, ref_file_num, ref_seq_num, 1,
                          gt_range_length(&gen_seq_bounds),
                          gen_seq_bounds.start,
                          gt_range_length(&ref_seq_bounds));
  gen_ranges = gt_array_new(sizeof (GtRange));
  gt_array_add_elem(gen_ranges, (void*) gen_range, sizeof (GtRange));
  rval = gth_align_dna(sa,
                       gen_ranges,
                       gen_strand_forward
                       ? gth_input_current_gen_seq_tran(input)
                       : gth_input_current_gen_seq_tran_rc(input),
                       gen_strand_forward
                       ? gth_input_current_gen_seq_orig(input)
                       : gth_input_current_gen_seq_orig_rc(input),
                       gth_input_current_ref_seq_tran(input),
                       gth_input_current_ref_seq_orig(input),
                       gt_range_length(&ref_seq_bounds),
                       gth_input_current_gen_alphabet(input),
                       gth_input_current_ref_alphabet(input),
                       false,
                       0,
                       false,
                       false,
                       false,
                       &gen_seq_bounds,
                       splice_site_model,
                       dp_options_core,
                       dp_options_est,
                       dp_options_postpro,
                       NULL,
                       NULL,
                       0,
                       stat,
                       NULL);
  if (rval) {
    gth_sa_delete(sa);
    sa = NULL;
  }
  gt_array_delete(gen_ranges);
  gth_stat_delete(stat);
  gth_dp_options_postpro_delete(dp_options_postpro);
  gth_dp_options_est_delete(dp_options_est);
  gth_dp_options_core_delete(dp_options_core);
  return sa;
}
