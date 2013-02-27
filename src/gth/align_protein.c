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

#include "core/codon_api.h"
#include "core/divmodmul.h"
#include "core/safearith.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "gth/array2dim_plain.h"
#include "gth/gthenum.h"
#include "gth/gtherror.h"
#include "gth/align_protein_imp.h"
#include "gth/compute_scores.h"

#define WSIZE_PROTEIN   20
#define WSIZE_DNA       60 /* (3 * WSIZE_PROTEIN) */

/* IMPORTANT: Definition has to be consistent with Retrace in
 * align_protein_imp.h. */
static const char *retracenames[]= {
  SHOWENUMTYPE(E_N3M),
  SHOWENUMTYPE(E_N2M),
  SHOWENUMTYPE(E_N1M),
  SHOWENUMTYPE(E_M),
  SHOWENUMTYPE(E_N3),
  SHOWENUMTYPE(E_N2),
  SHOWENUMTYPE(E_N1),
  SHOWENUMTYPE(IA_N3M),
  SHOWENUMTYPE(IB_N2M),
  SHOWENUMTYPE(IC_N1M),
  SHOWENUMTYPE(IA_N1),
  SHOWENUMTYPE(IB_N1),
  SHOWENUMTYPE(IC_N1),
  SHOWENUMTYPE(NUMOFRETRACE)
};

static void dp_table_core_free(DPtablecore *core)
{
  unsigned long t, n;

  /* freeing space for core->score and core->path */
  for (t =  E_STATE; t < PROTEIN_NUMOFSTATES; t++) {
    for (n = 0; n < PROTEIN_NUMOFSCORETABLES; n++)
      gt_free(core->score[t][n]);
  }
  gth_array2dim_plain_delete(core->path);
}

static GthPath path_e_state_read(GthDPtables *dpm, unsigned int n,
                                 unsigned int m)
{
  return dpm->core.path[n][m] & E_STATE_MASK;
}

static void path_e_state_write(GthDPtables *dpm, unsigned int n, unsigned int m,
                               GthPath e)
{
  dpm->core.path[n][m] = e;
}

static GthPath path_ia_state_read(GthDPtables *dpm, unsigned int n,
                                  unsigned int m)
{
  if ((dpm->core.path[n][m] & IA_STATE_MASK) >> 4)
    return E_N1;
  return IA_N1;
}

static void path_ia_state_write(GthDPtables *dpm, unsigned int n,
                                unsigned int m, GthPath e)
{
  if (e == E_N1)
    dpm->core.path[n][m] |= (1 << 4);
}

static GthPath path_ib_state_read(GthDPtables *dpm, unsigned int n,
                                  unsigned int m)
{
  if ((dpm->core.path[n][m] & IB_STATE_MASK) >> 5)
    return E_N2;
  return IB_N1;
}

static void path_ib_state_write(GthDPtables *dpm, unsigned int n,
                                unsigned int m, GthPath e)
{
  if (e == E_N2)
    dpm->core.path[n][m] |= (1 << 5);
}

static GthPath path_ic_state_read(GthDPtables *dpm, unsigned int n,
                                  unsigned int m)
{
  if ((dpm->core.path[n][m] & IC_STATE_MASK) >> 6)
    return E_N3;
  return IC_N1;
}

static void path_ic_state_write(GthDPtables *dpm, unsigned int n,
                                unsigned int m, GthPath e)
{
  if (e == E_N3)
    dpm->core.path[n][m] |= (1 << 6);
}

static int dp_table_core_init(DPtablecore *core, unsigned long gen_dp_length,
                              unsigned long ref_dp_length,
                              unsigned long autoicmaxmatrixsize,
                              bool introncutout, GthJumpTable *jump_table,
                              GthStat *stat)
{
  unsigned long matrixsize, t, n,
                sizeofpathtype = sizeof (GthPath);

  /* XXX: adjust this check for QUARTER_MATRIX case */
  if (PROTEIN_NUMOFSTATES * sizeofpathtype * (gen_dp_length + 1) >=
      (~0)/(ref_dp_length + 1)) {
    /* in this case the matrix would be larger than the addressable memory
       of this machine -> return ERROR_MATRIX_ALLOCATION_FAILED */
    return GTH_ERROR_MATRIX_ALLOCATION_FAILED;
  }

  matrixsize = gt_safe_mult_ulong(gen_dp_length + 1, ref_dp_length + 1);

  if (!introncutout && autoicmaxmatrixsize > 0) {
    /* in this case the automatic intron cutout technique is enabled
       check if allocated matrix would be larger as specified maximal
       matrix size. If so, return matrix allocation error */
    if (sizeofpathtype * matrixsize * PROTEIN_NUMOFSTATES >
        autoicmaxmatrixsize << 20) {
      return GTH_ERROR_MATRIX_ALLOCATION_FAILED;
    }
  }

  /* set everything to NULL */
  for  (t =  E_STATE; t < PROTEIN_NUMOFSTATES; t++) {
    for (n = 0; n < PROTEIN_NUMOFSCORETABLES; n++)
      core->score[t][n] = NULL;
  }
  core->path = NULL;

  /* allocating space for core->score and core->path */
  for (t = E_STATE; t < PROTEIN_NUMOFSTATES; t++) {
    for (n = 0; n < PROTEIN_NUMOFSCORETABLES; n++) {
      core->score[t][n] = gt_malloc(sizeof (GthFlt) * (ref_dp_length + 1));
    }
  }

  if (jump_table) {
    gth_array2dim_plain_calloc(core->path, gen_dp_length + 1,
                               ref_dp_length + 1);
  }
  else {
    gth_array2dim_plain_malloc(core->path, gen_dp_length + 1,
                               ref_dp_length + 1);
  }
  if (!core->path) {
    /* matrix allocation failed, return after free of allocated tables */
    dp_table_core_free(core);
    return GTH_ERROR_MATRIX_ALLOCATION_FAILED;
  }

  /* statistics */
  gth_stat_increment_numofbacktracematrixallocations(stat);
  gth_stat_increase_totalsizeofbacktracematricesinMB(stat,
                                           (sizeofpathtype * matrixsize) >> 20);

  return 0;
}

static const char* showretracenames(Retrace retrace)
{
  gt_assert(retrace <= NUMOFRETRACE);
  return retracenames[retrace];
}

/* the following function initalizes the input structure */
static void initinput(GthAlignInputProtein *input,
                      const unsigned char *ref_seq_orig,
                      GthInput *gth_input)
{
  input->ref_seq_orig       = ref_seq_orig;
  input->score_matrix       = gth_input_score_matrix(gth_input);
  input->score_matrix_alpha = gth_input_score_matrix_alpha(gth_input);
}

/* the following function allocates space for the DP tables for proteins */
static int dp_tables_alloc(GthDPtables *dpm, unsigned long gen_dp_length,
                           bool proteinexonpenal, unsigned long ref_dp_length,
                           unsigned long autoicmaxmatrixsize, bool introncutout,
                           GthJumpTable *jump_table, GthStat *stat)
{
  unsigned long n;
  int rval;

  /* allocate core */
  if ((rval = dp_table_core_init(&dpm->core, gen_dp_length, ref_dp_length,
                                 autoicmaxmatrixsize, introncutout, jump_table,
                                 stat))) {
    return rval;
  }

  /* allocating space for intronstart and splitcodon arrays */
  for (n = 0; n < PROTEIN_NUMOFSCORETABLES; n++) {
    dpm->intronstart_A[n] = gt_malloc(sizeof (unsigned long) *
                                      (ref_dp_length + 1));
    dpm->intronstart_B[n] = gt_malloc(sizeof (unsigned long) *
                                      (ref_dp_length + 1));
    dpm->intronstart_C[n] = gt_malloc(sizeof (unsigned long) *
                                      (ref_dp_length + 1));
    if (proteinexonpenal) {
      dpm->exonstart[n] = gt_malloc(sizeof (unsigned long) *
                                    (ref_dp_length + 1));
    }
    else
      dpm->exonstart[n] = NULL;
    dpm->splitcodon_B[n] = gt_malloc(sizeof (unsigned char) *
                                     (ref_dp_length + 1));
    dpm->splitcodon_C1[n] = gt_malloc(sizeof (unsigned char) *
                                      (ref_dp_length + 1));
    dpm->splitcodon_C2[n] = gt_malloc(sizeof (unsigned char) *
                                      (ref_dp_length + 1));
  }

  return 0;
}

/* the following function initializes the DP tables for proteins */
static void dp_tables_init(GthDPtables *dpm, bool proteinexonpenal,
                           unsigned long ref_dp_length)
{
  unsigned long n, m;

  /* initialize the DP matrices */
  for (n = 0; n <  PROTEIN_NUMOFSCORETABLES; n++) {
    SCORE(E_STATE,n,0)  = (GthFlt) 0.0;
    path_e_state_write(dpm, n, 0, E_N1);
    SCORE(IA_STATE,n,0) = (GthFlt) 0.0;
    path_ia_state_write(dpm, n, 0, IA_N1);
    SCORE(IB_STATE,n,0) = (GthFlt) 0.0;
    path_ib_state_write(dpm, n, 0, IB_N1);
    SCORE(IC_STATE,n,0) = (GthFlt) 0.0;
    path_ic_state_write(dpm, n, 0, IC_N1);

    for (m = 1; m <= ref_dp_length; m++) {
      SCORE(E_STATE,n,m) = (GthFlt) 0.0;
      path_e_state_write(dpm, n, m, E_M);
      /* disallow intron status for 5' non-matching cDNA letters: */
      SCORE(IA_STATE,n,m) = (GthFlt) GTH_MINUSINFINITY;
      path_ia_state_write(dpm, n, m, IA_N1);
      SCORE(IB_STATE,n,m) = (GthFlt) GTH_MINUSINFINITY;
      path_ib_state_write(dpm, n, m, IB_N1);
      SCORE(IC_STATE,n,m) = (GthFlt) GTH_MINUSINFINITY;
      path_ic_state_write(dpm, n, m, IC_N1);
    }
  }

  for (n = 0; n < PROTEIN_NUMOFSCORETABLES; n++) {
    for (m = 0; m <= ref_dp_length; m++) {
      /* XXX: replace by memset, in gthsahmtd, too */
      dpm->intronstart_A[n][m] = 0;
      dpm->intronstart_B[n][m] = 0;
      dpm->intronstart_C[n][m] = 0;
      if (proteinexonpenal)
        dpm->exonstart[n][m]   = 0;
    }

    /* set the splitcodon tables to "UNSET" values.
       this is needed for a check if a splitcodon has already been set.
       i.e., if an intron has already been introduced
       the splitcodon table "dpm->splitcodon_C2" needs not to be set, because
       it is always used in conjunction with table "dpm->splitcodon_C1". */
    memset(dpm->splitcodon_B[n], UNSET, sizeof (unsigned char) *
                                        (ref_dp_length + 1));
    memset(dpm->splitcodon_C1[n], UNSET, sizeof (unsigned char) *
                                         (ref_dp_length + 1));
  }

  path_e_state_write(dpm, 0, 0, E_N1M);
  path_ia_state_write(dpm, 0, 0, E_N1);
  path_ib_state_write(dpm, 0, 0, E_N2);
  path_ic_state_write(dpm, 0, 0, E_N3);
}

unsigned char gthgetcodon(unsigned char genomicchar1,
                          unsigned char genomicchar2,
                          unsigned char genomicchar3,
                          const GtUchar *gen_alphabet_characters,
                          const GtTransTable *transtable)
{
  char codon;
  GT_UNUSED int rval;

  /* translate dna into codon */
  rval = gt_trans_table_translate_codon(transtable,
                                        gen_alphabet_characters[genomicchar1],
                                        gen_alphabet_characters[genomicchar2],
                                        gen_alphabet_characters[genomicchar3],
                                        &codon, NULL);
  /* since the genomic sequence has been preprocessed before, the codon
     translation should not fail */
  gt_assert(!rval);

  return codon;
}

/* the following function evaluate the dynamic programming tables */
static void complete_path_matrix(GthDPtables *dpm, GthAlignInputProtein *input,
                                 bool proteinexonpenal,
                                 const unsigned char *gen_seq_tran,
                                 unsigned long gen_dp_length,
                                 unsigned long ref_dp_length,
                                 GthDPParam *dp_param,
                                 GthDPOptionsCore *dp_options_core,
                                 GthDPScoresProtein *dp_scores_protein)
{
  unsigned long n, m, modn, modnminus1, modnminus2, modnminus3;
  unsigned char origreferencechar;
  GthFlt value, maxvalue;
  GthPath retrace;

  /* stepping along the genomic sequence */
  for (n = GENOMICDPSTART; n <= gen_dp_length; n++) {
    modn       = GT_MOD4(n),
    modnminus1 = GT_MOD4(n-1),
    modnminus2 = GT_MOD4(n-2),
    modnminus3 = GT_MOD4(n-3);

    path_e_state_write(dpm, n, 0, E_N1);
    path_ia_state_write(dpm, n, 0, IA_N1);
    path_ib_state_write(dpm, n, 0, IB_N1);
    path_ic_state_write(dpm, n, 0, IC_N1);

    /* stepping along the protein sequence */
    for (m = REFERENCEDPSTART; m <= ref_dp_length; m++) {
      origreferencechar = input->ref_seq_orig[m-1];

      /* evaluate E_nm */
      /* 0. */
      maxvalue = SCORE(E_STATE, modnminus3, m-1) +
                 /* XXX: why is here no extra condition? */
                 (dp_param->log_1minusPdonor[n-3] +
                  GTHGETSCORE(dp_scores_protein, gen_seq_tran[n-3],
                              gen_seq_tran[n-2], gen_seq_tran[n-1],
                              origreferencechar));
      retrace  = (GthPath) E_N3M;

      /* 1. */
      value = SCORE(E_STATE, modnminus2, m-1);
      if (n < gen_dp_length || m < WSIZE_PROTEIN) {
        value += dp_param->log_1minusPdonor[n-2] +
                 GTHGETSCORE(dp_scores_protein, gen_seq_tran[n-2],
                             gen_seq_tran[n-1], DASH, origreferencechar);
      }
      UPDATEMAX(E_N2M);

      /* 2. */
      value = SCORE(E_STATE, modnminus1, m-1);
      if (n < gen_dp_length || m < WSIZE_PROTEIN) {
        value += dp_param->log_1minusPdonor[n-1] +
                 GTHGETSCORE(dp_scores_protein, gen_seq_tran[n-1], DASH, DASH,
                             origreferencechar);
      }
      UPDATEMAX(E_N1M);

      /* 3. */
      value = SCORE(E_STATE, modn, m-1);
      if (n < gen_dp_length || m < WSIZE_PROTEIN) {
        if (n == gen_dp_length) {
          /* in this case the value used in the 'else' branch below is not
             defined. */
          value += dp_param->log_1minusPdonor[n-1];
        }
        else {
          value += dp_param->log_1minusPdonor[n];
                                     /* XXX: ^^^  why n? */
        }
        value += GTHGETSCORE(dp_scores_protein, DASH, DASH, DASH,
                             origreferencechar);
      }
      UPDATEMAX(E_M);

      /* 4. */
      value = SCORE(E_STATE, modnminus3, m);
      if (m < ref_dp_length || n < WSIZE_DNA) {
        value += dp_param->log_1minusPdonor[n-3] +
                 GTHGETSCORE(dp_scores_protein, gen_seq_tran[n-3],
                             gen_seq_tran[n-2], gen_seq_tran[n-1], DASH);
      }
      UPDATEMAX(E_N3);

      /* 5. */
      value = SCORE(E_STATE, modnminus2, m);
      if (m < ref_dp_length || n < WSIZE_DNA) {
        value += dp_param->log_1minusPdonor[n-2] +
                 GTHGETSCORE(dp_scores_protein, gen_seq_tran[n-2],
                             gen_seq_tran[n-1], DASH, DASH);
      }
      UPDATEMAX(E_N2);

      /* 6. */
      value = SCORE(E_STATE, modnminus1, m);
      if (m < ref_dp_length || n < WSIZE_DNA) {
        value += dp_param->log_1minusPdonor[n-1] +
                 GTHGETSCORE(dp_scores_protein, gen_seq_tran[n-1], DASH, DASH,
                             DASH);
      }
      UPDATEMAX(E_N1);

      /* 7. */
      value = SCORE(IA_STATE, modnminus3, m-1);
      if (n > GENOMICDPSTART) /* the value below is only defined in this case */
        value += dp_param->log_Pacceptor[n-4];
      value += GTHGETSCORE(dp_scores_protein, gen_seq_tran[n-3],
                           gen_seq_tran[n-2], gen_seq_tran[n-1],
                           origreferencechar);
      if (n - 2 - dpm->intronstart_A[modnminus3][m-1] <
          dp_options_core->dpminintronlength) {
        value -= dp_options_core->shortintronpenalty;
      }
      UPDATEMAX(IA_N3M);

      /* 8.
         this recurrence is only used if an intron has already been introduced.
         (in this case "dpm->splitcodon_B[modnminus1][m-1]" is different from
         "UNSET". */
      if (dpm->splitcodon_B[modnminus1][m-1] != (unsigned char) UNSET) {
        value = SCORE(IB_STATE, modnminus2, m-1) +
                dp_param->log_Pacceptor[n-3] +
                GTHGETSCORE(dp_scores_protein,
                            dpm->splitcodon_B[modnminus2][m-1],
                            gen_seq_tran[n-2], gen_seq_tran[n-1],
                            origreferencechar);
        if (n - 1 - dpm->intronstart_B[modnminus2][m-1] <
            dp_options_core->dpminintronlength) {
          value -= dp_options_core->shortintronpenalty;
        }
        UPDATEMAX(IB_N2M);
      }

      /* 9.
         explanation for check see above.
         "dpm->splitcodon_C2[modnminus1][m-1]" needs not to be checked,
         because it is always set in conjunction with
         "dpm->splitcodon_C1[modnminus1][m-1]". */
      if (dpm->splitcodon_C1[modnminus1][m-1] != (unsigned char) UNSET) {
        value = SCORE(IC_STATE, modnminus1, m-1) +
                dp_param->log_Pacceptor[n-2] +
                GTHGETSCORE(dp_scores_protein,
                            dpm->splitcodon_C1[modnminus1][m-1],
                            dpm->splitcodon_C2[modnminus1][m-1],
                            gen_seq_tran[n-1], origreferencechar);
        if (n - dpm->intronstart_C[modnminus1][m-1] <
            dp_options_core->dpminintronlength) {
          value -= dp_options_core->shortintronpenalty;
        }
        UPDATEMAX(IC_N1M);
      }

      /* save maximum values */
      SCORE(E_STATE, modn, m) = maxvalue;
      path_e_state_write(dpm, n, m, retrace);

      if (proteinexonpenal) {
        switch (retrace) {
          case E_N3M:
            dpm->exonstart[modn][m] = dpm->exonstart[modnminus3][m-1];
            break;
          case E_N2M:
            dpm->exonstart[modn][m] = dpm->exonstart[modnminus2][m-1];
            break;
          case E_N1M:
            dpm->exonstart[modn][m] = dpm->exonstart[modnminus1][m-1];
            break;
          case E_M:
            dpm->exonstart[modn][m] = dpm->exonstart[modn][m-1];
            break;
          case E_N3:
            dpm->exonstart[modn][m] = dpm->exonstart[modnminus3][m];
            break;
          case E_N2:
            dpm->exonstart[modn][m] = dpm->exonstart[modnminus2][m];
            break;
          case E_N1:
            dpm->exonstart[modn][m] = dpm->exonstart[modnminus1][m];
            break;
          case IA_N3M:
          case IB_N2M:
          case IC_N1M:
            dpm->exonstart[modn][m] = n;
            break;
          case IC_N1:
            dpm->exonstart[modn][m] = n;
            break;
          default: gt_assert(0);
        }
      }

      /* evaluate IA_nm */
      maxvalue = SCORE(IA_STATE, modnminus1, m);
      if (!dp_options_core->freeintrontrans)
        maxvalue += dp_param->log_1minusPacceptor[n-2];
      retrace  = (GthPath) IA_N1;

      value = SCORE(E_STATE, modnminus1, m) + dp_param->log_Pdonor[n-1];
      if (proteinexonpenal) {
        if (n - dpm->exonstart[modnminus1][m] <
            dp_options_core->dpminexonlength) {
          value -= dp_options_core->shortexonpenalty;
        }
      }
      UPDATEMAX(E_N1);

      /* save maximum values */
      SCORE(IA_STATE, modn, m) = maxvalue;
      path_ia_state_write(dpm, n, m, retrace);

      switch (retrace) {
        case IA_N1:
          dpm->intronstart_A[modn][m] = dpm->intronstart_A[modnminus1][m];
          break;
        case E_N1:
          dpm->intronstart_A[modn][m] = n;
          break;
        default: gt_assert(0);
      }

      /* evaluate IB_nm */
      maxvalue = SCORE(IB_STATE, modnminus1, m);
      if (!dp_options_core->freeintrontrans)
        maxvalue += dp_param->log_1minusPacceptor[n-2];
      retrace  = (GthPath) IB_N1;

      value = SCORE(E_STATE, modnminus2, m) + dp_param->log_Pdonor[n-1];
      if (proteinexonpenal) {
        if (n - 1 - dpm->exonstart[modnminus2][m] <
            dp_options_core->dpminexonlength) {
          value -= dp_options_core->shortexonpenalty;
        }
      }
      UPDATEMAX(E_N2);

      /* save maximum values */
      SCORE(IB_STATE, modn, m) = maxvalue;
      path_ib_state_write(dpm, n, m, retrace);

      switch (retrace) {
        case(IB_N1):
          dpm->intronstart_B[modn][m] = dpm->intronstart_B[modnminus1][m];
          dpm->splitcodon_B[modn][m]  = dpm->splitcodon_B[modnminus1][m];
          break;
        case(E_N2):
          dpm->intronstart_B[modn][m] = n;
          dpm->splitcodon_B[modn][m]  = gen_seq_tran[n-2];
          break;
        default: gt_assert(0);
      }

      /* evaluate IC_nm */
      maxvalue = SCORE(IC_STATE, modnminus1, m);
      if (!dp_options_core->freeintrontrans)
        maxvalue += dp_param->log_1minusPacceptor[n-2];
      retrace = (GthPath) IC_N1;

      value = SCORE(E_STATE, modnminus3, m) + dp_param->log_Pdonor[n-1];
      if (proteinexonpenal) {
        if (n - 2 - dpm->exonstart[modnminus3][m] <
            dp_options_core->dpminexonlength) {
          value -= dp_options_core->shortexonpenalty;
        }
      }
      UPDATEMAX(E_N3);

      /* save maximum values */
      SCORE(IC_STATE, modn, m) = maxvalue;
      path_ic_state_write(dpm, n, m, retrace);

      switch (retrace) {
        case(IC_N1):
          dpm->intronstart_C[modn][m] = dpm->intronstart_C[modnminus1][m];
          dpm->splitcodon_C1[modn][m] = dpm->splitcodon_C1[modnminus1][m];
          dpm->splitcodon_C2[modn][m] = dpm->splitcodon_C2[modnminus1][m];
          break;
        case(E_N3):
          dpm->intronstart_C[modn][m] = n;
          dpm->splitcodon_C1[modn][m] = gen_seq_tran[n-3];
          dpm->splitcodon_C2[modn][m] = gen_seq_tran[n-2];
          break;
        default: gt_assert(0);
      }
    }
  }
}

static void include_exon(GthBacktracePath *backtrace_path,
                         unsigned long exonlength)
{
  unsigned long i,
                divresult = exonlength / GT_CODON_LENGTH,
                modresult = exonlength % GT_CODON_LENGTH;

  /* at least one editoperation already saved */
  gt_assert(gth_backtrace_path_length(backtrace_path));

  for (i = 0; i < divresult; i++)
    gth_backtrace_path_add_deletion(backtrace_path);

  switch (modresult) {
    case 0:
      /* nothing to do. */
      break;
    case 1:
      gth_backtrace_path_add_deletion_with_2_gaps(backtrace_path);
      break;
    case 2:
      gth_backtrace_path_add_deletion_with_1_gap(backtrace_path);
      break;
    default: gt_assert(0);
  }
}

static void include_intron(GthBacktracePath *backtrace_path, GthPath pathtype,
                           unsigned long intronlength)
{
  /* at least one editoperation already saved */
  gt_assert(gth_backtrace_path_length(backtrace_path));

  while (intronlength > 0) {
    switch (pathtype) {
      case IA_N1:
        gth_backtrace_path_add_intron(backtrace_path);
        break;
      case IB_N1:
        gth_backtrace_path_add_intron_with_1_base_left(backtrace_path);
        break;
      case IC_N1:
        gth_backtrace_path_add_intron_with_2_bases_left(backtrace_path);
        break;
      default: gt_assert(0);
    }
    intronlength--;
  }
}

/*
  This typedef is used in the backtracing for protein spliced alignments to
  ensure that introns which are not in phase are always flanked by edit
  operations of length 1. This is necessary for a correct funtioning of the
  STRICT cutoff functions.
*/
typedef enum {
  DUMMY_STATUS_UNDEFINED = 0,
  DUMMY_JUST_SET,
  FIRST_EOP_AFTER_DUMMY_SET,
  ENSURE_SINGLE_MATCH,
  NUMOFDUMMYSTATUSES
} Dummystatus;

static int evaltracepath(GthBacktracePath *backtrace_path, GthDPtables *dpm,
                         unsigned long ref_dp_length,
                         const unsigned char *gen_seq_tran,
                         unsigned long gen_dp_length, States actualstate,
                         bool introncutout, GthSplicedSeq *spliced_seq,
                         const GtTransTable *transtable, bool comments,
                         bool noicinintroncheck, GtAlphabet *gen_alphabet,
                         const unsigned char *ref_seq_orig,
                         GtFile *outfp)
{
  unsigned long genptr          = gen_dp_length, last_genptr = 0,
                genptr_tail     = gen_dp_length,
                refptr          = ref_dp_length,
                GT_UNUSED dummy_index     = GT_UNDEF_ULONG,
                dummy_b2_genptr = GT_UNDEF_ULONG,
                dummy_b3_genptr = GT_UNDEF_ULONG,
                dummy_c3_genptr = GT_UNDEF_ULONG;
#ifndef NDEBUG
  unsigned long numofincludedintrons = 0;
#endif
  GthPath pathtype = (GthPath) 0;
  Dummystatus dummystatus = DUMMY_STATUS_UNDEFINED;
  unsigned char codon;
  bool skipdummyprocessing = false;
  const GtUchar *gen_alphabet_characters = gt_alphabet_characters(gen_alphabet);

  gt_assert(!gth_backtrace_path_length(backtrace_path));

  /* skip the first processing of a dummy if starts with an intron state
     (to make sure the assertions don't fail) */
  if (actualstate == IB_STATE || actualstate == IC_STATE)
    skipdummyprocessing = true;

  while ((genptr > 0) || (refptr > 0)) {
    switch (actualstate) {
      case E_STATE:
        pathtype = path_e_state_read(dpm, genptr, refptr);
        break;
      case IA_STATE:
        pathtype = path_ia_state_read(dpm, genptr, refptr);
        break;
      case IB_STATE:
        pathtype = path_ib_state_read(dpm, genptr, refptr);
        break;
      case IC_STATE:
        pathtype = path_ic_state_read(dpm, genptr, refptr);
        break;
      default: gt_assert(0);
    }

    switch (dummystatus) {
      case DUMMY_STATUS_UNDEFINED:
        /* nothing to do */
        break;
      case DUMMY_JUST_SET:
        dummystatus = FIRST_EOP_AFTER_DUMMY_SET;
        break;
      case FIRST_EOP_AFTER_DUMMY_SET:
        dummystatus = ENSURE_SINGLE_MATCH;
        break;
      case ENSURE_SINGLE_MATCH:
        dummystatus = DUMMY_STATUS_UNDEFINED;
        break;
      default: gt_assert(0);
    }

    /* if we start with an intron with one or two bases left this ensures that
       the following match has length one */
    if (skipdummyprocessing)
      dummystatus = DUMMY_JUST_SET;

    if (introncutout) {
      while (genptr_tail >= genptr && genptr_tail > 0) {
        if (gth_spliced_seq_pos_is_border(spliced_seq, genptr_tail - 1)) {
          /* ensure that introns are only included into already existing introns
          */
          if ((pathtype != (GthPath) IA_N1 &&
               pathtype != (GthPath) IB_N1 &&
               pathtype != (GthPath) IC_N1) ||
               !gth_backtrace_path_last_is_intron(backtrace_path)) {
            if (noicinintroncheck) {
              /* intron cutout in intron check disabled,
                 include exon consisting of deletions */
              if (genptr != last_genptr) {
                last_genptr = genptr;
                include_exon(backtrace_path,
                             gth_spliced_seq_border_length(spliced_seq,
                                                           genptr_tail - 1));
              }
            }
            else {
              if (comments) {
                gt_file_xprintf(outfp, "%c abort backtracing, intron cutout "
                                   "at p=%s (genpos=%lu (actual strand!))\n",
                                   COMMENTCHAR,
                                   showretracenames((Retrace) pathtype),
                                   spliced_seq->positionmapping[genptr]);
              }
              return GTH_ERROR_CUTOUT_NOT_IN_INTRON;
            }
          }
          else {
            /* include intron */
            include_intron(backtrace_path, pathtype,
                           gth_spliced_seq_border_length(spliced_seq,
                                                         genptr_tail - 1));
#ifndef NDEBUG
            numofincludedintrons++;
#endif
          }
        }
        genptr_tail--;
      }
    }

    switch (actualstate) {
      case E_STATE:
        /* we are currently in an exon state, the possible pathtypes for exon
           states are handled here */
        switch (pathtype) {
          case E_N3M:
            codon = gthgetcodon(gen_seq_tran[genptr - 3],
                                gen_seq_tran[genptr - 2],
                                gen_seq_tran[genptr - 1],
                                gen_alphabet_characters,
                                transtable);

            if (codon == ref_seq_orig[refptr - 1]) {
              gth_backtrace_path_add_match(backtrace_path,
                                           dummystatus == ENSURE_SINGLE_MATCH);
            }
            else
              gth_backtrace_path_add_mismatch(backtrace_path);
            genptr-= 3;
            refptr--;
            break;
          case E_N2M:
            gth_backtrace_path_add_mismatch_with_1_gap(backtrace_path);
            genptr-= 2;
            refptr--;
            break;
          case E_N1M:
            gth_backtrace_path_add_mismatch_with_2_gaps(backtrace_path);
            if (dummystatus == FIRST_EOP_AFTER_DUMMY_SET)
              dummystatus = DUMMY_JUST_SET;
            genptr--;
            refptr--;
            break;
          case E_M:
            gth_backtrace_path_add_insertion(backtrace_path);
            if (dummystatus == FIRST_EOP_AFTER_DUMMY_SET)
              dummystatus = DUMMY_JUST_SET;
            refptr--;
            break;
          case E_N3:
            gth_backtrace_path_add_deletion(backtrace_path);
            if (dummystatus == FIRST_EOP_AFTER_DUMMY_SET)
              dummystatus = DUMMY_JUST_SET;
            genptr-= 3;
            break;
          case E_N2:
            gth_backtrace_path_add_deletion_with_1_gap(backtrace_path);
            if (dummystatus == FIRST_EOP_AFTER_DUMMY_SET)
              dummystatus = DUMMY_JUST_SET;
            genptr-= 2;
            break;
          case E_N1:
            gth_backtrace_path_add_deletion_with_2_gaps(backtrace_path);
            if (dummystatus == FIRST_EOP_AFTER_DUMMY_SET)
              dummystatus = DUMMY_JUST_SET;
            genptr--;
            break;
          case IA_N3M:
            codon = gthgetcodon(gen_seq_tran[genptr - 3],
                                gen_seq_tran[genptr - 2],
                                gen_seq_tran[genptr - 1],
                                gen_alphabet_characters,
                                transtable);
            if (codon == ref_seq_orig[refptr - 1]) {
              gth_backtrace_path_add_match(backtrace_path,
                                           dummystatus == ENSURE_SINGLE_MATCH);
            }
            else
              gth_backtrace_path_add_mismatch(backtrace_path);
            genptr-= 3;
            refptr--;
            actualstate = IA_STATE;
            break;
          case IB_N2M:
            /* add dummy */
            gth_backtrace_path_add_dummy(backtrace_path);
            gt_assert(dummy_b2_genptr == GT_UNDEF_ULONG);
            gt_assert(dummy_b3_genptr == GT_UNDEF_ULONG);

            dummystatus     = DUMMY_STATUS_UNDEFINED;
            dummy_b3_genptr = genptr - 1;
            dummy_b2_genptr = genptr - 2;

            genptr-= 2;
            refptr--;
            actualstate = IB_STATE;
            break;
          case IC_N1M:
            /* add dummy */
            gth_backtrace_path_add_dummy(backtrace_path);
            gt_assert(dummy_c3_genptr == GT_UNDEF_ULONG);

            dummystatus     = DUMMY_STATUS_UNDEFINED;
            dummy_c3_genptr = genptr - 1;

            genptr--;
            refptr--;
            actualstate = IC_STATE;
            break;
          default: gt_assert(0);
        }
        break;
      case IA_STATE:
        /* we are currently in an intron state ('A' type), the possible
           pathtypes for such intron states are handled here */
        switch (pathtype) {
          case IA_N1:
            gth_backtrace_path_add_intron(backtrace_path);
            genptr--;
            break;
          case E_N1:
            gth_backtrace_path_add_intron(backtrace_path);
            genptr--;
            actualstate = E_STATE;
            break;
          default: gt_assert(0);
        }
        break;
      case IB_STATE:
        /* we are currently in an intron state ('B' type), the possible
           pathtypes for such intron states are handled here */
        switch (pathtype) {
          case IB_N1:
            gth_backtrace_path_add_intron_with_1_base_left(backtrace_path);
            genptr--;
            break;
          case E_N2:
            gth_backtrace_path_add_intron_with_1_base_left(backtrace_path);
            if (skipdummyprocessing)
              skipdummyprocessing = false;
            else {
              /* resolve dummy */
              gt_assert(dummy_b2_genptr != GT_UNDEF_ULONG);
              gt_assert(dummy_b3_genptr != GT_UNDEF_ULONG);

              codon = gthgetcodon(gen_seq_tran[genptr - 2],
                                  gen_seq_tran[dummy_b2_genptr],
                                  gen_seq_tran[dummy_b3_genptr],
                                  gen_alphabet_characters, transtable);
              gth_backtrace_path_set_dummy(backtrace_path, codon
                                       == ref_seq_orig[refptr]);
              dummystatus     = DUMMY_JUST_SET;
              dummy_index     = GT_UNDEF_ULONG;
              dummy_b2_genptr = GT_UNDEF_ULONG;
              dummy_b3_genptr = GT_UNDEF_ULONG;
            }
            genptr-= 2;
            actualstate = E_STATE;
            break;
          default: gt_assert(0);
        }
        break;
      case IC_STATE:
        /* we are currently in an intron state ('C' type), the possible
           pathtypes for such intron states are handled here */
        switch (pathtype) {
          case IC_N1:
            gth_backtrace_path_add_intron_with_2_bases_left(backtrace_path);
            genptr--;
            break;
          case E_N3:
            gth_backtrace_path_add_intron_with_2_bases_left(backtrace_path);
            if (skipdummyprocessing)
              skipdummyprocessing = false;
            else {
              /* resolve dummy here */
              gt_assert(dummy_c3_genptr != GT_UNDEF_ULONG);

              codon = gthgetcodon(gen_seq_tran[genptr - 3],
                                  gen_seq_tran[genptr - 2],
                                  gen_seq_tran[dummy_c3_genptr],
                                  gen_alphabet_characters, transtable);
              gth_backtrace_path_set_dummy(backtrace_path, codon
                                       == ref_seq_orig[refptr]);
              dummystatus     = DUMMY_JUST_SET;
              dummy_index     = GT_UNDEF_ULONG;
              dummy_c3_genptr = GT_UNDEF_ULONG;
            }
            genptr-= 3;
            actualstate = E_STATE;
            break;
          default: gt_assert(0);
        }
        break;
      default: gt_assert(0);
    }
  }

  gt_assert(genptr == 0 && refptr == 0);
  gt_assert(!gth_backtrace_path_contain_dummy(backtrace_path));

#ifndef NDEBUG
  if (introncutout) {
    /* number of borders equals number of included introns */
    gt_assert(numofincludedintrons ==
              gth_spliced_seq_num_of_borders(spliced_seq));
  }
#endif

  return 0;
}

static int find_optimal_path(GthBacktracePath *backtrace_path, GthDPtables *dpm,
                             unsigned long ref_dp_length,
                             const unsigned char *gen_seq_tran,
                             unsigned long gen_dp_length, bool introncutout,
                             GthSplicedSeq *spliced_seq,
                             const GtTransTable *transtable, bool comments,
                             bool noicintroncheck, GtAlphabet *gen_alphabet,
                             const unsigned char *ref_seq_orig,
                             GtFile *outfp)
{
  int rval;
  GthFlt value, maxvalue;
  GthPath retrace;
  States state;

  maxvalue = SCORE(E_STATE, GT_MOD4(gen_dp_length), ref_dp_length);
  retrace  = (GthPath) E_STATE;

  for (state = (States) 1; state < PROTEIN_NUMOFSTATES; state++) {
    value = SCORE(state, GT_MOD4(gen_dp_length), ref_dp_length);
    UPDATEMAX(state);
  }

  if ((rval = evaltracepath(backtrace_path, dpm, ref_dp_length,
                            gen_seq_tran, gen_dp_length, (States) retrace,
                            introncutout, spliced_seq, transtable,
                            comments, noicintroncheck, gen_alphabet,
                            ref_seq_orig, outfp))) {
    return rval;
  }

  /* sum of backtrace_path equals ref_dp_length */
  gt_assert(gth_backtrace_path_is_valid(backtrace_path));

  return 0;
}

/* the following function frees the DP tables for proteins */
static void dp_tables_free(GthDPtables *dpm)
{
  unsigned long n;

  /* free core */
  dp_table_core_free(&dpm->core);

  /* freespace for intronstart and splitcodon arrays */
  for (n = 0; n < PROTEIN_NUMOFSCORETABLES; n++) {
    gt_free(dpm->intronstart_A[n]);
    gt_free(dpm->intronstart_B[n]);
    gt_free(dpm->intronstart_C[n]);
    gt_free(dpm->exonstart[n]);
    gt_free(dpm->splitcodon_B[n]);
    gt_free(dpm->splitcodon_C1[n]);
    gt_free(dpm->splitcodon_C2[n]);
  }
}

int gth_align_protein(GthSA *sa,
                      GtArray *gen_ranges,
                      const unsigned char *gen_seq_tran,
                      const unsigned char *ref_seq_tran,
                      const unsigned char *ref_seq_orig,
                      unsigned long ref_dp_length,
                      GtAlphabet *gen_alphabet,
                      GtAlphabet *ref_alphabet,
                      GthInput *gth_input,
                      bool introncutout,
                      unsigned long autoicmaxmatrixsize,
                      bool proteinexonpenal,
                      bool showeops,
                      bool comments,
                      bool gs2out,
                      unsigned long translationtable,
                      const GtRange *gen_seq_bounds,
                      GthSpliceSiteModel *splice_site_model,
                      GthDPOptionsCore *dp_options_core,
                      GthDPOptionsPostpro *dp_options_postpro,
                      GT_UNUSED GthProteinCompletePathMatrixJT
                      complete_path_matrix_jt,
                      GthJumpTable *jump_table,
                      GT_UNUSED unsigned long ref_offset,
                      GthStat *stat,
                      GtFile *outfp)
{
  unsigned long gen_dp_start, gen_dp_end, gen_dp_length;
  GthDPScoresProtein *dp_scores_protein;
  GthDPParam *dp_param;
  GthSplicedSeq *spliced_seq = NULL;
  GthAlignInputProtein input;
  GtTransTable *transtable;
  GthDPtables dpm;
  int rval;

  gt_assert(gen_ranges);

  dp_scores_protein =
    gth_dp_scores_protein_new(translationtable,
                              gth_input_score_matrix(gth_input),
                              gth_input_score_matrix_alpha(gth_input));

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
  initinput(&input, ref_seq_orig, gth_input);
  if ((rval = dp_tables_alloc(&dpm, introncutout ? spliced_seq->splicedseqlen
                                                 : gen_dp_length,
                              proteinexonpenal, ref_dp_length,
                              autoicmaxmatrixsize, introncutout, jump_table,
                              stat))) {
    gth_dp_param_delete(dp_param);
    gth_spliced_seq_delete(spliced_seq);
    gth_dp_scores_protein_delete(dp_scores_protein);
    return rval;
  }
  dp_tables_init(&dpm, proteinexonpenal, ref_dp_length);
  gth_sa_set(sa, PROTEIN_ALPHA, gen_dp_start, gen_dp_length);

  transtable = gt_trans_table_new(translationtable, NULL);
  /* XXX: the validity of the translation table has to be checked before */
  gt_assert(!rval);

  /* calculation */
  /* XXX */
#if 0
  if (jump_table) {
    gt_assert(complete_path_matrix_jt);
    complete_path_matrix_jt(&dpm, &input, proteinexonpenal,
                            introncutout ? spliced_seq->splicedseq
                                         : gen_seq_tran + gen_dp_start,
                            introncutout ? spliced_seq->splicedseqlen
                                         : gen_dp_length,
                            ref_dp_length, dp_param, dp_options_core,
                            dp_scores_protein, jump_table, gen_ranges,
                            ref_offset);
  }
  else {
#endif
    complete_path_matrix(&dpm, &input, proteinexonpenal,
                         introncutout ? spliced_seq->splicedseq
                                      : gen_seq_tran + gen_dp_start,
                         introncutout ? spliced_seq->splicedseqlen
                                      : gen_dp_length,
                         ref_dp_length, dp_param, dp_options_core,
                         dp_scores_protein);

  /* backtracing */
  if ((rval = find_optimal_path(gth_sa_backtrace_path(sa), &dpm, ref_dp_length,
                                introncutout ? spliced_seq->splicedseq
                                             : gen_seq_tran + gen_dp_start,
                                introncutout ? spliced_seq->splicedseqlen
                                             : gen_dp_length,
                                introncutout, spliced_seq, transtable,
                                comments, dp_options_core->noicinintroncheck,
                                gen_alphabet, input.ref_seq_orig, outfp))) {
    if (rval == GTH_ERROR_CUTOUT_NOT_IN_INTRON) {
      gt_trans_table_delete(transtable);
      dp_tables_free(&dpm);
      gth_dp_param_delete(dp_param);
      gth_spliced_seq_delete(spliced_seq);
      gth_dp_scores_protein_delete(dp_scores_protein);
    }
    return rval;
  }

  /* intron cutout is done after this point */
  if (showeops) {
    gt_file_xprintf(outfp, "showeops: ");
    gth_backtrace_path_show(gth_sa_backtrace_path(sa), false, 0, outfp);
  }

  /* determine cutoffs if not switched off by command line option */
  gth_sa_determine_cutoffs(sa, dp_options_postpro->leadcutoffsmode,
                               dp_options_postpro->termcutoffsmode,
                               dp_options_postpro->cutoffsminexonlen);

  /* remove zero base exons */
  gth_sa_remove_zero_base_exons(sa, stat);

  /* compute borders and scores */
  gth_compute_scores(sa, true, dp_param, NULL, gen_seq_tran + gen_dp_start,
                     ref_seq_tran, ref_seq_orig, transtable, gen_dp_start,
                     dp_options_postpro->scoreminexonlen, introncutout, gs2out,
                     spliced_seq, ref_dp_length, gen_alphabet, ref_alphabet,
                     dp_scores_protein);

  /* free */
  gt_trans_table_delete(transtable);
  dp_tables_free(&dpm);
  gth_dp_param_delete(dp_param);
  gth_spliced_seq_delete(spliced_seq);
  gth_dp_scores_protein_delete(dp_scores_protein);

  return 0;
}
