/*
  Copyright (c) 2003-2009 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2003-2005 Michael E Sparks <mespar1@iastate.edu>
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

#include "core/bioseq.h"
#include "core/fa.h"
#include "core/fileutils_api.h"
#include "core/ma_api.h"
#include "core/unused_api.h"
#include "core/xansi.h"
#include "gth/gtherror.h"
#include "gth/gthoutput.h"
#include "gth/gthprobdef.h"
#include "gth/bssm_param_hard_coded.h"
#include "gth/gthalignment.h"
#include "gth/gthspeciestab.h"
#include "gth/showbool.h"

/*
  This is a collection of functions associated with
  manipulating bssm_param objects for training purposes.

  Input data files are--strictly!--named as follows, using
  an obvious schema:
    F0_don  F1_don  F2_don  Fi_don  T0_don  T1_don  T2_don
    F0_acc  F1_acc  F2_acc  Fi_acc  T0_acc  T1_acc  T2_acc
  Phase is denoted as follows:
    1 -> C O D |
    2 -> C | O D
    0 -> C O | D
*/

#define PSEUDOPROB      0.05
#define NULLPROB        0.0
#define INITVAL_INT     0
#define MYVERSION       2  /* Version of BSSM param. Check with
                              BSSMPARAMVERSION in bssm_param.h */

#define NUMOFFILES (sizeof (filenames)/sizeof (filenames[0]))

#define BSSMENVNAME  "BSSMDIR"

static char *filenames[] =
{
  "T1",
  "T2",
  "T0",
  "F1",
  "F2",
  "F0",
  "Fi"
};

static int bssm_model_read(GthBssmModel *bssm_model, FILE *file, GtError *err)
{
  int had_err = 0;
  gt_error_check(err);
  gt_xfread(&bssm_model->hypothesisnum, sizeof (unsigned long), 1, file);
  if (bssm_model->hypothesisnum != HYPOTHESIS7 &&
      bssm_model->hypothesisnum != HYPOTHESIS2) {
    gt_error_set(err, "BSSM model contains unknown hypothesis number");
    had_err = -1;
  }
  if (!had_err) {
    gt_xfread(&bssm_model->windowsizeleft, sizeof (unsigned long), 1, file);
    gt_xfread(&bssm_model->windowsizeright, sizeof (unsigned long), 1, file);
    switch (bssm_model->hypothesisnum) {
      case HYPOTHESIS2:
        gt_xfread(&bssm_model->hypotables.hypo2table, sizeof (Hypo2table), 1,
                  file);
        break;
      case HYPOTHESIS7:
        gt_xfread(&bssm_model->hypotables.hypo7table, sizeof (Hypo7table), 1,
                  file);
        break;
      default: gt_assert(0);
    }
  }
  return had_err;
}

GthBssmParam* gth_bssm_param_new(void)
{
  return gt_calloc(1, sizeof (GthBssmParam));
}

GthBssmParam* gth_bssm_param_load(const char *filename, GtError *err)
{
  GthBssmParam *bssm_param = NULL;
  FILE *file = NULL;
  int had_err = 0;

  gt_error_check(err);

  if (gt_file_exists(filename))
    file = gt_xfopen(filename, "r");
  else {
    GtStr *path = gt_str_new();
    if (strchr(filename, '/')) {
      gt_error_set(err, "filename \"%s\" contains illegal symbol '/': the path "
                        "list specified by environment variable \"%s\" cannot "
                        "be searched for it", filename, BSSMENVNAME);
      had_err = -1;
    }
    if (!had_err)
      had_err = gt_file_find_in_env(path, filename, BSSMENVNAME, err);
    if (!had_err && !gt_str_length(path)) {
      gt_error_set(err, "file \"%s\" not found in directory list specified by "
                        "environment variable %s", filename, BSSMENVNAME);
      had_err = -1;
    }
    if (!had_err) {
      /* path found -> append filename */
      gt_str_append_char(path, '/');
      gt_str_append_cstr(path, filename);
      file = gt_xfopen(gt_str_get(path), "r");
    }
    gt_str_delete(path);
  }

  /* read version number and check if equals version number 2 */
  if (!had_err) {
    bssm_param = gt_malloc(sizeof *bssm_param);
    gt_xfread(&bssm_param->versionnum,  sizeof (unsigned char), 1, file);
    if (bssm_param->versionnum != (unsigned char) 2) {
      gt_error_set(err, "BSSM file %s has unrecognized version number %u\n",
                   filename, bssm_param->versionnum);
      had_err = -1;
    }
  }

  if (!had_err) {
    /* read in model variables */
    gt_xfread(&bssm_param->gt_donor_model_set,  sizeof (bool), 1, file);
    gt_xfread(&bssm_param->gc_donor_model_set,  sizeof (bool), 1, file);
    gt_xfread(&bssm_param->ag_acceptor_model_set,  sizeof (bool), 1, file);

    /* check if at least one model is set in file */
    if (!bssm_param->gt_donor_model_set &&
        !bssm_param->gc_donor_model_set &&
        !bssm_param->ag_acceptor_model_set) {
      gt_error_set(err, "BSSM file %s apparently contains no model", filename);
      had_err = -1;
    }
  }

  /* read GT donor site model */
  if (!had_err && bssm_param->gt_donor_model_set)
    had_err = bssm_model_read(&bssm_param->gt_donor_model, file, err);

  /* read GC donor site model */
  if (!had_err && bssm_param->gc_donor_model_set)
    had_err = bssm_model_read(&bssm_param->gc_donor_model, file, err);

  /* read AG acceptor site model */
  if (!had_err && bssm_param->ag_acceptor_model_set)
    had_err = bssm_model_read(&bssm_param->ag_acceptor_model, file, err);

  if (!had_err)
    gt_xfclose(file);

  if (had_err) {
    gth_bssm_param_delete(bssm_param);
    return NULL;
  }
  return bssm_param;
}

GthBssmParam* gth_bssm_param_extract(unsigned long speciesnum, GtError *err)
{
  GthBssmParam *bssm_param;
  unsigned long i, j, k, l;

  gt_error_check(err);

  bssm_param = gth_bssm_param_new();
  bssm_param->versionnum = (unsigned char) 2;
  bssm_param->gt_donor_model_set    = true;
  bssm_param->gc_donor_model_set    = false;
  bssm_param->ag_acceptor_model_set = true;

  if (speciesnum <= 1) {
    /* read in the human and mouse cases */
    bssm_param->gt_donor_model.hypothesisnum    = HYPOTHESIS2;
    bssm_param->ag_acceptor_model.hypothesisnum = HYPOTHESIS2;
    for (i = 0; i < HYPOTHESIS2; i++) {
      for (j = 0; j < WINSIZE + 2; j++) {
        for (k = 0; k < 4; k++) {
          for (l = 0; l < 4; l++) {
            bssm_param->gt_donor_model.hypotables.hypo2table[i][j][k][l] =
              (LOWPRECPROBTYPE) GU_2[speciesnum][i][j][k][l];
            bssm_param->ag_acceptor_model.hypotables.hypo2table[i][j][k][l] =
              (LOWPRECPROBTYPE) AG_2[speciesnum][i][j][k][l];
          }
        }
      }
    }
  }
  else if (speciesnum <  NUMOFSPECIES) {
    /* read in all others */
    bssm_param->gt_donor_model.hypothesisnum    = HYPOTHESIS7;
    bssm_param->ag_acceptor_model.hypothesisnum = HYPOTHESIS7;
    for (i = 0; i < HYPOTHESIS7; i++) {
      for (j = 0; j < WINSIZE + 2; j++) {
        for (k = 0; k < 4; k++) {
          for (l = 0; l < 4; l++) {
            bssm_param->gt_donor_model.hypotables.hypo7table[i][j][k][l] =
              (LOWPRECPROBTYPE) GU_7[speciesnum - 2][i][j][k][l];
            bssm_param->ag_acceptor_model.hypotables.hypo7table[i][j][k][l] =
              (LOWPRECPROBTYPE) AG_7[speciesnum - 2][i][j][k][l];
          }
        }
      }
    }
  }
  else {
    gt_error_set(err, "illegal speciesnum (speciesnum == %lu)", speciesnum);
    gth_bssm_param_delete(bssm_param);
    return NULL;
  }

  /* setting the window sizes */
  bssm_param->gt_donor_model.windowsizeleft     =  wsize[speciesnum][0][0];
  bssm_param->gt_donor_model.windowsizeright    =  wsize[speciesnum][0][1];
  bssm_param->ag_acceptor_model.windowsizeleft  =  wsize[speciesnum][1][0];
  bssm_param->ag_acceptor_model.windowsizeright =  wsize[speciesnum][1][1];

  return bssm_param;
}

void gth_bssm_param_delete(GthBssmParam *bssm_param)
{
  if (!bssm_param) return;
  gt_free(bssm_param);
}

#ifndef NDEBUG
static bool bssm_models_are_equal(GthBssmModel *checkmodel,
                                 GthBssmModel *testmodel)
{
  unsigned long i, j, k, l;

  if (checkmodel->hypothesisnum != testmodel->hypothesisnum)
    return  false;
  if (checkmodel->windowsizeleft != testmodel->windowsizeleft)
    return  false;
  if (checkmodel->windowsizeright != testmodel->windowsizeright)
    return  false;

  if (checkmodel->hypothesisnum == HYPOTHESIS2) {
    for (i = 0; i < HYPOTHESIS2; i++) {
        for (j = 0; j < WINSIZE + 2; j++) {
          for (k = 0; k < 4; k++) {
            for (l = 0; l < 4; l++) {
              if (checkmodel->hypotables.hypo2table[i][j][k][l] !=
                  testmodel->hypotables.hypo2table[i][j][k][l]) {
                return false;
              }
            }
          }
        }
      }
  }
  else if (checkmodel->hypothesisnum == HYPOTHESIS7) {
    for (i = 0; i < HYPOTHESIS7; i++) {
      for (j = 0; j < WINSIZE + 2; j++) {
        for (k = 0; k < 4; k++) {
          for (l = 0; l < 4; l++) {
            if (checkmodel->hypotables.hypo7table[i][j][k][l] !=
                testmodel->hypotables.hypo7table[i][j][k][l]) {
              return false;
            }
          }
        }
      }
    }
  }
  else
    return false;

  return true;
}
#endif

#ifndef NDEBUG
static bool bssmfile_equals_param(const char *filename,
                                  GthBssmParam *checkparam)
{
  GthBssmParam *testparam;

  /* reading test parameters from file */
  testparam = gth_bssm_param_load(filename, NULL);
  gt_assert(testparam);

  if (checkparam->versionnum != testparam->versionnum) {
    gth_bssm_param_delete(testparam);
    return  false;
  }
  if (checkparam->gt_donor_model_set != testparam->gt_donor_model_set) {
    gth_bssm_param_delete(testparam);
    return  false;
  }
  if (checkparam->gc_donor_model_set != testparam->gc_donor_model_set) {
    gth_bssm_param_delete(testparam);
    return  false;
  }
  if (checkparam->ag_acceptor_model_set != testparam->ag_acceptor_model_set) {
    gth_bssm_param_delete(testparam);
    return  false;
  }

  if (checkparam->gt_donor_model_set) {
    if (!bssm_models_are_equal(&checkparam->gt_donor_model,
                              &testparam->gt_donor_model)) {
      gth_bssm_param_delete(testparam);
      return false;
    }
  }
  if (checkparam->gc_donor_model_set) {
    if (!bssm_models_are_equal(&checkparam->gc_donor_model,
                              &testparam->gc_donor_model)) {
      gth_bssm_param_delete(testparam);
      return false;
    }
  }
  if (checkparam->ag_acceptor_model_set) {
    if (!bssm_models_are_equal(&checkparam->ag_acceptor_model,
                              &testparam->ag_acceptor_model)) {
      gth_bssm_param_delete(testparam);
      return false;
    }
  }

    gth_bssm_param_delete(testparam);
  return true;
}
#endif

static int writeBssmmodeltofile(FILE *file, GthBssmModel *bssm_model,
                                GtError *err)
{
  int had_err = 0;
  gt_error_check(err);
  gt_xfwrite(&bssm_model->hypothesisnum, sizeof (unsigned long), 1, file);
  if (bssm_model->hypothesisnum != HYPOTHESIS7 &&
      bssm_model->hypothesisnum != HYPOTHESIS2) {
    gt_error_set(err, "BSSM model contains unknown hypothesis number");
    had_err = -1;
  }
  if (!had_err) {
    gt_xfwrite(&bssm_model->windowsizeleft, sizeof (unsigned long), 1, file);
    gt_xfwrite(&bssm_model->windowsizeright, sizeof (unsigned long), 1, file);
    switch (bssm_model->hypothesisnum) {
      case HYPOTHESIS2:
        gt_xfwrite(&bssm_model->hypotables.hypo2table,  sizeof (Hypo2table), 1,
                   file);
        break;
      case HYPOTHESIS7:
        gt_xfwrite(&bssm_model->hypotables.hypo7table,  sizeof (Hypo7table), 1,
                   file);
        break;
      default: gt_assert(0);
    }
  }
  return had_err;
}

int gth_bssm_param_save(GthBssmParam *bssm_param, const char *filename,
                        GtError *err)
{
  FILE *file;
  int had_err = 0;

  gt_error_check(err);

  file = gt_fa_xfopen(filename, "w");
  gt_assert(file);

  /* write version number */
  gt_xfwrite(&bssm_param->versionnum,  sizeof (unsigned char), 1, file);

  /* write in model variables */
  gt_xfwrite(&bssm_param->gt_donor_model_set,  sizeof (bool), 1, file);
  gt_xfwrite(&bssm_param->gc_donor_model_set,  sizeof (bool), 1, file);
  gt_xfwrite(&bssm_param->ag_acceptor_model_set,  sizeof (bool), 1, file);

  /* check if at least one model is set */
  if (!bssm_param->gt_donor_model_set &&
      !bssm_param->gc_donor_model_set &&
      !bssm_param->ag_acceptor_model_set) {
    gt_error_set(err, "BSSM parameter to write contain no model");
    had_err = -1;
  }

  /* write GT donor site model */
  if (!had_err && bssm_param->gt_donor_model_set)
    had_err = writeBssmmodeltofile(file, &bssm_param->gt_donor_model, err);

  /* write GC donor site model */
  if (!had_err && bssm_param->gc_donor_model_set)
    had_err = writeBssmmodeltofile(file, &bssm_param->gc_donor_model, err);

  /* write AG acceptor site model */
  if (!had_err && bssm_param->ag_acceptor_model_set)
    had_err = writeBssmmodeltofile(file, &bssm_param->ag_acceptor_model, err);

  gt_fa_xfclose(file);

#ifndef NDEBUG
  /* make sure written model equals the input model */
  if (!had_err) {
    gt_assert(bssmfile_equals_param(filename, bssm_param));
  }
#endif

  return had_err;
}

/* The following function outouts <bssm_model>.
   It is assumed that the Hypo7table is used. */
static void bssm_model_echo(const GthBssmModel *bssm_model, FILE *outfp)
{
  unsigned long i, j, k, l;

  for (i = 0; i < HYPOTHESIS7; i++) {
    fprintf(outfp,"\n\nHypothesis: %lu", i);
    for (j = 0; j < STRINGSIZE; j++) {
      fprintf(outfp,"\n");
      for (k = 0; k < ALPHSIZE; k++) {
        fprintf(outfp,"\n");
        for (l = 0; l < ALPHSIZE; l++) {
          fprintf(outfp,"%.4f ", bssm_model->hypotables.hypo7table[i][j][k][l]);
        }
      }
    }
  }
  fprintf(outfp,"\n\n");
}

void gth_bssm_param_echo(const GthBssmParam *bssm_param, FILE *outfp)
{
  gt_assert(bssm_param && outfp);
  fprintf(outfp,"BSSMPARAMVERSION is %u\n\n", bssm_param->versionnum);
  fprintf(outfp,"Is the GT donor model set? -> %s\n",
          GTH_SHOWBOOL(bssm_param->gt_donor_model_set));
  fprintf(outfp,"Is the GC donor model set? -> %s\n\n",
          GTH_SHOWBOOL(bssm_param->gc_donor_model_set));
  fprintf(outfp,"Is the AG acceptor model set? -> %s\n\n",
          GTH_SHOWBOOL(bssm_param->ag_acceptor_model_set));

  if (bssm_param->gt_donor_model_set) {
    fprintf(outfp,"reporting GT donor model parameterization");
    bssm_model_echo(&bssm_param->gt_donor_model, outfp);
  }

  if (bssm_param->gc_donor_model_set) {
    fprintf(outfp,"reporting GC donor model parameterization");
    bssm_model_echo(&bssm_param->gc_donor_model, outfp);
  }

  if (bssm_param->ag_acceptor_model_set) {
    fprintf(outfp,"reporting AG acceptor model parameterization");
    bssm_model_echo(&bssm_param->ag_acceptor_model, outfp);
  }
}

void gth_bssm_param_show_info(const GthBssmParam *bssm_param, GtFile *outfp)
{
#define SEVENCLASSSTRING        "seven-class"
#define TWOCLASSSTRING          "two-class"

#define PRINT_CLASS_STRING(MODEL) \
  if (bssm_param->MODEL##_model_set) \
  { \
    gt_file_xprintf(outfp, " (%s)", \
                    bssm_param->MODEL##_model.hypothesisnum == HYPOTHESIS7 \
                    ? SEVENCLASSSTRING : TWOCLASSSTRING); \
  } \
  gt_file_xfputc('\n', outfp);

  gt_file_xprintf(outfp,
                  "%c the specified BSSM parameter file contains the following "
                  "models:\n", COMMENTCHAR);
  gt_file_xprintf(outfp, "%c GT donor sites   = %s", COMMENTCHAR,
                  GTH_SHOWBOOL(bssm_param->gt_donor_model_set));
  PRINT_CLASS_STRING(gt_donor);

  gt_file_xprintf(outfp, "%c GC donor sites   = %s", COMMENTCHAR,
                  GTH_SHOWBOOL(bssm_param->gc_donor_model_set));
  PRINT_CLASS_STRING(gc_donor);

  gt_file_xprintf(outfp, "%c AG acceptor sites= %s", COMMENTCHAR,
                  GTH_SHOWBOOL(bssm_param->ag_acceptor_model_set));
  PRINT_CLASS_STRING(ag_acceptor);
}

static void set_window_sizes_in_Bssmmodel(GthBssmModel *bssm_model)
{
  bssm_model->hypothesisnum      = HYPOTHESIS7;

  /* We have decided to leave maximum splice signal
     extent window at...a maximal value!            */
  bssm_model->windowsizeleft  = MAXSPLICESIG;
  bssm_model->windowsizeright = MAXSPLICESIG;
}

/* updates the BSSM parameterization file */
static void build_bssm(GtBioseq *bioseq, GthBssmModel *bssm_model,
                       unsigned int hypothesisnum)
{
  unsigned long mono_ct[STRINGSIZE-1][ALPHSIZE],         /* Mononuc freq */
                di_ct[STRINGSIZE-1][ALPHSIZE][ALPHSIZE]; /* Dinuc freq */
  double mono_freq,      /* Mononuc relative freq */
         di_freq;        /* Dinuc relative freq */
  unsigned long i, j, k, /* Iterator variables */
                num_entries = gt_bioseq_number_of_sequences(bioseq);
  const GtUchar *encoded_seq;
  GtSeq *seq;

  /* Inits of local variables */
  for (i = 0; i < (STRINGSIZE-1); i++) {
    for (j = 0; j < ALPHSIZE; j++) {
      mono_ct[i][j] = INITVAL_INT;
      for (k = 0; k < ALPHSIZE; k++)
        di_ct[i][j][k] = INITVAL_INT;
    }
  }

  /* mononucleotides */
  for (i = 0; i < (STRINGSIZE-1); i++) {
    for (j = 0; j < num_entries; j++) {
      seq = gt_bioseq_get_seq(bioseq, j);
      encoded_seq = gt_seq_get_encoded(seq);
      mono_ct[i][encoded_seq[i]]++;
    }
  }

  /* dinucleotides */
  for (i = 0; i < (STRINGSIZE-1); i++) {
    for (j = 0; j < num_entries; j++) {
      seq = gt_bioseq_get_seq(bioseq, j);
      encoded_seq = gt_seq_get_encoded(seq);
      di_ct[i][encoded_seq[i]]
              [encoded_seq[i + 1]]++;
    }
  }

  /* Record equilibrium frequencies (1st ``slot" in transition freqs) */
  for (i = 0; i < ALPHSIZE; i++) {
    for (j = 0; j < ALPHSIZE; j++) {
      bssm_model->hypotables
      .hypo7table[hypothesisnum][0][i][j] = (LOWPRECPROBTYPE)
                                            mono_ct[0][i] / num_entries;
    }
  }

  /* Populate the remaining transition frequencies */
  for (k = 1; k < STRINGSIZE; k++) {
    for (i = 0; i < ALPHSIZE; i++) {
      mono_freq = (double) mono_ct[k-1][i] / num_entries;
      for (j = 0; j < ALPHSIZE; j++) {
        di_freq = (double) di_ct[k-1][i][j] / num_entries;
        if (mono_freq == 0.0) {
          bssm_model->hypotables
          .hypo7table[hypothesisnum][k][i][j] = (LOWPRECPROBTYPE) NULLPROB;
        }
        else {
          bssm_model->hypotables
          .hypo7table[hypothesisnum][k][i][j] = (LOWPRECPROBTYPE)
                                                (di_freq / mono_freq);
        }
      }

      /* Remove non-zero transition probabilities:
         Briefly, 0.0 entries (dinucleotide absent in training corpus) are
         replaced arbitrarily by PSEUDOPROB, and non-0.0 entries p are replaced
         by p = p * (1 - 4 * PSEUDOPROB) + PSEUDOPROB */
      for (j = 0; j < ALPHSIZE; ++j) {
        /* If any entry is NULLPROB, ALL elements in the row need fixed */
        if (bssm_model->hypotables
            .hypo7table[hypothesisnum][k][i][j] == NULLPROB) {
          /* Fix all elements in the row, then break */
          for (j = 0; j < ALPHSIZE; j++) {
            if (bssm_model->hypotables
                .hypo7table[hypothesisnum][k][i][j] == NULLPROB) {
               bssm_model->hypotables
                .hypo7table[hypothesisnum][k][i][j] = (LOWPRECPROBTYPE)
                                                      PSEUDOPROB;
            }
            else {
              /* Adjust non-zero transition prob */
              bssm_model->hypotables.hypo7table[hypothesisnum][k][i][j] =
                (LOWPRECPROBTYPE)
                (bssm_model->hypotables.hypo7table[hypothesisnum][k][i][j] *
                 (1 - (4 * PSEUDOPROB)) + PSEUDOPROB);
            }
          }
          break;
        }
      }
    }
  }
}

int gth_bssm_param_parameterize(GthBssmParam *bssm_param, const char *path,
                                Termtype termtype, bool gzip, GtError *err)
{
  GtAlphabet *alphabet = NULL;
  GtBioseq *bioseq;
  char file2proc[PATH_MAX+1];
  unsigned long i, j;
  int had_err = 0;

  gt_error_check(err);

  /* set version number */
  bssm_param->versionnum = (unsigned char) MYVERSION;

  /* set model to true and set window sizes */
  switch (termtype) {
    case GT_DONOR_TYPE:
      bssm_param->gt_donor_model_set = true;
      set_window_sizes_in_Bssmmodel(&bssm_param->gt_donor_model);
      break;
    case GC_DONOR_TYPE:
      bssm_param->gc_donor_model_set = true;
      set_window_sizes_in_Bssmmodel(&bssm_param->gc_donor_model);
      break;
    case AG_ACCEPTOR_TYPE:
      bssm_param->ag_acceptor_model_set = true;
      set_window_sizes_in_Bssmmodel(&bssm_param->ag_acceptor_model);
      break;
    default: gt_assert(0);
  }

  for (i = 0; !had_err && i < NUMOFFILES; i++) {
    /* process datafile */
    strcpy(file2proc, path);
    switch (termtype) {
      case GT_DONOR_TYPE:
        strcat(file2proc, "/GT_donor/");
        strcat(file2proc, filenames[i]);
        break;
      case GC_DONOR_TYPE:
        strcat(file2proc, "/GC_donor/");
        strcat(file2proc, filenames[i]);
        break;
      case AG_ACCEPTOR_TYPE:
        strcat(file2proc, "/AG_acceptor/");
        strcat(file2proc, filenames[i]);
        break;
      default: gt_assert(0);
    }

    if (gzip)
      strcat(file2proc, ".gz");

    if (!(bioseq = gt_bioseq_new(file2proc, err)))
      had_err = -1;

    if (!had_err)
      alphabet = gt_bioseq_get_alphabet(bioseq);

    /* check here if all sequences have the length 102 and correct bases at
       positions 51 and 52 (i.e., GT, GC, or AG) */
    for (j = 0; !had_err && j < gt_bioseq_number_of_sequences(bioseq); j++) {
      const GtUchar *encoded_seq;
      GtSeq *seq;
      /* check length */
      if (gt_bioseq_get_sequence_length(bioseq, j) != STRINGSIZE) {
        gt_error_set(err, "sequence %lu in file \"%s\" does not have length %u",
                     j, file2proc, STRINGSIZE);
        had_err = -1;
      }
      seq = gt_bioseq_get_seq(bioseq, j);
      encoded_seq = gt_seq_get_encoded(seq);
      if (!had_err) {
        /* check base correctness */
        switch (termtype) {
          case GT_DONOR_TYPE:
            if (encoded_seq[50] != gt_alphabet_encode(alphabet, 'G') ||
                encoded_seq[51] != gt_alphabet_encode(alphabet, 'T')) {
              gt_error_set(err, "sequence %lu in file \"%s\" is not a GT "
                                "sequence", j, file2proc);
              had_err = -1;
            }
            break;
          case GC_DONOR_TYPE:
            if (encoded_seq[50] != gt_alphabet_encode(alphabet, 'G') ||
                encoded_seq[51] != gt_alphabet_encode(alphabet, 'C')) {
              gt_error_set(err, "sequence %lu in file \"%s\" is not a GC "
                                "sequence", j, file2proc);
              had_err = -1;
            }
            break;
          case AG_ACCEPTOR_TYPE:
            if (encoded_seq[50] != gt_alphabet_encode(alphabet, 'A') ||
                encoded_seq[51] != gt_alphabet_encode(alphabet, 'G')) {
              gt_error_set(err, "sequence %lu in file \"%s\" is not a AG "
                                "sequence", j, file2proc);
              had_err = -1;
            }
            break;
          default: gt_assert(0);
        }
      }
    }

    if (!had_err) {
      switch (termtype) {
        case GT_DONOR_TYPE:
          build_bssm(bioseq, &bssm_param->gt_donor_model, i);
          break;
        case GC_DONOR_TYPE:
          build_bssm(bioseq, &bssm_param->gc_donor_model, i);
          break;
        case AG_ACCEPTOR_TYPE:
          build_bssm(bioseq, &bssm_param->ag_acceptor_model, i);
          break;
        default: gt_assert(0);
      }
    }

    /* free space */
    gt_bioseq_delete(bioseq);
  }

  return had_err;
}
