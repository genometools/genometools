/*
  Copyright (c) 2004-2009 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2004      Michael E Sparks <mespar1@iastate.edu>
  Copyright (c) 2004-2008 Center for Bioinformatics, University of Hamburg

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

/*
  Main implementation file for a program to build a 1st order Markov model for
  GTH/GSQ splice sites.

  Please note: Transliteration scheme
    A   -> 0
    C   -> 1
    G   -> 2
    T/U -> 3

  Input data files are--strictly!--named as follows, using
  an obvious schema:
    F0_don  F1_don  F2_don  Fi_don  T0_don  T1_don  T2_don
    F0_acc  F1_acc  F2_acc  Fi_acc  T0_acc  T1_acc  T2_acc
  Phase is denoted as follows (Brendel conventions) :
    1 -> C O D |
    2 -> C | O D
    0 -> C O | D
*/

#include "core/option_api.h"
#include "gth/bssm_param.h"
#include "gth/gt_gthbssmbuild.h"

#define GTDONORMODEL_OPT_CSTR     "gtdonor"
#define GCDONORMODEL_OPT_CSTR     "gcdonor"
#define AGACCEPTORMODEL_OPT_CSTR  "agacceptor"

typedef struct {
  GtStr *bssmfile,
        *datapath;
  bool calcGTdonormodel,
       calcGCdonormodel,
       calcAGacceptormodel,
       gzip;
} Commandlineopts;

static void initcommandlineopts(Commandlineopts *commandlineopts)
{
  commandlineopts->bssmfile = gt_str_new();
  commandlineopts->datapath = gt_str_new();
}

static GtOPrval gthbssmbuild_parse_options(int *parsed_args,
                                           Commandlineopts *commandlineopts,
                                           int argc, const char **argv,
                                           GtShowVersionFunc gth_version_func,
                                           GtError *err)
{
  GtOptionParser *op;
  GtOption *optbssmfile, *optdatapath, *optgtdonormodel, *optgcdonormodel,
         *optagacceptormodel, *optgzip;
  GtOPrval oprval;
  gt_error_check(err);

  op = gt_option_parser_new("[option ...] -datapath dir -bssmfile file ",
                            "Build a BSSM file from a directory tree "
                            "containing the training data.");

  /* specify all options with the corresponding help-text */
  optbssmfile = gt_option_new_filename("bssmfile", "specify name of BSSM file "
                                       "to store parameters in",
                                       commandlineopts->bssmfile);
  gt_option_parser_add_option(op, optbssmfile);

  optdatapath = gt_option_new_string("datapath", "specify root of "
                                     "species-specific training data directory "
                                     "tree", commandlineopts->datapath, NULL);
  gt_option_parser_add_option(op, optdatapath);

  optgtdonormodel = gt_option_new_bool(GTDONORMODEL_OPT_CSTR, "train GT donor "
                                       "model",
                                       &commandlineopts->calcGTdonormodel,
                                       false);
  gt_option_parser_add_option(op, optgtdonormodel);

  optgcdonormodel = gt_option_new_bool(GCDONORMODEL_OPT_CSTR, "train GC donor "
                                       "model",
                                       &commandlineopts->calcGCdonormodel,
                                       false);
  gt_option_parser_add_option(op, optgcdonormodel);

  optagacceptormodel = gt_option_new_bool(AGACCEPTORMODEL_OPT_CSTR, "train AG "
                                          "acceptor model",
                                          &commandlineopts->calcAGacceptormodel,
                                          false);
  gt_option_parser_add_option(op, optagacceptormodel);

  optgzip = gt_option_new_bool("gzip", "use gzip'ed input files",
                               &commandlineopts->gzip, false);
  gt_option_parser_add_option(op, optgzip);

  /* mandatory options */
  gt_option_is_mandatory(optbssmfile);
  gt_option_is_mandatory(optdatapath);

  gt_option_parser_set_max_args(op, 0);
  gt_option_parser_set_mail_address(op, "<gremme@zbh.uni-hamburg.de>");
  oprval = gt_option_parser_parse(op, parsed_args, argc, argv, gth_version_func,
                                  err);

  /* some checks */
  if (oprval == GT_OPTION_PARSER_OK &&
      !gt_option_is_set(optgtdonormodel) &&
      !gt_option_is_set(optgcdonormodel) &&
      !gt_option_is_set(optagacceptormodel)) {
    gt_error_set(err, "at least one of the options -%s, -%s, or -%s have to "
                      "be used", GTDONORMODEL_OPT_CSTR, GCDONORMODEL_OPT_CSTR,
              AGACCEPTORMODEL_OPT_CSTR);
    oprval = GT_OPTION_PARSER_ERROR;
  }

  gt_option_parser_delete(op);

  return oprval;
}

static void freecommandlineopts(Commandlineopts *commandlineopts)
{
  gt_str_delete(commandlineopts->bssmfile);
  gt_str_delete(commandlineopts->datapath);
}

int gt_gthbssmbuild(int argc, const char **argv,
                    GtShowVersionFunc gth_version_func, GtError *err)
{
  GthBSSMParam *bssm_param; /* stores model parameterization  */
  Commandlineopts commandlineopts;
  int parsed_args, had_err = 0;

  /* process command line args */
  initcommandlineopts(&commandlineopts);
  switch (gthbssmbuild_parse_options(&parsed_args, &commandlineopts, argc, argv,
                                     gth_version_func, err)) {
    case GT_OPTION_PARSER_OK: break;
    case GT_OPTION_PARSER_ERROR: return -1;
    case GT_OPTION_PARSER_REQUESTS_EXIT: return 0;
  }
  gt_assert(parsed_args = argc);

  /* initialize our parameter object */
  bssm_param = gth_bssm_param_new();

  /* process each model type */
  if (commandlineopts.calcGTdonormodel) {
    had_err = gth_bssm_param_parameterize(bssm_param,
                                          gt_str_get(commandlineopts.datapath),
                                          GT_DONOR_TYPE,
                                          commandlineopts.gzip, err);
  }
  if (!had_err && commandlineopts.calcGCdonormodel) {
    had_err = gth_bssm_param_parameterize(bssm_param,
                                          gt_str_get(commandlineopts.datapath),
                                          GC_DONOR_TYPE,
                                          commandlineopts.gzip, err);
  }
  if (!had_err && commandlineopts.calcAGacceptormodel) {
    had_err = gth_bssm_param_parameterize(bssm_param,
                                          gt_str_get(commandlineopts.datapath),
                                          AG_ACCEPTOR_TYPE,
                                          commandlineopts.gzip, err);
  }

  /* overwrite the binary file with our new object */
  if (!had_err) {
    had_err = gth_bssm_param_save(bssm_param,
                                  gt_str_get(commandlineopts.bssmfile), err);
  }

  /* free memory */
  freecommandlineopts(&commandlineopts);
  gth_bssm_param_delete(bssm_param);

  return had_err;
}
