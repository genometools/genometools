/*
  Copyright (c) 2005-2010 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2005-2008 Center for Bioinformatics, University of Hamburg

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
#include "core/ma_api.h"
#include "gth/default.h"
#include "gth/bssm_param.h"
#include "gth/splice_site_model_rep.h"

#define SET_PROB_VALUE(V, PROB)                           \
        ssm->V             = PROB;                        \
        ssm->log_##V       = (GthFlt) log(PROB); \
        ssm->log1minus_##V = (GthFlt) log(1.0 - (PROB))

GthSpliceSiteModel* gth_splice_site_model_new(void)
{
  GthSpliceSiteModel *ssm = gt_malloc(sizeof *ssm);

  ssm->useU12intronmodel = !GTH_DEFAULT_DISABLEU12INTRONMODEL;

  SET_PROB_VALUE(genericGTdonorprob, GTH_DEFAULT_GENERIC_GT_DONORPROB);
  SET_PROB_VALUE(nongenericGTdonorprob, GTH_DEFAULT_NONGENERIC_GT_DONORPROB);
  SET_PROB_VALUE(genericGCdonorprob, GTH_DEFAULT_GENERIC_GC_DONORPROB);
  SET_PROB_VALUE(nongenericGCdonorprob, GTH_DEFAULT_NONGENERIC_GC_DONORPROB);
  SET_PROB_VALUE(genericATdonorprob, GTH_DEFAULT_GENERIC_AT_DONORPROB);
  SET_PROB_VALUE(nongenericATdonorprob, GTH_DEFAULT_NONGENERIC_AT_DONORPROB);
  SET_PROB_VALUE(genericAGacceptorprob, GTH_DEFAULT_GENERIC_AG_ACCEPTORPROB);
  SET_PROB_VALUE(nongenericAGacceptorprob,
                 GTH_DEFAULT_NONGENERIC_AG_ACCEPTORPROB);
  SET_PROB_VALUE(genericACacceptorprob, GTH_DEFAULT_GENERIC_AC_ACCEPTORPROB);
  SET_PROB_VALUE(nongenericACacceptorprob,
                 GTH_DEFAULT_NONGENERIC_AC_ACCEPTORPROB);
  SET_PROB_VALUE(genericothersplicesitep,
                 GTH_DEFAULT_GENERIC_OTHERSPLICESITEPROB);
  SET_PROB_VALUE(nongenericothersplicesitep,
                 GTH_DEFAULT_NONGENERIC_OTHERSPLICESITEPROB);
  gth_splice_site_model_set_U12typedonorprob(ssm,
                                             GTH_DEFAULT_U12_TYPEDONORPROB);
  gth_splice_site_model_set_U12typedonorprob_one_mismatch(ssm,
                                      GTH_DEFAULT_U12_TYPEDONORPROBONEMISMATCH);

  ssm->bssm_param = NULL;

  return ssm;
}

void gth_splice_site_model_delete(GthSpliceSiteModel *ssm)
{
  if (!ssm) return;
  gth_bssm_param_delete(ssm->bssm_param);
  gt_free(ssm);
}

int gth_splice_site_model_load_bssm(GthSpliceSiteModel *ssm,
                                    const char *bssmfile, GtError *err)
{
  GtStr *bssmfilename;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(ssm && bssmfile && strlen(bssmfile));
  if (ssm->bssm_param) {
    gth_bssm_param_delete(ssm->bssm_param);
    ssm->bssm_param = NULL;
  }
  bssmfilename = gt_str_new();
  gt_str_append_cstr(bssmfilename, bssmfile);
  gt_str_append_char(bssmfilename, '.');
  gt_str_append_cstr(bssmfilename, BSSMFILEENDING);
  if (!(ssm->bssm_param = gth_bssm_param_load(gt_str_get(bssmfilename), err)))
    had_err = -1;
  gt_str_delete(bssmfilename);
  return had_err;
}

void gth_splice_site_model_U12intronmodel_set_usage(GthSpliceSiteModel *ssm,
                                                    bool useU12intronmodel)
{
  ssm->useU12intronmodel = useU12intronmodel;
}

void gth_splice_site_model_set_U12typedonorprob(GthSpliceSiteModel *ssm,
                                                GthFlt prob)
{
  SET_PROB_VALUE(U12typedonorprob, prob);
}

void gth_splice_site_model_set_U12typedonorprob_one_mismatch(GthSpliceSiteModel
                                                             *ssm,
                                                             GthFlt prob)
{
  SET_PROB_VALUE(U12typedonorprobonemismatch, prob);
}
