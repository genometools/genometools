/*
  Copyright (c) 2010 Gordon Gremme <gordon@gremme.org>

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
#include "gth/gthprobdef.h"
#include "gth/bssm_param_rep.h"
#include "gth/bssm_param_rmsd.h"

static void compute_rmsd(GthDbl *rmsd, GtUword *n,
                         GthBSSMModel *bssm_model_1, GthBSSMModel *bssm_model_2)
{
  GtUword i, j, k, l;
  GthDbl x, y;
  gt_assert(rmsd && n && bssm_model_1 && bssm_model_2);
  *rmsd = 0.0;
  *n = 0;
  for (i = 0; i < HYPOTHESIS7; i++) {
    for (j = 0; j < STRINGSIZE; j++) {
      for (k = 0; k < ALPHSIZE; k++) {
        for (l = 0; l < ALPHSIZE; l++) {
          x = bssm_model_1->hypotables.hypo7table[i][j][k][l];
          y = bssm_model_2->hypotables.hypo7table[i][j][k][l];
          *rmsd = (x - y) * (x - y);
          (*n)++;
        }
      }
    }
  }
  *rmsd /= *n;
  *rmsd = sqrt(*rmsd);
}

static int show_rmsd(GthBSSMParam *bssm_1, GthBSSMParam *bssm_2, GtError *err)
{
  unsigned int current = 0;
  GtUword ns[3];
  GthDbl rmsds[3];
  gt_error_check(err);
  gt_assert(bssm_1 && bssm_2);
  gt_assert(gth_bssm_param_is_seven_class(bssm_1));
  gt_assert(gth_bssm_param_is_seven_class(bssm_2));
  if (bssm_1->gt_donor_model_set && bssm_2->gt_donor_model_set) {
    compute_rmsd(&rmsds[current], &ns[current], &bssm_1->gt_donor_model,
                 &bssm_2->gt_donor_model);
    printf("RMSD for GT donor site model:    %f\n", rmsds[current]);
    current++;
  }
  if (bssm_1->gc_donor_model_set && bssm_2->gc_donor_model_set) {
    compute_rmsd(&rmsds[current], &ns[current], &bssm_1->gc_donor_model,
                 &bssm_2->gc_donor_model);
    printf("RMSD for GC donor site model:    %f\n", rmsds[current]);
    current++;
  }
  if (bssm_1->ag_acceptor_model_set && bssm_2->ag_acceptor_model_set) {
    compute_rmsd(&rmsds[current], &ns[current], &bssm_1->ag_acceptor_model,
                 &bssm_2->ag_acceptor_model);
    printf("RMSD for AG acceptor site model: %f\n", rmsds[current]);
    current++;
  }
  if (current) {
    GthDbl overall_rmsd = 0.0;
    GtUword i, n = 0;
    for (i = 0; i < current; i++) {
      overall_rmsd += rmsds[i] * rmsds[i] * ns[i];
      n += ns[i];
    }
    overall_rmsd /= n;
    overall_rmsd = sqrt(overall_rmsd);
    printf("overall RMSD:                    %f\n", overall_rmsd);
  }
  else {
    gt_error_set(err, "given BSSM files have no common site models");
    return -1;
  }
  return 0;
}

int gth_bssm_param_rmsd_show(const char *bssm_file_1, const char *bssm_file_2,
                             GtError *err)
{
  GthBSSMParam *bssm_1, *bssm_2 = NULL;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(bssm_file_1 && bssm_file_2);
  if (!(bssm_1 = gth_bssm_param_load(bssm_file_1, err)))
    had_err = -1;
  if (!had_err) {
    if (!gth_bssm_param_is_seven_class(bssm_1)) {
      gt_error_set(err, "BSSM file '%s' is not seven-class", bssm_file_1);
      had_err = -1;
    }
  }
  if (!had_err) {
    if (!(bssm_2 = gth_bssm_param_load(bssm_file_2, err)))
      had_err = -1;
  }
  if (!had_err) {
    if (!gth_bssm_param_is_seven_class(bssm_2)) {
      gt_error_set(err, "BSSM file '%s' is not seven-class", bssm_file_2);
      had_err = -1;
    }
  }
  if (!had_err)
    had_err = show_rmsd(bssm_1, bssm_2, err);
  gth_bssm_param_delete(bssm_2);
  gth_bssm_param_delete(bssm_1);
  return had_err;
}
