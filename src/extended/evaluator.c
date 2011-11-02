/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg

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

#include <string.h>
#include "core/assert_api.h"
#include "core/ensure.h"
#include "core/ma.h"
#include "extended/evaluator.h"

struct GtEvaluator {
  unsigned long T, /* true */
                A, /* actual */
                P; /* predicted */
};

GtEvaluator* gt_evaluator_new(void)
{
  return gt_calloc(1, sizeof (GtEvaluator));
}

void gt_evaluator_add_true(GtEvaluator *evaluator)
{
  gt_assert(evaluator);
  gt_assert(evaluator->T < evaluator->A && evaluator->T < evaluator->P);
  evaluator->T++;
}

void gt_evaluator_add_actual(GtEvaluator *evaluator, unsigned long inc)
{
  gt_assert(evaluator);
  evaluator->A += inc;
}

void gt_evaluator_add_predicted(GtEvaluator *evaluator, unsigned long inc)
{
  gt_assert(evaluator);
  evaluator->P += inc;
}

double gt_evaluator_get_sensitivity(const GtEvaluator *evaluator)
{
  double sensitivity = 1.0;
  gt_assert(evaluator);
  gt_assert(evaluator->T <= evaluator->A);
  if (evaluator->A)
    sensitivity = (double) evaluator->T / evaluator->A;
  gt_assert(sensitivity >= 0.0 && sensitivity <= 1.0);
  return sensitivity;
}

double gt_evaluator_get_specificity(const GtEvaluator *evaluator)
{
  double specificity = 1.0;
  gt_assert(evaluator);
  gt_assert(evaluator->T <= evaluator->P);
  if (evaluator->P)
    specificity = (double) evaluator->T / evaluator->P;
  gt_assert(specificity >= 0.0 && specificity <= 1.0);
  return specificity;
}

void gt_evaluator_show_sensitivity(const GtEvaluator *evaluator, GtFile *outfp)
{
  gt_assert(evaluator);
  gt_file_xprintf(outfp, "%6.2f%% (%lu/%lu)",
                  gt_evaluator_get_sensitivity(evaluator) * 100.0,
                  evaluator->T, evaluator->A);
}

void gt_evaluator_show_specificity(const GtEvaluator *evaluator, GtFile *outfp)
{
  gt_assert(evaluator);
  gt_file_xprintf(outfp, "%6.2f%% (%lu/%lu)",
                  gt_evaluator_get_specificity(evaluator) * 100.0,
                  evaluator->T, evaluator->P);
}
void gt_evaluator_reset(GtEvaluator *evaluator)
{
  gt_assert(evaluator);
  memset(evaluator, 0, sizeof *evaluator);
}

int gt_evaluator_unit_test(GtError *err)
{
  GtEvaluator *evaluator = gt_evaluator_new();
  int had_err = 0;
  gt_error_check(err);

  gt_ensure(had_err, gt_evaluator_get_sensitivity(evaluator) == 1.0);
  gt_ensure(had_err, gt_evaluator_get_specificity(evaluator) == 1.0);

  gt_evaluator_add_actual(evaluator, 1);
  gt_ensure(had_err, gt_evaluator_get_sensitivity(evaluator) == 0.0);
  gt_ensure(had_err, gt_evaluator_get_specificity(evaluator) == 1.0);

  gt_evaluator_add_predicted(evaluator, 1);
  gt_ensure(had_err, gt_evaluator_get_sensitivity(evaluator) == 0.0);
  gt_ensure(had_err, gt_evaluator_get_specificity(evaluator) == 0.0);

  gt_evaluator_add_true(evaluator);
  gt_ensure(had_err, gt_evaluator_get_sensitivity(evaluator) == 1.0);
  gt_ensure(had_err, gt_evaluator_get_specificity(evaluator) == 1.0);

  gt_evaluator_reset(evaluator);
  gt_ensure(had_err, gt_evaluator_get_sensitivity(evaluator) == 1.0);
  gt_ensure(had_err, gt_evaluator_get_specificity(evaluator) == 1.0);

  gt_evaluator_add_predicted(evaluator, 1);
  gt_ensure(had_err, gt_evaluator_get_sensitivity(evaluator) == 1.0);
  gt_ensure(had_err, gt_evaluator_get_specificity(evaluator) == 0.0);

  gt_evaluator_reset(evaluator);
  gt_ensure(had_err, gt_evaluator_get_sensitivity(evaluator) == 1.0);
  gt_ensure(had_err, gt_evaluator_get_specificity(evaluator) == 1.0);

  gt_evaluator_add_actual(evaluator, 2);
  gt_evaluator_add_predicted(evaluator, 2);
  gt_evaluator_add_true(evaluator);
  gt_ensure(had_err, gt_evaluator_get_sensitivity(evaluator) == 0.5);
  gt_ensure(had_err, gt_evaluator_get_specificity(evaluator) == 0.5);

  gt_evaluator_reset(evaluator);
  gt_ensure(had_err, gt_evaluator_get_sensitivity(evaluator) == 1.0);
  gt_ensure(had_err, gt_evaluator_get_specificity(evaluator) == 1.0);

  gt_evaluator_add_actual(evaluator, 4);
  gt_evaluator_add_predicted(evaluator, 4);
  gt_evaluator_add_true(evaluator);
  gt_ensure(had_err, gt_evaluator_get_sensitivity(evaluator) == 0.25);
  gt_ensure(had_err, gt_evaluator_get_specificity(evaluator) == 0.25);
  gt_evaluator_add_true(evaluator);
  gt_ensure(had_err, gt_evaluator_get_sensitivity(evaluator) == 0.5);
  gt_ensure(had_err, gt_evaluator_get_specificity(evaluator) == 0.5);

  gt_evaluator_add_true(evaluator);
  gt_ensure(had_err, gt_evaluator_get_sensitivity(evaluator) == 0.75);
  gt_ensure(had_err, gt_evaluator_get_specificity(evaluator) == 0.75);

  gt_evaluator_add_true(evaluator);
  gt_ensure(had_err, gt_evaluator_get_sensitivity(evaluator) == 1.0);
  gt_ensure(had_err, gt_evaluator_get_specificity(evaluator) == 1.0);

  gt_evaluator_delete(evaluator);

  return had_err;
}

void gt_evaluator_delete(GtEvaluator *evaluator)
{
  if (!evaluator) return;
  gt_free(evaluator);
}
