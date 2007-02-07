/*
  Copyright (c) 2006 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include <string.h>
#include "evaluator.h"
#include "xansi.h"

struct Evaluator {
  unsigned long T, /* true */
                A, /* actual */
                P; /* predicted */
};

Evaluator* evaluator_new(void)
{
  return xcalloc(1, sizeof(Evaluator));
}

void evaluator_add_true(Evaluator *e)
{
  assert(e);
  assert(e->T < e->A && e->T < e->P);
  e->T++;
}

void evaluator_add_actual(Evaluator *e, unsigned long inc)
{
  assert(e);
  e->A += inc;
}

void evaluator_add_predicted(Evaluator *e, unsigned long inc)
{
  assert(e);
  e->P += inc;
}

double evaluator_get_sensitivity(const Evaluator *e)
{
  double sensitivity = 1.0;
  assert(e);
  assert(e->T <= e->A);
  if (e->A)
    sensitivity = (double) e->T / e->A;
  assert(sensitivity >= 0.0 && sensitivity <= 1.0);
  return sensitivity;
}

double evaluator_get_specificity(const Evaluator *e)
{
  double specificity = 1.0;
  assert(e);
  assert(e->T <= e->P);
  if (e->P)
    specificity = (double) e->T / e->P;
  assert(specificity >= 0.0 && specificity <= 1.0);
  return specificity;
}

void evaluator_show_sensitivity(const Evaluator *e, FILE *outfp)
{
  assert(e && outfp);
  fprintf(outfp, "%6.2f%% (%lu/%lu)", evaluator_get_sensitivity(e) * 100.0,
          e->T, e->A);
}

void evaluator_show_specificity(const Evaluator *e, FILE *outfp)
{
  assert(e && outfp);
  fprintf(outfp, "%6.2f%% (%lu/%lu)", evaluator_get_specificity(e) * 100.0,
          e->T, e->P);
}
void evaluator_reset(Evaluator *e)
{
  assert(e);
  memset(e, 0, sizeof(Evaluator));
}

int evaluator_unit_test(void)
{
  Evaluator *e = evaluator_new();

  assert(evaluator_get_sensitivity(e) == 1.0);
  assert(evaluator_get_specificity(e) == 1.0);

  evaluator_add_actual(e, 1);
  assert(evaluator_get_sensitivity(e) == 0.0);
  assert(evaluator_get_specificity(e) == 1.0);

  evaluator_add_predicted(e, 1);
  assert(evaluator_get_sensitivity(e) == 0.0);
  assert(evaluator_get_specificity(e) == 0.0);

  evaluator_add_true(e);
  assert(evaluator_get_sensitivity(e) == 1.0);
  assert(evaluator_get_specificity(e) == 1.0);

  evaluator_reset(e);
  assert(evaluator_get_sensitivity(e) == 1.0);
  assert(evaluator_get_specificity(e) == 1.0);

  evaluator_add_predicted(e, 1);
  assert(evaluator_get_sensitivity(e) == 1.0);
  assert(evaluator_get_specificity(e) == 0.0);

  evaluator_reset(e);
  assert(evaluator_get_sensitivity(e) == 1.0);
  assert(evaluator_get_specificity(e) == 1.0);

  evaluator_add_actual(e, 2);
  evaluator_add_predicted(e, 2);
  evaluator_add_true(e);
  assert(evaluator_get_sensitivity(e) == 0.5);
  assert(evaluator_get_specificity(e) == 0.5);

  evaluator_reset(e);
  assert(evaluator_get_sensitivity(e) == 1.0);
  assert(evaluator_get_specificity(e) == 1.0);

  evaluator_add_actual(e, 4);
  evaluator_add_predicted(e, 4);
  evaluator_add_true(e);
  assert(evaluator_get_sensitivity(e) == 0.25);
  assert(evaluator_get_specificity(e) == 0.25);
  evaluator_add_true(e);
  assert(evaluator_get_sensitivity(e) == 0.5);
  assert(evaluator_get_specificity(e) == 0.5);

  evaluator_add_true(e);
  assert(evaluator_get_sensitivity(e) == 0.75);
  assert(evaluator_get_specificity(e) == 0.75);

  evaluator_add_true(e);
  assert(evaluator_get_sensitivity(e) == 1.0);
  assert(evaluator_get_specificity(e) == 1.0);

  evaluator_free(e);

  return EXIT_SUCCESS;
}

void evaluator_free(Evaluator *e)
{
  if (!e) return;
  free(e);
}
