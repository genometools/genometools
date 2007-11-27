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

#include <assert.h>
#include <string.h>
#include "libgtcore/ensure.h"
#include "libgtcore/ma.h"
#include "libgtext/evaluator.h"

struct Evaluator {
  unsigned long T, /* true */
                A, /* actual */
                P; /* predicted */
};

Evaluator* evaluator_new(Env *env)
{
  return ma_calloc(1, sizeof (Evaluator));
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
  memset(e, 0, sizeof (Evaluator));
}

int evaluator_unit_test(Env *env)
{
  Evaluator *e = evaluator_new(env);
  int had_err = 0;
  env_error_check(env);

  ensure(had_err, evaluator_get_sensitivity(e) == 1.0);
  ensure(had_err, evaluator_get_specificity(e) == 1.0);

  evaluator_add_actual(e, 1);
  ensure(had_err, evaluator_get_sensitivity(e) == 0.0);
  ensure(had_err, evaluator_get_specificity(e) == 1.0);

  evaluator_add_predicted(e, 1);
  ensure(had_err, evaluator_get_sensitivity(e) == 0.0);
  ensure(had_err, evaluator_get_specificity(e) == 0.0);

  evaluator_add_true(e);
  ensure(had_err, evaluator_get_sensitivity(e) == 1.0);
  ensure(had_err, evaluator_get_specificity(e) == 1.0);

  evaluator_reset(e);
  ensure(had_err, evaluator_get_sensitivity(e) == 1.0);
  ensure(had_err, evaluator_get_specificity(e) == 1.0);

  evaluator_add_predicted(e, 1);
  ensure(had_err, evaluator_get_sensitivity(e) == 1.0);
  ensure(had_err, evaluator_get_specificity(e) == 0.0);

  evaluator_reset(e);
  ensure(had_err, evaluator_get_sensitivity(e) == 1.0);
  ensure(had_err, evaluator_get_specificity(e) == 1.0);

  evaluator_add_actual(e, 2);
  evaluator_add_predicted(e, 2);
  evaluator_add_true(e);
  ensure(had_err, evaluator_get_sensitivity(e) == 0.5);
  ensure(had_err, evaluator_get_specificity(e) == 0.5);

  evaluator_reset(e);
  ensure(had_err, evaluator_get_sensitivity(e) == 1.0);
  ensure(had_err, evaluator_get_specificity(e) == 1.0);

  evaluator_add_actual(e, 4);
  evaluator_add_predicted(e, 4);
  evaluator_add_true(e);
  ensure(had_err, evaluator_get_sensitivity(e) == 0.25);
  ensure(had_err, evaluator_get_specificity(e) == 0.25);
  evaluator_add_true(e);
  ensure(had_err, evaluator_get_sensitivity(e) == 0.5);
  ensure(had_err, evaluator_get_specificity(e) == 0.5);

  evaluator_add_true(e);
  ensure(had_err, evaluator_get_sensitivity(e) == 0.75);
  ensure(had_err, evaluator_get_specificity(e) == 0.75);

  evaluator_add_true(e);
  ensure(had_err, evaluator_get_sensitivity(e) == 1.0);
  ensure(had_err, evaluator_get_specificity(e) == 1.0);

  evaluator_delete(e, env);

  return had_err;
}

void evaluator_delete(Evaluator *e, Env *env)
{
  if (!e) return;
  ma_free(e);
}
