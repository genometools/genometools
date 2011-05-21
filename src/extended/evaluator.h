/*
  Copyright (c) 2006-2007, 2011 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007       Center for Bioinformatics, University of Hamburg

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

#ifndef EVALUATOR_H
#define EVALUATOR_H

#include "core/error_api.h"
#include "core/file.h"

typedef struct GtEvaluator GtEvaluator;

GtEvaluator* gt_evaluator_new(void);
void         gt_evaluator_add_true(GtEvaluator*);
void         gt_evaluator_add_actual(GtEvaluator*, unsigned long);
void         gt_evaluator_add_predicted(GtEvaluator*, unsigned long);
double       gt_evaluator_get_sensitivity(const GtEvaluator*);
double       gt_evaluator_get_specificity(const GtEvaluator*);
void         gt_evaluator_show_sensitivity(const GtEvaluator*, GtFile*);
void         gt_evaluator_show_specificity(const GtEvaluator*, GtFile*);
void         gt_evaluator_reset(GtEvaluator*);
int          gt_evaluator_unit_test(GtError*);
void         gt_evaluator_delete(GtEvaluator*);

#endif
