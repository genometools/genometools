/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef EVALUATOR_H
#define EVALUATOR_H

#include <stdio.h>
#include <libgt/env.h>

typedef struct Evaluator Evaluator;

Evaluator* evaluator_new(Env*);
void       evaluator_add_true(Evaluator*);
void       evaluator_add_actual(Evaluator*, unsigned long);
void       evaluator_add_predicted(Evaluator*, unsigned long);
double     evaluator_get_sensitivity(const Evaluator*);
double     evaluator_get_specificity(const Evaluator*);
void       evaluator_show_sensitivity(const Evaluator*, FILE*);
void       evaluator_show_specificity(const Evaluator*, FILE*);
void       evaluator_reset(Evaluator*);
int        evaluator_unit_test(Env*);
void       evaluator_delete(Evaluator*, Env*);

#endif
