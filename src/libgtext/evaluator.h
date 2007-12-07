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

#ifndef EVALUATOR_H
#define EVALUATOR_H

#include <stdio.h>
#include "libgtcore/error.h"

typedef struct Evaluator Evaluator;

Evaluator* evaluator_new(void);
void       evaluator_add_true(Evaluator*);
void       evaluator_add_actual(Evaluator*, unsigned long);
void       evaluator_add_predicted(Evaluator*, unsigned long);
double     evaluator_get_sensitivity(const Evaluator*);
double     evaluator_get_specificity(const Evaluator*);
void       evaluator_show_sensitivity(const Evaluator*, FILE*);
void       evaluator_show_specificity(const Evaluator*, FILE*);
void       evaluator_reset(Evaluator*);
int        evaluator_unit_test(Error*);
void       evaluator_delete(Evaluator*);

#endif
