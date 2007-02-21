/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include <stdio.h>
#include "array.h"
#include "disc_distri.h"
#include "xansi.h"

struct DiscDistri {
  Array *values;
  unsigned long num_of_occurrences;
};

DiscDistri* disc_distri_new(void)
{
  return xcalloc(1, sizeof (DiscDistri));
}

void disc_distri_add(DiscDistri *d, unsigned long value)
{
  unsigned long *distri, zero = 0;
  assert(d);

  if (!d->values)
    d->values = array_new(sizeof (unsigned long));

  while (array_size(d->values) <= value)
    array_add(d->values, zero);

  distri = array_get_space(d->values);
  distri[value]++;
  d->num_of_occurrences++;
}

void disc_distri_show(const DiscDistri *d)
{
  unsigned long value, occurrences;
  double probability, cumulative_probability = 0.0;
  assert(d);

  for (value = 0; value < array_size(d->values); value++) {
    occurrences = *(unsigned long*) array_get(d->values, value);
    probability = (double) occurrences / d->num_of_occurrences;
    cumulative_probability += probability;
    if (occurrences)
      printf("%lu: %lu (prob=%.4f,cumulative=%.4f)\n", value, occurrences,
             probability, cumulative_probability);
  }
}

void disc_distri_delete(DiscDistri *d)
{
  if (!d) return;
  array_delete(d->values);
  free(d);
}
