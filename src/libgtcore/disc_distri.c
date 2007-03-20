/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include <stdio.h>
#include <libgtcore/array.h>
#include <libgtcore/disc_distri.h>
#include <libgtcore/xansi.h>

struct DiscDistri {
  Array *values;
  unsigned long num_of_occurrences;
};

DiscDistri* disc_distri_new(Env *env)
{
  return env_ma_calloc(env, 1, sizeof (DiscDistri));
}

void disc_distri_add(DiscDistri *d, unsigned long value, Env *env)
{
  unsigned long *distri, zero = 0;
  assert(d);

  if (!d->values)
    d->values = array_new(sizeof (unsigned long), env);

  while (array_size(d->values) <= value)
    array_add(d->values, zero, env);

  distri = array_get_space(d->values);
  distri[value]++;
  d->num_of_occurrences++;
}

unsigned long disc_distri_get(const DiscDistri *d, unsigned long value)
{
  assert(d);
  if (value < array_size(d->values))
    return 0;
  return *(unsigned long*) array_get(d->values, value);
}

void disc_distri_show(const DiscDistri *d)
{
  assert(d);
  disc_distri_show_generic(d, NULL);
}

void disc_distri_show_generic(const DiscDistri *d, GenFile *genfile)
{
  unsigned long value, occurrences;
  double probability, cumulative_probability = 0.0;
  assert(d);

  for (value = 0; value < array_size(d->values); value++) {
    occurrences = *(unsigned long*) array_get(d->values, value);
    probability = (double) occurrences / d->num_of_occurrences;
    cumulative_probability += probability;
    if (occurrences)
      genfile_xprintf(genfile, "%lu: %lu (prob=%.4f,cumulative=%.4f)\n", value,
                      occurrences, probability, cumulative_probability);
  }
}

void disc_distri_delete(DiscDistri *d, Env *env)
{
  if (!d) return;
  array_delete(d->values, env);
  env_ma_free(d, env);
}
