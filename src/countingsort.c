/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include <string.h>
#include "countingsort.h"
#include "ensure.h"
#include "xansi.h"

void countingsort(void *out, const void *in, size_t elem_size,
                  unsigned long size, unsigned long max_elemvalue, void *data,
                  unsigned long (*get_elemvalue)(const void *elem, void *data))
{
  unsigned long i, k, *c;
  assert(out && in && elem_size && size && max_elemvalue && get_elemvalue);

  /* allocate count array */
  c = xcalloc(sizeof(unsigned long), max_elemvalue + 1);

  /* count number of elements of a given value */
  for (i = 0; i < size; i++) {
    k = get_elemvalue(in + elem_size * i, data);
    assert(k <= max_elemvalue);
    c[k]++;
  }

  /* compute running sum of the count array */
  for (i = 1; i <= max_elemvalue; i++)
    c[i] += c[i-1];

  /* sorting (stable) */
  for (i = size; i > 0; i--) {
    k = get_elemvalue(in + elem_size * (i-1), data);
    memcpy(out + elem_size * (c[k] - 1), in + elem_size * (i-1), elem_size);
    c[k]--;
  }

  free(c);
}

unsigned long countingsort_get_max(const void *in, size_t elem_size,
                                   unsigned long size, void *data,
                                   unsigned long (*get_elemvalue)
                                                 (const void *elem, void *data))
{
  unsigned long i, value, max_value = 0;
  for (i = 0; i < size; i++) {
    value = get_elemvalue(in + elem_size * i, data);
    if (value > max_value)
      max_value = value;
  }
  return max_value;
}

static unsigned long get_int(const void *elem, void *data)
{
  assert(elem);
  return *(unsigned int*) elem;
}

int countingsort_unit_test(void)
{
  unsigned int numbers[]        = { 1, 2, 1, 2, 0 }, numbers_out[5],
               sorted_numbers[] = { 0, 1, 1, 2, 2 };
  countingsort(numbers_out, numbers, sizeof(unsigned int), 5,
               countingsort_get_max(numbers, sizeof(unsigned int), 5, NULL,
                                    get_int), NULL,  get_int);
  ensure(!memcmp(sorted_numbers, numbers_out, sizeof(unsigned int) * 5));
  return EXIT_SUCCESS;
}
