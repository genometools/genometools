/*
  Copyright (c) 2006-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2008 Center for Bioinformatics, University of Hamburg

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
#include <math.h>
#include <stdlib.h>
#include "libgtcore/mathsupport.h"
#include "libgtcore/yarandom.h"

#define EPSILON  0.0000000000001

double logsum(double p1, double p2)
{
  if (p1 > p2)
    return (p1-p2 > 50.0) ? p1 : p1 + log(1.0 + exp(p2-p1));
  return (p2-p1 > 50.0) ? p2 : p2 + log(1.0 + exp(p1-p2));
}

bool double_equals_one(double d)
{
  if (fabs(1.0 - d) <= EPSILON)
    return true;
  return false;
}

bool double_equals_double(double d1, double d2)
{
  if (fabs(d1 - d2) <= EPSILON)
    return true;
  return false;
}

unsigned long rand_max(unsigned long maximal_value)
{
  unsigned long r;
  assert(maximal_value);
  r = ((double) random() / ((double) RAND_MAX + 1) * (maximal_value + 1));
  assert(r <= maximal_value);
  return r;
}

double rand_max_double(double maximal_value)
{
  double r;
  assert(maximal_value);
  r = ((double) random() / RAND_MAX) * maximal_value; /* XXX */
  assert(r >= 0.0 && r <= maximal_value);
  return r;
}

double rand_0_to_1(void)
{
  double r;
  r = (double) random() / RAND_MAX;
  assert(r >= 0.0 && r <= 1.0);
  return r;
}

char rand_char(void)
{
  int offset;
  char c;
  offset = rand_max(25);
  c = 97 + offset;
  assert(c >= 'a' && c <= 'z');
  return c;
}
