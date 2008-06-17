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

#ifndef MATHSUPPORT_H
#define MATHSUPPORT_H

#include <stdbool.h>

/* Returns the log of the sum of two log probabilities. */
double        logsum(double p1, double p2);
bool          double_equals_one(double);
bool          double_equals_double(double, double);

/* Returns a random number between 0 and maximal_value. */
unsigned long rand_max(unsigned long maximal_value);
/* Returns a random double between 0.0 and maximal_value. */
double        rand_max_double(double maximal_value);
/* Returns a random double between 0.0 and 1.0. */
double        rand_0_to_1(void);
/* Returns a random character from 'a' to 'z'. */
char          rand_char(void);

#endif
