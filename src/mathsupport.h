/*
  Copyright (c) 2006 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef MATHSUPPORT_H
#define MATHSUPPORT_H

/* returns the log of the sum of two log probabilities */
double       logsum(double p1, double p2);
unsigned int double_equals_one(double);
unsigned int double_equals_double(double, double);

/* returns a random number between 0 and maximal_value (employs rand(3)) */
unsigned long rand_max(unsigned long maximal_value);
/* returns a random double between 0.0 and maximal_value (employs rand(3)) */
double rand_max_double(double maximal_value);
/* returns a random double between 0.0 and 1.0 (employs rand(3)) */
double rand_0_to_1(void);

#endif
