/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef LINEAREDIST_H
#define LINEAREDIST_H

#include <gtcore.h>

/* Compute the edit distance of sequences u and v in O(max{|u|,|v|}) space */
unsigned long linearedist(const char *u, unsigned long n,
                          const char *v, unsigned long m, Env*);

#endif
