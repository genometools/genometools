/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef CSTR_H
#define CSTR_H

#include <stdio.h>
#include "env.h"

void          cstr_show(const char*, unsigned long length, FILE*);

/* use p and a blank as prefix for cstr_array[0] and return the result */
char**        cstr_array_prefix_first(char **cstr_array, const char *p);
unsigned long cstr_array_size(char **cstr_array); /* O(n) */
void          cstr_array_delete(char **cstr_array, Env*);

#endif
