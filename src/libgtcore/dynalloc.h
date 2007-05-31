/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef DYNALLOC_H
#define DYNALLOC_H

#include <libgtcore/env.h>

/*
  Do not use this function directly! It is only used internally to implement
  dynamic arrays, dynamic strings, and the like.

  It dynamically re-/allocates memory. Thereby, usually more memory
  is allocated than what was asked for (to avoid frequent realloc calls).

  ptr: the previously allocated memory
  allocated: allocated memory size before and after the call
  size: requested memory size
  env: the environment object
*/

void* dynalloc(void *ptr, size_t *allocated, size_t size, Env *env);

#endif
