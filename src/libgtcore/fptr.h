/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef FPTR_H
#define FPTR_H

#include <libgtcore/env.h>

/* the generic function pointers */
typedef int  (*Compare)(const void*, const void*);
typedef int  (*CompareWithData)(const void*, const void*, const void*);
typedef void (*FreeFunc)(void*, Env*);

#endif
