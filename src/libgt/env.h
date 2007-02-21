/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef ENV_H
#define ENV_H

#include "error.h"
#include "ma.h"

/* the enviroment class (creates and holds all singular objects) */
typedef struct Env Env;

Env*   env_new(void);
MA*    env_ma(const Env*);    /* return the memory allocator */
Error* env_error(const Env*); /* return the error object */
void   env_delete(Env*);

#endif
