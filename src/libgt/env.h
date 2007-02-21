/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef ENV_H
#define ENV_H

#include "error.h"
#include "log.h"
#include "ma.h"

/* the enviroment class (creates and holds all singular objects) */
typedef struct Env Env;

Env*   env_new(void);
MA*    env_ma(const Env*);    /* return the memory allocator */
Error* env_error(const Env*); /* return the error object */
Log*   env_log(const Env*);   /* return the log object or NULL */
void   env_set_log(Env*, Log*);
void   env_delete(Env*);

/* wrapper for error functions (to make many env_error() calls unnecessary */
void    env_error_set(Env*, const char *format, ...)
          __attribute__ ((format (printf, 2, 3)));
#define env_error_is_set(env)\
        error_is_set(env_error(env))
#define env_error_unset(env)\
        error_unset(env_error(env))
#define env_error_get(env)\
        error_get(env_error(env))

/* make sure that the error is not set, has to be used at the beginning of
   every routine which has an Env* argument */
#define env_error_check(env)\
        assert(!env || !error_is_set(env_error(env)))

#endif
