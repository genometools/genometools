/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

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

#ifndef ENV_H
#define ENV_H

/* the enviroment class (creates and holds all singular objects) */
typedef struct Env Env;

#include "libgtcore/error.h"

Env*   env_new(void);
Error* env_error(const Env*); /* return the error object */
void   env_set_spacepeak(Env*, bool);
/* returns 0 if no memory map, file pointer, or memory has been leaked and a
   value != 0 otherwise */
int    env_delete(Env*);

/* wrapper for error functions */
void    env_error_set(Env*, const char *format, ...)
          __attribute__ ((format (printf, 2, 3)));
#define env_error_is_set(env)\
        error_is_set(env_error(env))
#define env_error_unset(env)\
        error_unset(env_error(env))
#define env_error_get(env)\
        error_get(env_error(env))
#define env_error_set_progname(env, progname)\
        error_set_progname(env_error(env), progname, env)
#define env_error_get_progname(env)\
        error_get_progname(env_error(env))
/* make sure that the error is not set, should be used at the beginning of
   every routine which has an Env* argument, except for destructors! */
#define env_error_check(env)\
        assert(!env || !env_error(env) || !error_is_set(env_error(env)))

#endif
