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

#include "libgtcore/ma.h"
#include "libgtcore/error.h"
#include "libgtcore/fa.h"
#include "libgtcore/log.h"

Env*   env_new(void);
MA*    env_ma(const Env*);    /* return the memory allocator */
FA*    env_fa(const Env*);    /* return the file allocator */
Error* env_error(const Env*); /* return the error object */
void   env_set_spacepeak(Env*, bool);
/* returns 0 if no memory map, file pointer, or memory has been leaked and a
   value != 0 otherwise */
int    env_delete(Env*);

/* wrapper for memory functions */
#define env_ma_malloc(env, size)\
        ma_malloc_mem(env_ma(env), size, __FILE__, __LINE__)
#define env_ma_calloc(env, nmemb, size)\
        ma_calloc_mem(env_ma(env), nmemb, size, __FILE__, __LINE__)
#define env_ma_realloc(env, ptr, size)\
        ma_realloc_mem(env_ma(env), ptr, size, __FILE__, __LINE__)
#define env_ma_free(ptr, env)\
        ma_free_mem(ptr, env_ma(env), __FILE__, __LINE__)
/* free functions get the data object (here the env object) _always_ as the
   last object */
void    env_ma_free_func(void *ptr, Env*);
#define env_ma_get_space_peak(env)\
        ma_get_space_peak(env_ma(env))

/* wrapper for file functions */
#define env_fa_fopen(env, path, mode)\
        fa_fopen(env_fa(env), path, mode, __FILE__, __LINE__)
#define env_fa_xfopen(env, path, mode)\
        fa_xfopen(env_fa(env), path, mode, __FILE__, __LINE__)
void    env_fa_fclose(FILE *stream, Env*);
void    env_fa_xfclose(FILE *stream, Env*);

#define env_fa_gzopen(env, path, mode)\
        fa_gzopen(env_fa(env), path, mode, __FILE__, __LINE__)
#define env_fa_xgzopen(env, path, mode)\
        fa_xgzopen(env_fa(env), path, mode, __FILE__, __LINE__)
void    env_fa_gzclose(gzFile stream, Env*);
void    env_fa_xgzclose(gzFile stream, Env*);

#define env_fa_bzopen(env, path, mode)\
        fa_bzopen(env_fa(env), path, mode, __FILE__, __LINE__)
#define env_fa_xbzopen(env, path, mode)\
        fa_xbzopen(env_fa(env), path, mode, __FILE__, __LINE__)
void    env_fa_bzclose(BZFILE *stream, Env*);
void    env_fa_xbzclose(BZFILE *stream, Env*);

#define env_fa_xtmpfile(env, template)\
        fa_xtmpfile(env_fa(env), template, __FILE__, __LINE__)

#define env_fa_mmap_read(env, filename, len)\
        fa_mmap_read(env_fa(env), filename, len, __FILE__, __LINE__)
#define env_fa_mmap_write(env, filename, len)\
        fa_mmap_write(env_fa(env), filename, len, __FILE__, __LINE__)
#define env_fa_xmmap_read(env, filename, len)\
        fa_xmmap_read(env_fa(env), filename, len, __FILE__, __LINE__)
#define env_fa_xmmap_write(env, filename, len)\
        fa_xmmap_write(env_fa(env), filename, len, __FILE__, __LINE__)
void    env_fa_xmunmap(void *addr, Env*);

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

/* wrapper for log functions */
void    env_log_log(Env*, const char *format, ...)
          __attribute__ ((format (printf, 2, 3)));

#endif
