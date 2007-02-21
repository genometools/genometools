/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef ERROR_H
#define ERROR_H

#include <assert.h>
#include <stdarg.h>
#include <stdbool.h>

/* the error class */
typedef struct Error Error;

Error*      error_new(void);
void        error_set(Error*, const char *format, ...)
              __attribute__ ((format (printf, 2, 3)));
void        error_vset(Error*, const char *format, va_list);
bool        error_is_set(const Error*);
/* make sure that the error is not set, has to be used at the beginning of
   every routine which has an Error* argument */
#define     error_check(e)\
            assert(!e || !error_is_set(e))
void        error_unset(Error*);
/* get the error string (the error must be set) */
const char* error_get(const Error*);
void        error_delete(Error*);

#endif
