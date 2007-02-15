/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef ERROR_H
#define ERROR_H

#include <assert.h>
#include <stdbool.h>

/* simple error interface: print an error and exit */
void error(const char *format, ...)
  __attribute__ ((format (printf, 1, 2)));

/* sophisticated error interface: full-fledged error class */
typedef struct Error Error;

Error*      error_new(void);
void        error_set(Error *, const char *format, ...)
              __attribute__ ((format (printf, 2, 3)));
bool        error_is_set(const Error*);
/* make sure that the error is not set, has to be used at the beginning of
   every routine which has an Error* argument */
#define     error_check(e)\
            assert(!e || !error_is_set(e))
void        error_unset(Error*);
/* get the error string (the error must be set) */
const char* error_get(const Error*);
/* aborts and shows error string, if the error is set */
void        error_abort(const Error*);
void        error_free(Error*);

#endif
