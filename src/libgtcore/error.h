/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg

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
void        error_unset(Error*);
/* get the error string (the error must be set) */
const char* error_get(const Error*);
void        error_set_progname(Error*, const char *progname);
const char* error_get_progname(const Error*);
void        error_delete(Error*);

/* make sure that the error is not set, should be used at the beginning of
   every routine which has an Error* argument */

#define error_check(e)\
        assert(!e|| !error_is_set(e))

#endif
