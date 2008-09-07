/*
  Copyright (c) 2006-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2008 Center for Bioinformatics, University of Hamburg

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

#ifndef ERROR_API_H
#define ERROR_API_H

#include <assert.h>
#include <stdarg.h>
#include <stdbool.h>

/* the error class */
typedef struct GT_Error GT_Error;

GT_Error*   gt_error_new(void);
void        gt_error_set(GT_Error*, const char *format, ...)
              __attribute__ ((format (printf, 2, 3)));
void        gt_error_vset(GT_Error*, const char *format, va_list);
bool        gt_error_is_set(const GT_Error*);
void        gt_error_unset(GT_Error*);
/* Get the error string (the error must be set). */
const char* gt_error_get(const GT_Error*);
void        gt_error_delete(GT_Error*);

/* Make sure that the error is not set, should be used at the beginning of
   every routine which has an GT_Error* argument. */
#define gt_error_check(err)\
        assert(!err || !gt_error_is_set(err))

#endif
