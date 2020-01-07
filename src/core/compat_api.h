/*
  Copyright (c) 2013 Gordon Gremme <gordon@gremme.org>

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

#ifndef COMPAT_API_H
#define COMPAT_API_H

#include "core/types_api.h"

/* Compat module */

#ifndef _WIN32
/* Path separator, as char. */
#define GT_PATH_SEPARATOR     '/'
/* Path separator, as string. */
#define GT_PATH_SEPARATOR_STR "/"
/* Path component separator, as char. */
#define GT_PATH_VAR_SEPARATOR ':'
#else
/* Path separator, as char. */
#define GT_PATH_SEPARATOR     '\\'
/* Path separator, as string. */
#define GT_PATH_SEPARATOR_STR "\\"
/* Path component separator, as char. */
#define GT_PATH_VAR_SEPARATOR ';'
#endif

/* Return (read-write) handle of temporary file, with template <templ>. */
int     gt_mkstemp(char *templ);
/* Returns the page size of the current platform. */
GtUword gt_pagesize(void);

#endif
