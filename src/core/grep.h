/*
  Copyright (c) 2014 Sascha Steinbiss <ss34@sanger.ac.uk>
  Copyright (c) 2014 Genome Research Ltd.

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

#ifndef GREP_H
#define GREP_H

#include <stdbool.h>
#include "core/error_api.h"
#include "core/grep_api.h"
#include "core/str_api.h"

/* Escapes the POSIX extended regular expression special characters
   (.^$*+?()[{\|) in <str> of length <len> and writes the escaped string
   to <dest>, which is truncated first.
   Any special characters appearing in <str> will be treated as literals when
   passing <dest> as a pattern to <gt_grep_*()>. */
void gt_grep_escape_extended(GtStr *dest, const char *str, size_t len);

#endif
