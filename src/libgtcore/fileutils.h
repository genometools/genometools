/*
  Copyright (c) 2006-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c)      2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
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

#ifndef FILEUTILS_H
#define FILEUTILS_H

#include <stdbool.h>
#include <sys/types.h>
#include <stdio.h>
#include "libgtcore/str.h"
#include "libgtcore/strarray.h"

bool           file_exists(const char*);
/* Returns true if the file with path <a> has a later modification time than the
   file with path <b>, false otherwise. */
bool           file_is_newer(const char *a, const char *b);
unsigned long  file_number_of_lines(const char*);
/* Returns the suffix of <path>, if there is any. Returns "" otherwise
   The suffix is the part after and including the last '.' but after the last
   '/'. Except if <path> ends with ".gz" or ".bz2", then the suffix is the part
   after and including the second last '.'. */
const char*    file_suffix(const char *path);

/* Set <path> to the dirname of <file>, if it has one, to "" otherwise. */
void           file_dirname(Str *path, const char *file);

/* Find <file> in $PATH, if it has no dirname; set <path> to dirname otherwise.
   Sets <path> to the empty string if <file> could not be found in $PATH. */
int            file_find_in_path(Str *path, const char *file, Error*);

off_t          files_estimate_total_size(const StrArray *filenames);

#endif
