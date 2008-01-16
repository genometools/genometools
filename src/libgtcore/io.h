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

#ifndef IO_H
#define IO_H

#include <stdbool.h>
#include <stdio.h>

/* the I/O class */
typedef struct IO IO;

IO*           io_new(const char *path, const char *mode);
int           io_get_char(IO*, char*); /* returns -1 if no char is left,
                                          0 otherwise */
void          io_unget_char(IO*, char); /* can only be used once at a time */
bool          io_line_start(const IO*);
bool          io_has_char(IO*);
unsigned long io_get_line_number(const IO*);
const char*   io_get_filename(const IO*);
void          io_delete(IO*);

#endif
