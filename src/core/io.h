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
#include "core/str.h"

/* the I/O class */
typedef struct GT_IO GT_IO;

GT_IO*        gt_io_new(const char *path, const char *mode);
/* Returns -1 if no char is left, 0 otherwise. */
int           gt_io_get_char(GT_IO*, char*);
/* Can only be used once at a time.*/
void          gt_io_unget_char(GT_IO*, char);
bool          gt_io_line_start(const GT_IO*);
bool          gt_io_has_char(GT_IO*);
char          gt_io_peek(GT_IO*);
char          gt_io_next(GT_IO*);
unsigned long gt_io_get_line_number(const GT_IO*);
const char*   gt_io_get_filename(const GT_IO*);
GtStr*       gt_io_get_filename_str(const GT_IO*);
void          gt_io_delete(GT_IO*);

#endif
