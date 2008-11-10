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

#define GT_CARRIAGE_RETURN  '\r'
#define GT_END_OF_LINE      '\n'
#define GT_END_OF_FILE      EOF

/* the I/O class */
typedef struct GtIO GtIO;

GtIO*         gt_io_new(const char *path, const char *mode);
void          gt_io_delete(GtIO*);
/* Returns -1 if no char is left, 0 otherwise. */
int           gt_io_get_char(GtIO*, char*);
/* Can only be used once at a time.*/
void          gt_io_unget_char(GtIO*, char);
bool          gt_io_line_start(const GtIO*);
bool          gt_io_has_char(GtIO*);
char          gt_io_peek(GtIO*);
char          gt_io_next(GtIO*);
unsigned long gt_io_get_line_number(const GtIO*);
const char*   gt_io_get_filename(const GtIO*);
GtStr*        gt_io_get_filename_str(const GtIO*);

int           gt_io_expect(GtIO*, char expected_char, GtError*);

#endif
