/*
  Copyright (c) 2011 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2011 Center for Bioinformatics, University of Hamburg

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

#ifndef DESC_BUFFER_H
#define DESC_BUFFER_H

#include <stdio.h>
#include "core/error_api.h"

typedef struct GtDescBuffer GtDescBuffer;

/* Return an empty <GtDescBuffer*> object. */
GtDescBuffer* gt_desc_buffer_new(void);
/* Increase the reference count for <db> and return it. */
GtDescBuffer* gt_desc_buffer_ref(GtDescBuffer *db);
unsigned long gt_desc_buffer_length(const GtDescBuffer *s);
const char*   gt_desc_buffer_get_next(GtDescBuffer *db);
/* Append character <c> to <db>. */
void          gt_desc_buffer_append_char(GtDescBuffer *db, char c);
void          gt_desc_buffer_finish(GtDescBuffer *db);
/* Reset <db> to length 0. */
void          gt_desc_buffer_reset(GtDescBuffer *db);
/* Decrease the reference count for <db> or delete it, if this was the last
   reference. */
void          gt_desc_buffer_delete(GtDescBuffer *db);

int           gt_desc_buffer_unit_test(GtError *err);

#endif
