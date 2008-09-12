/*
  Copyright (c) 2007-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007-2008 Center for Bioinformatics, University of Hamburg

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

#ifndef STRARRAY_H
#define STRARRAY_H

#include "core/str.h"

/* the string array class  */
typedef struct GtStrArray GtStrArray;

GtStrArray*  gt_strarray_new(void);
GtStrArray*  gt_strarray_new_file(const char *path);
void          gt_strarray_add_cstr(GtStrArray*, const char*);
void          gt_strarray_add_cstr_nt(GtStrArray*, const char*,unsigned long);
void          gt_strarray_add(GtStrArray*, const GtStr*);
const char*   gt_strarray_get(const GtStrArray*, unsigned long strnum);
/* Returns an internal GtStr pointer (i.e., _not_ a new reference!). */
GtStr*       gt_strarray_get_str(const GtStrArray*, unsigned long strnum);
void          gt_strarray_set_size(GtStrArray*, unsigned long);
/* Returns number of strings. */
unsigned long gt_strarray_size(const GtStrArray*);
void          gt_strarray_delete(GtStrArray*);

#endif
