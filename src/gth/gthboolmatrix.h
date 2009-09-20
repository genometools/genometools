/*
  Copyright (c) 2005-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2005-2007 Center for Bioinformatics, University of Hamburg

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

#ifndef GTHBOOLMATRIX_H
#define GTHBOOLMATRIX_H

#include <stdbool.h>
#include "core/error.h"

/* a two-dimensional matrix containing boolean values */
typedef struct GthBoolMatrix GthBoolMatrix;

GthBoolMatrix* gthboolmatrix_new(void);
bool           gthboolmatrix_get(GthBoolMatrix*, unsigned long firstdim,
                                 unsigned long seconddim);
void           gthboolmatrix_set(GthBoolMatrix*, unsigned long firstdim,
                                 unsigned long seconddim, bool);
void           gthboolmatrix_delete(GthBoolMatrix*);

#endif
