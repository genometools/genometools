/*
  Copyright (c) 2005-2007, 2013 Gordon Gremme <gordon@gremme.org>
  Copyright (c) 2005-2007       Center for Bioinformatics, University of Hamburg

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

#ifndef BOOL_MATRIX_H
#define BOOL_MATRIX_H

#include <stdbool.h>
#include "core/error.h"
#include "core/types_api.h"

/* a two-dimensional matrix containing boolean values */
typedef struct GtBoolMatrix GtBoolMatrix;

GtBoolMatrix* gt_bool_matrix_new(void);
bool          gt_bool_matrix_get(GtBoolMatrix*, GtUword firstdim,
                                 GtUword seconddim);
void          gt_bool_matrix_set(GtBoolMatrix*, GtUword firstdim,
                                 GtUword seconddim, bool);
GtUword       gt_bool_matrix_get_first_column(const GtBoolMatrix*,
                                              GtUword firstdim);
GtUword       gt_bool_matrix_get_last_column(const GtBoolMatrix*,
                                             GtUword firstdim);
GtUword       gt_bool_matrix_get_next_column(const GtBoolMatrix*,
                                             GtUword firstdim,
                                             GtUword i);
void          gt_bool_matrix_delete(GtBoolMatrix*);

#endif
