/*
  Copyright (c) 2011 Gordon Gremme <gremme@zbh.uni-hamburg.de>

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

#ifndef PATH_MATRIX_H
#define PATH_MATRIX_H

#include "core/array2dim_sparse.h"

typedef struct GthPathMatrix GthPathMatrix;

GthPathMatrix* gth_path_matrix_new(GthPath **path,
                                   unsigned long gen_dp_length,
                                   unsigned long ref_dp_length,
                                   const GtRange *btmatrixgenrange,
                                   const GtRange *btmatrixrefrange,
                                   GtRowInfo *ri);
void           gth_path_matrix_show(GthPathMatrix*);
void           gth_path_matrix_set_max_path(GthPathMatrix*,
                                            unsigned long genptr,
                                            unsigned long refptr,
                                            bool e_state);
void           gth_path_matrix_delete(GthPathMatrix*);

#endif
