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

#include <stdbool.h>
#include "core/array2dim_api.h"
#include "core/divmodmul.h"
#include "core/minmax.h"
#include "gth/align_common.h"
#include "gth/align_dna_imp.h"
#include "gth/path_matrix.h"

typedef struct {
  bool used;
  PATHTYPE e_entry,
           i_entry;
} PMEntry;

struct GthPathMatrix {
  GtRange gen_range,
          ref_range;
  PMEntry **entries;
};

static void path_matrix_fill(GthPathMatrix *pm, PATHTYPE **path)
{
  unsigned long genptr, refptr;
  PATHTYPE e_path, i_path;
  bool lower;

  for (genptr = pm->gen_range.start; genptr <= pm->gen_range.end; genptr++) {
    for (refptr = pm->ref_range.start; refptr <= pm->ref_range.end; refptr++) {
      e_path = path[GT_DIV2(genptr)][refptr];
      i_path = path[GT_DIV2(genptr)][refptr];
      lower = (bool) !GT_MOD2(genptr);
      if (lower) {
        e_path &= LOWER_E_STATE_MASK;
        i_path &= LOWER_I_STATE_MASK;
        i_path >>= 3;
      }
      else {
        e_path &= UPPER_E_STATE_MASK;
        e_path >>= 4;
        i_path &= UPPER_I_STATE_MASK;
        i_path >>= 7;
      }

    }
  }

}

GthPathMatrix* gth_path_matrix_new(PATHTYPE **path,
                                   unsigned long gen_dp_length,
                                   unsigned long ref_dp_length,
                                   const GtRange *btmatrixgenrange,
                                   const GtRange *btmatrixrefrange)
{
  GthPathMatrix *pm = gt_malloc(sizeof *pm);
  pm->gen_range.start = btmatrixgenrange->start;
  pm->gen_range.end   = MIN(btmatrixgenrange->end, gen_dp_length);
  pm->ref_range.start = btmatrixrefrange->start;
  pm->ref_range.end   = MIN(btmatrixrefrange->end, ref_dp_length);
  gt_array2dim_calloc(pm->entries, gt_range_length(&pm->ref_range),
                      gt_range_length(&pm->gen_range));
  path_matrix_fill(pm, path);
  return pm;
}

void gth_path_matrix_show(GthPathMatrix *pm)
{
  gt_assert(pm);
}

void gth_path_matrix_delete(GthPathMatrix *pm)
{
  if (!pm) return;
  gt_array2dim_delete(pm->entries);
  gt_free(pm);
}
