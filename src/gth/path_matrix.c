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
  bool used, on_max_path_e, on_max_path_i;
  GthPath e_path,
          i_path;
} PMEntry;

struct GthPathMatrix {
  GtRange gen_range,
          ref_range;
  PMEntry **entries;
};

static void path_matrix_fill(GthPathMatrix *pm, GthPath **path, GtRowInfo *ri)
{
  unsigned long genptr, refptr, genidx, refidx;
  GthPath e_path, i_path;
  bool lower;

  for (genptr = pm->gen_range.start; genptr <= pm->gen_range.end; genptr++) {
    for (refptr = pm->ref_range.start; refptr <= pm->ref_range.end; refptr++) {
      if (!ri ||
          (refptr >= ri[genptr].offset &&
           refptr < ri[genptr].offset + ri[genptr].length)) {
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
        gt_assert(i_path == DNA_E_N || i_path == DNA_I_N);
        genidx = genptr - pm->gen_range.start;
        refidx = refptr - pm->ref_range.start;
        pm->entries[genidx][refidx].used = true;
        pm->entries[genidx][refidx].e_path = e_path;
        pm->entries[genidx][refidx].i_path = i_path;
      }
    }
  }

}

GthPathMatrix* gth_path_matrix_new(GthPath **path,
                                   unsigned long gen_dp_length,
                                   unsigned long ref_dp_length,
                                   const GtRange *btmatrixgenrange,
                                   const GtRange *btmatrixrefrange,
                                   GtRowInfo *ri)
{
  GthPathMatrix *pm = gt_malloc(sizeof *pm);
  pm->gen_range.start = btmatrixgenrange->start;
  pm->gen_range.end   = MIN(btmatrixgenrange->end, gen_dp_length);
  pm->ref_range.start = btmatrixrefrange->start;
  pm->ref_range.end   = MIN(btmatrixrefrange->end, ref_dp_length);
  gt_array2dim_calloc(pm->entries, gt_range_length(&pm->gen_range),
                      gt_range_length(&pm->ref_range));
  path_matrix_fill(pm, path, ri);
  return pm;
}

static char on_max_path_char(bool on_max_path)
{
  return on_max_path ? '*' : ' ';
}

static char matrix_char(GthPath path)
{
  switch (path) {
    case DNA_E_N:
    case DNA_E_NM:
    case DNA_E_M:
      return 'E';
    case DNA_I_N:
    case DNA_I_NM:
    case DNA_I_M:
      return 'I';
    default:
      gt_assert(0);
      return 'X'; /* cannot happen */
  }
}

static char direction_char(GthPath path)
{
  switch (path) {
    case DNA_E_N:
    case DNA_I_N:
      return '-';
    case DNA_E_NM:
    case DNA_I_NM:
      return '\\';
    case DNA_E_M:
    case DNA_I_M:
      return '|';
    default:
      gt_assert(0);
      return 'X'; /* cannot happen */
  }
}

void gth_path_matrix_show(GthPathMatrix *pm)
{
  unsigned long genptr, refptr, genidx, refidx;
  PMEntry entry;
  gt_assert(pm);
  printf("    ");
  for (genptr = pm->gen_range.start; genptr <= pm->gen_range.end; genptr++)
    printf("%4lu", genptr);
  printf("\n\n");

  for (refptr = pm->ref_range.start; refptr <= pm->ref_range.end; refptr++) {
    refidx = refptr - pm->ref_range.start;
    printf("%4lu", refptr);
    for (genptr = pm->gen_range.start; genptr <= pm->gen_range.end; genptr++) {
      genidx = genptr - pm->gen_range.start;
      entry = pm->entries[genidx][refidx];
      if (entry.used) {
        printf(" %c%c%c", on_max_path_char(entry.on_max_path_e),
               matrix_char(entry.e_path), direction_char(entry.e_path));
      }
      else
        printf("    ");
    }
    printf("\n    ");
    for (genptr = pm->gen_range.start; genptr <= pm->gen_range.end; genptr++) {
      genidx = genptr - pm->gen_range.start;
      entry = pm->entries[genidx][refidx];
      if (entry.used) {
        printf(" %c%c%c", on_max_path_char(entry.on_max_path_i),
               matrix_char(entry.i_path), direction_char(entry.i_path));
      }
      else
        printf("    ");
    }
    putchar('\n');
  }
}

void gth_path_matrix_set_max_path(GthPathMatrix *pm, unsigned long genptr,
                                  unsigned long refptr, bool e_state)
{
  unsigned long genidx, refidx;
  gt_assert(pm);
  if (genptr >= pm->gen_range.start && genptr <= pm->gen_range.end &&
      refptr >= pm->ref_range.start && refptr <= pm->ref_range.end) {
    genidx = genptr - pm->gen_range.start;
    refidx = refptr - pm->ref_range.start;
    if (e_state)
      pm->entries[genidx][refidx].on_max_path_e = true;
    else
      pm->entries[genidx][refidx].on_max_path_i = true;
  }
}

void gth_path_matrix_delete(GthPathMatrix *pm)
{
  if (!pm) return;
  gt_array2dim_delete(pm->entries);
  gt_free(pm);
}
