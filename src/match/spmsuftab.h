/*
  Copyright (c) 2011 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
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

#ifndef SPMSUFTAB_H
#define SPMSUFTAB_H

#include "core/unused_api.h"
#include "core/logger_api.h"
#include "core/compact_ulong_store.h"

typedef struct
{
  unsigned long partoffset, numofentries, maxvalue;
  bool usebitsforpositions;
  GtCompactUlongStore *bitpackarray;
} GtSpmsuftab;

GT_UNUSED static inline void gt_spmsuftab_set(GtSpmsuftab *spmsuftab,
                                              unsigned long idx,
                                              unsigned long value)
{
  gt_assert(idx >= spmsuftab->partoffset);
  idx -= spmsuftab->partoffset;
  gt_assert(idx < spmsuftab->numofentries && value <= spmsuftab->maxvalue);
  gt_assert(spmsuftab->bitpackarray != NULL);
  gt_compact_ulong_store_update(spmsuftab->bitpackarray,idx,value);
}

GT_UNUSED static inline unsigned long gt_spmsuftab_get(
                                      const GtSpmsuftab *spmsuftab,
                                      unsigned long idx)
{
  gt_assert(idx >= spmsuftab->partoffset);
  idx -= spmsuftab->partoffset;
  gt_assert(idx < spmsuftab->numofentries);
  return gt_compact_ulong_store_get(spmsuftab->bitpackarray,idx);
}

GtSpmsuftab *gt_spmsuftab_new(unsigned long numofentries,
                              unsigned long maxvalue,
                              unsigned int bitsforseqnumrelpos,
                              GtLogger *logger);

bool gt_spmsuftab_usebitsforpositions(const GtSpmsuftab *spmsuftab);

void gt_spmsuftab_delete(GtSpmsuftab *spmsuftab);

size_t gt_spmsuftab_requiredspace(unsigned long numofentries,
                                  unsigned long maxvalue,
                                  unsigned int bitsforseqnumrelpos);

void gt_spmsuftab_partoffset(GtSpmsuftab *spmsuftab,
                             unsigned long offset);

#endif
