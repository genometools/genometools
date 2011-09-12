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

#define SPMSUFTABBITPACK
#ifdef SPMSUFTABBITPACK
#include "core/bitpackarray.h"
#endif

typedef struct
{
#ifdef SPMSUFTABBITPACK
  BitPackArray *bitpackarray;
#else
  unsigned long *suftab;
#endif
} GtSpmsuftab;

GT_UNUSED static inline void gt_spmsuftab_set(GtSpmsuftab *spmsuftab,
                                              unsigned long idx,
                                              unsigned long value)
{
#ifdef SPMSUFTABBITPACK
  gt_assert(spmsuftab->bitpackarray != NULL);
#ifdef _LP64
  bitpackarray_store_uint64(spmsuftab->bitpackarray,(BitOffset) idx,
                            (uint64_t) value);
#else
  bitpackarray_store_uint32(spmsuftab->bitpackarray,(BitOffset) idx,
                            (uint32_t) value);
#endif
#else
  spmsuftab->suftab[idx] = value;
#endif
}

GT_UNUSED static inline unsigned long gt_spmsuftab_get(
                                      const GtSpmsuftab *spmsuftab,
                                      unsigned long idx)
{
#ifdef SPMSUFTABBITPACK
  gt_assert(spmsuftab->bitpackarray != NULL);
#ifdef _LP64
  return bitpackarray_get_uint64(spmsuftab->bitpackarray,(BitOffset) idx);
#else
  return (unsigned long)
         bitpackarray_get_uint32(spmsuftab->bitpackarray,(BitOffset) idx);
#endif
#else
  gt_assert(spmsuftab->suftab != NULL);
  return spmsuftab->suftab[idx];
#endif
}

GtSpmsuftab *gt_spmsuftab_new(unsigned long numofentries,
                              unsigned long maxvalue);

void gt_spmsuftab_delete(GtSpmsuftab *spmsuftab);

size_t gt_spmsuftab_requiredspace(unsigned long numofentries,
                                  unsigned long maxvalue);

#endif
