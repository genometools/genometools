/*
  Copyright (c) 2011 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2012 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
  Copyright (c) 2011-2012 Center for Bioinformatics, University of Hamburg

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

#ifndef COMPACT_ULONG_STORE_H
#define COMPACT_ULONG_STORE_H

#include "core/error_api.h"

/* Class <GtCompactUlongStore> stores fixed bit witdh unsigned integers. Maximum
   bitwidth is sizeof (unsigned long) */
typedef struct GtCompactUlongStore GtCompactUlongStore;

/* Return a new <GtCompactUlongStore> object with <numofentries> elements of
   <bitsperentry> bit width */
GtCompactUlongStore* gt_compact_ulong_store_new(unsigned long numofentries,
                                                unsigned int bitsperentry);

/* Deletes <cus> object and frees all associated memory */
void                 gt_compact_ulong_store_delete(GtCompactUlongStore *cus);

/* Calculates the size in bytes for a <GtCompactUlongStore> with <numofentries>
   elements of <bitsperentry> bit width */
size_t               gt_compact_ulong_store_size(unsigned long numofentries,
                                                 unsigned int bitsperentry);

/* Return element stored in <cus> at position <idx> cast to unsigned long */
unsigned long        gt_compact_ulong_store_get(const GtCompactUlongStore *cus,
                                                unsigned long idx);

/* Set element at position <idx> in <cus> to <value> */
void                 gt_compact_ulong_store_update(GtCompactUlongStore *cus,
                                                   unsigned long idx,
                                                   unsigned long value);

int                  gt_compact_ulong_store_unit_test(GtError *err);

#endif
