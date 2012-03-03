/*
  Copyright (c) 2011 Giorgio Gonnella <gonnella@zbh.uni-hamburg.de>
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

#ifndef UINT64HASHTABLE_H
#define UINT64HASHTABLE_H

#include <inttypes.h>
#include "core/timer_api.h"
#include "core/error.h"

/* An hash table of uint64_t values without any associated information */
typedef struct GtUint64hashtable GtUint64hashtable;

/* Create a new <GtUint64hashtable> with space for at least
   <nof_elements> values */
GtUint64hashtable* gt_uint64hashtable_new(size_t nof_elements);

/* Deletes a <GtUint64hashtable> and frees all associated memory */
void gt_uint64hashtable_delete(GtUint64hashtable *table);

/* Searches key in table; returns true if found, false otherwise;
 * if insert_if_not_found is true, key is added to the table if not already
 * contained */
bool gt_uint64hashtable_search(GtUint64hashtable *table, uint64_t key,
                               bool insert_if_not_found);

unsigned long gt_uint64hashtable_countsum_get(const GtUint64hashtable *table);

unsigned long gt_uint64hashtable_partialsums(GtUint64hashtable *table,
                                             GtTimer *timer);

unsigned long gt_uint64hashtable_insertionindex(GtUint64hashtable *table,
                                                uint64_t key);

int gt_uint64hashtable_unit_test(GtError *err);

#endif
