/*
  Copyright (c) 2014 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
  Copyright (c) 2014 Center for Bioinformatics, University of Hamburg

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
#ifndef INTSET_H
#define INTSET_H

#include <inttypes.h>
#include <stdbool.h>
#include <stdlib.h>

#include "core/error_api.h"
#include "core/types_api.h"

/* The GtIntset interface used to store a fixed set of sorted 64bit integers
   with reduced space. (Mathematical set, no duplicates allowed) */
typedef struct GtIntset GtIntset;

/* Type and size of elements used to represent the intset. */
typedef enum {
  GT_INTSET_8  = sizeof (uint8_t),
  GT_INTSET_16 = sizeof (uint16_t),
  GT_INTSET_32 = sizeof (uint32_t)
} GtIntsetType;

/* Increases the reference count of the <GtIntset>. */
GtIntset*    gt_intset_ref(GtIntset *intset);

/* Add <elem> to <intset>, <elem> has to be larger than the previous added
   element. */
void         gt_intset_add(GtIntset *intset, GtUword elem);

/* Returns true if <elem> is a member of set <intset>. */
bool         gt_intset_is_member(GtIntset *intset, GtUword elem);

/* Returns the number of the element in <intset> that is the smallest element
   larger than <pos>.
   This is used for sets representing the separator positions in a set of
   sequences, to determine the sequence number corresponding to any position in
   the concatenated string of the sequence set. */
GtUword      gt_intset_pos2seqnum(GtIntset *intset, GtUword pos);

/* Return the optimal type of intset for your set. */
GtIntsetType gt_intset_best_type(GtUword maxelement, GtUword num_of_elems);

/* Free the memory of <intset>. */
void         gt_intset_delete(GtIntset *intset);

/* Function for unit tests within implementations of this class. Fails if
   <gt_intset_is_member()> called with any number between and including <start>
   and <end> returns true.
   */
int gt_intset_unit_test_notinset(GtIntset *intset, GtUword start,
                                 GtUword end, GtError *err);

/* Function for unit tests within implementations of this class. Fails if
   <gt_intset_pos2seqnum()> called with any number between and including <start>
   and <end> returns any number different than <num>. */
int gt_intset_unit_test_check_seqnum(GtIntset *intset, GtUword start,
                                     GtUword end, GtUword num, GtError *err);

#endif
