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

#include <stdlib.h>
#include <stdbool.h>
#include <inttypes.h>

#include "core/types_api.h"

/* Type and size of elements used to represent the intset. */
typedef enum {
  GT_INTSET_8  = sizeof (uint8_t),
  GT_INTSET_16 = sizeof (uint16_t),
  GT_INTSET_32 = sizeof (uint32_t)
} GtIntsetType;

/* Return the optimal type of intset for your set. */
GtIntsetType gt_intset_best_type(GtUword maxelement, GtUword num_of_elems);

/* The <GtIntset8> class. Used to store a fixed number of sorted 64bit integers
   with reduced space using 8bit elements. */
typedef struct GtIntset8 GtIntset8;
/* The <GtIntset16> class. Used to store a fixed number of sorted 64bit integers
   with reduced space using 16bit elements. */
typedef struct GtIntset16 GtIntset16;
/* The <GtIntset32> class. Used to store a fixed number of sorted 64bit integers
   with reduced space using 32bit elements. */
typedef struct GtIntset32 GtIntset32;

/* Return a new <GtIntset8> object with space for <num_of_elems> <GtUword>
   integers, where the largest number is <max_elem>. */
static GtIntset8  *gt_intset_8_new(GtUword max_elem, GtUword num_of_elems);
/* Return a new <GtIntset16> object with space for <num_of_elems> <GtUword>
   integers, where the largest number is <max_elem>. */
static GtIntset16 *gt_intset_16_new(GtUword max_elem, GtUword num_of_elems);
/* Return a new <GtIntset32> object with space for <num_of_elems> <GtUword>
   integers, where the largest number is <max_elem>. */
static GtIntset32 *gt_intset_32_new(GtUword max_elem, GtUword num_of_elems);

/* Free the memory of <intset>. */
static void        gt_intset_8_delete(GtIntset8 *intset);
/* Free the memory of <intset>. */
static void        gt_intset_16_delete(GtIntset16 *intset);
/* Free the memory of <intset>. */
static void        gt_intset_32_delete(GtIntset32 *intset);

/* Add <elem> to <intset>. <elem> has to be larger than the previous <elem>
   added. */
static void        gt_intset_8_add(GtIntset8 *intset, GtUword elem);
/* Add <elem> to <intset>. <elem> has to be larger than the previous <elem>
   added. */
static void        gt_intset_16_add(GtIntset16 *intset, GtUword elem);
/* Add <elem> to <intset>. <elem> has to be larger than the previous <elem>
   added. */
static void        gt_intset_32_add(GtIntset32 *intset, GtUword elem);

/* Returns <true> if <elem> is a member of the set <intset>. */
static bool        gt_intset_8_is_member(const GtIntset8 *intset,
                                         GtUword elem);
/* Returns <true> if <elem> is a member of the set <intset>. */
static bool        gt_intset_16_is_member(const GtIntset16 *intset,
                                          GtUword elem);
/* Returns <true> if <elem> is a member of the set <intset>. */
static bool        gt_intset_32_is_member(const GtIntset32 *intset,
                                          GtUword elem);

/* Returns the number of the element in <intset> that is the smallest element
   larger than <pos>.
   This is used for sets representing the separator positions in a set of
   sequences, to determine the sequence number corresponding to any position in
   the concatenated string of the sequence set. */
static GtUword     gt_intset_8_pos2seqnum(const GtIntset8 *intset,
                                          GtUword pos);
/* Returns the number of the element in <intset> that is the smallest element
   larger than <pos>.
   This is used for sets representing the separator positions in a set of
   sequences, to determine the sequence number corresponding to any position in
   the concatenated string of the sequence set. */
static GtUword     gt_intset_16_pos2seqnum(const GtIntset16 *intset,
                                           GtUword pos);
/* Returns the number of the element in <intset> that is the smallest element
   larger than <pos>.
   This is used for sets representing the separator positions in a set of
   sequences, to determine the sequence number corresponding to any position in
   the concatenated string of the sequence set. */
static GtUword     gt_intset_32_pos2seqnum(const GtIntset32 *intset,
                                           GtUword pos);

/* Returns the size of an intset with given number of elements
   <num_of_elems> and maximum value <maxelement>.
 */
static size_t      gt_intset_8_size(GtUword maxelement, GtUword num_of_elems);
/* Returns the size of an intset with given number of elements
   <num_of_elems> and maximum value <maxelement>.
 */
static size_t      gt_intset_16_size(GtUword maxelement, GtUword num_of_elems);
/* Returns the size of an intset with given number of elements
   <num_of_elems> and maximum value <maxelement>.
 */
static size_t      gt_intset_32_size(GtUword maxelement, GtUword num_of_elems);

#include "extended/intset_impl.h"

int gt_intset_unit_test(GtError *err);

#endif
