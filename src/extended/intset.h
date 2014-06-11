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
#include "extended/io_function_pointers.h"

/* The <GtIntset> interface used to store a fixed set of sorted 64bit integers
   with reduced space. (Mathematical set, no duplicates allowed) */
typedef struct GtIntset GtIntset;

/* Return <GtIntset> of optimal size, choosing one of the available
   implementations of this class.
   <maxelement> has to be larger or equal <num_of_elems>. */
GtIntset* gt_intset_best_new(GtUword maxelement, GtUword num_of_elems);

/* Allocates new memory for <GtIntset> and fills it with data read from file
   <fp>. Exits on error reading the file. */
GtIntset* gt_intset_new_from_file(FILE *fp, GtError *err);

/* Increases the reference count of the <GtIntset>. */
GtIntset* gt_intset_ref(GtIntset *intset);

/* Add <elem> to <intset>, <elem> has to be larger than the previous added
   element. */
void      gt_intset_add(GtIntset *intset, GtUword elem);

/* Returns the element at index <idx> in the sorted set <intset>. */
GtUword   gt_intset_get(GtIntset *intset, GtUword idx);

/* Returns actual number of stored elements */
GtUword   gt_intset_size(GtIntset *intset);

/* Returns true if <elem> is a member of set <intset>. */
bool      gt_intset_is_member(GtIntset *intset, GtUword elem);

/* Returns the number of the element in <intset> that is the smallest element
   larger than or equal <pos> or <num_of_elems> if there is no such <element>.
   This can be used for sets representing the separator positions in a set of
   sequences, to determine the sequence number corresponding to any position in
   the concatenated string of the sequence set.
   Fails for <pos> > <maxelement>! */
GtUword   gt_intset_get_idx_smaller_geq(GtIntset *intset, GtUword pos);

/* Returns the size in bytes of the <GtIntset>-structure. */
size_t    gt_intset_size_of_struct(GtIntset *intset);

/* Returns the size of the representation of an <intset> with its given number
   of elements <num_of_elems> and maximum value <maxelement>, in bytes. This
   does not include the size of the structure. */
size_t    gt_intset_size_of_rep(GtIntset *intset);

/* Write <intset> to file <fp>. Fails with exit on IO-error. Returns NULL if
   data error occures and writes it to <err>, <intset> will be deleted at that
   point. */
GtIntset* gt_intset_write(GtIntset *intset, FILE *fp, GtError *err);

/* Read or write to/from File, depending on <intset>. If <NULL>, it allocates
   memory for a new <GtIntset> object and tries to fill it from file <fp>.
   If not <NULL> it writs the content of <intset> to <fp>.
   Returns <NULL> on error, in which case <intset> will be deleted and <err>
   will be set. */
GtIntset* gt_intset_io(GtIntset *intset, FILE *fp, GtError *err);

/* Free the memory of <intset>. */
void      gt_intset_delete(GtIntset *intset);

/* Runs unit tests for all implementations of the <GtIntset> class. */
int       gt_intset_unit_test(GtError *err);
#endif
