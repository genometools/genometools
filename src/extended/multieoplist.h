/*
  Copyright (c) 2013 Ole Eigenbrod <ole.eigenbrod@gmx.de>
  Copyright (c) 2013 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
  Copyright (c) 2013 Center for Bioinformatics, University of Hamburg

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

#ifndef MULTIEOPLIST_H
#define MULTIEOPLIST_H

#include "core/array.h"
#include "core/encseq_api.h"
#include "core/error_api.h"
#include "core/types_api.h"
#include "extended/io_function_pointers.h"

/* Class <GtMultieoplist> stores a list of edit operations that define an
   alignment between two sequences. */
typedef struct GtMultieoplist GtMultieoplist;

/* <AlignmentEoptype> represents the types of operations allowed in
   <GtMultieoplist>. <Replacement> will be stored as a <Match>. Its use is
   deprecated and should be replaced by <Match> or <Mismatch>. */
typedef enum {
  Match,
  Mismatch,
  Replacement,
  Deletion,
  Insertion
} AlignmentEoptype;

/* <GtMultieop> one edit operation of type <type> repeated <steps> times in the
   alignment. */
typedef struct {
  AlignmentEoptype type;
  GtUword steps;
} GtMultieop;

/* Returns new empty <GtMultieoplist> object */
GtMultieoplist* gt_multieoplist_new(void);

/* Returns new <GtMultieoplist> object of size size*/
GtMultieoplist* gt_multieoplist_new_with_size(GtUword size);

/* Returns reference to <multieops> */
GtMultieoplist* gt_multieoplist_ref(GtMultieoplist *multieops);
void            gt_multieoplist_delete(GtMultieoplist *multieops);

/* Add one <Replacement> operation to the list of edit operations. Multiple
   additions will be combined by increasing <steps>. */
void            gt_multieoplist_add_replacement(GtMultieoplist *multieops);

/* Add one <Insertion> operation to the list of edit operations. Multiple
   additions will be combined by increasing <steps>. */
void            gt_multieoplist_add_insertion(GtMultieoplist *multieops);

/* Add one <Deletion> operation to the list of edit operations. Multiple
   additions will be combined by increasing <steps>. */
void            gt_multieoplist_add_deletion(GtMultieoplist *multieops);

/* Add one <Mismatch> operation to the list of edit operations. Multiple
   additions will be combined by increasing <steps>. */
void            gt_multieoplist_add_mismatch(GtMultieoplist *multieops);

/* Add one <Match> operation to the list of edit operations. Multiple
   additions will be combined by increasing <steps>. */
void            gt_multieoplist_add_match(GtMultieoplist *multieops);

/* Remove all operations from <multieops> */
void            gt_multieoplist_reset(GtMultieoplist *multieops);

/* Remove the last added edit operation. Will remove exactly one! (decrease
   <steps>) */
void            gt_multieoplist_remove_last(GtMultieoplist *multieops);

/* TODO add function to add a <GtMultieop> */

/* Returns pointer to a copy of contents from <GtMultieoplist> <source> to
   <GtMultieoplist> <copy>. <copy> may be NULL, returns new <GtMultieoplist>
   object in that case. <source> may not be NULL!
 */
GtMultieoplist* gt_multieoplist_clone(GtMultieoplist *copy,
                                      GtMultieoplist *source);

/* Returns the number of <GtMultieop> elements in <multieops>, each of which
   can have <steps> >= 1, therefor this represents not the length of the
   alignment. */
GtUword         gt_multieoplist_get_num_entries(GtMultieoplist *multieops);

/* Returns <GtMultieop> number <index>. */
GtMultieop      gt_multieoplist_get_entry(GtMultieoplist *multieops,
                                          GtUword index);

/* Returns sum of <Replacement>, <Match>, <Mismatch> and <Deletion> including
   their <steps>. This corresponds to the length of the first sequence. */
GtUword         gt_multieoplist_get_repdel_length(GtMultieoplist *multieops);

/* Returns sum of <Replacement>, <Match>, <Mismatch> and <Insertion> including
   their <steps>. This corresponds to the length of the second sequence. */
GtUword         gt_multieoplist_get_repins_length(GtMultieoplist *multieops);

/* Returns sum of <Replacement>, <Match>, <Mismatch>, <Deletion> and
   <Insertion> including their <steps>. This corresponds to the length of the
   whole alignment. */
GtUword         gt_multieoplist_get_length(GtMultieoplist *multieops);

/* Print a string representation of <multieops> to <fp> ending with a newline.
   For example: [M5,R2,M3,I1,M6] R is equivalent for <Mismatch> and
   <Replacement>.
   Assumes reverse ordered <GtMultieoplist>. */
void            gt_multieoplist_show(GtMultieoplist *multieops, FILE *fp);

/* combine two <GtMultieoplist>s. Adds <multieops_to_add> to the end of
   <multieops> (which is usualy the start of the alignment). <forward> defines
   if <multieops_to_add> should be read in forward or reverse direction.
   If <multieops_to_add> and <multieops> are from consequtive alignments:
   [...<-...][...<-...] or [...->...][...->...]
   use true
   if they ar from alignments with the same startpoint:
   [...->...][...<-...]
   use false
   and add the one that is at the end of the other.
   if the layout is like this:
   [...<-...][...->...]
   none of them can be added to the end of the other without reversing one
   beforehand. */
void            gt_multieoplist_combine(GtMultieoplist *multieops,
                             GtMultieoplist *multieops_to_add,
                             bool forward);

/* Read or write to/from File, depending on <multieops>. If <NULL>, it allocates
   memory for a new <GtMultieoplist> object and tries to fill it from file <fp>.
   If not <NULL> it writs the content of <multieops> to <fp>.
   Returns <NULL> on error, in which case <multieops> will be deleted and <err>
   will be set. */
GtMultieoplist* gt_meoplist_io(GtMultieoplist *multieops, FILE *fp,
                               GtError *err);

#endif
