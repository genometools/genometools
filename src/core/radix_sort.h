/*
  Copyright (c) 2011-2012 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2011-2012 Center for Bioinformatics, University of Hamburg

  Permission to use, copy, modify, and distribute this software for any
  purpose with or without fee is hereby granted, provided that the above
  copyright notice and this permission notice appear in all copies.

  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WA(RR)ANTIES
  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WA(RR)ANTIES OF
  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/

#ifndef RADIX_SORT_H
#define RADIX_SORT_H

#include <stdbool.h>
#include "core/types_api.h"

/* This file describes the interface of a radixsort implementation
   which allows for keys of type <unsigned long>, optionally associated
   with values of type <unsigned long>. Thus the implementation can
   sorted arrays referenced by pointers of type <unsigned long *>,
   or arrays referenced by pointers of type <GtUlongPair *>.
   In the latter case, the component <a> is the key and <b> is the
   value.

   The implementation of radixsort allows to split the array
   to be sorted into a given number of at least two parts, which
   are sorted independently. Let <n> be the length of the array and
   <p> be the number of parts.
   - When running radix sort with <p> threads, each part is sorted by
     its own thread, requiring an array of size <n> as auxiliary space.
     In the ideal case, the speedup for sorting all parts is <1/p>.
   - When running radixsort with one thread, the parts are sorted
     one after the other, requiring an array of size <n/p> as
     auxiliary space.

   The time and space improvements cannot be combined.

   The disadvantage of the splitting method is the fact, that,
   in order to obtain a completely sorted array, the different
   parts have to be merged. As the methods to merge the
   arrays in-place require extra space and thus would sacrifice the
   space advantage, we have not added appropriate code and suppose that
   the results of the final sorting are immediately processed anyway.
   Thus we only need a mechanism to simultaneously read the <p> sorted parts.
   This mechanism is implemented by the class GtRadixreader which, for
   efficiency reasons, is not opaque and mainly accessed via macros.
*/

/* Sascha: move radixreader code to own header file */

/*
   The following macro reads the next values from the different parts
   represented by <RR> and stores it in <VALUE>. <RR> refers to a
   set of independently sorted arrays of type <unsigned long>.
   If there are no more values left in the different arrays,
   then <STOPSTATEMENT> is executed.
*/

/*
   The following macro reads the next values from the different parts
   represented by <RR> and stores it in <VALUE>. <RR> refers to a
   set of independently sorted arrays of type <GtUlongPair>.
   If there are no more values left in the different arrays,
   then <STOPSTATEMENT> is executed.
*/

/* If the numbers of parts is larger than 2, then for all parts, we
   store the smallest non-processed element in a priority
   queue according to the given <sortkey>. For the case of sorting
   pairs, we have an additional value <additionalvalue>.
   <part> stores the ordinal number of the part, from which the values
   stem.
*/

/* If <currentptr> is smaller than <endptr>, the <currentptr> refers
   to the next element to be processed. The analogue holds for
   <currentptr_pair> and <endptr_pair> in case we sort pairs.
*/

/* the following function stores an element with from the array number <part>
   with its <sortkey>  and <additionalvalue> value. */

/* The following type is opaque and implemented in <radix_sort.c>. */

/* The following function creates an object of class <GtRadixsort> and return
   a pointer to it. The object can be used to sort arrays of
   <unsigned long>-integers.
   If <smalltables> is <true>, then the radix-keys consists of units of 1 byte.
   If <smalltables> is <false>, then the radix-keys consists of units of
   2 bytes. We recommend to set <smalltables> to <true>. <maxlen> is the
   maximum size of the array to be sorted. <rparts> is the number of parts
   in which the arrays are sorted. If <withthreads> is <true>, then
   <rparts> threads are used to sort <rparts> subarrays of similar size.
   If <withthreads> is <false>, then one thread is used to
   sort <rparts> subarrays of similar size one after the other.
   The memory area, in which the elements to be sorted are found can
   be supplied by the argument <arr>. If <arr> is <NULL>, then the
   an an array of the appropiate size is created.
   */

/* The following function is like the previous, except that the
   created object can be used to sort arrays of <GtUlongPair>-values. */

/* This is the function to perform the sorting task for the first
   <len> elements of the array stored in the object <radixsort>.
   If at least two parts have been used in the sorting, then
   a <GtRadixreader> object is return. Otherwise, <NULL> is returned.
*/

/* Delete a <GtRadixsortinfo>-object. */

/* Determine of maximum number of entries in an array such that
   the given memory limit <memlimit> (in bytes) for the array itself and
   the auxiliary array is not exceeded.
   <rparts> and <withthreads> have the same meaning as the corresponding
   parameters in <gt_radixsort_new_ulong>.
*/

/* The analogue function as before, but for the case that
   arrays over type <GtUlongPair> are to be sorted.
*/

/* Return the size of the <GtRadixsortinfo>-object. */

/* Return a pointer to the memory area in which the elements to
   be sorted can be stored.
*/

/* Same as before, but for the case that pairs are to be sorted.
*/

/* Verify the corrected of the sorting for the given <GtRadixreader>-object.
   Only works for the case that integers are sorted.
*/

void gt_radixsort_lsb_linear(unsigned long *source,unsigned long len);

/* The following function implements a radixsort which does not require
   extra workspace, i.e. it is inplace. The main idea is adapoted from
   http://drdobbs.com/architecture-and-design/221600153
   Instead of a recursive approach we use an iterative approach.
*/

void gt_radixsort_inplace_GtUlong(unsigned long *source, unsigned long len);

typedef struct GtRadixsortIPinfo GtRadixsortIPinfo;

GtRadixsortIPinfo *gt_radixsortinfo2_new(bool pairs,unsigned long maxlen);

size_t gt_radixsortinfo2_size(const GtRadixsortIPinfo *radixsortinfo);

unsigned long gt_radixsortinfo2_max_num_of_entries_ulong(size_t memlimit);

unsigned long gt_radixsortinfo2_max_num_of_entries_ulongpair(size_t memlimit);

void gt_radixsort_inplace_sort(GtRadixsortIPinfo *radixsortinfo,
                               unsigned long len);

void gt_radixsortinfo2_delete(GtRadixsortIPinfo *radixsortinfo);

unsigned long *gt_radixsortinfo2_space_ulong(GtRadixsortIPinfo *radixsortinfo);

GtUlongPair *gt_radixsortinfo2_space_ulongpair(
                                            GtRadixsortIPinfo *radixsortinfo);

#endif
