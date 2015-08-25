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
#include <stdint.h>
#include "core/types_api.h"

/* This file describes the interface of a radixsort implementation
   which allows for keys of type <GtUword>, optionally associated
   with values of type <GtUword>. Thus the implementation can
   sort arrays referenced by pointers of type <GtUword *>,
   or arrays referenced by pointers of type <GtUwordPair *>.
   In the latter case, the component <a> is the key and <b> is the
   value.

   The implementation of radixsort allows one to split the array
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
   efficiency reasons, is not opaque and mainly accessed via macros. */

typedef struct
{
  uint64_t uint64_a, uint64_b;
} Gtuint64keyPair;

/* The following type represents the workspace for the radixsort
   implementation. We use it in applications with on the order of hundred
   calls to the sorting function. Each such call is supplied with the
   same working space thus saving many workspace creations and deletions. */
typedef struct GtRadixsortinfo GtRadixsortinfo;

/* The following function creates an object of class <GtRadixsortinfo> and
   returns a pointer to it. The object can be used to sort arrays of
   <GtUword>-integers. <maxlen> is the
   maximum size of the array to be sorted. */
GtRadixsortinfo* gt_radixsort_new_ulong(GtUword maxlen);

/* The following function is like the previous, except that the
   created object can be used to sort arrays of <GtUwordPair> values. */
GtRadixsortinfo* gt_radixsort_new_ulongpair(GtUword maxlen);

/* The following function is like the previous, except that the
   created object can be used to sort arrays of <Gtuint64keypair> values. */

GtRadixsortinfo* gt_radixsort_new_uint64keypair(GtUword maxlen);

/* Return the size of the <GtRadixsortinfo> object. */
size_t           gt_radixsort_size(const GtRadixsortinfo *radixsortinfo);

/* This is the function to perform the sorting task for the first
   <len> elements of the array stored in the object <radixsortinfo> */
void             gt_radixsort_inplace_sort(GtRadixsortinfo *radixsortinfo,
                                           GtUword len);

/* sort an array of values of type <GtUword> */
void             gt_radixsort_inplace_ulong(GtUword *source, GtUword len);

/* Determine the maximum number of entries in an array of <GtUword> such
   that the given memory limit <memlimit> (in bytes) for the array itself and
   the auxiliary array is not exceeded. */
GtUword          gt_radixsort_max_num_of_entries_ulong(size_t memlimit);

/* Return a pointer to the memory area in which the elements to
   be sorted can be stored. */
GtUword*         gt_radixsort_space_ulong(GtRadixsortinfo *radixsortinfo);

/* Special version of function to sort Pairs of GtUwords, where the
   first component is the key */
void             gt_radixsort_inplace_GtUwordPair(GtUwordPair *source,
                                                  GtUword len);

/* Determine the maximum number of entries in an array of <GtUwordPair>s such
   that the given memory limit <memlimit> (in bytes) for the array itself and
   the auxiliary array is not exceeded. */
GtUword          gt_radixsort_max_num_of_entries_ulongpair(size_t memlimit);

/* The analogue function as before, but for the case that
   arrays over type <GtUwordPair> are to be sorted. */
GtUwordPair*     gt_radixsort_space_ulongpair(GtRadixsortinfo *radixsortinfo);

/* For testing purposes we also provide a radixsort which start with
   the least signifcant bits. */
void             gt_radixsort_lsb_linear(GtUword *source,GtUword len);

/* Special version of function to sort Pairs of uint64_t-values, where
   both components serve as keys. */
void            gt_radixsort_inplace_Gtuint64keyPair(Gtuint64keyPair *source,
                                                     GtUword len);

/* Return the space storing the uint64_t key pairs. */
Gtuint64keyPair* gt_radixsort_space_uint64keypair(
                                                GtRadixsortinfo *radixsortinfo);

/* Delete a <GtRadixsortinfo> object. */
void             gt_radixsort_delete(GtRadixsortinfo *radixsortinfo);

#endif
