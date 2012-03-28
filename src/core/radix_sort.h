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
#define GT_RADIXREADER_NEXT(VALUE,RR,STOPSTATEMENT)\
        if ((RR)->ptrtab == NULL)\
        {\
          if ((RR)->ptr1 < (RR)->end1)\
          {\
            if ((RR)->ptr2 < (RR)->end2)\
            {\
              if (*(RR)->ptr1 <= *(RR)->ptr2)\
              {\
                VALUE = *(RR)->ptr1++;\
              } else\
              {\
                VALUE = *(RR)->ptr2++;\
              }\
            } else\
            {\
              VALUE = *(RR)->ptr1++;\
            }\
          } else\
          {\
            if ((RR)->ptr2 < (RR)->end2)\
            {\
              VALUE = *(RR)->ptr2++;\
            } else\
            {\
              STOPSTATEMENT;\
            }\
          }\
        } else\
        {\
          if ((RR)->pq_numofelements > 0)\
          {\
            GtRadixreaderPQelemtype *minelem = (RR)->pq_values +\
                                               (--(RR)->pq_numofelements);\
            GtRadixreaderPointerpair *ptrtabptr = (RR)->ptrtab + minelem->part;\
            VALUE = minelem->sortkey;\
            if (ptrtabptr->currentptr < ptrtabptr->endptr)\
            {\
              gt_radixreaderPQadd((RR),\
                                  *(ptrtabptr->currentptr++),\
                                  minelem->part,0);\
            }\
          } else\
          {\
            STOPSTATEMENT;\
          }\
        }

/*
   The following macro reads the next values from the different parts
   represented by <RR> and stores it in <VALUE>. <RR> refers to a
   set of independently sorted arrays of type <GtUlongPair>.
   If there are no more values left in the different arrays,
   then <STOPSTATEMENT> is executed.
*/
#define GT_RADIXREADER_NEXT_PAIR(VALUE,RR,STOPSTATEMENT)\
        if ((RR)->ptrtab == NULL)\
        {\
          if ((RR)->ptr1_pair < (RR)->end1_pair)\
          {\
            if ((RR)->ptr2_pair < (RR)->end2_pair)\
            {\
              if ((RR)->ptr1_pair->a <= (RR)->ptr2_pair->a)\
              {\
                VALUE = *(RR)->ptr1_pair++;\
              } else\
              {\
                VALUE = *(RR)->ptr2_pair++;\
              }\
            } else\
            {\
              VALUE = *(RR)->ptr1_pair++;\
            }\
          } else\
          {\
            if ((RR)->ptr2_pair < (RR)->end2_pair)\
            {\
              VALUE = *(RR)->ptr2_pair++;\
            } else\
            {\
              STOPSTATEMENT;\
            }\
          }\
        } else\
        {\
          if ((RR)->pq_numofelements > 0)\
          {\
            GtRadixreaderPQelemtype *minelem = (RR)->pq_values +\
                                               (--(RR)->pq_numofelements);\
            GtRadixreaderPointerpair *ptrtabptr = (RR)->ptrtab + minelem->part;\
            VALUE.a = minelem->sortkey;\
            VALUE.b = minelem->additionalvalue;\
            if (ptrtabptr->currentptr_pair < ptrtabptr->endptr_pair)\
            {\
              gt_radixreaderPQadd(RR,\
                                  ptrtabptr->currentptr_pair->a,\
                                  minelem->part,\
                                  ptrtabptr->currentptr_pair->b);\
              ptrtabptr->currentptr_pair++;\
            }\
          } else\
          {\
            STOPSTATEMENT;\
          }\
        }

/* If the numbers of parts is larger than 2, then for all parts, we
   store the smallest non-processed element in a priority
   queue according to the given <sortkey>. For the case of sorting
   pairs, we have an additional value <additionalvalue>.
   <part> stores the ordinal number of the part, from which the values
   stem.
*/
typedef struct
{
  unsigned long sortkey,
                additionalvalue;
  unsigned int part;
} GtRadixreaderPQelemtype;

/* If <currentptr> is smaller than <endptr>, the <currentptr> refers
   to the next element to be processed. The analogue holds for
   <currentptr_pair> and <endptr_pair> in case we sort pairs.
*/
typedef struct
{
  unsigned long *currentptr, *endptr;
  GtUlongPair *currentptr_pair, *endptr_pair;
} GtRadixreaderPointerpair;

typedef struct
{
  unsigned long *ptr1, *ptr2, *end1, *end2; /* for 2 parts */
  GtUlongPair *ptr1_pair, *ptr2_pair, *end1_pair, *end2_pair; /* for 2 parts */
  /* the remaining components are for more than 2 parts */
  GtRadixreaderPointerpair *ptrtab; /* for each part store <currentptr>
                                       and <endptr> or
                                       <currentptr_pair> and <endptr_pair>. */
  unsigned long pq_numofelements; /* for > 2 parts; number of elements
                                     in priority queue */
  GtRadixreaderPQelemtype *pq_values; /* priority queue */
} GtRadixreader;

/* the following function stores an element with from the array number <part>
   with its <sortkey>  and <additionalvalue> value. */

/*@unused@*/ static inline void gt_radixreaderPQadd(GtRadixreader *rr,
                                                    unsigned long sortkey,
                                                    unsigned int part,
                                                    unsigned long
                                                      additionalvalue)
{
  GtRadixreaderPQelemtype *ptr;

    /* store elements in reverse order, i.e.\ with the minimum element
       at the last index; move elements to the right until an element
       larger or equal than the key is found. */
  for (ptr = rr->pq_values + rr->pq_numofelements; ptr > rr->pq_values; ptr--)
  {
    if ((ptr-1)->sortkey < sortkey)
    {
      *ptr = *(ptr-1);
    } else
    {
      break;
    }
  }
  ptr->sortkey = sortkey;
  ptr->part = part;
  ptr->additionalvalue = additionalvalue;
  rr->pq_numofelements++;
}

/* The following type is opaque and implemented in <radix_sort.c>. */

typedef struct GtRadixsortinfo GtRadixsortinfo;

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

GtRadixsortinfo *gt_radixsort_new_ulong(bool smalltables,
                                        unsigned long maxlen,
                                        unsigned int rparts,
                                        bool withthreads,
                                        unsigned long *arr);

/* The following function is like the previous, except that the
   created object can be used to sort arrays of <GtUlongPair>-values. */

GtRadixsortinfo *gt_radixsort_new_ulongpair(bool smalltables,
                                            unsigned long maxlen,
                                            unsigned int rparts,
                                            bool withthreads,
                                            GtUlongPair *arr);

/* This is the function to perform the sorting task for the first
   <len> elements of the array stored in the object <radixsort>.
   If at least two parts have been used in the sorting, then
   a <GtRadixreader> object is return. Otherwise, <NULL> is returned.
*/
GtRadixreader *gt_radixsort_sort(GtRadixsortinfo *radixsort,
                                 unsigned long len);

/* Delete a <GtRadixsortinfo>-object. */

void gt_radixsort_delete(GtRadixsortinfo *radixsort);

/* Determine of maximum number of entries in an array such that
   the given memory limit <memlimit> (in bytes) for the array itself and
   the auxiliary array is not exceeded.
   <rparts> and <withthreads> have the same meaning as the corresponding
   parameters in <gt_radixsort_new_ulong>.
*/
unsigned long gt_radixsort_max_num_of_entries_ulong(unsigned int rparts,
                                                   bool withthreads,
                                                   size_t memlimit);

/* The analogue function as before, but for the case that
   arrays over type <GtUlongPair> are to be sorted.
*/
unsigned long gt_radixsort_max_num_of_entries_ulongpair(unsigned int rparts,
                                                        bool withthreads,
                                                        size_t memlimit);

/* Return the size of the <GtRadixsortinfo>-object. */

size_t gt_radixsort_size(const GtRadixsortinfo *radixsort);

/* Return a pointer to the memory area in which the elements to
   be sorted can be stored.
*/
unsigned long *gt_radixsort_space_ulong(GtRadixsortinfo *radixsort);

/* Same as before, but for the case that pairs are to be sorted.
*/
GtUlongPair *gt_radixsort_space_ulongpair(GtRadixsortinfo *radixsort);

/* Verify the corrected of the sorting for the given <GtRadixreader>-object.
   Only works for the case that integers are sorted.
*/
void gt_radixsort_verify(GtRadixreader *rr);

/* The following function implements a radixsort which does not require
   extra workspace, i.e. it is inplace. The main idea is adapoted from
   http://drdobbs.com/architecture-and-design/221600153
   Instead of a recursive approach we use an iterative approach.
*/

void gt_radixsort_inplace_GtUlong(unsigned long *source, unsigned long len);

#endif
