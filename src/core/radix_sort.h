/*
  Copyright (c) 2011 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2011 Center for Bioinformatics, University of Hamburg

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

/* move radixreader code to own header file */

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
            VALUE.b = minelem->suffixref;\
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

typedef struct
{
  unsigned long sortkey,
                suffixref;
  unsigned int part;
} GtRadixreaderPQelemtype;

typedef struct
{
  unsigned long *currentptr, *endptr;
  GtUlongPair *currentptr_pair, *endptr_pair;
} GtRadixreaderPointerpair;

typedef struct
{
  unsigned long *ptr1, *ptr2, *end1, *end2;
  GtUlongPair *ptr1_pair, *ptr2_pair, *end1_pair, *end2_pair;
  GtRadixreaderPointerpair *ptrtab;
  unsigned long pq_numofelements;
  GtRadixreaderPQelemtype *pq_values;
} GtRadixreader;

/*@unused@*/ static inline void gt_radixreaderPQadd(GtRadixreader *rr,
                                                    unsigned long sortkey,
                                                    unsigned int part,
                                                    unsigned long suffixref)
{
  GtRadixreaderPQelemtype *ptr;

    /* store elements in reverse order, i.e.\ with the minimum element
       at the last index */
    /* move elements to the right until an element larger or equal than
       the key is found. */
  for (ptr = rr->pq_values + rr->pq_numofelements; ptr > rr->pq_values;
       ptr--)
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
  ptr->suffixref = suffixref;
  rr->pq_numofelements++;
}

typedef struct GtRadixsortinfo GtRadixsortinfo;

GtRadixsortinfo *gt_radixsort_new_ulong(bool smalltables,
                                        unsigned long maxlen,
                                        unsigned int rparts,
                                        bool withthreads,
                                        unsigned long *arr);

GtRadixsortinfo *gt_radixsort_new_ulongpair(bool smalltables,
                                            unsigned long maxlen,
                                            unsigned int rparts,
                                            bool withthreads,
                                            GtUlongPair *arr);

GtRadixreader *gt_radixsort_sort(GtRadixsortinfo *radixsort,
                                 unsigned long len);

void gt_radixsort_delete(GtRadixsortinfo *radixsort);

unsigned long gt_radixsort_estimate_num_of_entries(bool pair,
                                                   unsigned int rparts,
                                                   size_t memlimit,
                                                   bool withthreads);

size_t gt_radixsort_size(const GtRadixsortinfo *radixsort);

unsigned long *gt_radixsort_space_ulong(GtRadixsortinfo *radixsort);

GtUlongPair *gt_radixsort_space_ulongpair(GtRadixsortinfo *radixsort);

void gt_radixsort_verify(GtRadixreader *rr);

void gt_radixsort_divide(unsigned long *source,
                         unsigned long *dest,
                         unsigned long len);

void gt_radixsort_recursive(unsigned long *source,
                            unsigned long *dest,
                            unsigned long len);

#endif
