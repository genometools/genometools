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

#ifndef RADIX_INTSORT_H
#define RADIX_INTSORT_H

#include "core/types_api.h"
#include "extended/priorityqueue.h"

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
          if (gt_priorityqueue_is_empty((RR)->priorityqueue))\
          {\
            STOPSTATEMENT;\
          } else\
          {\
            GtPQelementtype *minelem \
              = gt_priorityqueue_delete_min((RR)->priorityqueue);\
            VALUE = minelem->sortkey;\
            if ((RR)->ptrtab[minelem->value].currentptr <\
                (RR)->ptrtab[minelem->value].endptr)\
            {\
              gt_priorityqueue_add((RR)->priorityqueue,\
                                   *(RR)->ptrtab[minelem->value].currentptr,\
                                   minelem->value);\
              (RR)->ptrtab[minelem->value].currentptr++;\
            }\
          }\
        }

#define GT_RADIXREADER_NEXT_PAIR(VALUE,RR,STOPSTATEMENT)\
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
        }

typedef struct
{
  GtUlong *currentptr, *endptr;
} GtRadixreaderPointerpair;

typedef struct
{
  GtUlong *ptr1, *ptr2, *end1, *end2;
  GtUlongPair *ptr1_pair, *ptr2_pair, *end1_pair, *end2_pair;
  GtRadixreaderPointerpair *ptrtab;
  GtPriorityQueue *priorityqueue;
} GtRadixreader;

typedef struct GtRadixsortinfo GtRadixsortinfo;

GtRadixsortinfo *gt_radixsort_new(bool pair,
                                  bool smalltables,
                                  unsigned long maxlen,
                                  unsigned int parts,
                                  void *arr);

unsigned long gt_radixsort_entries(bool pair,unsigned int parts,
                                   size_t memlimit);

void gt_radixsort_verify(GtRadixreader *rr);

GtUlong *gt_radixsort_arr(GtRadixsortinfo *radixsort);

GtUlongPair *gt_radixsort_arrpair(GtRadixsortinfo *radixsort);

size_t gt_radixsort_size(const GtRadixsortinfo *radixsort);

void gt_radixsort_delete(GtRadixsortinfo *radixsort);

void gt_radixsort_linear(GtRadixsortinfo *radixsort,unsigned long len);

GtRadixreader *gt_radixsort_linear_rr(GtRadixsortinfo *radixsort,
                                      unsigned long len);

void gt_radixsort_GtUlong_divide(GtUlong *source,
                                 GtUlong *dest,
                                 unsigned long len);

void gt_radixsort_GtUlong_recursive(GtUlong *source,
                                    GtUlong *dest,
                                    unsigned long len);

#endif
