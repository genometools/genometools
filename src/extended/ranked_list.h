/*
  Copyright (c) 2012 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2012 Center for Bioinformatics, University of Hamburg

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

#ifndef RANKED_LIST_H
#define RANKED_LIST_H

#include "extended/rbtree.h"

typedef GtRBTreeCompareFunc GtRankedlistCompareFunc;
typedef GtRBTreeIter GtRankedListIter;
typedef struct GtRankedlist GtRankedlist;

GtRankedlist *gt_ranked_list_new(unsigned long maxsize,
                                 GtRBTreeCompareFunc comparefunction,
                                 void *cmpinfo);

void gt_ranked_list_insert(GtRankedlist *ranked_list,void *elemin);

void *gt_ranked_list_minimum_key(const GtRankedlist *ranked_list);

void *gt_ranked_list_maximum_key(const GtRankedlist *ranked_list);

unsigned long gt_ranked_list_currentsize(const GtRankedlist *ranked_list);

GtRankedListIter *gt_ranked_list_iter_new_from_first(GtRankedlist *ranked_list);

GtRankedListIter *gt_ranked_list_iter_new_from_last(GtRankedlist *ranked_list);

void *gt_ranked_list_iter_next(GtRankedListIter *trav);

void *gt_ranked_list_iter_prev(GtRankedListIter *trav);

void gt_ranked_list_delete(GtRankedlist *ranked_list);

#endif
