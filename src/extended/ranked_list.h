#ifndef RANKED_LIST_H
#define RANKED_LIST_H

#include "extended/rbtree.h"

typedef struct GtRankedlist GtRankedlist;

GtRankedlist *gt_ranked_list_new(unsigned long maxsize,
                                 GtRBTreeCompareFunc comparefunction,
                                 void *cmpinfo);

void gt_ranked_list_insert(GtRankedlist *ranked_list,void *elemin);

void *gt_ranked_list_minimum_key(const GtRankedlist *ranked_list);

void gt_ranked_list_delete(GtRankedlist *ranked_list);

#endif
