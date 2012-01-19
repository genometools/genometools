/*
  Copyright (c) 2010 Joachim Bonnet <joachim.bonnet@studium.uni-hamburg.de>

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

#ifndef CLUSTERED_SET_H
#define CLUSTERED_SET_H

#include "core/error_api.h"

typedef struct GtClusteredSet GtClusteredSet;
typedef struct GtClusteredSetClass GtClusteredSetClass;
typedef struct GtClusteredSetMembers GtClusteredSetMembers;
typedef struct GtClusteredSetIterator GtClusteredSetIterator;

typedef enum {
  GT_CLUSTERED_SET_ITERATOR_STATUS_OK    =  0,
  GT_CLUSTERED_SET_ITERATOR_STATUS_END   = -1,
  GT_CLUSTERED_SET_ITERATOR_STATUS_ERROR = -2
} GtClusteredSetIteratorStatus;

int gt_clustered_set_unit_test(GtError*);

void gt_clustered_set_delete(GtClusteredSet*, GtError*);

int gt_clustered_set_merge_clusters(GtClusteredSet*,
                                    unsigned long,
                                    unsigned long,
                                    GtError*);

unsigned long gt_clustered_set_num_of_clusters(GtClusteredSet*, GtError*);

unsigned long gt_clustered_set_num_of_elements(GtClusteredSet*, GtError*);

GtClusteredSetIterator* gt_clustered_set_get_iterator(GtClusteredSet*,
                                                      unsigned long,
                                                      GtError*);

GtClusteredSetIteratorStatus
gt_clustered_set_iterator_next(GtClusteredSetIterator*,
                               unsigned long*,
                               GtError*);

void gt_clustered_set_iterator_delete(GtClusteredSetIterator*, GtError*);
#endif
