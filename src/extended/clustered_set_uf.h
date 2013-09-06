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

#ifndef CLUSTERED_SET_UF_H
#define CLUSTERED_SET_UF_H

#include "extended/clustered_set.h"

typedef struct GtClusteredSetUFElem GtClusteredSetUFElem;
typedef struct GtClusteredSetUFClusterInfo GtClusteredSetUFClusterInfo;
typedef struct GtClusteredSetUF GtClusteredSetUF;

const GtClusteredSetClass* gt_clustered_set_union_find_class(void);

GtClusteredSet*  gt_clustered_set_union_find_new(GtUword, GtError*);

int gt_clustered_set_union_find_merge_clusters(GtClusteredSet*,
                                               GtUword c1,
                                               GtUword c2,
                                               GtError*);

GtUword gt_clustered_set_union_find_num_of_clusters(GtClusteredSet*,
                                                             GtError*);

GtUword gt_clustered_set_union_find_num_of_elements(GtClusteredSet*,
                                                          GtError*);

void gt_clustered_set_union_find_delete(GtClusteredSet*, GtError*);

GtClusteredSetIterator*
gt_clustered_set_union_find_iterator_new(GtClusteredSet*,
                                         GtUword,
                                         GtError*);
int gt_clustered_set_union_find_unit_test(GtError *err);

#endif
