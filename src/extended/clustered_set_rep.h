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

#ifndef CLUSTERED_SET_REP_H
#define CLUSTERED_SET_REP_H

#include "extended/clustered_set.h"
#include <string.h>
#include "core/array.h"

typedef unsigned long
  (*GtClusteredSetNumberOfClustersFunc)(GtClusteredSet*, GtError*);
typedef int
  (*GtClusteredSetMergeClustersFunc)(GtClusteredSet*, unsigned long, unsigned
   long, GtError*);
typedef void
  (*GtClusteredSetFreeFunc)(GtClusteredSet*, GtError*);
typedef GtClusteredSetIterator*
  (*GtClusteredSetIteratorFunc)(GtClusteredSet*, unsigned long, GtError*);
typedef unsigned long
  (*GtClusteredSetNumberOfElementsFunc)(GtClusteredSet*, GtError*);
typedef unsigned long
  (*GtClusteredSetClusterNumFunc)(GtClusteredSet*, unsigned long, GtError*);

struct GtClusteredSetClass {
  size_t size;
  GtClusteredSetMergeClustersFunc merge_clusters;
  GtClusteredSetNumberOfClustersFunc number_of_clusters;
  GtClusteredSetFreeFunc free;
  GtClusteredSetIteratorFunc iterator;
  GtClusteredSetNumberOfElementsFunc number_of_elements;
  GtClusteredSetClusterNumFunc cluster_num;
};

struct GtClusteredSet {
  const GtClusteredSetClass *c_class;
  GtClusteredSetMembers *pvt;
};

struct GtClusteredSetIterator{
  unsigned long curpos, length, *elems;
};

const GtClusteredSetClass* gt_clustered_set_class_new(
  size_t size,
  GtClusteredSetNumberOfClustersFunc number_of_clusters,
  GtClusteredSetMergeClustersFunc merge_clusters,
  GtClusteredSetFreeFunc free,
  GtClusteredSetIteratorFunc iterator,
  GtClusteredSetNumberOfElementsFunc number_of_elements,
  GtClusteredSetClusterNumFunc cluster_num
);

GtClusteredSet* gt_clustered_set_create(const GtClusteredSetClass*);
void* gt_clustered_set_cast(const GtClusteredSetClass*, GtClusteredSet*);

#endif
