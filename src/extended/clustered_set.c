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

#include "clustered_set_rep.h"
#include "clustered_set_uf.h"
#include "core/ma.h"
#include "core/unused_api.h"
#include "core/class_alloc.h"

const GtClusteredSetClass* gt_clustered_set_class_new(
  size_t size,
  GtClusteredSetNumberOfClustersFunc number_of_clusters,
  GtClusteredSetMergeClustersFunc merge_clusters,
  GtClusteredSetFreeFunc free,
  GtClusteredSetIteratorFunc iterator,
  GtClusteredSetNumberOfElementsFunc number_of_elements,
  GtClusteredSetClusterNumFunc cluster_num)
{
  GtClusteredSetClass  *c_class = gt_class_alloc(sizeof *c_class);
  c_class->size = size;
  c_class->number_of_clusters = number_of_clusters;
  c_class->merge_clusters = merge_clusters;
  c_class->free = free;
  c_class->iterator = iterator;
  c_class->number_of_elements = number_of_elements;
  c_class->cluster_num = cluster_num;
  return c_class;
}

GtClusteredSet* gt_clustered_set_create(const GtClusteredSetClass *csc)
{
  GtClusteredSet *cs;
  gt_assert(csc && csc->size);
  cs = gt_calloc(1, csc->size) ;
  cs->c_class = csc;
  return cs;
}

void gt_clustered_set_delete(GtClusteredSet *cs, GtError *err)
{
 if (!cs) return;
 if (cs->c_class->free)
  cs->c_class->free(cs, err);
 gt_free(cs);
}

int gt_clustered_set_merge_clusters(GtClusteredSet *cs,
                                    unsigned long c1,
                                    unsigned long c2,
                                    GtError *err)
{
  return cs->c_class->merge_clusters(cs, c1, c2, err);
}

unsigned long gt_clustered_set_num_of_clusters(GtClusteredSet *cs, GtError *err)
{
  return cs->c_class->number_of_clusters(cs, err);
}

unsigned long gt_clustered_set_num_of_elements(GtClusteredSet *cs,
                                                  GtError *err)
{
  return cs->c_class->number_of_elements(cs, err);
}

unsigned long gt_clustered_set_cluster_num(GtClusteredSet *cs,
                                           GtError *err,
                                           unsigned long e)
{
  return cs->c_class->cluster_num(cs, e, err);
}

GtClusteredSetIterator* gt_clustered_set_get_iterator(GtClusteredSet *cs,
                                                      unsigned long c,
                                                      GtError *err)
{
    return cs->c_class->iterator(cs, c, err);
}

GtClusteredSetIteratorStatus
gt_clustered_set_iterator_next(GtClusteredSetIterator *cs_i,
                               unsigned long *elm,
                               GT_UNUSED GtError *err)
{
  if (cs_i->curpos >= cs_i->length)
    return GT_CLUSTERED_SET_ITERATOR_STATUS_END;
  *elm = cs_i->elems[cs_i->curpos];
  cs_i->curpos++;
  return GT_CLUSTERED_SET_ITERATOR_STATUS_OK;
}

void gt_clustered_set_iterator_delete(GtClusteredSetIterator *cs_i,
                                      GT_UNUSED GtError *err)
{
  if (!cs_i)
    return;
  gt_free(cs_i->elems);
  gt_free(cs_i);
  return;
}
