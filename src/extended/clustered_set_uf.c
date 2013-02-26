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

#include <stdio.h>
#include "core/array_api.h"
#include "core/bittab_api.h"
#include "core/class_alloc_lock.h"
#include "core/ensure.h"
#include "core/ma.h"
#include "core/unused_api.h"
#include "extended/clustered_set.h"
#include "extended/clustered_set_uf.h"
#include "extended/clustered_set_rep.h"

#define CLUSTERED_SET_UNION_FIND_TEST_SIZE  1024

#define CLUSTERNIL (cs_uf->num_of_elems)

#define SINGLETON(E) !gt_bittab_bit_is_set(cs_uf->in_cluster, E)

#define CINFO(C) ((GtClusteredSetUFClusterInfo*)\
                 gt_array_get(cs_uf->cluster_info, C))

#define CHECKCLUSTER(C)\
        if ((C) >= cs_uf->next_free_cluster_info)\
        {\
         gt_error_set(err, "cluster %lu is too large, must be smaller than"\
                      "%lu", C, cs_uf->next_free_cluster_info);\
          had_err = -1;\
        }\
        if (CINFO(C)->cluster_size == 0)\
        {\
         gt_error_set(err, "cluster %lu is empty", C);\
         had_err = -1;\
        }

struct GtClusteredSetUFElem {
  unsigned long cluster_num, next_elem;
};

struct GtClusteredSetUFClusterInfo {
  unsigned long first_elem, last_elem, cluster_size;
};

struct GtClusteredSetUF {
  const GtClusteredSet parent_instance;
  GtClusteredSetUFElem *cluster_elems;
  GtArray *cluster_info;
  GtBittab *in_cluster;
  unsigned long num_of_elems, next_free_cluster_info;
};

static void
gt_clustered_set_union_find_change_cluster_number(GtClusteredSetUF *cs_uf,
                                                  unsigned long first_elem,
                                                  unsigned long target)
{

  unsigned long i;
  i = first_elem;
  while (1)
  {
    cs_uf->cluster_elems[i].cluster_num = target;
    if ((i = cs_uf->cluster_elems[i].next_elem) == cs_uf->num_of_elems)
      break;
  }
}

static void
gt_clustered_set_union_find_make_new_cluster(GtClusteredSetUF *cs_uf,
                                             unsigned long elem1,
                                             unsigned long elem2,
                                             GT_UNUSED GtError *err)
{
  unsigned long cluster_num;
  cluster_num = cs_uf->next_free_cluster_info++;

  cs_uf->cluster_elems[elem1].cluster_num = cluster_num;
  cs_uf->cluster_elems[elem2].cluster_num = cluster_num;
  cs_uf->cluster_elems[elem1].next_elem = elem2;
  cs_uf->cluster_elems[elem2].next_elem = cs_uf->num_of_elems;

  GtClusteredSetUFClusterInfo cluster_info;
  cluster_info.cluster_size = (unsigned long) 2;
  cluster_info.first_elem = elem1;
  cluster_info.last_elem = elem2;

  gt_array_add(cs_uf->cluster_info, cluster_info);
}

static void
gt_clustered_set_union_find_append_elem(GtClusteredSetUF *cs_uf,
                                        unsigned long c,
                                        unsigned long elem,
                                        GT_UNUSED GtError *err)
{
  GtClusteredSetUFClusterInfo *cluster_info;
  cs_uf->cluster_elems[elem].cluster_num = c;
  cs_uf->cluster_elems[elem].next_elem = cs_uf->num_of_elems;
  cluster_info = CINFO(c);
  cs_uf->cluster_elems[cluster_info->last_elem].next_elem = elem;
  cluster_info->last_elem = elem;
  cluster_info->cluster_size++;
}

static void
gt_clustered_set_union_find_join_clusters(GtClusteredSetUF *cs_uf,
                                           unsigned long c1,
                                           unsigned long c2,
                                           GT_UNUSED GtError *err)
{
  GtClusteredSetUFClusterInfo *cluster_info_c1;
  GtClusteredSetUFClusterInfo *cluster_info_c2;

  cluster_info_c1 = CINFO(c1);
  cluster_info_c2 = CINFO(c2);

  gt_clustered_set_union_find_change_cluster_number(cs_uf,
                                                    cluster_info_c2->first_elem,
                                                    c1);

  cs_uf->cluster_elems[cluster_info_c1->last_elem].next_elem =
    cluster_info_c2->first_elem;

  cluster_info_c2->first_elem = cs_uf->num_of_elems;
  cluster_info_c1->cluster_size += cluster_info_c2->cluster_size;
  cluster_info_c1->last_elem = cluster_info_c2->last_elem;
  cluster_info_c2->cluster_size = (unsigned long) 0;
}

GtClusteredSetIterator*
gt_clustered_set_union_find_iterator_new(GtClusteredSet *cs,
                                         unsigned long c,
                                         GT_UNUSED GtError *err)
{
  gt_assert(cs);
  GtClusteredSetUF *cs_uf = (GtClusteredSetUF*) cs;
  GtClusteredSetIterator *cs_i = gt_calloc(1, sizeof (GtClusteredSetIterator));
  unsigned long i = 0, j = 0;
  if (SINGLETON(c)) {
    cs_i->length = (unsigned long) 1;
    cs_i->curpos = (unsigned long) 0;
    cs_i->elems = gt_calloc(1, sizeof (unsigned long));
    cs_i->elems[j] = c;
  }
  else {
    if (CINFO(c)->cluster_size > 0) {
      GtClusteredSetUFClusterInfo *cluster_info;
      cluster_info = CINFO(c);
      cs_i->length = cluster_info->cluster_size;
      cs_i->curpos = (unsigned long) 0;
      cs_i->elems = gt_calloc(cluster_info->cluster_size,
                              sizeof (unsigned long));
      i = cluster_info->first_elem;

      do {
        cs_i->elems[j] = i;
        j++;
      } while ((i = cs_uf->cluster_elems[i].next_elem) < cs_uf->num_of_elems);

    } else {
      gt_free(cs_i);
      return NULL;
    }
  }
  return cs_i;
}

void gt_clustered_set_union_find_delete(GtClusteredSet *cs,
                                        GT_UNUSED GtError *err)
{
  if (!cs) return;
  GtClusteredSetUF *cs_uf = (GtClusteredSetUF*) cs;
  gt_free(cs_uf->cluster_elems);
  gt_array_delete(cs_uf->cluster_info);
  gt_bittab_delete(cs_uf->in_cluster);
}

unsigned long gt_clustered_set_union_find_num_of_clusters(GtClusteredSet *cs,
                                                       GT_UNUSED GtError *err)
{
  gt_assert(cs);
  GtClusteredSetUF *cs_uf = (GtClusteredSetUF*) cs;
  return cs_uf->next_free_cluster_info;
}

unsigned long gt_clustered_set_union_find_num_of_elements(GtClusteredSet *cs,
                                                       GT_UNUSED GtError *err)
{
  gt_assert(cs);
  GtClusteredSetUF *cs_uf = (GtClusteredSetUF*) cs;
  return cs_uf->num_of_elems;
}

unsigned long gt_clustered_set_union_find_cluster_num(GtClusteredSet *cs,
                                                      unsigned long e,
                                                      GT_UNUSED GtError *err
                                                      )
{
  gt_assert(cs);
  GtClusteredSetUF *cs_uf = (GtClusteredSetUF*) cs;
  if (SINGLETON(e))
    return CLUSTERNIL;
  else
    return cs_uf->cluster_elems[e].cluster_num;
}

int gt_clustered_set_union_find_merge_clusters(GtClusteredSet *cs,
                                               unsigned long e1,
                                               unsigned long e2,
                                               GtError *err)
{
  gt_assert(cs);
  int had_err = 0;
  GtClusteredSetUFClusterInfo *cluster_info_c1 = NULL;
  GtClusteredSetUFClusterInfo *cluster_info_c2 = NULL;
  unsigned long target = 0, source = 0, c1 = 0, c2 = 0;
  GtClusteredSetUF *cs_uf = (GtClusteredSetUF*) cs;
  if (e1 == e2) {
    gt_error_set(err, "expected %lu to be unequal %lu", e1, e2 );
    had_err = -1;
  }

  if (e1 >= cs_uf->num_of_elems || e2 >= cs_uf->num_of_elems) {
    gt_error_set(err, "%lu and %lu must not be larger than %lu",
                 e1, e2, cs_uf->num_of_elems);
    had_err = -1;
  }

  if (!had_err) {
    if (SINGLETON(e1)) {
       /* printf("%lu is singleton\n", e1); */
      if (SINGLETON(e2)) {
        /* printf("%lu is singleton\n", e2);*/
        gt_clustered_set_union_find_make_new_cluster(cs_uf, e1, e2, err);
        gt_bittab_set_bit(cs_uf->in_cluster, e2);
      }
      else {
        c2 = cs_uf->cluster_elems[e2].cluster_num;
        CHECKCLUSTER(c2);
        gt_clustered_set_union_find_append_elem(cs_uf, c2, e1, err);
      }
      gt_bittab_set_bit(cs_uf->in_cluster, e1);
    }
    else {
      c1 = cs_uf->cluster_elems[e1].cluster_num;
      CHECKCLUSTER(c1);
      if (SINGLETON(e2)) {
        gt_clustered_set_union_find_append_elem(cs_uf, c1, e2, err);
        gt_bittab_set_bit(cs_uf->in_cluster, e2);
      }
      else {
        c2 = cs_uf->cluster_elems[e2].cluster_num;
        CHECKCLUSTER(c2);
        cluster_info_c1 = CINFO(c1);
        cluster_info_c2 = CINFO(c2);
        if (cluster_info_c1->cluster_size > cluster_info_c2->cluster_size) {
          target = c1;
          source = c2;
        }
        else {
          target = c1;
          source = c2;
        }
        if (target != source)
          gt_clustered_set_union_find_join_clusters(cs_uf, target, source, err);
      }
    }
  }
  return had_err;
}

const GtClusteredSetClass* gt_clustered_set_union_find_class(void)
{

  static const GtClusteredSetClass *csc = NULL;
  gt_class_alloc_lock_enter();
  if (!csc) {
    csc = gt_clustered_set_class_new(sizeof (GtClusteredSetUF),
      gt_clustered_set_union_find_num_of_clusters,
      gt_clustered_set_union_find_merge_clusters,
      gt_clustered_set_union_find_delete,
      gt_clustered_set_union_find_iterator_new,
      gt_clustered_set_union_find_num_of_elements,
      gt_clustered_set_union_find_cluster_num);
  }
  gt_class_alloc_lock_leave();
  return csc;
}

GtClusteredSet*
gt_clustered_set_union_find_new(unsigned long num_of_elems,
                                GT_UNUSED GtError *err)
{
  GtClusteredSet *cs;
  GtClusteredSetUF *cs_uf;
  cs = gt_clustered_set_create(gt_clustered_set_union_find_class());
  cs_uf = (GtClusteredSetUF*) cs;
  cs_uf->num_of_elems = num_of_elems;
  cs_uf->cluster_elems = gt_calloc(num_of_elems, sizeof (GtClusteredSetUFElem));
  cs_uf->cluster_info = gt_array_new(sizeof (GtClusteredSetUFClusterInfo));
  cs_uf->in_cluster = gt_bittab_new(num_of_elems);
  cs_uf->next_free_cluster_info = (unsigned long) 0;

  return cs;
}

int gt_clustered_set_union_find_unit_test(GtError *err)
{
  int had_err = 0, i = 0, j = 0;
  GtClusteredSet *cs = NULL;
  cs = gt_clustered_set_union_find_new(1, err);
  gt_ensure(had_err,
            gt_clustered_set_union_find_num_of_elements(cs, err) ==
            gt_clustered_set_union_find_cluster_num(cs, 0, err));
  gt_ensure(had_err,
            gt_clustered_set_union_find_num_of_clusters(cs, err) == 0);
  gt_clustered_set_union_find_delete(cs, err);

  if (!had_err) {
    cs = gt_clustered_set_union_find_new(2, err);
    gt_clustered_set_union_find_merge_clusters(cs, 0, 1, err);

    gt_ensure(had_err,
              gt_clustered_set_union_find_num_of_clusters(cs, err) == 1);
    gt_ensure(had_err,
              gt_clustered_set_union_find_cluster_num(cs, 0, err) ==
              gt_clustered_set_union_find_cluster_num(cs, 1, err));
    gt_ensure(had_err,
              gt_clustered_set_union_find_num_of_clusters(cs, err) == 1);

    gt_clustered_set_union_find_delete(cs, err);
  }

  if (!had_err) {
    cs = gt_clustered_set_union_find_new(3, err);
    gt_clustered_set_union_find_merge_clusters(cs, 0, 1, err);
    gt_clustered_set_union_find_merge_clusters(cs, 1, 2, err);

    gt_ensure(had_err,
               gt_clustered_set_union_find_num_of_clusters(cs, err) == 1);
    gt_ensure(had_err,
               gt_clustered_set_union_find_cluster_num(cs, 0, err) ==
               gt_clustered_set_union_find_cluster_num(cs, 1, err));

    gt_ensure(had_err,
               gt_clustered_set_union_find_cluster_num(cs, 0, err) ==
               gt_clustered_set_union_find_cluster_num(cs, 2, err));

    gt_ensure(had_err,
              gt_clustered_set_union_find_cluster_num(cs, 1, err) ==
              gt_clustered_set_union_find_cluster_num(cs, 2, err));

    gt_clustered_set_union_find_delete(cs, err);
  }

  if (!had_err) {
    cs = gt_clustered_set_union_find_new(4, err);
    gt_clustered_set_union_find_merge_clusters(cs, 0, 1, err);
    gt_clustered_set_union_find_merge_clusters(cs, 2, 3, err);

    gt_ensure(had_err,
              gt_clustered_set_union_find_num_of_clusters(cs, err) == 2);

    gt_clustered_set_union_find_merge_clusters(cs, 0, 2, err);

    for (i = 0; i < 4 - 1; i++) {
      for (j = i + 1; j < 4; j++) {
        gt_ensure(had_err,
                  gt_clustered_set_union_find_cluster_num(cs, i, err) ==
                  gt_clustered_set_union_find_cluster_num(cs, j, err));
      }
    }
    gt_clustered_set_union_find_delete(cs, err);
  }

  if (!had_err) {
    cs = gt_clustered_set_union_find_new(CLUSTERED_SET_UNION_FIND_TEST_SIZE,
                                         err);

    for (i = 0; !had_err && i < CLUSTERED_SET_UNION_FIND_TEST_SIZE; i++) {
      gt_ensure(had_err,
                gt_clustered_set_union_find_num_of_elements(cs, err) ==
                gt_clustered_set_union_find_cluster_num(cs, i, err));
    }
    for (i = 1;!had_err && i < CLUSTERED_SET_UNION_FIND_TEST_SIZE; i++) {
      gt_clustered_set_union_find_merge_clusters(cs, 0, i, err);
    }

    for (i = 0; !had_err && i < CLUSTERED_SET_UNION_FIND_TEST_SIZE; i++) {
      gt_ensure(had_err,
                gt_clustered_set_union_find_cluster_num(cs, i, err) == 0);
    }
    gt_clustered_set_union_find_delete(cs, err);
  }

  return had_err;
}
