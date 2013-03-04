/*
  Copyright (c) 2008-2010 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2008-2010 Center for Bioinformatics, University of Hamburg

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

#include <string.h>
#include "core/assert_api.h"
#include "core/class_alloc.h"
#include "core/array.h"
#include "core/ensure.h"
#include "core/ma.h"
#include "core/thread_api.h"
#include "core/unused_api.h"
#include "extended/feature_index_rep.h"
#include "extended/feature_node.h"
#ifdef GT_THREADS_ENABLED
#include "extended/genome_node.h"
#endif
#include "extended/feature_visitor.h"
#include "extended/genome_node.h"
#include "extended/gff3_in_stream.h"

struct GtFeatureIndexClass {
  size_t size;
  GtFeatureIndexAddRegionNodeFunc add_region_node;
  GtFeatureIndexAddFeatureNodeFunc add_feature_node;
  GtFeatureIndexRemoveNodeFunc remove_node;
  GtFeatureIndexGetFeatsForSeqidFunc get_features_for_seqid;
  GtFeatureIndexGetFeatsForRangeFunc get_features_for_range;
  GtFeatureIndexGetFirstSeqidFunc get_first_seqid;
  GtFeatureIndexSaveFunc save_func;
  GtFeatureIndexGetSeqidsFunc get_seqids;
  GtFeatureIndexGetRangeForSeqidFunc get_range_for_seqid;
  GtFeatureIndexHasSeqidFunc has_seqid;
  GtFeatureIndexFreeFunc free;
};

struct GtFeatureIndexMembers {
  unsigned int reference_count;
  GtRWLock *lock;
};

const GtFeatureIndexClass* gt_feature_index_class_new(size_t size,
                                         GtFeatureIndexAddRegionNodeFunc
                                                 add_region_node,
                                         GtFeatureIndexAddFeatureNodeFunc
                                                 add_feature_node,
                                         GtFeatureIndexRemoveNodeFunc
                                                 remove_node,
                                         GtFeatureIndexGetFeatsForSeqidFunc
                                                 get_features_for_seqid,
                                         GtFeatureIndexGetFeatsForRangeFunc
                                                 get_features_for_range,
                                         GtFeatureIndexGetFirstSeqidFunc
                                                 get_first_seqid,
                                         GtFeatureIndexSaveFunc
                                                 save_func,
                                         GtFeatureIndexGetSeqidsFunc
                                                 get_seqids,
                                         GtFeatureIndexGetRangeForSeqidFunc
                                                 get_range_for_seqid,
                                         GtFeatureIndexHasSeqidFunc
                                                 has_seqid,
                                         GtFeatureIndexFreeFunc
                                                 free)
{
  GtFeatureIndexClass *c_class = gt_class_alloc(sizeof *c_class);
  c_class->size = size;
  c_class->add_region_node = add_region_node;
  c_class->add_feature_node = add_feature_node;
  c_class->remove_node = remove_node;
  c_class->get_features_for_seqid = get_features_for_seqid;
  c_class->get_features_for_range = get_features_for_range;
  c_class->get_first_seqid = get_first_seqid;
  c_class->save_func = save_func;
  c_class->get_seqids = get_seqids;
  c_class->get_range_for_seqid = get_range_for_seqid;
  c_class->has_seqid = has_seqid;
  c_class->free = free;
  return c_class;
}

GtFeatureIndex* gt_feature_index_create(const GtFeatureIndexClass *fic)
{
  GtFeatureIndex *fi;
  gt_assert(fic && fic->size);
  fi = gt_calloc(1, fic->size);
  fi->c_class = fic;
  fi->pvt = gt_calloc(1, sizeof (GtFeatureIndexMembers));
  fi->pvt->lock = gt_rwlock_new();
  return fi;
}

GtFeatureIndex* gt_feature_index_ref(GtFeatureIndex *fi)
{
  gt_assert(fi);
  gt_rwlock_wrlock(fi->pvt->lock);
  fi->pvt->reference_count++;
  gt_rwlock_unlock(fi->pvt->lock);
  return fi;
}

void gt_feature_index_delete(GtFeatureIndex *fi)
{
  if (!fi) return;
  gt_rwlock_wrlock(fi->pvt->lock);
  if (fi->pvt->reference_count) {
    fi->pvt->reference_count--;
    gt_rwlock_unlock(fi->pvt->lock);
    return;
  }
  gt_assert(fi->c_class);
  if (fi->c_class->free)
    fi->c_class->free(fi);
  gt_rwlock_unlock(fi->pvt->lock);
  gt_rwlock_delete(fi->pvt->lock);
  gt_free(fi->pvt);
  gt_free(fi);
}

int  gt_feature_index_add_region_node(GtFeatureIndex *feature_index,
                                      GtRegionNode *region_node,
                                      GtError *err)
{
  int ret;
  gt_assert(feature_index && feature_index->c_class && region_node);
  gt_rwlock_wrlock(feature_index->pvt->lock);
  ret = feature_index->c_class->add_region_node(feature_index, region_node,
                                                err);
  gt_rwlock_unlock(feature_index->pvt->lock);
  return ret;
}

int  gt_feature_index_add_feature_node(GtFeatureIndex *feature_index,
                                       GtFeatureNode *feature_node,
                                       GtError *err)
{
  int ret;
  gt_assert(feature_index && feature_index->c_class && feature_node);
  gt_rwlock_wrlock(feature_index->pvt->lock);
  ret = feature_index->c_class->add_feature_node(feature_index, feature_node,
                                                 err);
  gt_rwlock_unlock(feature_index->pvt->lock);
  return ret;
}

int  gt_feature_index_remove_node(GtFeatureIndex *feature_index,
                                  GtFeatureNode *node, GtError *err)
{
  int ret;
  gt_assert(feature_index && feature_index->c_class && node);
  gt_rwlock_wrlock(feature_index->pvt->lock);
  ret = feature_index->c_class->remove_node(feature_index, node, err);
  gt_rwlock_unlock(feature_index->pvt->lock);
  return ret;
}

int gt_feature_index_add_gff3file(GtFeatureIndex *feature_index,
                                  const char *gff3file, GtError *err)
{
  GtNodeStream *gff3_in_stream;
  GtGenomeNode *gn;
  GtArray *tmp;
  int had_err = 0;
  unsigned long i;
  gt_error_check(err);
  gt_assert(feature_index && gff3file);
  tmp = gt_array_new(sizeof (GtGenomeNode*));
  gff3_in_stream = gt_gff3_in_stream_new_unsorted(1, &gff3file);
  while (!(had_err = gt_node_stream_next(gff3_in_stream, &gn, err)) && gn)
    gt_array_add(tmp, gn);
  if (!had_err) {
    GtNodeVisitor *feature_visitor = gt_feature_visitor_new(feature_index);
    for (i=0;i<gt_array_size(tmp);i++) {
      gn = *(GtGenomeNode**) gt_array_get(tmp, i);
      /* no need to lock, add_*_node() is synchronized */
      had_err = gt_genome_node_accept(gn, feature_visitor, NULL);
      gt_assert(!had_err); /* cannot happen */
    }
    gt_node_visitor_delete(feature_visitor);
  }
  gt_node_stream_delete(gff3_in_stream);
  for (i=0;i<gt_array_size(tmp);i++)
    gt_genome_node_delete(*(GtGenomeNode**) gt_array_get(tmp, i));
  gt_array_delete(tmp);
  return had_err;
}

GtArray* gt_feature_index_get_features_for_seqid(GtFeatureIndex *fi,
                                                 const char *seqid,
                                                 GtError *err)
{
  GtArray *arr;
  gt_assert(fi && fi->c_class && seqid);
  gt_rwlock_rdlock(fi->pvt->lock);
  arr = fi->c_class->get_features_for_seqid(fi, seqid, err);
  gt_rwlock_unlock(fi->pvt->lock);
  return arr;
}

int gt_feature_index_get_features_for_range(GtFeatureIndex *feature_index,
                                            GtArray *results,
                                            const char *seqid,
                                            const GtRange *range, GtError *err)
{
  int ret;
  gt_assert(feature_index && feature_index->c_class && results && seqid &&
            range);
  gt_assert(gt_range_length(range) > 0);
  gt_rwlock_rdlock(feature_index->pvt->lock);
  ret = feature_index->c_class->get_features_for_range(feature_index, results,
                                                       seqid, range, err);
  gt_rwlock_unlock(feature_index->pvt->lock);
  return ret;
}

char* gt_feature_index_get_first_seqid(const GtFeatureIndex
                                             *feature_index,
                                              GtError *err)
{
  const char *str;
  gt_assert(feature_index && feature_index->c_class);
  gt_rwlock_rdlock(feature_index->pvt->lock);
  str = feature_index->c_class->get_first_seqid(feature_index, err);
  gt_rwlock_unlock(feature_index->pvt->lock);
  return (char*) str;
}

int gt_feature_index_save(GtFeatureIndex *feature_index, GtError *err)
{
  int ret;
  gt_assert(feature_index && feature_index->c_class);
  gt_rwlock_wrlock(feature_index->pvt->lock);
  ret = feature_index->c_class->save_func(feature_index, err);
  gt_rwlock_unlock(feature_index->pvt->lock);
  return ret;
}

GtStrArray* gt_feature_index_get_seqids(const GtFeatureIndex *feature_index,
                                        GtError *err)
{
  GtStrArray *strarr;
  gt_assert(feature_index && feature_index->c_class);
  gt_rwlock_rdlock(feature_index->pvt->lock);
  strarr = feature_index->c_class->get_seqids(feature_index, err);
  gt_rwlock_unlock(feature_index->pvt->lock);
  return strarr;
}

int  gt_feature_index_get_range_for_seqid(GtFeatureIndex *feature_index,
                                          GtRange *range,
                                          const char *seqid,
                                          GtError *err)
{
  int ret;
  gt_assert(feature_index && feature_index->c_class && range && seqid);
  gt_rwlock_rdlock(feature_index->pvt->lock);
  ret = feature_index->c_class->get_range_for_seqid(feature_index, range,
                                                    seqid, err);
  gt_rwlock_unlock(feature_index->pvt->lock);
  return ret;
}

int  gt_feature_index_has_seqid(const GtFeatureIndex *feature_index,
                                bool *has_seqid,
                                const char *seqid,
                                GtError *err)
{
  int ret;
  gt_assert(feature_index && feature_index->c_class && seqid);
  gt_rwlock_rdlock(feature_index->pvt->lock);
  ret = feature_index->c_class->has_seqid(feature_index, has_seqid,
                                          seqid, err);
  gt_rwlock_unlock(feature_index->pvt->lock);
  return ret;
}

void* gt_feature_index_cast(GT_UNUSED const GtFeatureIndexClass *fic,
                            GtFeatureIndex *fi)
{
  gt_assert(fic && fi && fi->c_class == fic);
  return fi;
}

#define GT_FI_TEST_FEATURES_PER_THREAD 1000
#define GT_FI_TEST_START 1
#define GT_FI_TEST_END 10000000
#define GT_FI_TEST_FEATURE_WIDTH 2000
#define GT_FI_TEST_QUERY_WIDTH 50000
#define GT_FI_TEST_SEQID "testseqid"

typedef struct {
  GtFeatureIndex *fi;
  GtError *err;
  GtArray *nodes;
  GtMutex *mutex;
  unsigned long next_node_idx,
                error_count;
} GtFeatureIndexTestShared;

static void* gt_feature_index_unit_test_add(void *data)
{
  GtFeatureNode *feature;
  GtFeatureIndexTestShared *shm = (GtFeatureIndexTestShared*) data;

  while (true) {
    gt_mutex_lock(shm->mutex);
    if (shm->next_node_idx == GT_FI_TEST_FEATURES_PER_THREAD * gt_jobs) {
      gt_mutex_unlock(shm->mutex);
      return NULL;
    }
    gt_assert(shm->next_node_idx < gt_array_size(shm->nodes));
    feature = *(GtFeatureNode**) gt_array_get(shm->nodes, shm->next_node_idx++);
    gt_mutex_unlock(shm->mutex);
    gt_feature_index_add_feature_node(shm->fi, feature, shm->err);
    gt_genome_node_delete((GtGenomeNode*) feature);
  }
  return NULL;
}

static int cmp_range_start(const void *v1, const void *v2)
{
  return gt_genome_node_compare((GtGenomeNode**) v1, (GtGenomeNode**) v2);
}

static void* gt_feature_index_unit_test_query(void *data)
{
  GtFeatureIndexTestShared *shm = (GtFeatureIndexTestShared*) data;
  GtRange rng;
  GtError *err = shm->err;
  unsigned long i;
  int had_err = 0;
  GtArray *arr, *arr_ref;

  gt_mutex_lock(shm->mutex);
  if (gt_error_is_set(shm->err)) {
    gt_mutex_unlock(shm->mutex);
    return NULL;
  }
  gt_mutex_unlock(shm->mutex);

  arr = gt_array_new(sizeof (GtFeatureNode*));
  arr_ref = gt_array_new(sizeof (GtFeatureNode*));
  rng.start = random() % (GT_FI_TEST_END - GT_FI_TEST_QUERY_WIDTH);
  rng.end = rng.start + random() % (GT_FI_TEST_QUERY_WIDTH);

  /* get reference set by linear search */
  gt_mutex_lock(shm->mutex);
  for (i=0; i<GT_FI_TEST_FEATURES_PER_THREAD * gt_jobs; i++) {
    GtRange rng2;
    GtFeatureNode *fn;
    fn = *(GtFeatureNode**) gt_array_get(shm->nodes, i);
    rng2 = gt_genome_node_get_range((GtGenomeNode*) fn);
    if (gt_range_overlap(&rng, &rng2)) {
      gt_array_add(arr_ref, fn);
    }
  }
  gt_mutex_unlock(shm->mutex);

  /* query feature index */
  gt_feature_index_get_features_for_range(shm->fi, arr, GT_FI_TEST_SEQID, &rng,
                                          err);

  /* result size must be equal */
  if (gt_array_size(arr) != gt_array_size(arr_ref))
    had_err = -1;

  /* nodes must be the same (note that we should not rely on ptr equality) */
  if (!had_err) {
    gt_array_sort(arr_ref, cmp_range_start);
    gt_array_sort(arr    , cmp_range_start);

    for (i=0;i<gt_array_size(arr);i++) {
      if (had_err)
        break;
      if (!gt_feature_node_is_similar(*(GtFeatureNode**) gt_array_get(arr, i),
                                      *(GtFeatureNode**)
                                      gt_array_get(arr_ref, i))) {
        had_err = -1;
      }
    }
  }

  if (had_err) {
    gt_mutex_lock(shm->mutex);
    shm->error_count++;
    gt_mutex_unlock(shm->mutex);
  }

  gt_array_delete(arr);
  gt_array_delete(arr_ref);
  return NULL;
}

/* to be called from implementing class! */
int gt_feature_index_unit_test(GtFeatureIndex *fi, GtError *err)
{
  int had_err = 0, i, rval;
  GtFeatureIndexTestShared sh;
  GtStrArray *seqids;
  GtStr *seqid;
  GtRange check_range;
  GtRegionNode *rn;
  bool has_seqid;
  gt_error_check(err);

  sh.mutex = gt_mutex_new();
  sh.nodes = gt_array_new(sizeof (GtFeatureNode*));
  sh.error_count = 0;
  sh.next_node_idx = 0;
  sh.fi = fi;
  sh.err = gt_error_new();

  /* create region */
  seqid = gt_str_new_cstr(GT_FI_TEST_SEQID);
  rn = (GtRegionNode*) gt_region_node_new(seqid, GT_FI_TEST_START,
                                          GT_FI_TEST_END);

  /* test seqid is not supposed to exist */
  gt_ensure(had_err, gt_feature_index_has_seqid(sh.fi, &has_seqid,
                                                 GT_FI_TEST_SEQID, err) == 0);
  gt_ensure(had_err, !has_seqid);

  /* add a sequence region directly and check if it has been added */
  rval = gt_feature_index_add_region_node(sh.fi, rn, err);
  gt_ensure(had_err, rval == 0);
  gt_genome_node_delete((GtGenomeNode*) rn);
  gt_ensure(had_err, gt_feature_index_has_seqid(sh.fi, &has_seqid,
                                                GT_FI_TEST_SEQID, err) == 0);
  gt_ensure(had_err, has_seqid);

  gt_feature_index_get_range_for_seqid(sh.fi, &check_range, GT_FI_TEST_SEQID,
                                       err);
  gt_ensure(had_err, check_range.start == GT_FI_TEST_START
                    && check_range.end == GT_FI_TEST_END);

  /* set up nodes to store */
  for (i=0;i<GT_FI_TEST_FEATURES_PER_THREAD*gt_jobs;i++) {
    unsigned long start, end;
    GtFeatureNode *fn;
    start = random() % (GT_FI_TEST_END - GT_FI_TEST_FEATURE_WIDTH);
    end = start + random() % (GT_FI_TEST_FEATURE_WIDTH);
    fn = gt_feature_node_cast(gt_feature_node_new(seqid, "gene", start, end,
                                                  GT_STRAND_FORWARD));
    gt_array_add(sh.nodes, fn);
  }
  /* test parallel addition */
  gt_multithread(gt_feature_index_unit_test_add, &sh, err);
  seqids = gt_feature_index_get_seqids(fi, err);
  gt_ensure(had_err, seqids);
  gt_ensure(had_err, gt_feature_index_has_seqid(fi, &has_seqid,GT_FI_TEST_SEQID,
                                                err) == 0);
  gt_ensure(had_err, has_seqid);
  gt_ensure(had_err, gt_str_array_size(seqids) == 1);

  /* test parallel query */
  if (!had_err)
    gt_multithread(gt_feature_index_unit_test_query, &sh, err);
  gt_ensure(had_err, sh.error_count == 0);

  gt_mutex_delete(sh.mutex);
  gt_error_delete(sh.err);
  gt_str_array_delete(seqids);
  gt_array_delete(sh.nodes);
  gt_str_delete(seqid);
  return had_err;
}
