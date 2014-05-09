/*
  Copyright (c) 2006-2011 Gordon Gremme <gordon@gremme.org>
  Copyright (c) 2006-2008 Center for Bioinformatics, University of Hamburg

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

#include <stdarg.h>
#include "core/assert_api.h"
#include "core/class_alloc.h"
#include "core/cstr_api.h"
#include "core/ensure.h"
#include "core/hashtable.h"
#include "core/ma.h"
#include "core/md5_seqid.h"
#include "core/msort.h"
#include "core/queue_api.h"
#include "core/unused_api.h"
#include "extended/eof_node_api.h"
#include "extended/genome_node_rep.h"
#include "extended/gff3_defines.h"
#include "extended/gff3_visitor_api.h"
#include "extended/region_node_api.h"

typedef struct {
  void *ptr;
  GtFree free_func;
} GtGenomeNodeUserData;

GtGenomeNode* gt_genome_node_ref(GtGenomeNode *gn)
{
  gt_assert(gn);
  gt_rwlock_wrlock(gn->lock);
  gn->reference_count++;
  gt_rwlock_unlock(gn->lock);
  return gn;
}

const GtGenomeNodeClass*
gt_genome_node_class_new(size_t size,
                         GtGenomeNodeFreeFunc free,
                         GtGenomeNodeSetSeqidFunc get_seqid,
                         GtGenomeNodeGetIdstrFunc get_idstr,
                         GtGenomeNodeGetRangeFunc get_range,
                         GtGenomeNodeSetRangeFunc set_range,
                         GtGenomeNodeChangeSeqidFunc change_seqid,
                         GtGenomeNodeAcceptFunc accept)
{
  GtGenomeNodeClass *c_class;
  gt_assert(size);
  c_class = gt_class_alloc(sizeof *c_class);
  c_class->size = size;
  c_class->free = free;
  c_class->get_seqid = get_seqid;
  c_class->get_idstr = get_idstr;
  c_class->get_range = get_range;
  c_class->set_range = set_range;
  c_class->change_seqid = change_seqid;
  c_class->accept = accept;
  return c_class;
}

static void userdata_delete(void *data)
{
  GtGenomeNodeUserData *ud;
  if (!data) return;
  ud = (GtGenomeNodeUserData*) data;
  if (ud->free_func)
    ud->free_func(ud->ptr);
  ud->free_func = ud->ptr = NULL;
  gt_free(ud);
}

static int compare_genome_node_type(GtGenomeNode *gn_a, GtGenomeNode *gn_b)
{
  void *rn_a, *rn_b, *sn_a, *sn_b, *en_a, *en_b;
  GtMetaNode *mn_a, *mn_b;

  /* meta nodes first */
  mn_a = gt_meta_node_try_cast(gn_a);
  mn_b = gt_meta_node_try_cast(gn_b);

  if (mn_a && !mn_b)
    return -1;
  if (!mn_a && mn_b)
    return 1;
  if (mn_a && mn_b) {
    int ret = 0;
    if (!ret && strcmp(gt_meta_node_get_directive(mn_a),
                       GT_GFF_VERSION_DIRECTIVE) == 0)
      ret = -1;
    if (!ret && strcmp(gt_meta_node_get_directive(mn_b),
                       GT_GFF_VERSION_DIRECTIVE) == 0)
      ret = 1;
    if (!ret && strcmp(gt_meta_node_get_directive(mn_a),
                       GT_GVF_VERSION_DIRECTIVE) == 0)
      ret = -1;
    if (!ret && strcmp(gt_meta_node_get_directive(mn_b),
                       GT_GVF_VERSION_DIRECTIVE) == 0)
      ret = 1;
    return ret;
  }

  /* then region nodes */
  rn_a = gt_region_node_try_cast(gn_a);
  rn_b = gt_region_node_try_cast(gn_b);

  if (rn_a && !rn_b)
    return -1;
  if (!rn_a && rn_b)
    return 1;

  /* sequence nodes last, but before eof nodes */
  sn_a = gt_sequence_node_try_cast(gn_a);
  sn_b = gt_sequence_node_try_cast(gn_b);
  en_a = gt_eof_node_try_cast(gn_a);
  en_b = gt_eof_node_try_cast(gn_b);

  if (sn_a && !sn_b)
    return en_b ? -1 : 1;
  if (!sn_a && sn_b)
    return en_a ? 1 : -1;
  if (en_a && !en_b)
    return 1;
  if (!en_a && en_b)
    return -1;

  return 0;
}

int gt_genome_node_cmp(GtGenomeNode *gn_a, GtGenomeNode *gn_b)
{
  GtRange range_a, range_b;
  int rval;
  const char *id_a, *id_b;
  gt_assert(gn_a && gn_b);
  /* ensure that region nodes come first and sequence nodes come last,
     otherwise we don't get a valid GFF3 stream */
  if ((rval = compare_genome_node_type(gn_a, gn_b)))
    return rval;

  id_a = gt_str_get(gt_genome_node_get_idstr(gn_a));
  id_b = gt_str_get(gt_genome_node_get_idstr(gn_b));
  if ((rval = gt_md5_seqid_cmp_seqids(id_a, id_b))) {
    return rval;
  }
  range_a = gt_genome_node_get_range(gn_a),
  range_b = gt_genome_node_get_range(gn_b);
  return gt_range_compare(&range_a, &range_b);
}

static int compare_genome_nodes_with_delta(GtGenomeNode *gn_a,
                                           GtGenomeNode *gn_b,
                                           GtUword delta)
{
  GtRange range_a, range_b;
  int rval;
  const char *id_a, *id_b;
  gt_assert(gn_a && gn_b);
  /* ensure that sequence regions come first, otherwise we don't get a valid
     gff3 stream */
  if ((rval = compare_genome_node_type(gn_a, gn_b)))
    return rval;

  id_a = gt_str_get(gt_genome_node_get_idstr(gn_a));
  id_b = gt_str_get(gt_genome_node_get_idstr(gn_b));
  if ((rval = gt_md5_seqid_cmp_seqids(id_a, id_b))) {
    return rval;
  }
  range_a = gt_genome_node_get_range(gn_a);
  range_b = gt_genome_node_get_range(gn_b);
  return gt_range_compare_with_delta(&range_a, &range_b, delta);
}

GtGenomeNode* gt_genome_node_create(const GtGenomeNodeClass *gnc)
{
  GtGenomeNode *gn;
  gt_assert(gnc && gnc->size);
  gn                     = gt_malloc(gnc->size);
  gn->c_class            = gnc;
  gn->filename           = NULL; /* means the node is generated */
  gn->line_number        = 0;
  gn->reference_count    = 0;
  gn->userdata           = NULL;
  gn->userdata_nof_items = 0;
#ifdef GT_THREADS_ENABLED
  gn->lock              = gt_rwlock_new();
#endif
  return gn;
}

void gt_genome_node_set_origin(GtGenomeNode *gn, GtStr *filename,
                               unsigned int line_number)
{
  gt_assert(gn && filename && line_number);
  gt_str_delete(gn->filename);
  gn->filename = gt_str_ref(filename);
  gn->line_number = line_number;
}

void* gt_genome_node_cast(GT_UNUSED const GtGenomeNodeClass *gnc,
                          GtGenomeNode *gn)
{
  gt_assert(gnc && gn && gn->c_class == gnc);
  return gn;
}

void* gt_genome_node_try_cast(const GtGenomeNodeClass *gnc, GtGenomeNode *gn)
{
  gt_assert(gnc && gn);
  if (gn->c_class == gnc)
    return gn;
  return NULL;
}

const char* gt_genome_node_get_filename(const GtGenomeNode *gn)
{
  gt_assert(gn);
  if (gn->filename)
    return gt_str_get(gn->filename);
  return "generated";
}

unsigned int gt_genome_node_get_line_number(const GtGenomeNode *gn)
{
  gt_assert(gn);
  return gn->line_number;
}

GtStr* gt_genome_node_get_seqid(GtGenomeNode *gn)
{
  gt_assert(gn && gn->c_class);
  if (gn->c_class->get_seqid)
    return gn->c_class->get_seqid(gn);
  return NULL;
}

GtStr* gt_genome_node_get_idstr(GtGenomeNode *gn)
{
  gt_assert(gn && gn->c_class && gn->c_class->get_idstr);
  return gn->c_class->get_idstr(gn);
}

GtUword gt_genome_node_get_start(GtGenomeNode *gn)
{
  return gt_genome_node_get_range(gn).start;
}

GtUword gt_genome_node_get_end(GtGenomeNode *gn)
{
  return gt_genome_node_get_range(gn).end;
}

GtRange gt_genome_node_get_range(GtGenomeNode *gn)
{
  gt_assert(gn && gn->c_class && gn->c_class->get_range);
  return gn->c_class->get_range(gn);
}

GtUword gt_genome_node_get_length(GtGenomeNode *gn)
{
  GtRange range;
  gt_assert(gn && gn->c_class && gn->c_class->get_range);
  range = gt_genome_node_get_range(gn);
  return gt_range_length(&range);
}

void gt_genome_node_set_range(GtGenomeNode *gn, const GtRange *range)
{
  gt_assert(gn && gn->c_class && gn->c_class->set_range);
  gt_assert(range->start <= range->end);
  gn->c_class->set_range(gn, range);
}

void gt_genome_node_change_seqid(GtGenomeNode *gn, GtStr *seqid)
{
  gt_assert(gn && gn->c_class && gn->c_class->change_seqid && seqid);
  gn->c_class->change_seqid(gn, seqid);
}

int gt_genome_node_accept(GtGenomeNode *gn, GtNodeVisitor *nv, GtError *err)
{
  gt_error_check(err);
  gt_assert(gn && nv && gn->c_class && gn->c_class->accept);
  return gn->c_class->accept(gn, nv, err);
}

int gt_genome_node_compare(GtGenomeNode **gn_a, GtGenomeNode **gn_b)
{
  return gt_genome_node_cmp(*gn_a, *gn_b);
}

int gt_genome_node_compare_with_data(GtGenomeNode **gn_a, GtGenomeNode **gn_b,
                                     GT_UNUSED void *unused)
{
  return gt_genome_node_cmp(*gn_a, *gn_b);
}

int gt_genome_node_compare_delta(GtGenomeNode **gn_a, GtGenomeNode **gn_b,
                                 void *delta)
{
  GtUword *deltaptr = delta;
  gt_assert(delta);
  return compare_genome_nodes_with_delta(*gn_a, *gn_b, *deltaptr);
}

void gt_genome_nodes_sort(GtArray *nodes)
{
  qsort(gt_array_get_space(nodes), gt_array_size(nodes),
        sizeof (GtGenomeNode*), (GtCompare) gt_genome_node_compare);
}

void gt_genome_nodes_sort_stable(GtArray *nodes)
{
  gt_msort(gt_array_get_space(nodes), gt_array_size(nodes),
           sizeof (GtGenomeNode*), (GtCompare) gt_genome_node_compare);
}

void gt_genome_nodes_show(GtArray *nodes, GtFile *outfp)
{
  GtNodeVisitor *gff3_visitor;
  GtUword i;
  gt_assert(nodes);
  gff3_visitor = gt_gff3_visitor_new(outfp);
  for (i = 0; i < gt_array_size(nodes); i++) {
    GT_UNUSED int had_err;
    had_err = gt_genome_node_accept(*(GtGenomeNode**) gt_array_get(nodes, i),
                                    gff3_visitor, NULL);
    gt_assert(!had_err); /* should not happen */
  }
  gt_node_visitor_delete(gff3_visitor);
}

bool gt_genome_nodes_are_equal_region_nodes(GtGenomeNode *gn_a,
                                            GtGenomeNode *gn_b)
{
  void *sr_a, *sr_b;

  sr_a = gn_a ? gt_region_node_try_cast(gn_a) : NULL;
  sr_b = gn_b ? gt_region_node_try_cast(gn_b) : NULL;

  if (sr_a && sr_b && !gt_str_cmp(gt_genome_node_get_seqid(gn_a),
                                  gt_genome_node_get_seqid(gn_b))) {
    return true;
  }
  return false;
}

bool gt_genome_nodes_are_sorted(const GtArray *nodes)
{
  GtUword i;
  gt_assert(nodes);
  for (i = 1; i < gt_array_size(nodes); i++) {
    if (gt_genome_node_compare(gt_array_get(nodes, i-1),
                               gt_array_get(nodes, i)) > 0) {
      return false;
    }
  }
  return true;
}

void gt_genome_node_add_user_data(GtGenomeNode *gn, const char *key, void *data,
                                  GtFree free_func)
{
  GtGenomeNodeUserData *ud, *myud;
  gt_assert(gn && key);
  ud = gt_malloc(sizeof (GtGenomeNodeUserData));
  ud->ptr = data;
  ud->free_func = free_func;
  if (!gn->userdata) {
    gn->userdata = gt_hashmap_new(GT_HASH_STRING, gt_free_func,
                                  userdata_delete);
  }
  gt_assert(gn->userdata != NULL);
  /* remove old data entry first if there is one */
  if ((myud = gt_hashmap_get(gn->userdata, key)) != NULL) {
    gt_hashmap_remove(gn->userdata, key);
  } else {
    gn->userdata_nof_items++;
  }
  gt_hashmap_add(gn->userdata, gt_cstr_dup((char*) key), ud);
}

void* gt_genome_node_get_user_data(const GtGenomeNode *gn, const char *key)
{
  GtGenomeNodeUserData *ud;
  gt_assert(gn && key);
  if (!gn->userdata)
    return NULL;
  ud = (GtGenomeNodeUserData*) gt_hashmap_get(gn->userdata, key);
  return (ud ? ud->ptr : NULL);
}

void gt_genome_node_release_user_data(GtGenomeNode *gn, const char *key)
{
  GtGenomeNodeUserData *ud;
  gt_assert(gn && key);
  if (!gn->userdata)
    return;
  if ((ud = (GtGenomeNodeUserData*) gt_hashmap_get(gn->userdata, key))) {
    gt_hashmap_remove(gn->userdata, (char*) key);
    if (--gn->userdata_nof_items == 0) {
      gt_hashmap_delete(gn->userdata);
      gn->userdata = NULL;
    }
  }
}

int gt_genome_node_unit_test(GtError *err)
{
  int had_err = 0;
  GtGenomeNodeClass *gnc;
  GtGenomeNode *gn;
  void *testptr1 = "foo bar",
       *testptr2 = NULL;
  const char *testkey1 = "key1",
             *testkey2 = "key2";
  gt_error_check(err);

  testptr2 = gt_malloc(sizeof (char)*4);

  /* this is a very simple GtGenomeNodeClass without any callbacks */
  gnc = gt_calloc(1, sizeof (GtGenomeNodeClass));
  gnc->size = sizeof (GtGenomeNode);

  gn = gt_genome_node_create(gnc);
  gt_ensure(gt_genome_node_get_user_data(gn, testkey1) == NULL);
  gt_ensure(gt_genome_node_get_user_data(gn, testkey2) == NULL);
  gt_ensure(gn->userdata_nof_items == 0);
  gt_ensure(gn->userdata == NULL);

  gt_genome_node_add_user_data(gn, testkey1, testptr1, NULL);
  gt_ensure(gt_genome_node_get_user_data(gn, testkey1) != NULL);
  gt_ensure(gt_genome_node_get_user_data(gn, testkey1) == testptr1);
  gt_ensure(gn->userdata_nof_items == 1);
  gt_ensure(gn->userdata != NULL);

  gt_genome_node_add_user_data(gn, testkey2, testptr2, gt_free_func);
  gt_ensure(gt_genome_node_get_user_data(gn, testkey2) != NULL);
  gt_ensure(gt_genome_node_get_user_data(gn, testkey2) == testptr2);
  gt_ensure(gn->userdata_nof_items == 2);

  gt_genome_node_release_user_data(gn, testkey1);
  gt_ensure(gt_genome_node_get_user_data(gn, testkey1) == NULL);
  gt_ensure(gn->userdata_nof_items == 1);

  gt_genome_node_release_user_data(gn, testkey2);
  gt_ensure(gt_genome_node_get_user_data(gn, testkey2) == NULL);
  gt_ensure(gn->userdata_nof_items == 0);
  gt_ensure(gn->userdata == NULL);

  testptr2 = gt_malloc(sizeof (char)*4);
  gt_genome_node_add_user_data(gn, testkey1, testptr1, NULL);
  gt_genome_node_add_user_data(gn, testkey2, testptr2, gt_free_func);
  gt_ensure(gn->userdata != NULL);
  gt_genome_node_delete(gn);

  gt_free(gnc);

  return had_err;
}

void gt_genome_node_delete(GtGenomeNode *gn)
{
  if (!gn) return;
  gt_rwlock_wrlock(gn->lock);
  if (gn->reference_count) {
    gn->reference_count--;
    gt_rwlock_unlock(gn->lock);
    return;
  }
  gt_assert(gn->c_class);
  if (gn->c_class->free)
    gn->c_class->free(gn);
  gt_str_delete(gn->filename);
  if (gn->userdata)
    gt_hashmap_delete(gn->userdata);
  gt_rwlock_unlock(gn->lock);
#ifdef GT_THREADS_ENABLED
  gt_rwlock_delete(gn->lock);
#endif
  gt_free(gn);
}
