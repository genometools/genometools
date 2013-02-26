/*
  Copyright (c) 2010-2012 Gordon Gremme <gremme@zbh.uni-hamburg.de>

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

#include "core/assert_api.h"
#include "core/class_alloc_lock.h"
#include "core/md5_seqid.h"
#include "core/str_array.h"
#include "core/undef_api.h"
#include "extended/feature_node.h"
#include "extended/genome_node.h"
#include "extended/gff3_defines.h"
#include "extended/gff3_parser.h"
#include "extended/md5_to_id_visitor.h"
#include "extended/node_visitor_api.h"
#include "extended/regular_seqid.h"

struct GtMD5ToSeqidsVisitor {
  const GtNodeVisitor parent_instance;
  GtRegionMapping *region_mapping;
};

#define  md5_to_id_visitor_cast(GV)\
         gt_node_visitor_cast(gt_md5_to_id_visitor_class(), GV)

static void md5_to_id_visitor_free(GtNodeVisitor *nv)
{
  GtMD5ToSeqidsVisitor *md5_to_id_visitor =
    md5_to_id_visitor_cast(nv);
  gt_assert(md5_to_id_visitor);
  gt_region_mapping_delete(md5_to_id_visitor->region_mapping);
}

typedef struct {
  GtStr *new_seqid;
  GtRegionMapping *region_mapping;
} M2IChangeSeqidInfo;

static int m2i_change_target_seqids(GtFeatureNode *fn, const char *target,
                                    GtRegionMapping *region_mapping,
                                    GtError *err)
{
  GtStrArray *target_ids;
  GtArray *target_ranges, *target_strands;
  GtStr *desc, *new_seqid;
  unsigned long i;
  int had_err;
  gt_error_check(err);
  gt_assert(fn && target);
  target_ids = gt_str_array_new();
  target_ranges = gt_array_new(sizeof (GtRange));
  target_strands = gt_array_new(sizeof (GtStrand));
  desc = gt_str_new();
  new_seqid = gt_str_new();
  had_err = gt_gff3_parser_parse_all_target_attributes(target, false,
                                                       target_ids,
                                                       target_ranges,
                                                       target_strands, "", 0,
                                                       err);
  for (i = 0; !had_err && i < gt_str_array_size(target_ids); i++) {
    GtStr *seqid;
    gt_str_reset(desc);
    gt_str_reset(new_seqid);
    seqid = gt_str_array_get_str(target_ids, i);
    if (gt_str_length(seqid) > GT_MD5_SEQID_TOTAL_LEN)
      gt_str_append_cstr(new_seqid, gt_str_get(seqid) + GT_MD5_SEQID_TOTAL_LEN);
    else {
      /* for backward compatibility */
      if (!region_mapping) {
        gt_error_set(err, "no region mapping defined");
        had_err = -1;
      }
      if (!had_err) {
        had_err = gt_region_mapping_get_description(region_mapping, desc, seqid,
                                                    err);
      }
      if (!had_err)
        gt_regular_seqid_save(new_seqid, desc);
    }
    if (!had_err)
      gt_str_array_set(target_ids, i, new_seqid);
  }
  if (!had_err) {
    GtStr *new_target = gt_str_new();
    gt_gff3_parser_build_target_str(new_target, target_ids, target_ranges,
                                    target_strands);
    gt_feature_node_set_attribute(fn, GT_GFF_TARGET, gt_str_get(new_target));
    gt_str_delete(new_target);
  }
  gt_str_delete(new_seqid);
  gt_str_delete(desc);
  gt_array_delete(target_strands);
  gt_array_delete(target_ranges);
  gt_str_array_delete(target_ids);
  return had_err;
}

static int m2i_change_seqid(GtFeatureNode *fn, void *data, GtError *err)
{
  const char *target;
  M2IChangeSeqidInfo *info = (M2IChangeSeqidInfo*) data;
  gt_error_check(err);
  gt_assert(fn && info);
  gt_genome_node_change_seqid((GtGenomeNode*) fn, info->new_seqid);
  if ((target = gt_feature_node_get_attribute(fn, GT_GFF_TARGET)))
    return m2i_change_target_seqids(fn, target, info->region_mapping, err);
  return 0;
}

static int md5_to_seqid(GtGenomeNode *gn, GtRegionMapping *region_mapping,
                        GtError *err)
{
  GtStr *seqid;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(gn);
  seqid = gt_genome_node_get_seqid(gn);
  if (gt_md5_seqid_has_prefix(gt_str_get(seqid))) {
    GtStr *new_seqid = gt_str_new();
    /* seqid is a MD5 seqid -> change id */
    if (gt_str_length(seqid) > GT_MD5_SEQID_TOTAL_LEN)
      gt_str_append_cstr(new_seqid, gt_str_get(seqid) + GT_MD5_SEQID_TOTAL_LEN);
    else {
      /* for backward compatibility */
      GtStr *desc = gt_str_new();
      if (!region_mapping) {
        gt_error_set(err, "no region mapping defined");
        had_err = -1;
      }
      if (!had_err) {
        had_err = gt_region_mapping_get_description(region_mapping, desc, seqid,
                                                    err);
      }
      if (!had_err)
        gt_regular_seqid_save(new_seqid, desc);
      gt_str_delete(desc);
    }
    if (!had_err) {
      if (gt_feature_node_try_cast(gn)) {
        M2IChangeSeqidInfo info;
        info.new_seqid = new_seqid;
        info.region_mapping = region_mapping;
        had_err = gt_feature_node_traverse_children((GtFeatureNode*) gn, &info,
                                                    m2i_change_seqid, true,
                                                    err);
      }
      else
        gt_genome_node_change_seqid(gn, new_seqid);
    }
    gt_str_delete(new_seqid);
  }
  return had_err;
}

static int md5_to_id_visitor_feature_node(GtNodeVisitor *nv,
                                              GtFeatureNode *fn,
                                              GtError *err)
{
  GtMD5ToSeqidsVisitor *v = md5_to_id_visitor_cast(nv);
  gt_error_check(err);
  return md5_to_seqid((GtGenomeNode*) fn, v->region_mapping, err);
}

static int md5_to_id_visitor_region_node(GtNodeVisitor *nv,
                                             GtRegionNode *rn,
                                             GtError *err)
{
  GtMD5ToSeqidsVisitor *v = md5_to_id_visitor_cast(nv);
  gt_error_check(err);
  return md5_to_seqid((GtGenomeNode*) rn, v->region_mapping, err);
}

const GtNodeVisitorClass* gt_md5_to_id_visitor_class()
{
  static const GtNodeVisitorClass *nvc = NULL;
  gt_class_alloc_lock_enter();
  if (!nvc) {
    nvc = gt_node_visitor_class_new(sizeof (GtMD5ToSeqidsVisitor),
                                    md5_to_id_visitor_free,
                                    NULL,
                                    md5_to_id_visitor_feature_node,
                                    md5_to_id_visitor_region_node,
                                    NULL,
                                    NULL);
  }
  gt_class_alloc_lock_leave();
  return nvc;
}

GtNodeVisitor* gt_md5_to_id_visitor_new(GtRegionMapping *region_mapping)
{
  GtNodeVisitor *nv;
  GtMD5ToSeqidsVisitor *md5_to_id_visitor;
  nv = gt_node_visitor_create(gt_md5_to_id_visitor_class());
  md5_to_id_visitor = md5_to_id_visitor_cast(nv);
  md5_to_id_visitor->region_mapping = region_mapping;
  return nv;
}
