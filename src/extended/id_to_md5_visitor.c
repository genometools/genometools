/*
  Copyright (c) 2010 Gordon Gremme <gremme@zbh.uni-hamburg.de>

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
#include "core/md5_seqid.h"
#include "core/undef.h"
#include "extended/gff3_parser.h"
#include "extended/node_visitor_rep.h"
#include "extended/id_to_md5_visitor.h"

struct GtSeqidsToMD5Visitor {
  const GtNodeVisitor parent_instance;
  GtRegionMapping *region_mapping;
  bool substitute_target_ids;
};

#define  id_to_md5_visitor_cast(GV)\
         gt_node_visitor_cast(gt_id_to_md5_visitor_class(), GV)

static void id_to_md5_visitor_free(GtNodeVisitor *nv)
{
  GtSeqidsToMD5Visitor *id_to_md5_visitor =
    id_to_md5_visitor_cast(nv);
  gt_assert(id_to_md5_visitor);
  gt_region_mapping_delete(id_to_md5_visitor->region_mapping);
}

typedef struct {
  GtStr *new_seqid;
  GtRegionMapping *region_mapping;
  bool substitute_target_ids;
  unsigned long offset;
} I2MChangeSeqidInfo;

static void i2m_build_new_target(GtStr *target, GtStrArray *target_ids,
                                 GtArray *target_ranges,
                                 GtArray *target_strands)
{
  unsigned long i;
  gt_assert(target && target_ids && target_ranges && target_strands);
  for (i = 0; i < gt_str_array_size(target_ids); i++) {
    GtRange *range;
    GtStrand *strand;
    range = gt_array_get(target_ranges, i);
    strand = gt_array_get(target_strands, i);
    if (i)
      gt_str_append_char(target, ',');
    gt_str_append_cstr(target, GT_MD5_SEQID_PREFIX);
    gt_str_append_cstr(target, gt_str_array_get(target_ids, i));
    gt_str_append_char(target, ' ');
    gt_str_append_ulong(target, range->start);
    gt_str_append_char(target, ' ');
    gt_str_append_ulong(target, range->end);
    if (*strand != GT_NUM_OF_STRAND_TYPES) {
      gt_str_append_char(target, ' ');
      gt_str_append_char(target, GT_STRAND_CHARS[*strand]);
    }
  }
}

static int i2m_change_target_seqids(GtFeatureNode *fn, const char *target,
                                    GtRegionMapping *region_mapping,
                                    GtError *err)
{
  GtStrArray *target_ids;
  GtArray *target_ranges, *target_strands;
  unsigned long i;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(fn && target && region_mapping);
  target_ids = gt_str_array_new();
  target_ranges = gt_array_new(sizeof (GtRange));
  target_strands = gt_array_new(sizeof (GtStrand));
  gt_gff3_parser_parse_all_target_attributes(target, target_ids, target_ranges,
                                             target_strands);
  for (i = 0; !had_err && i < gt_str_array_size(target_ids); i++) {
    GtStr *seqid;
    unsigned long offset;
    const char *md5;
    GtRange *range;
    seqid = gt_str_array_get_str(target_ids, i);
    range = gt_array_get(target_ranges, i);
    if (!(md5 = gt_region_mapping_get_md5_fingerprint(region_mapping, seqid,
                                                      range, &offset, err))) {
      had_err = -1;
    }
    if (!had_err) {
      GtRange transformed_range;
      gt_str_array_set_cstr(target_ids, i, md5);
      gt_assert(offset);
      transformed_range = gt_range_offset(range, -(offset - 1));
      range->start = transformed_range.start;
      range->end = transformed_range.end;
    }
  }
  if (!had_err) {
    GtStr *new_target = gt_str_new();
    i2m_build_new_target(new_target, target_ids, target_ranges, target_strands);
    gt_feature_node_set_attribute(fn, GT_GFF_TARGET, gt_str_get(new_target));
    gt_str_delete(new_target);
  }
  gt_array_delete(target_strands);
  gt_array_delete(target_ranges);
  gt_str_array_delete(target_ids);
  return had_err;
}

static int i2m_change_seqid(GtFeatureNode *fn, void *data, GtError *err)
{
  I2MChangeSeqidInfo *info = (I2MChangeSeqidInfo*) data;
  const char *target;
  gt_error_check(err);
  gt_assert(fn && info);
  gt_genome_node_change_seqid((GtGenomeNode*) fn, info->new_seqid);
  if (info->offset) {
    GtRange old_range, new_range;
    old_range = gt_genome_node_get_range((GtGenomeNode*) fn);
    new_range = gt_range_offset(&old_range, -info->offset);
    gt_genome_node_set_range((GtGenomeNode*) fn, &new_range);
  }
  if ((target = gt_feature_node_get_attribute(fn, GT_GFF_TARGET)) &&
      info->substitute_target_ids) {
    return i2m_change_target_seqids(fn, target, info->region_mapping, err);
  }
  return 0;
}

static int seqid_to_md5(GtGenomeNode *gn, GtRegionMapping *region_mapping,
                        bool substitute_target_ids, GtError *err)
{
  GtStr *seqid;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(gn && region_mapping);
  seqid = gt_genome_node_get_seqid(gn);
  if (!gt_md5_seqid_has_prefix(gt_str_get(seqid))) {
    /* seqid is not already a MD5 seqid -> change id */
    unsigned long offset;
    const char *md5;
    GtRange range = gt_genome_node_get_range(gn);
    if (!(md5 = gt_region_mapping_get_md5_fingerprint(region_mapping, seqid,
                                                      &range, &offset, err))) {
      had_err = -1;
    }
    if (!had_err) {
      GtStr *new_seqid = gt_str_new_cstr(GT_MD5_SEQID_PREFIX);
      gt_str_append_cstr(new_seqid, md5);
      if (gt_feature_node_try_cast(gn)) {
        I2MChangeSeqidInfo info;
        info.new_seqid = new_seqid;
        info.region_mapping = region_mapping;
        info.substitute_target_ids = substitute_target_ids;
        gt_assert(offset);
        info.offset = offset - 1;
        had_err = gt_feature_node_traverse_children((GtFeatureNode*) gn, &info,
                                                    i2m_change_seqid, true,
                                                    err);
      }
      else
        gt_genome_node_change_seqid(gn, new_seqid);
      gt_str_delete(new_seqid);
    }
  }
  return had_err;
}

static int id_to_md5_visitor_feature_node(GtNodeVisitor *nv,
                                              GtFeatureNode *fn,
                                              GtError *err)
{
  GtSeqidsToMD5Visitor *v = id_to_md5_visitor_cast(nv);
  gt_error_check(err);
  return seqid_to_md5((GtGenomeNode*) fn, v->region_mapping,
                      v->substitute_target_ids, err);
}

static int id_to_md5_visitor_region_node(GtNodeVisitor *nv,
                                             GtRegionNode *rn,
                                             GtError *err)
{
  GtSeqidsToMD5Visitor *v = id_to_md5_visitor_cast(nv);
  gt_error_check(err);
  return seqid_to_md5((GtGenomeNode*) rn, v->region_mapping,
                      v->substitute_target_ids, err);
}

const GtNodeVisitorClass* gt_id_to_md5_visitor_class()
{
  static const GtNodeVisitorClass *nvc = NULL;
  if (!nvc) {
    nvc = gt_node_visitor_class_new(sizeof (GtSeqidsToMD5Visitor),
                                    id_to_md5_visitor_free,
                                    NULL,
                                    id_to_md5_visitor_feature_node,
                                    id_to_md5_visitor_region_node,
                                    NULL,
                                    NULL);
  }
  return nvc;
}

GtNodeVisitor* gt_id_to_md5_visitor_new(GtRegionMapping *region_mapping,
                                        bool substitute_target_ids)
{
  GtNodeVisitor *nv;
  GtSeqidsToMD5Visitor *id_to_md5_visitor;
  nv = gt_node_visitor_create(gt_id_to_md5_visitor_class());
  id_to_md5_visitor = id_to_md5_visitor_cast(nv);
  id_to_md5_visitor->region_mapping = region_mapping;
  id_to_md5_visitor->substitute_target_ids = substitute_target_ids;
  return nv;
}
