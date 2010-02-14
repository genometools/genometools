/*
  Copyright (c) 2006-2010 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#include "extended/feature_node_iterator_api.h"
#include "extended/genome_node.h"
#include "extended/region_mapping.h"
#include "extended/reverse.h"

static int extract_join_feature(GtGenomeNode *gn, const char *type,
                                GtRegionMapping *region_mapping,
                                GtStr *sequence, bool *reverse_strand,
                                GtError *err)
{
  const char *raw_sequence;
  unsigned long raw_sequence_length;
  GtFeatureNode *gf;
  GtRange range;
  int had_err = 0;

  gt_error_check(err);
  gf = gt_genome_node_cast(gt_feature_node_class(), gn);
  gt_assert(gf);

  if (gt_feature_node_has_type(gf, type)) {
    had_err = gt_region_mapping_get_raw_sequence(region_mapping, &raw_sequence,
                                                 &raw_sequence_length,
                                                 gt_genome_node_get_seqid(gn),
                                                 err);
    if (!had_err) {
      range = gt_genome_node_get_range(gn);
      gt_assert(range.start); /* 1-based coordinates */
      raw_sequence += range.start - 1;
      gt_assert(range.end <= raw_sequence_length);
      gt_str_append_cstr_nt(sequence, raw_sequence, gt_range_length(&range));
      if (gt_feature_node_get_strand(gf) == GT_STRAND_REVERSE)
        *reverse_strand = true;
    }
  }
  return had_err;
}

int gt_extract_feat_sequence(GtStr *sequence, GtGenomeNode *gn,
                             const char *type,
                             bool join, GtRegionMapping *region_mapping,
                             GtError *err)
{
  GtFeatureNode *gf;
  GtRange range;
  const char *raw_sequence;
  unsigned long raw_sequence_length;
  int had_err = 0;

  gt_error_check(err);
  gf = gt_genome_node_cast(gt_feature_node_class(), gn);
  gt_assert(gf);

  if (join) {
    GtFeatureNodeIterator *fni;
    GtFeatureNode *child;
    bool reverse_strand = false;
    /* in this case we have to traverse the children */
    fni = gt_feature_node_iterator_new_direct(gt_feature_node_cast(gn));
    while (!had_err && (child = gt_feature_node_iterator_next(fni))) {
      if (extract_join_feature((GtGenomeNode*) child, type, region_mapping,
                               sequence, &reverse_strand, err)) {
        had_err = -1;
      }
    }
    gt_feature_node_iterator_delete(fni);
    if (!had_err && gt_str_length(sequence)) {
      if (reverse_strand) {
        had_err = gt_reverse_complement(gt_str_get(sequence),
                                     gt_str_length(sequence), err);
      }
    }
  }
  else if (gt_feature_node_get_type(gf) == type) {
    gt_assert(!had_err);
    /* otherwise we only have to look this feature */
    range = gt_genome_node_get_range(gn);
    gt_assert(range.start); /* 1-based coordinates */
    had_err = gt_region_mapping_get_raw_sequence(region_mapping,
                                                 &raw_sequence,
                                                 &raw_sequence_length,
                                                 gt_genome_node_get_seqid(gn),
                                                 err);
    if (!had_err) {
      gt_assert(range.end <= raw_sequence_length);
      gt_str_append_cstr_nt(sequence, raw_sequence + range.start - 1,
                            gt_range_length(&range));
      if (gt_feature_node_get_strand(gf) == GT_STRAND_REVERSE) {
        had_err = gt_reverse_complement(gt_str_get(sequence),
                                     gt_str_length(sequence), err);
      }
    }
  }
  return had_err;
}
