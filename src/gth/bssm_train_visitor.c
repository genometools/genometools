/*
  Copyright (c) 2010-2013 Gordon Gremme <gordon@gremme.org>

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

#include "core/alphabet.h"
#include "core/cstr_api.h"
#include "core/ma_api.h"
#include "core/warning_api.h"
#include "core/xansi_api.h"
#include "extended/feature_node_iterator_api.h"
#include "extended/feature_type.h"
#include "extended/node_visitor_api.h"
#include "extended/reverse_api.h"
#include "gth/bssm_helper.h"
#include "gth/bssm_seq_processor.h"
#include "gth/bssm_train_visitor.h"

struct GthBSSMTrainVisitor {
  const GtNodeVisitor parent_instance;
  GtRegionMapping *region_mapping;
  char  *filter_type,
        *extract_type;
  GthBSSMSeqProcessor *bsp;
  unsigned int good_exon_count;
  double cutoff;
};

#define bssm_train_visitor_cast(GV)\
        gt_node_visitor_cast(gth_bssm_train_visitor_class(), GV)

static void bssm_train_visitor_free(GtNodeVisitor *nv)
{
  GthBSSMTrainVisitor *bssm_train_visitor = bssm_train_visitor_cast(nv);
  gth_bssm_seq_processor_delete(bssm_train_visitor->bsp);
  gt_free(bssm_train_visitor->extract_type);
  gt_free(bssm_train_visitor->filter_type);
  gt_region_mapping_delete(bssm_train_visitor->region_mapping);
}

static int get_seq(GtStr *seq, bool *process, const GtRange *range,
                   GtStr *seqid, bool reverse, GtRegionMapping *region_mapping,
                   GtError *err)
{
  GtUword sequence_length;
  char *sequence = NULL;
  int had_err;
  gt_error_check(err);
  gt_assert(seq && process && range && seqid && region_mapping);
  had_err = gt_region_mapping_get_sequence_length(region_mapping,
                                                  &sequence_length, seqid, err);
  if (!had_err) {
    had_err = gt_region_mapping_get_sequence(region_mapping, &sequence, seqid,
                                             range->start, range->end, err);
  }
  if (!had_err) {
    gt_assert(range->start && range->end); /* 1-based coordinates */
    gt_assert(range->end <= sequence_length);
    gt_str_reset(seq);
    gt_str_append_cstr_nt(seq, sequence, gt_range_length(range));
    if (!gth_seq_contains_wildcard(seq)) {
    /* no ambiguous bases found -> process */
      if (reverse) {
        had_err = gt_reverse_complement(gt_str_get(seq), gt_str_length(seq),
                                        err);
      }
      *process = true;
    }
    else
      *process = false;
  }
  gt_free(sequence);
  return had_err;
}

static int process_ranges(GtArray *ranges, GtStr *seqid, bool reverse,
                          GthBSSMTrainVisitor *visitor, GtError *err)
{
  GtUword i, phase = 0;
  GtRange range;
  GtStr *seq;
  bool process;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(ranges && seqid && visitor);
  seq = gt_str_new();
  if (reverse)
    gt_array_reverse(ranges);
  for (i = 0; !had_err && i < gt_array_size(ranges); i++) {
    if (i) {
      /* get intron */
      if (reverse) {
        range.start = ((GtRange*) gt_array_get(ranges, i))->end + 1;
        range.end   = ((GtRange*) gt_array_get(ranges, i-1))->start - 1;
      }
      else{
        range.start = ((GtRange*) gt_array_get(ranges, i-1))->end + 1;
        range.end   = ((GtRange*) gt_array_get(ranges, i))->start - 1;
      }
      gt_assert(range.start <= range.end);
      if (gt_range_length(&range) >= 2) {
        had_err = get_seq(seq, &process, &range, seqid, reverse,
                          visitor->region_mapping, err);
        if (!had_err && process) {
          /* no ambiguous bases found -> process intron */
          gth_bssm_seq_processor_proc_intron(visitor->bsp, phase, seqid, &range,
                                             reverse, seq);
        }
      }
      else {
        gt_warning("ignoring intron of length < 2 for sequence ID '%s'",
                   gt_str_get(seqid));
      }
    }
    /* get exon */
    if (!had_err) {
      range = *(GtRange*) gt_array_get(ranges, i);
      had_err = get_seq(seq, &process, &range, seqid, reverse,
                        visitor->region_mapping, err);
    }
    if (!had_err && process) {
      /* no ambiguous bases found -> process exon */
      gth_bssm_seq_processor_proc_exon(visitor->bsp, phase, seqid, &range,
                                       reverse, seq);
    }
    /* update phase */
    if (!had_err)
      phase = (phase + gt_range_length(&range)) % 3;
  }
  gt_str_delete(seq);
  return had_err;
}

static int traverse_direct(GtFeatureNode *fn, GthBSSMTrainVisitor *visitor,
                           GtError *err)
{
  GtFeatureNode *node, *first_node = NULL;
  GtFeatureNodeIterator *fni;
  bool type_found = false;
  GtStrand found_strand = GT_NUM_OF_STRAND_TYPES;
  unsigned int count = 0;
  GtArray *ranges;
  GtRange range;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(fn && visitor);
  ranges = gt_array_new(sizeof (GtRange));
  fni = gt_feature_node_iterator_new_direct(fn);
  while ((node = gt_feature_node_iterator_next(fni))) {
    if (gt_feature_node_has_type(node, visitor->filter_type)) {
      if (!gt_feature_node_score_is_defined(node) ||
          gt_feature_node_get_score(node) >= visitor->cutoff) {
        count++;
      }
    }
    if (gt_feature_node_has_type(node, visitor->extract_type)) {
      if (!type_found) {
        type_found = true;
        found_strand = gt_feature_node_get_strand(node);
        if (!(found_strand == GT_STRAND_FORWARD ||
              found_strand == GT_STRAND_REVERSE)) {
          gt_error_set(err, "strand (%c) of feature on line %u in file \"%s\" "
                       "is neither forward (%c) nor reverse (%c)",
                       GT_STRAND_CHARS[found_strand],
                       gt_genome_node_get_line_number((GtGenomeNode*) node),
                       gt_genome_node_get_filename((GtGenomeNode*) node),
                       GT_STRAND_CHARS[GT_STRAND_FORWARD],
                       GT_STRAND_CHARS[GT_STRAND_REVERSE]);
          had_err = -1;
          break;
        }
        first_node = node;
      }
      else if (gt_feature_node_get_strand(node) != found_strand) {
        gt_error_set(err, "strand (%c) of feature on line %u in file \"%s\" is "
                     "different from strand (%c) of feature on line %u in file "
                     "\"%s\"",
                     GT_STRAND_CHARS[gt_feature_node_get_strand(node)],
                     gt_genome_node_get_line_number((GtGenomeNode*) node),
                     gt_genome_node_get_filename((GtGenomeNode*) node),
                     GT_STRAND_CHARS[found_strand],
                     gt_genome_node_get_line_number((GtGenomeNode*) first_node),
                     gt_genome_node_get_filename((GtGenomeNode*) first_node));
        had_err = -1;
        break;
      }
      /* everything is fine with feature of desired type -> store range */
      range = gt_genome_node_get_range((GtGenomeNode*) node);
      gt_array_add(ranges, range);
    }
  }
  gt_feature_node_iterator_delete(fni);
  if (!had_err && gt_array_size(ranges) && count >= visitor->good_exon_count) {
    gt_assert(type_found && first_node);
    gt_assert(found_strand == GT_STRAND_FORWARD ||
              found_strand == GT_STRAND_REVERSE);
    had_err = process_ranges(ranges, gt_genome_node_get_seqid((GtGenomeNode*)
                                                              first_node),
                             found_strand == GT_STRAND_REVERSE, visitor, err);
  }
  gt_array_delete(ranges);
  return had_err;
}

static int bssm_train_visitor_feature_node(GtNodeVisitor *nv, GtFeatureNode *fn,
                                           GtError *err)
{
  GthBSSMTrainVisitor *v = bssm_train_visitor_cast(nv);
  GtFeatureNodeIterator *fni;
  GtFeatureNode *node;
  int had_err = 0;
  gt_error_check(err);

  if (!had_err) {
    fni = gt_feature_node_iterator_new(fn);
    while (!had_err && (node = gt_feature_node_iterator_next(fni)))
      had_err = traverse_direct(node, v, err);
    gt_feature_node_iterator_delete(fni);
  }

  return had_err;
}

const GtNodeVisitorClass* gth_bssm_train_visitor_class()
{
  static const GtNodeVisitorClass *nvc = NULL;
  if (!nvc) {
    nvc = gt_node_visitor_class_new(sizeof (GthBSSMTrainVisitor),
                                    bssm_train_visitor_free,
                                    NULL,
                                    bssm_train_visitor_feature_node,
                                    NULL,
                                    NULL,
                                    NULL);
  }
  return nvc;
}

GtNodeVisitor* gth_bssm_train_visitor_new(GtRegionMapping *region_mapping,
                                          GthBSSMSeqProcessor *bsp,
                                          const char *filter_type,
                                          const char *extract_type,
                                          unsigned int good_exon_count,
                                          double cutoff)
{
  GtNodeVisitor *nv;
  GthBSSMTrainVisitor *btv;
  gt_assert(region_mapping);
  nv = gt_node_visitor_create(gth_bssm_train_visitor_class());
  btv = bssm_train_visitor_cast(nv);
  btv->region_mapping = region_mapping;
  btv->filter_type = gt_cstr_dup(filter_type);
  btv->extract_type = gt_cstr_dup(extract_type);
  btv->bsp = bsp;
  btv->good_exon_count = good_exon_count;
  btv->cutoff = cutoff;
  return nv;
}
