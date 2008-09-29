/*
  Copyright (c) 2007-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007-2008 Center for Bioinformatics, University of Hamburg

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
#include "core/minmax.h"
#include "core/queue.h"
#include "core/undef.h"
#include "core/unused_api.h"
#include "extended/filter_visitor.h"
#include "extended/gff3_parser.h"
#include "extended/node_visitor_rep.h"

struct GtFilterVisitor {
  const GtNodeVisitor parent_instance;
  GtQueue *gt_genome_node_buffer;
  GtStr *seqid,
      *typefilter;
  GtRange contain_range,
        overlap_range;
  GtStrand strand,
         targetstrand;
  bool has_CDS;
  unsigned long max_gene_length,
                gene_num,     /* the number of passed genes */
                max_gene_num, /* the maximal number of genes which can pass */
                current_feature,
                feature_num;
  double min_gene_score,
         max_gene_score,
         min_average_splice_site_prob;
};

#define filter_visitor_cast(GV)\
        gt_node_visitor_cast(gt_filter_visitor_class(), GV)

static void filter_visitor_free(GtNodeVisitor *gv)
{
  GtFilterVisitor *filter_visitor = filter_visitor_cast(gv);
  gt_queue_delete(filter_visitor->gt_genome_node_buffer);
  gt_str_delete(filter_visitor->seqid);
  gt_str_delete(filter_visitor->typefilter);
}

static int filter_visitor_comment(GtNodeVisitor *gv, GtCommentNode *c,
                                  GT_UNUSED GtError *err)
{
  GtFilterVisitor *filter_visitor;
  gt_error_check(err);
  filter_visitor = filter_visitor_cast(gv);
  gt_queue_add(filter_visitor->gt_genome_node_buffer, c);
  return 0;
}

static bool filter_contain_range(GtFeatureNode *gf, GtRange contain_range)
{
  GtRange range;
  gt_assert(gf);
  range = gt_genome_node_get_range((GtGenomeNode*) gf);
  if (contain_range.start != UNDEF_ULONG &&
      !gt_range_contains(&contain_range, &range)) {
    return true;
  }
  return false;
}

static bool filter_overlap_range(GtFeatureNode *gf, GtRange overlap_range)
{
  GtRange feature_range;
  gt_assert(gf);
  feature_range = gt_genome_node_get_range((GtGenomeNode*) gf);
  if (overlap_range.start != UNDEF_ULONG &&
      !gt_range_overlap(&overlap_range, &feature_range))
    return true;
  return false;
}

static bool filter_strand(GtFeatureNode *gf, GtStrand strand)
{
  gt_assert(gf);
  if (strand != GT_NUM_OF_STRAND_TYPES &&
      gt_feature_node_get_strand(gf) != strand)
    return true;
  return false;
}

static bool filter_targetstrand(GtFeatureNode *fn, GtStrand targetstrand)
{
  const char *target;
  gt_assert(fn);
  if (targetstrand != GT_NUM_OF_STRAND_TYPES &&
      (target = gt_feature_node_get_attribute(fn, TARGET_STRING))) {
    unsigned long num_of_targets;
    GtStrand parsed_strand;
    int had_err;
    had_err = gt_gff3_parser_parse_target_attributes(target, &num_of_targets,
                                                     NULL, NULL, &parsed_strand,
                                                     "", 0, NULL);
    gt_assert(!had_err);
    if (num_of_targets == 1 && parsed_strand != GT_NUM_OF_STRAND_TYPES &&
        parsed_strand != targetstrand) {
      return true;
    }
  }
  return false;
}

static bool filter_has_CDS(GtFeatureNode *gf, bool has_CDS)
{
  gt_assert(gf);
  if (has_CDS && !gt_feature_node_has_CDS(gf))
    return true;
  return false;
}

static bool filter_min_average_ssp(GtFeatureNode *gf, double minaveragessp)
{
  gt_assert(gf);
  if (minaveragessp != UNDEF_DOUBLE &&
      gt_feature_node_has_splice_site(gf) &&
      gt_feature_node_average_splice_site_prob(gf) < minaveragessp) {
    return true;
  }
  return false;
}

static int filter_visitor_genome_feature(GtNodeVisitor *gv,
                                         GtFeatureNode *gf,
                                         GT_UNUSED GtError *err)
{
  GtFilterVisitor *fv;
  bool filter_node = false;
  gt_error_check(err);
  fv = filter_visitor_cast(gv);
  fv->current_feature++;
  if (!gt_str_length(fv->seqid) || /* no seqid was specified or seqids are
                                      equal */
      !gt_str_cmp(fv->seqid, gt_genome_node_get_seqid((GtGenomeNode*) gf))) {
    GtRange range = gt_genome_node_get_range((GtGenomeNode*) gf);
    /* enforce maximum gene length */
    /* XXX: we (spuriously) assume that genes are always root nodes */
    if (gf && gt_feature_node_has_type(gf, gft_gene)) {
      if (fv->max_gene_length != UNDEF_ULONG &&
          gt_range_length(&range) > fv->max_gene_length) {
        filter_node = true;
      }
      else if (fv->max_gene_num != UNDEF_ULONG &&
               fv->gene_num >= fv->max_gene_num) {
        filter_node = true;
      }
      else if (fv->min_gene_score != UNDEF_DOUBLE &&
               gt_feature_node_get_score(gf) < fv->min_gene_score) {
        filter_node = true;
      }
      else if (fv->max_gene_score != UNDEF_DOUBLE &&
               gt_feature_node_get_score(gf) > fv->max_gene_score) {
        filter_node = true;
      }
      else if (fv->feature_num != UNDEF_ULONG &&
               fv->feature_num != fv->current_feature) {
        filter_node = true;
      }
      if (!filter_node)
        fv->gene_num++; /* gene passed filter */
    }
  }
  else
    filter_node = true;

  if (!filter_node)
    filter_node = filter_contain_range(gf, fv->contain_range);

  if (!filter_node)
    filter_node = filter_overlap_range(gf, fv->overlap_range);

  if (!filter_node)
    filter_node = filter_strand(gf, fv->strand);

  if (!filter_node)
    filter_node = filter_targetstrand(gf, fv->targetstrand);

  if (!filter_node)
    filter_node = filter_has_CDS(gf, fv->has_CDS);

  if (!filter_node)
    filter_node = filter_min_average_ssp(gf, fv->min_average_splice_site_prob);

  if (filter_node)
    gt_genome_node_rec_delete((GtGenomeNode*) gf);
  else
    gt_queue_add(fv->gt_genome_node_buffer, gf);

  return 0;
}

static int filter_visitor_region_node(GtNodeVisitor *gv, GtRegionNode *rn,
                                      GT_UNUSED GtError *err)
{
  GtFilterVisitor *filter_visitor;
  gt_error_check(err);
  filter_visitor = filter_visitor_cast(gv);
  if (!gt_str_length(filter_visitor->seqid) || /* no seqid was specified */
      !gt_str_cmp(filter_visitor->seqid,       /* or seqids are equal */
               gt_genome_node_get_seqid((GtGenomeNode*) rn))) {
    if (filter_visitor->contain_range.start != UNDEF_ULONG) {
      GtRange range = gt_genome_node_get_range((GtGenomeNode*) rn);
      if (gt_range_overlap(&range, &filter_visitor->contain_range)) {
        /* an overlapping contain range was defined -> update range  */
        range.start = MAX(range.start, filter_visitor->contain_range.start);
        range.end = MIN(range.end, filter_visitor->contain_range.end);
        gt_genome_node_set_range((GtGenomeNode*) rn, range);
        gt_queue_add(filter_visitor->gt_genome_node_buffer, rn);
      }
      else /* contain range does not overlap with <rn> range -> delete <rn> */
        gt_genome_node_delete((GtGenomeNode*) rn);
    }
    else
      gt_queue_add(filter_visitor->gt_genome_node_buffer, rn);
  }
  else
    gt_genome_node_rec_delete((GtGenomeNode*) rn);
  return 0;
}

static int filter_visitor_sequence_node(GtNodeVisitor *gv, GtSequenceNode *sn,
                                        GT_UNUSED GtError *err)
{
  GtFilterVisitor *filter_visitor;
  gt_error_check(err);
  filter_visitor = filter_visitor_cast(gv);
  if (!gt_str_length(filter_visitor->seqid) || /* no seqid was specified */
      !gt_str_cmp(filter_visitor->seqid,       /* or seqids are equal */
               gt_genome_node_get_seqid((GtGenomeNode*) sn))) {
    gt_queue_add(filter_visitor->gt_genome_node_buffer, sn);
  }
  else
    gt_genome_node_rec_delete((GtGenomeNode*) sn);
  return 0;
}

const GtNodeVisitorClass* gt_filter_visitor_class()
{
  static const GtNodeVisitorClass *gvc = NULL;
  if (!gvc) {
    gvc = gt_node_visitor_class_new(sizeof (GtFilterVisitor),
                                    filter_visitor_free,
                                    filter_visitor_comment,
                                    filter_visitor_genome_feature,
                                    filter_visitor_region_node,
                                    filter_visitor_sequence_node);
  }
  return gvc;
}

GtNodeVisitor* gt_filter_visitor_new(GtStr *seqid, GtStr *typefilter,
                                     GtRange contain_range,
                                     GtRange overlap_range,
                                     GtStrand strand, GtStrand targetstrand,
                                     bool has_CDS,
                                     unsigned long max_gene_length,
                                     unsigned long max_gene_num,
                                     double min_gene_score,
                                     double max_gene_score,
                                     double min_average_splice_site_prob,
                                     unsigned long feature_num)
{
  GtNodeVisitor *gv = gt_node_visitor_create(gt_filter_visitor_class());
  GtFilterVisitor *filter_visitor = filter_visitor_cast(gv);
  filter_visitor->gt_genome_node_buffer = gt_queue_new();
  filter_visitor->seqid = gt_str_ref(seqid);
  filter_visitor->typefilter = gt_str_ref(typefilter);
  filter_visitor->contain_range = contain_range;
  filter_visitor->overlap_range = overlap_range;
  filter_visitor->strand = strand;
  filter_visitor->targetstrand = targetstrand;
  filter_visitor->has_CDS = has_CDS;
  filter_visitor->max_gene_length = max_gene_length;
  filter_visitor->gene_num = 0;
  filter_visitor->max_gene_num = max_gene_num;
  filter_visitor->min_gene_score = min_gene_score;
  filter_visitor->max_gene_score = max_gene_score;
  filter_visitor->min_average_splice_site_prob = min_average_splice_site_prob;
  filter_visitor->feature_num = feature_num;
  return gv;
}

unsigned long gt_filter_visitor_node_buffer_size(GtNodeVisitor *gv)
{
  GtFilterVisitor *filter_visitor = filter_visitor_cast(gv);
  return gt_queue_size(filter_visitor->gt_genome_node_buffer);
}

GtGenomeNode* gt_filter_visitor_get_node(GtNodeVisitor *gv)
{
  GtFilterVisitor *filter_visitor;
  filter_visitor = filter_visitor_cast(gv);
  return gt_queue_get(filter_visitor->gt_genome_node_buffer);
}
