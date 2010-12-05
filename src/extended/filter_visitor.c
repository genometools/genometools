/*
  Copyright (c) 2007-2010 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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
  GtQueue *node_buffer;
  GtStr *seqid,
        *source,
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

static void filter_visitor_free(GtNodeVisitor *nv)
{
  GtFilterVisitor *filter_visitor = filter_visitor_cast(nv);
  gt_str_delete(filter_visitor->typefilter);
  gt_str_delete(filter_visitor->source);
  gt_str_delete(filter_visitor->seqid);
  gt_queue_delete(filter_visitor->node_buffer);
}

static int filter_visitor_comment_node(GtNodeVisitor *nv, GtCommentNode *c,
                                       GT_UNUSED GtError *err)
{
  GtFilterVisitor *filter_visitor;
  gt_error_check(err);
  filter_visitor = filter_visitor_cast(nv);
  gt_queue_add(filter_visitor->node_buffer, c);
  return 0;
}

static bool filter_contain_range(GtFeatureNode *fn, GtRange contain_range)
{
  GtRange range;
  gt_assert(fn);
  range = gt_genome_node_get_range((GtGenomeNode*) fn);
  if (contain_range.start != GT_UNDEF_ULONG &&
      !gt_range_contains(&contain_range, &range)) {
    return true;
  }
  return false;
}

static bool filter_overlap_range(GtFeatureNode *fn, GtRange overlap_range)
{
  GtRange feature_range;
  gt_assert(fn);
  feature_range = gt_genome_node_get_range((GtGenomeNode*) fn);
  if (overlap_range.start != GT_UNDEF_ULONG &&
      !gt_range_overlap(&overlap_range, &feature_range))
    return true;
  return false;
}

static bool filter_strand(GtFeatureNode *fn, GtStrand strand)
{
  gt_assert(fn);
  if (strand != GT_NUM_OF_STRAND_TYPES &&
      gt_feature_node_get_strand(fn) != strand)
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

static bool filter_has_CDS(GtFeatureNode *fn, bool has_CDS)
{
  gt_assert(fn);
  if (has_CDS && !gt_feature_node_has_CDS(fn))
    return true;
  return false;
}

static bool filter_min_average_ssp(GtFeatureNode *fn, double minaveragessp)
{
  gt_assert(fn);
  if (minaveragessp != GT_UNDEF_DOUBLE &&
      gt_feature_node_has_splice_site(fn) &&
      gt_feature_node_average_splice_site_prob(fn) < minaveragessp) {
    return true;
  }
  return false;
}

static int filter_visitor_feature_node(GtNodeVisitor *nv,
                                       GtFeatureNode *fn,
                                       GT_UNUSED GtError *err)
{
  GtFilterVisitor *fv;
  bool filter_node = false;
  gt_error_check(err);
  fv = filter_visitor_cast(nv);
  fv->current_feature++;
  if ((!gt_str_length(fv->seqid) || /* no seqid was specified or seqids are
                                       equal */
       !gt_str_cmp(fv->seqid, gt_genome_node_get_seqid((GtGenomeNode*) fn))) &&
      (!gt_str_length(fv->source) || /* no source was specified or sources are
                                        equal */
       !strcmp(gt_str_get(fv->source), gt_feature_node_get_source(fn)))) {
    GtRange range = gt_genome_node_get_range((GtGenomeNode*) fn);
    /* enforce maximum gene length */
    /* XXX: we (spuriously) assume that genes are always root nodes */
    if (fn && gt_feature_node_has_type(fn, gt_ft_gene)) {
      if (fv->max_gene_length != GT_UNDEF_ULONG &&
          gt_range_length(&range) > fv->max_gene_length) {
        filter_node = true;
      }
      else if (fv->max_gene_num != GT_UNDEF_ULONG &&
               fv->gene_num >= fv->max_gene_num) {
        filter_node = true;
      }
      else if (fv->min_gene_score != GT_UNDEF_DOUBLE &&
               gt_feature_node_get_score(fn) < fv->min_gene_score) {
        filter_node = true;
      }
      else if (fv->max_gene_score != GT_UNDEF_DOUBLE &&
               gt_feature_node_get_score(fn) > fv->max_gene_score) {
        filter_node = true;
      }
      else if (fv->feature_num != GT_UNDEF_ULONG &&
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
    filter_node = filter_contain_range(fn, fv->contain_range);

  if (!filter_node)
    filter_node = filter_overlap_range(fn, fv->overlap_range);

  if (!filter_node)
    filter_node = filter_strand(fn, fv->strand);

  if (!filter_node)
    filter_node = filter_targetstrand(fn, fv->targetstrand);

  if (!filter_node)
    filter_node = filter_has_CDS(fn, fv->has_CDS);

  if (!filter_node)
    filter_node = filter_min_average_ssp(fn, fv->min_average_splice_site_prob);

  if (filter_node)
    gt_genome_node_delete((GtGenomeNode*) fn);
  else
    gt_queue_add(fv->node_buffer, fn);

  return 0;
}

static int filter_visitor_region_node(GtNodeVisitor *nv, GtRegionNode *rn,
                                      GT_UNUSED GtError *err)
{
  GtFilterVisitor *filter_visitor;
  gt_error_check(err);
  filter_visitor = filter_visitor_cast(nv);
  if (!gt_str_length(filter_visitor->seqid) || /* no seqid was specified */
      !gt_str_cmp(filter_visitor->seqid,       /* or seqids are equal */
               gt_genome_node_get_seqid((GtGenomeNode*) rn))) {
    if (filter_visitor->contain_range.start != GT_UNDEF_ULONG) {
      GtRange range = gt_genome_node_get_range((GtGenomeNode*) rn);
      if (gt_range_overlap(&range, &filter_visitor->contain_range)) {
        /* an overlapping contain range was defined -> update range  */
        range.start = MAX(range.start, filter_visitor->contain_range.start);
        range.end = MIN(range.end, filter_visitor->contain_range.end);
        gt_genome_node_set_range((GtGenomeNode*) rn, &range);
        gt_queue_add(filter_visitor->node_buffer, rn);
      }
      else /* contain range does not overlap with <rn> range -> delete <rn> */
        gt_genome_node_delete((GtGenomeNode*) rn);
    }
    else
      gt_queue_add(filter_visitor->node_buffer, rn);
  }
  else
    gt_genome_node_delete((GtGenomeNode*) rn);
  return 0;
}

static int filter_visitor_sequence_node(GtNodeVisitor *nv, GtSequenceNode *sn,
                                        GT_UNUSED GtError *err)
{
  GtFilterVisitor *filter_visitor;
  gt_error_check(err);
  filter_visitor = filter_visitor_cast(nv);
  if (!gt_str_length(filter_visitor->seqid) || /* no seqid was specified */
      !gt_str_cmp(filter_visitor->seqid,       /* or seqids are equal */
                  gt_genome_node_get_seqid((GtGenomeNode*) sn))) {
    gt_queue_add(filter_visitor->node_buffer, sn);
  }
  else
    gt_genome_node_delete((GtGenomeNode*) sn);
  return 0;
}

static int filter_visitor_eof_node(GtNodeVisitor *nv, GtEOFNode *eofn,
                                   GT_UNUSED GtError *err)
{
  GtFilterVisitor *filter_visitor;
  gt_error_check(err);
  filter_visitor = filter_visitor_cast(nv);
  gt_queue_add(filter_visitor->node_buffer, eofn);
  return 0;
}

const GtNodeVisitorClass* gt_filter_visitor_class()
{
  static const GtNodeVisitorClass *nvc = NULL;
  if (!nvc) {
    nvc = gt_node_visitor_class_new(sizeof (GtFilterVisitor),
                                    filter_visitor_free,
                                    filter_visitor_comment_node,
                                    filter_visitor_feature_node,
                                    filter_visitor_region_node,
                                    filter_visitor_sequence_node,
                                    filter_visitor_eof_node);
  }
  return nvc;
}

GtNodeVisitor* gt_filter_visitor_new(GtStr *seqid,
                                     GtStr *source,
                                     GtStr *typefilter,
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
  GtNodeVisitor *nv = gt_node_visitor_create(gt_filter_visitor_class());
  GtFilterVisitor *filter_visitor = filter_visitor_cast(nv);
  filter_visitor->node_buffer = gt_queue_new();
  filter_visitor->seqid = gt_str_ref(seqid);
  filter_visitor->source = gt_str_ref(source);
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
  return nv;
}

unsigned long gt_filter_visitor_node_buffer_size(GtNodeVisitor *nv)
{
  GtFilterVisitor *filter_visitor = filter_visitor_cast(nv);
  return gt_queue_size(filter_visitor->node_buffer);
}

GtGenomeNode* gt_filter_visitor_get_node(GtNodeVisitor *nv)
{
  GtFilterVisitor *filter_visitor = filter_visitor_cast(nv);
  return gt_queue_get(filter_visitor->node_buffer);
}
