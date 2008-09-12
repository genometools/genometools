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
#include "core/minmax.h"
#include "core/queue.h"
#include "core/undef.h"
#include "core/unused_api.h"
#include "extended/filter_visitor.h"
#include "extended/genome_visitor_rep.h"
#include "extended/gff3_parser.h"

struct FilterVisitor {
  const GenomeVisitor parent_instance;
  GT_Queue *gt_genome_node_buffer;
  GtStr *seqid,
      *typefilter;
  GT_Range contain_range,
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
        genome_visitor_cast(filter_visitor_class(), GV)

static void filter_visitor_free(GenomeVisitor *gv)
{
  FilterVisitor *filter_visitor = filter_visitor_cast(gv);
  gt_queue_delete(filter_visitor->gt_genome_node_buffer);
  gt_str_delete(filter_visitor->seqid);
  gt_str_delete(filter_visitor->typefilter);
}

static int filter_visitor_comment(GenomeVisitor *gv, GT_Comment *c,
                                  GT_UNUSED GT_Error *err)
{
  FilterVisitor *filter_visitor;
  gt_error_check(err);
  filter_visitor = filter_visitor_cast(gv);
  gt_queue_add(filter_visitor->gt_genome_node_buffer, c);
  return 0;
}

static bool filter_contain_range(GT_GenomeFeature *gf, GT_Range contain_range)
{
  assert(gf);
  if (contain_range.start != UNDEF_ULONG &&
      !gt_range_contains(contain_range,
                         gt_genome_node_get_range((GT_GenomeNode*) gf))) {
    return true;
  }
  return false;
}

static bool filter_overlap_range(GT_GenomeFeature *gf, GT_Range overlap_range)
{
  assert(gf);
  if (overlap_range.start != UNDEF_ULONG &&
      !gt_range_overlap(overlap_range,
                        gt_genome_node_get_range((GT_GenomeNode*) gf))) {
    return true;
  }
  return false;
}

static bool filter_strand(GT_GenomeFeature *gf, GtStrand strand)
{
  assert(gf);
  if (strand != GT_NUM_OF_STRAND_TYPES &&
      gt_genome_feature_get_strand(gf) != strand)
    return true;
  return false;
}

static bool filter_targetstrand(GT_GenomeFeature *gf, GtStrand targetstrand)
{
  const char *target;
  assert(gf);
  if (targetstrand != GT_NUM_OF_STRAND_TYPES &&
      (target = gt_genome_feature_get_attribute((GT_GenomeNode*) gf,
                                             TARGET_STRING))) {
    unsigned long num_of_targets;
    GtStrand parsed_strand;
    int had_err;
    had_err = gt_gff3_parser_parse_target_attributes(target, &num_of_targets,
                                                     NULL, NULL, &parsed_strand,
                                                     "", 0, NULL);
    assert(!had_err);
    if (num_of_targets == 1 && parsed_strand != GT_NUM_OF_STRAND_TYPES &&
        parsed_strand != targetstrand) {
      return true;
    }
  }
  return false;
}

static bool filter_has_CDS(GT_GenomeFeature *gf, bool has_CDS)
{
  assert(gf);
  if (has_CDS && !gt_genome_feature_has_CDS(gf))
    return true;
  return false;
}

static bool filter_min_average_ssp(GT_GenomeFeature *gf, double minaveragessp)
{
  assert(gf);
  if (minaveragessp != UNDEF_DOUBLE &&
      gt_genome_feature_has_splice_site(gf) &&
      gt_genome_feature_average_splice_site_prob(gf) < minaveragessp) {
    return true;
  }
  return false;
}

static int filter_visitor_genome_feature(GenomeVisitor *gv,
                                         GT_GenomeFeature *gf,
                                         GT_UNUSED GT_Error *err)
{
  FilterVisitor *fv;
  bool filter_node = false;
  gt_error_check(err);
  fv = filter_visitor_cast(gv);
  fv->current_feature++;
  if (!gt_str_length(fv->seqid) || /* no seqid was specified or seqids are
                                      equal */
      !gt_str_cmp(fv->seqid, gt_genome_node_get_seqid((GT_GenomeNode*) gf))) {
    /* enforce maximum gene length */
    /* XXX: we (spuriously) assume that genes are always root nodes */
    if (gf && gt_genome_feature_has_type(gf, gft_gene)) {
      if (fv->max_gene_length != UNDEF_ULONG &&
          gt_range_length(gt_genome_node_get_range((GT_GenomeNode*) gf)) >
          fv->max_gene_length) {
        filter_node = true;
      }
      else if (fv->max_gene_num != UNDEF_ULONG &&
               fv->gene_num >= fv->max_gene_num) {
        filter_node = true;
      }
      else if (fv->min_gene_score != UNDEF_DOUBLE &&
               gt_genome_feature_get_score(gf) < fv->min_gene_score) {
        filter_node = true;
      }
      else if (fv->max_gene_score != UNDEF_DOUBLE &&
               gt_genome_feature_get_score(gf) > fv->max_gene_score) {
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
    gt_genome_node_rec_delete((GT_GenomeNode*) gf);
  else
    gt_queue_add(fv->gt_genome_node_buffer, gf);

  return 0;
}

static int filter_visitor_sequence_region(GenomeVisitor *gv,
                                          GT_SequenceRegion *sr,
                                          GT_UNUSED GT_Error *err)
{
  FilterVisitor *filter_visitor;
  gt_error_check(err);
  filter_visitor = filter_visitor_cast(gv);
  if (!gt_str_length(filter_visitor->seqid) || /* no seqid was specified */
      !gt_str_cmp(filter_visitor->seqid,       /* or seqids are equal */
               gt_genome_node_get_seqid((GT_GenomeNode*) sr))) {
    if (filter_visitor->contain_range.start != UNDEF_ULONG) {
      GT_Range range = gt_genome_node_get_range((GT_GenomeNode*) sr);
      if (gt_range_overlap(range, filter_visitor->contain_range)) {
        /* an overlapping contain range was defined -> update range  */
        range.start = MAX(range.start, filter_visitor->contain_range.start);
        range.end = MIN(range.end, filter_visitor->contain_range.end);
        gt_genome_node_set_range((GT_GenomeNode*) sr, range);
        gt_queue_add(filter_visitor->gt_genome_node_buffer, sr);
      }
      else /* contain range does not overlap with <sr> range -> delete <sr> */
        gt_genome_node_delete((GT_GenomeNode*) sr);
    }
    else
      gt_queue_add(filter_visitor->gt_genome_node_buffer, sr);
  }
  else
    gt_genome_node_rec_delete((GT_GenomeNode*) sr);
  return 0;
}

static int filter_visitor_sequence_node(GenomeVisitor *gv, GT_SequenceNode *sn,
                                        GT_UNUSED GT_Error *err)
{
  FilterVisitor *filter_visitor;
  gt_error_check(err);
  filter_visitor = filter_visitor_cast(gv);
  if (!gt_str_length(filter_visitor->seqid) || /* no seqid was specified */
      !gt_str_cmp(filter_visitor->seqid,       /* or seqids are equal */
               gt_genome_node_get_seqid((GT_GenomeNode*) sn))) {
    gt_queue_add(filter_visitor->gt_genome_node_buffer, sn);
  }
  else
    gt_genome_node_rec_delete((GT_GenomeNode*) sn);
  return 0;
}

const GenomeVisitorClass* filter_visitor_class()
{
  static const GenomeVisitorClass gvc = { sizeof (FilterVisitor),
                                          filter_visitor_free,
                                          filter_visitor_comment,
                                          filter_visitor_genome_feature,
                                          filter_visitor_sequence_region,
                                          filter_visitor_sequence_node };
  return &gvc;
}

GenomeVisitor* filter_visitor_new(GtStr *seqid, GtStr *typefilter,
                                  GT_Range contain_range,
                                  GT_Range overlap_range,
                                  GtStrand strand, GtStrand targetstrand,
                                  bool has_CDS, unsigned long max_gene_length,
                                  unsigned long max_gene_num,
                                  double min_gene_score, double max_gene_score,
                                  double min_average_splice_site_prob,
                                  unsigned long feature_num)
{
  GenomeVisitor *gv = genome_visitor_create(filter_visitor_class());
  FilterVisitor *filter_visitor = filter_visitor_cast(gv);
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

unsigned long filter_visitor_node_buffer_size(GenomeVisitor *gv)
{
  FilterVisitor *filter_visitor = filter_visitor_cast(gv);
  return gt_queue_size(filter_visitor->gt_genome_node_buffer);
}

GT_GenomeNode* filter_visitor_get_node(GenomeVisitor *gv)
{
  FilterVisitor *filter_visitor;
  filter_visitor = filter_visitor_cast(gv);
  return gt_queue_get(filter_visitor->gt_genome_node_buffer);
}
