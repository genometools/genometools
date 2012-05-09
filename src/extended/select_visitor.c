/*
  Copyright (c) 2007-2012 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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
#include "core/ma.h"
#include "core/minmax.h"
#include "core/queue_api.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "extended/feature_node.h"
#include "extended/gff3_defines.h"
#include "extended/gff3_parser.h"
#include "extended/node_visitor_api.h"
#include "extended/script_filter.h"
#include "extended/select_visitor.h"

typedef enum {
  GT_SELECT_AND,
  GT_SELECT_OR
} GtSelectLogic;

struct GtSelectVisitor {
  const GtNodeVisitor parent_instance;
  GtQueue *node_buffer;
  GtStr *seqid,
        *source;
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
         min_average_splice_site_prob,
         single_intron_factor;
  GtStrArray *select_files;
  GtSelectLogic select_logic;
  bool is_lua;
  GtArray *script_filters;
  GtSelectNodeFunc drophandler;
  void *data;
};

#define select_visitor_cast(GV)\
        gt_node_visitor_cast(gt_select_visitor_class(), GV)

static void select_visitor_free(GtNodeVisitor *nv)
{
  int i;
  GtSelectVisitor *select_visitor = select_visitor_cast(nv);
  gt_str_delete(select_visitor->source);
  gt_str_delete(select_visitor->seqid);
  if (gt_array_size(select_visitor->script_filters) > 0) {
    for (i = 0; i < gt_array_size(select_visitor->script_filters); i++) {
      GtScriptFilter *todelete =
            *(GtScriptFilter**) gt_array_get(select_visitor->script_filters, i);
      gt_script_filter_delete(todelete);
    }
  }
  gt_array_delete(select_visitor->script_filters);
  gt_queue_delete(select_visitor->node_buffer);
}

static int select_visitor_comment_node(GtNodeVisitor *nv, GtCommentNode *c,
                                       GT_UNUSED GtError *err)
{
  GtSelectVisitor *select_visitor;
  gt_error_check(err);
  select_visitor = select_visitor_cast(nv);
  gt_queue_add(select_visitor->node_buffer, c);
  return 0;
}

static int select_visitor_meta_node(GtNodeVisitor *nv, GtMetaNode *mn,
                                    GT_UNUSED GtError *err)
{
  GtSelectVisitor *select_visitor;
  gt_error_check(err);
  select_visitor = select_visitor_cast(nv);
  gt_queue_add(select_visitor->node_buffer, mn);
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
      (target = gt_feature_node_get_attribute(fn, GT_GFF_TARGET))) {
    unsigned long num_of_targets;
    GtStrand parsed_strand;
    GT_UNUSED int had_err;
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

static bool filter_min_average_ssp(GtFeatureNode *fn, double minaveragessp,
                                   double singleintronfactor)
{
  gt_assert(fn);
  if (minaveragessp != GT_UNDEF_DOUBLE && gt_feature_node_has_splice_site(fn)) {
    unsigned long num_ss;
    double avg_ssp = gt_feature_node_average_splice_site_prob(fn, &num_ss);
    if (num_ss <= 2 && avg_ssp < singleintronfactor * minaveragessp)
      return true;
    else if (avg_ssp < minaveragessp)
      return true;
  }
  return false;
}

static int filter_lua(GtArray *filters, GtFeatureNode *fn, GtSelectLogic logic,
                      bool *select_node, GtError *err)
{
  int had_err = 0;
  unsigned long i;
  gt_assert(filters && fn && select_node);
  gt_error_check(err);

  for (i = 0; !had_err && i < gt_array_size(filters); i++) {
    bool result;
    GtScriptFilter *sf = *(GtScriptFilter**) gt_array_get(filters, i);
    had_err = gt_script_filter_run(sf, fn, &result, err);

    if (!had_err && i == 0) {
      *select_node = result;
    } else if (!had_err && logic == GT_SELECT_AND) {
      *select_node = *select_node || result;
      if (*select_node) {
        break;
      }
    } else if (!had_err && logic == GT_SELECT_OR) {
      *select_node = *select_node && result;
      if (!*select_node) {
        break;
      }
    }
  }
  return had_err;
}

static int select_visitor_feature_node(GtNodeVisitor *nv,
                                       GtFeatureNode *fn,
                                       GtError *err)
{
  GtSelectVisitor *fv;
  bool select_node = false;
  int had_err = 0;
  gt_error_check(err);
  fv = select_visitor_cast(nv);
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
        select_node = true;
      }
      else if (fv->max_gene_num != GT_UNDEF_ULONG &&
               fv->gene_num >= fv->max_gene_num) {
        select_node = true;
      }
      else if (fv->min_gene_score != GT_UNDEF_DOUBLE &&
               gt_feature_node_get_score(fn) < fv->min_gene_score) {
        select_node = true;
      }
      else if (fv->max_gene_score != GT_UNDEF_DOUBLE &&
               gt_feature_node_get_score(fn) > fv->max_gene_score) {
        select_node = true;
      }
      else if (fv->feature_num != GT_UNDEF_ULONG &&
               fv->feature_num != fv->current_feature) {
        select_node = true;
      }
      if (!select_node)
        fv->gene_num++; /* gene passed filter */
    }
  }
  else
    select_node = true;

  if (!select_node)
    select_node = filter_contain_range(fn, fv->contain_range);

  if (!select_node)
    select_node = filter_overlap_range(fn, fv->overlap_range);

  if (!select_node)
    select_node = filter_strand(fn, fv->strand);

  if (!select_node)
    select_node = filter_targetstrand(fn, fv->targetstrand);

  if (!select_node)
    select_node = filter_has_CDS(fn, fv->has_CDS);

  if (!select_node) {
    select_node = filter_min_average_ssp(fn, fv->min_average_splice_site_prob,
                                         fv->single_intron_factor);
  }

  if (fv->is_lua && !select_node)
    had_err = filter_lua(fv->script_filters, fn, fv->select_logic,
                         &select_node, err);

  if (select_node && !had_err)
    fv->drophandler((GtGenomeNode*) fn, fv->data, err);

  if (select_node)
    gt_genome_node_delete((GtGenomeNode*) fn);
  else
    gt_queue_add(fv->node_buffer, fn);

  return had_err;
}

static int select_visitor_region_node(GtNodeVisitor *nv, GtRegionNode *rn,
                                      GT_UNUSED GtError *err)
{
  GtSelectVisitor *select_visitor;
  gt_error_check(err);
  select_visitor = select_visitor_cast(nv);
  if (!gt_str_length(select_visitor->seqid) || /* no seqid was specified */
      !gt_str_cmp(select_visitor->seqid,       /* or seqids are equal */
               gt_genome_node_get_seqid((GtGenomeNode*) rn))) {
    if (select_visitor->contain_range.start != GT_UNDEF_ULONG) {
      GtRange range = gt_genome_node_get_range((GtGenomeNode*) rn);
      if (gt_range_overlap(&range, &select_visitor->contain_range)) {
        /* an overlapping contain range was defined -> update range  */
        range.start = MAX(range.start, select_visitor->contain_range.start);
        range.end = MIN(range.end, select_visitor->contain_range.end);
        gt_genome_node_set_range((GtGenomeNode*) rn, &range);
        gt_queue_add(select_visitor->node_buffer, rn);
      }
      else {
        /* contain range does not overlap with <rn> range -> handle <rn> */
        select_visitor->drophandler((GtGenomeNode*) rn, select_visitor->data,
                                    err);
        gt_genome_node_delete((GtGenomeNode*) rn);
      }
    }
    else
      gt_queue_add(select_visitor->node_buffer, rn);
  }
  else {
    select_visitor->drophandler((GtGenomeNode*) rn, select_visitor->data, err);
    gt_genome_node_delete((GtGenomeNode*) rn);
  }
  return 0;
}

static int select_visitor_sequence_node(GtNodeVisitor *nv, GtSequenceNode *sn,
                                        GT_UNUSED GtError *err)
{
  GtSelectVisitor *select_visitor;
  gt_error_check(err);
  select_visitor = select_visitor_cast(nv);
  if (!gt_str_length(select_visitor->seqid) || /* no seqid was specified */
      !gt_str_cmp(select_visitor->seqid,       /* or seqids are equal */
                  gt_genome_node_get_seqid((GtGenomeNode*) sn))) {
    gt_queue_add(select_visitor->node_buffer, sn);
  }
  else {
    select_visitor->drophandler((GtGenomeNode*) sn, select_visitor->data, err);
    gt_genome_node_delete((GtGenomeNode*) sn);
  }
  return 0;
}

static int select_visitor_eof_node(GtNodeVisitor *nv, GtEOFNode *eofn,
                                   GT_UNUSED GtError *err)
{
  GtSelectVisitor *select_visitor;
  gt_error_check(err);
  select_visitor = select_visitor_cast(nv);
  gt_queue_add(select_visitor->node_buffer, eofn);
  return 0;
}

const GtNodeVisitorClass* gt_select_visitor_class()
{
  static GtNodeVisitorClass *nvc = NULL;
  if (!nvc) {
    nvc = gt_node_visitor_class_new(sizeof (GtSelectVisitor),
                                    select_visitor_free,
                                    select_visitor_comment_node,
                                    select_visitor_feature_node,
                                    select_visitor_region_node,
                                    select_visitor_sequence_node,
                                    select_visitor_eof_node);
    gt_node_visitor_class_set_meta_node_func(nvc, select_visitor_meta_node);
  }
  return nvc;
}

GtNodeVisitor* gt_select_visitor_new(GtStr *seqid,
                                     GtStr *source,
                                     const GtRange *contain_range,
                                     const GtRange *overlap_range,
                                     GtStrand strand, GtStrand targetstrand,
                                     bool has_CDS,
                                     unsigned long max_gene_length,
                                     unsigned long max_gene_num,
                                     double min_gene_score,
                                     double max_gene_score,
                                     double min_average_splice_site_prob,
                                     unsigned long feature_num,
                                     GtStrArray *select_files,
                                     GtStr *select_logic,
                                     GtError *err)
{
  GtNodeVisitor *nv = gt_node_visitor_create(gt_select_visitor_class());
  GtSelectVisitor *select_visitor = select_visitor_cast(nv);
  select_visitor->node_buffer = gt_queue_new();
  select_visitor->seqid = gt_str_ref(seqid);
  select_visitor->source = gt_str_ref(source);
  if (contain_range)
    select_visitor->contain_range = *contain_range;
  else {
    select_visitor->contain_range.start = GT_UNDEF_ULONG;
    select_visitor->contain_range.end   = GT_UNDEF_ULONG;
  }
  if (overlap_range)
    select_visitor->overlap_range = *overlap_range;
  else {
    select_visitor->overlap_range.start = GT_UNDEF_ULONG;
    select_visitor->overlap_range.end   = GT_UNDEF_ULONG;
  }
  select_visitor->strand = strand;
  select_visitor->targetstrand = targetstrand;
  select_visitor->has_CDS = has_CDS;
  select_visitor->max_gene_length = max_gene_length;
  select_visitor->gene_num = 0;
  select_visitor->max_gene_num = max_gene_num;
  select_visitor->min_gene_score = min_gene_score;
  select_visitor->max_gene_score = max_gene_score;
  select_visitor->min_average_splice_site_prob = min_average_splice_site_prob;
  select_visitor->single_intron_factor = 1.0;
  select_visitor->feature_num = feature_num;
  select_visitor->select_files = select_files;
  select_visitor->is_lua = false;

  if (gt_str_array_size(select_visitor->select_files) > 0) {
    int i;
    select_visitor->is_lua = true;
    select_visitor->script_filters = gt_array_new(sizeof (GtScriptFilter*));
    for (i = 0; i < gt_str_array_size(select_visitor->select_files); i++) {
      GtScriptFilter *sf;
      sf = gt_script_filter_new(gt_str_array_get(select_visitor->select_files,
                                                 i),
                                err);
      if (!sf) {
        gt_node_visitor_delete(nv);
        return NULL;
        break;
      } else {
        gt_array_add(select_visitor->script_filters, sf);
      }
    }
  }
  if (strcmp(gt_str_get(select_logic), "AND") == 0) {
    select_visitor->select_logic = GT_SELECT_AND;
  } else {
    select_visitor->select_logic = GT_SELECT_OR;
  }
  return nv;
}

void gt_select_visitor_set_single_intron_factor(GtNodeVisitor *nv,
                                                double single_intron_factor)
{
  GtSelectVisitor *select_visitor = select_visitor_cast(nv);
  gt_assert(single_intron_factor >= 1.0);
  select_visitor->single_intron_factor = single_intron_factor;
}

unsigned long gt_select_visitor_node_buffer_size(GtNodeVisitor *nv)
{
  GtSelectVisitor *select_visitor = select_visitor_cast(nv);
  return gt_queue_size(select_visitor->node_buffer);
}

GtGenomeNode* gt_select_visitor_get_node(GtNodeVisitor *nv)
{
  GtSelectVisitor *select_visitor = select_visitor_cast(nv);
  return gt_queue_get(select_visitor->node_buffer);
}

void gt_select_visitor_set_drophandler(GtSelectVisitor *fv,
                                       GtSelectNodeFunc fp,
                                       void *data)
{
  gt_assert(fv && fp != NULL);
  fv->drophandler = fp;
  fv->data = data;
}
