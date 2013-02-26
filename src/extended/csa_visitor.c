/*
  Copyright (c) 2006-2012 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#include "core/assert_api.h"
#include "core/class_alloc_lock.h"
#include "core/log.h"
#include "core/queue_api.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "extended/csa_variable_strands.h"
#include "extended/csa_visitor.h"
#include "extended/feature_node.h"
#include "extended/feature_type.h"
#include "extended/gff3_defines.h"
#include "extended/node_visitor_api.h"

#define GT_CSA_SOURCE_TAG "gt csa"

struct CSAVisitor {
  const GtNodeVisitor parent_instance;
  GtQueue *gt_genome_node_buffer;
  unsigned long join_length;
  GtArray *cluster;
  GtFeatureNode *buffered_feature;
  GtRange first_range,
          second_range;
  GtStr *first_str,
        *second_str,
        *gt_csa_source_str;
};

#define csa_visitor_cast(GV)\
        gt_node_visitor_cast(gt_csa_visitor_class(), GV)

static void csa_visitor_free(GtNodeVisitor *nv)
{
  CSAVisitor *csa_visitor = csa_visitor_cast(nv);
  gt_queue_delete(csa_visitor->gt_genome_node_buffer);
  gt_array_delete(csa_visitor->cluster);
  gt_str_delete(csa_visitor->gt_csa_source_str);
}

static int csa_visitor_feature_node(GtNodeVisitor *nv, GtFeatureNode *fn,
                                    GT_UNUSED GtError *err)
{
  CSAVisitor *csa_visitor;
  gt_error_check(err);
  csa_visitor = csa_visitor_cast(nv);

  /* determine the first range if necessary */
  if (csa_visitor->buffered_feature) {
    csa_visitor->first_range =
      gt_genome_node_get_range((GtGenomeNode*) csa_visitor->buffered_feature);
    csa_visitor->first_str =
      gt_genome_node_get_seqid((GtGenomeNode*) csa_visitor->buffered_feature);
    gt_assert(!gt_array_size(csa_visitor->cluster));
    gt_array_add(csa_visitor->cluster, csa_visitor->buffered_feature);
    csa_visitor->buffered_feature = NULL;
  }
  else if (!gt_array_size(csa_visitor->cluster)) {
    csa_visitor->first_range = gt_genome_node_get_range((GtGenomeNode*) fn);
    csa_visitor->first_str = gt_genome_node_get_seqid((GtGenomeNode*) fn);
    gt_array_add(csa_visitor->cluster, fn);
    return 0;
  }

  gt_assert(!csa_visitor->buffered_feature);
  csa_visitor->second_range = gt_genome_node_get_range((GtGenomeNode*) fn);
  csa_visitor->second_str = gt_genome_node_get_seqid((GtGenomeNode*) fn);

  if ((gt_str_cmp(csa_visitor->first_str, csa_visitor->second_str) == 0) &&
      (csa_visitor->first_range.end + csa_visitor->join_length >=
       csa_visitor->second_range.start)) {
      /* we are still in the cluster */
      gt_array_add(csa_visitor->cluster, fn);
      /* update first range */
      gt_assert(csa_visitor->second_range.start >=
                csa_visitor->first_range.start);
      if (csa_visitor->second_range.end > csa_visitor->first_range.end)
        csa_visitor->first_range.end = csa_visitor->second_range.end;
  }
  else {
    /* end of cluster -> process it */
    gt_log_log("process cluster");
    csa_visitor->buffered_feature = fn;
    gt_csa_visitor_process_cluster(nv, false);
    csa_visitor->first_range = csa_visitor->second_range;
    csa_visitor->first_str = csa_visitor->second_str;
  }
  return 0;
}

static int csa_visitor_default_func(GtNodeVisitor *nv, GtGenomeNode *gn,
                                    GT_UNUSED GtError *err)
{
  CSAVisitor *csa_visitor;
  gt_error_check(err);
  csa_visitor = csa_visitor_cast(nv);
  gt_queue_add(csa_visitor->gt_genome_node_buffer, gn);
  return 0;
}

static int csa_visitor_comment_node(GtNodeVisitor *nv, GtCommentNode *c,
                                    GtError *err)
{
  return csa_visitor_default_func(nv, (GtGenomeNode*) c, err);
}

static int csa_visitor_meta_node(GtNodeVisitor *nv, GtMetaNode *mn,
                                 GtError *err)
{
  return csa_visitor_default_func(nv, (GtGenomeNode*) mn, err);
}

static int csa_visitor_region_node(GtNodeVisitor *nv, GtRegionNode *rn,
                                   GtError *err)
{
  return csa_visitor_default_func(nv, (GtGenomeNode*) rn, err);
}

static int csa_visitor_sequence_node(GtNodeVisitor *nv, GtSequenceNode *sn,
                                     GtError *err)
{
  gt_csa_visitor_process_cluster(nv, true);
  return csa_visitor_default_func(nv, (GtGenomeNode*) sn, err);
}

static int csa_visitor_eof_node(GtNodeVisitor *nv, GtEOFNode *eofn,
                                GtError *err)
{
  gt_csa_visitor_process_cluster(nv, true);
  return csa_visitor_default_func(nv, (GtGenomeNode*) eofn, err);
}

const GtNodeVisitorClass* gt_csa_visitor_class()
{
  static GtNodeVisitorClass *nvc = NULL;
  gt_class_alloc_lock_enter();
  if (!nvc) {
    nvc = gt_node_visitor_class_new(sizeof (CSAVisitor),
                                    csa_visitor_free,
                                    csa_visitor_comment_node,
                                    csa_visitor_feature_node,
                                    csa_visitor_region_node,
                                    csa_visitor_sequence_node,
                                    csa_visitor_eof_node);
    gt_node_visitor_class_set_meta_node_func(nvc, csa_visitor_meta_node);
  }
  gt_class_alloc_lock_leave();
  return nvc;
}

GtNodeVisitor* gt_csa_visitor_new(unsigned long join_length)
{
  GtNodeVisitor *nv = gt_node_visitor_create(gt_csa_visitor_class());
  CSAVisitor *csa_visitor = csa_visitor_cast(nv);
  csa_visitor->gt_genome_node_buffer = gt_queue_new();
  csa_visitor->join_length = join_length;
  csa_visitor->cluster = gt_array_new(sizeof (GtFeatureNode*));
  csa_visitor->buffered_feature = NULL;
  csa_visitor->gt_csa_source_str = gt_str_new_cstr(GT_CSA_SOURCE_TAG);
  return nv;
}

unsigned long gt_csa_visitor_node_buffer_size(GtNodeVisitor *nv)
{
  CSAVisitor *csa_visitor = csa_visitor_cast(nv);
  return gt_queue_size(csa_visitor->gt_genome_node_buffer);
}

GtGenomeNode* gt_csa_visitor_get_node(GtNodeVisitor *nv)
{
  CSAVisitor *csa_visitor;
  csa_visitor = csa_visitor_cast(nv);
  return gt_queue_get(csa_visitor->gt_genome_node_buffer);
}

static GtRange get_genomic_range(const void *sa)
{
  GtFeatureNode *fn = *(GtFeatureNode**) sa;
  gt_assert(fn && gt_feature_node_has_type(fn, gt_ft_gene));
  return gt_genome_node_get_range((GtGenomeNode*) fn);
}

static GtStrand get_strand(const void *sa)
{
  GtFeatureNode *fn = *(GtFeatureNode**) sa;
  gt_assert(fn && gt_feature_node_has_type(fn, gt_ft_gene));
  return gt_feature_node_get_strand(fn);
}

static int csa_visitor_save_exon(GtFeatureNode *fn, void *data,
                                 GT_UNUSED GtError *err)
{
  GtArray *exon_ranges = (GtArray*) data;
  GtRange range;
  gt_error_check(err);
  gt_assert(fn && exon_ranges);
  if (gt_feature_node_has_type(fn, gt_ft_exon)) {
    range = gt_genome_node_get_range((GtGenomeNode*) fn);
    gt_array_add(exon_ranges, range);
  }
  return 0;
}

static void get_exons(GtArray *exon_ranges, const void *sa)
{
  GtFeatureNode *fn = *(GtFeatureNode**) sa;
  GT_UNUSED int had_err;
  gt_assert(exon_ranges && fn && gt_feature_node_has_type(fn, gt_ft_gene));
  had_err = gt_feature_node_traverse_children(fn, exon_ranges,
                                              csa_visitor_save_exon, false,
                                              NULL);
  gt_assert(!had_err); /* csa_visitor_save_exon() is sane */
  /* we got at least one exon */
  gt_assert(gt_array_size(exon_ranges));
  gt_assert(gt_ranges_are_sorted_and_do_not_overlap(exon_ranges));
}

static void add_sa_to_exon_feature_array(GtArray *exon_nodes,
                                         GtFeatureNode *sa,
                                         GtStr *seqid,
                                         GtStr *gt_csa_source_str,
                                         GtStrand gene_strand)
{
  GtArray *exons_from_sa;
  unsigned long i,
                exon_feature_index = 0,
                exons_from_sa_index = 0;
  GtFeatureNode *exon_feature, *exons_from_sa_feature;
  GtGenomeNode *new_feature;
  GtRange exon_feature_range, exons_from_sa_range;

  gt_assert(exon_nodes && sa);
  gt_assert(gene_strand != GT_STRAND_BOTH); /* is defined */

  exons_from_sa = gt_array_new(sizeof (GtFeatureNode*));
  gt_feature_node_get_exons(sa, exons_from_sa);
  gt_genome_nodes_sort(exons_from_sa);

  while (exon_feature_index < gt_array_size(exon_nodes) &&
         exons_from_sa_index < gt_array_size(exons_from_sa)) {
    exon_feature = *(GtFeatureNode**)
                   gt_array_get(exon_nodes, exon_feature_index);
    exons_from_sa_feature = *(GtFeatureNode**)
                            gt_array_get(exons_from_sa, exons_from_sa_index);

    exon_feature_range = gt_genome_node_get_range((GtGenomeNode*)
                                                  exon_feature);
    exons_from_sa_range =
      gt_genome_node_get_range((GtGenomeNode*) exons_from_sa_feature);

    switch (gt_range_compare(&exon_feature_range, &exons_from_sa_range)) {
      case -1:
        if (gt_range_overlap(&exon_feature_range, &exons_from_sa_range)) {
          if (!gt_range_contains(&exon_feature_range, &exons_from_sa_range)) {
            gt_assert(gt_genome_node_get_start((GtGenomeNode*) exon_feature) <=
              gt_genome_node_get_start((GtGenomeNode*) exons_from_sa_feature));
            gt_assert(gt_genome_node_get_end((GtGenomeNode*) exon_feature) <
                gt_genome_node_get_end((GtGenomeNode*) exons_from_sa_feature));
            /* update right border and score */
            gt_feature_node_set_end(exon_feature,
                                   gt_genome_node_get_end((GtGenomeNode*)
                                                       exons_from_sa_feature));
            if (gt_feature_node_score_is_defined(exons_from_sa_feature)) {
              gt_feature_node_set_score(exon_feature,
                            gt_feature_node_get_score(exons_from_sa_feature));
            }
          }
          exons_from_sa_index++;
        }
        exon_feature_index++;
        break;
      case 0:
        gt_assert(gt_range_overlap(&exon_feature_range, &exons_from_sa_range));
        /* update score if necessary */
        if ((gt_feature_node_score_is_defined(exon_feature) &&
             gt_feature_node_score_is_defined(exons_from_sa_feature) &&
             gt_feature_node_get_score(exon_feature) <
             gt_feature_node_get_score(exons_from_sa_feature)) ||
            (!gt_feature_node_score_is_defined(exon_feature) &&
             gt_feature_node_score_is_defined(exons_from_sa_feature))) {
          gt_feature_node_set_score(exon_feature,
                            gt_feature_node_get_score(exons_from_sa_feature));
        }
        exon_feature_index++;
        exons_from_sa_index++;
        break;
      case 1:
        gt_assert(gt_range_overlap(&exon_feature_range, &exons_from_sa_range));
        gt_assert(gt_genome_node_get_start((GtGenomeNode*) exon_feature) <=
              gt_genome_node_get_start((GtGenomeNode*) exons_from_sa_feature));
        /* update right border and score, if necessary */
        if (gt_genome_node_get_end((GtGenomeNode*) exons_from_sa_feature) >
            gt_genome_node_get_end((GtGenomeNode*) exon_feature)) {
          gt_feature_node_set_end(exon_feature,
                                 gt_genome_node_get_end((GtGenomeNode*)
                                                     exons_from_sa_feature));
          if (gt_feature_node_score_is_defined(exons_from_sa_feature)) {
            gt_feature_node_set_score(exon_feature,
                            gt_feature_node_get_score(exons_from_sa_feature));
          }
        }
        exon_feature_index++;
        exons_from_sa_index++;
        break;
      default: gt_assert(0);
    }
  }

  /* add remaining exons */
  for (i = exons_from_sa_index; i < gt_array_size(exons_from_sa); i++) {
    GtRange range;
    exons_from_sa_feature = *(GtFeatureNode**)
                            gt_array_get(exons_from_sa, i);
    range = gt_genome_node_get_range((GtGenomeNode*) exons_from_sa_feature),
    new_feature = gt_feature_node_new(seqid, gt_ft_exon, range.start, range.end,
                                      gene_strand);
    if (gt_feature_node_score_is_defined(exons_from_sa_feature)) {
      gt_feature_node_set_score((GtFeatureNode*) new_feature,
                            gt_feature_node_get_score(exons_from_sa_feature));
    }
    gt_feature_node_set_source((GtFeatureNode*) new_feature, gt_csa_source_str);
    gt_array_add(exon_nodes, new_feature);
  }

  gt_array_delete(exons_from_sa);
}

#ifndef NDEBUG
static bool genome_nodes_are_sorted_and_do_not_overlap(GtArray *exon_nodes)
{
  GtArray *ranges = gt_array_new(sizeof (GtRange));
  unsigned long i;
  GtRange range;
  bool rval;
  gt_assert(exon_nodes);
  for (i = 0; i < gt_array_size(exon_nodes); i++) {
    range = gt_genome_node_get_range(*(GtGenomeNode**)
                                     gt_array_get(exon_nodes, i));
    gt_array_add(ranges, range);
  }
  rval = gt_ranges_are_sorted_and_do_not_overlap(ranges);
  gt_array_delete(ranges);
  return rval;
}
#endif

static void mRNA_set_target_attribute(GtFeatureNode *mRNA_feature,
                                      const GtCSASpliceForm *csa_splice_form)
{
  unsigned long i;
  GtStr *targets;
  gt_assert(mRNA_feature && csa_splice_form);
  targets = gt_str_new();
  for (i = 0; i < gt_csa_splice_form_num_of_sas(csa_splice_form); i++) {
    GtFeatureNode *sa = *(GtFeatureNode**)
                        gt_csa_splice_form_get_sa(csa_splice_form, i);
    if (gt_feature_node_get_attribute(sa, GT_GFF_TARGET)) {
      if (gt_str_length(targets))
        gt_str_append_char(targets, ',');
      gt_str_append_cstr(targets,
                         gt_feature_node_get_attribute(sa, GT_GFF_TARGET));
    }
  }
  if (gt_str_length(targets)) {
    gt_feature_node_add_attribute(mRNA_feature, GT_GFF_TARGET,
                                  gt_str_get(targets));
  }
  gt_str_delete(targets);
}

static GtFeatureNode* create_mRNA_feature(GtCSASpliceForm *csa_splice_form,
                                          GtStr *gt_csa_source_str)
{
  GtFeatureNode *mRNA_feature;
  GtArray *exon_nodes;
  unsigned long i;
  GtRange range;
  GtStrand strand;
  GtStr *seqid;
  gt_assert(csa_splice_form && gt_csa_source_str);

  /* create the mRNA feature itself */
  seqid = gt_genome_node_get_seqid(*(GtGenomeNode**)
                               gt_csa_splice_form_get_representative(
                                                              csa_splice_form));
  range = gt_csa_splice_form_genomic_range(csa_splice_form),
  strand = gt_csa_splice_form_strand(csa_splice_form),
  mRNA_feature = (GtFeatureNode*)
                 gt_feature_node_new(seqid, gt_ft_mRNA, range.start, range.end,
                                     strand);
  gt_feature_node_set_source(mRNA_feature, gt_csa_source_str);
  mRNA_set_target_attribute(mRNA_feature, csa_splice_form);

  /* create exon features */
  exon_nodes = gt_array_new(sizeof (GtGenomeNode*));
  for (i = 0; i < gt_csa_splice_form_num_of_sas(csa_splice_form); i++) {
    add_sa_to_exon_feature_array(exon_nodes,
                                 *(GtFeatureNode**)
                                 gt_csa_splice_form_get_sa(csa_splice_form, i),
                                 seqid, gt_csa_source_str, strand);
  }
  gt_assert(genome_nodes_are_sorted_and_do_not_overlap(exon_nodes));

  /* add exon features to mRNA feature */
  for (i = 0; i < gt_array_size(exon_nodes); i++) {
    gt_feature_node_add_child(mRNA_feature,
                             *(GtFeatureNode**) gt_array_get(exon_nodes, i));
  }

  gt_array_delete(exon_nodes);

  return mRNA_feature;
}

static GtFeatureNode* create_gene_feature(GtCSAGene *csa_gene,
                                          GtStr *gt_csa_source_str)
{
  GtFeatureNode *gene_feature, *mRNA_feature;
  GtRange range;
  unsigned long i;
  gt_assert(csa_gene && gt_csa_source_str);

  /* create top-level gene feature */
  range = gt_csa_gene_genomic_range(csa_gene),
  gene_feature = (GtFeatureNode*)
    gt_feature_node_new(gt_genome_node_get_seqid(*(GtGenomeNode**)
                                      gt_csa_gene_get_representative(csa_gene)),
                        gt_ft_gene, range.start, range.end,
                        gt_csa_gene_strand(csa_gene));
  gt_feature_node_set_source(gene_feature, gt_csa_source_str);

  /* create mRNA features representing the splice forms */
  for (i = 0; i < gt_csa_gene_num_of_splice_forms(csa_gene); i++) {
    GtCSASpliceForm *csa_splice_form = gt_csa_gene_get_splice_form(csa_gene, i);
    mRNA_feature = create_mRNA_feature(csa_splice_form, gt_csa_source_str);
    gt_feature_node_add_child(gene_feature, mRNA_feature);
  }

  return gene_feature;
}

static void process_csa_genes(GtQueue *gt_genome_node_buffer,
                              GtArray *csa_genes,
                              GtStr *gt_csa_source_str)
{
  unsigned long i;
  gt_assert(csa_genes);
  for (i = 0; i < gt_array_size(csa_genes); i++) {
    GtFeatureNode *gene_feature = create_gene_feature(*(GtCSAGene**)
                                                     gt_array_get(csa_genes, i),
                                                      gt_csa_source_str);
    gt_queue_add(gt_genome_node_buffer, gene_feature);
  }
}

void gt_csa_visitor_process_cluster(GtNodeVisitor *nv, bool final_cluster)
{
  CSAVisitor *csa_visitor = csa_visitor_cast(nv);
  GT_UNUSED GtFeatureNode *first_feature;
  GtArray *csa_genes;
  unsigned long i;

  if (final_cluster) {
    gt_assert(!gt_array_size(csa_visitor->cluster) ||
              !csa_visitor->buffered_feature);
    if (csa_visitor->buffered_feature) {
      gt_array_add(csa_visitor->cluster, csa_visitor->buffered_feature);
      csa_visitor->buffered_feature = NULL;
    }
  }

  if (!gt_array_size(csa_visitor->cluster)) {
    gt_assert(final_cluster);
    return;
  }

  /* compute the consensus spliced alignments */
  first_feature = *(GtFeatureNode**)
                  gt_array_get_first(csa_visitor->cluster);
  csa_genes = gt_csa_variable_strands(gt_array_get_space(csa_visitor->cluster),
                                      gt_array_size(csa_visitor->cluster),
                                      sizeof (GtFeatureNode*),
                                      get_genomic_range,
                                      get_strand, get_exons);

  process_csa_genes(csa_visitor->gt_genome_node_buffer, csa_genes,
                    csa_visitor->gt_csa_source_str);

  for (i = 0; i < gt_array_size(csa_genes); i++)
    gt_csa_gene_delete(*(GtCSAGene**) gt_array_get(csa_genes, i));
  gt_array_delete(csa_genes);

  /* remove the cluster genome nodes */
  for (i = 0; i < gt_array_size(csa_visitor->cluster); i++) {
    gt_genome_node_delete(*(GtGenomeNode**)
                              gt_array_get(csa_visitor->cluster, i));
  }
  gt_array_reset(csa_visitor->cluster);
}
