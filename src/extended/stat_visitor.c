/*
  Copyright (c) 2006-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#include <assert.h>
#include "core/disc_distri.h"
#include "core/unused_api.h"
#include "extended/node_visitor_rep.h"
#include "extended/stat_visitor.h"

struct StatVisitor {
  const GtNodeVisitor parent_instance;
  unsigned long number_of_sequence_regions,
                number_of_genes,
                number_of_protein_coding_genes,
                number_of_mRNAs,
                number_of_exons,
                number_of_CDSs,
                number_of_LTR_retrotransposons,
                exon_number_for_distri;
  unsigned long long total_length_of_sequence_regions;
  DiscDistri *gene_length_distribution,
             *gene_score_distribution,
             *exon_length_distribution,
             *exon_number_distribution,
             *intron_length_distribution;
};

#define stat_visitor_cast(GV)\
        gt_node_visitor_cast(stat_visitor_class(), GV)

static void stat_visitor_free(GtNodeVisitor *gv)
{
  StatVisitor *stat_visitor = stat_visitor_cast(gv);
  disc_distri_delete(stat_visitor->gene_length_distribution);
  disc_distri_delete(stat_visitor->gene_score_distribution);
  disc_distri_delete(stat_visitor->exon_length_distribution);
  disc_distri_delete(stat_visitor->exon_number_distribution);
  disc_distri_delete(stat_visitor->intron_length_distribution);
}

static int add_exon_number(GtGenomeNode *gn, void *data,
                           GT_UNUSED GtError *err)
{
  StatVisitor *stat_visitor = (StatVisitor*) data;
  GtFeatureNode *gf = (GtFeatureNode*) gn;
  gt_error_check(err);
  assert(stat_visitor && gf);
  if (gt_feature_node_has_type(gf, gft_exon))
    stat_visitor->exon_number_for_distri++;
  return 0;
}

static int compute_statistics(GtGenomeNode *gn, void *data, GtError *err)
{
  StatVisitor *stat_visitor;
  GtFeatureNode *gf;
  int rval;
  gt_error_check(err);
  assert(data);
  stat_visitor = (StatVisitor*) data;
  gf = (GtFeatureNode*) gn;
  if (gt_feature_node_has_type(gf, gft_gene)) {
    stat_visitor->number_of_genes++;
    if (gt_feature_node_has_CDS(gf))
      stat_visitor->number_of_protein_coding_genes++;
    if (stat_visitor->gene_length_distribution) {
      disc_distri_add(stat_visitor->gene_length_distribution,
                      gt_range_length(gt_genome_node_get_range((GtGenomeNode*)
                                                               gf)));
    }
    if (stat_visitor->gene_score_distribution) {
      disc_distri_add(stat_visitor->gene_score_distribution,
                     gt_feature_node_get_score(gf) * 100.0);
    }
  }
  else if (gt_feature_node_has_type(gf, gft_mRNA)) {
    stat_visitor->number_of_mRNAs++;
  }
  else if (gt_feature_node_has_type(gf, gft_exon)) {
    stat_visitor->number_of_exons++;
    if (stat_visitor->exon_length_distribution) {
      disc_distri_add(stat_visitor->exon_length_distribution,
                      gt_range_length(gt_genome_node_get_range((GtGenomeNode*)
                                                               gf)));
    }
  }
  else if (gt_feature_node_has_type(gf, gft_CDS)) {
    stat_visitor->number_of_CDSs++;
  }
  else if (gt_feature_node_has_type(gf, gft_intron)) {
    if (stat_visitor->intron_length_distribution) {
      disc_distri_add(stat_visitor->intron_length_distribution,
                      gt_range_length(gt_genome_node_get_range((GtGenomeNode*)
                                                               gf)));
    }
  }
  else if (gt_feature_node_has_type(gf, gft_LTR_retrotransposon)) {
    stat_visitor->number_of_LTR_retrotransposons++;
  }
  if (stat_visitor->exon_number_distribution) {
    stat_visitor->exon_number_for_distri = 0;
    rval = gt_genome_node_traverse_direct_children(gn, stat_visitor,
                                                add_exon_number, err);
    assert(!rval); /* add_exon_number() is sane */
    if (stat_visitor->exon_number_for_distri) {
      disc_distri_add(stat_visitor->exon_number_distribution,
                     stat_visitor->exon_number_for_distri);
    }
  }
  return 0;
}

static int stat_visitor_genome_feature(GtNodeVisitor *gv, GtFeatureNode *gf,
                                       GtError *err)
{
  StatVisitor *stat_visitor;
  gt_error_check(err);
  stat_visitor = stat_visitor_cast(gv);
  return gt_genome_node_traverse_children((GtGenomeNode*) gf, stat_visitor,
                                       compute_statistics, false, err);
}

static int stat_visitor_region_node(GtNodeVisitor *gv, GtRegionNode *rn,
                                    GT_UNUSED GtError *err)
{
  StatVisitor *stat_visitor;
  gt_error_check(err);
  stat_visitor = stat_visitor_cast(gv);
  stat_visitor->number_of_sequence_regions++;
  stat_visitor->total_length_of_sequence_regions +=
    gt_range_length(gt_genome_node_get_range((GtGenomeNode*) rn));
  return 0;
}

const GtNodeVisitorClass* stat_visitor_class()
{
  static const GtNodeVisitorClass gvc = { sizeof (StatVisitor),
                                          stat_visitor_free,
                                          NULL,
                                          stat_visitor_genome_feature,
                                          stat_visitor_region_node,
                                          NULL };
  return &gvc;
}

GtNodeVisitor* stat_visitor_new(bool gene_length_distri,
                                bool gene_score_distri,
                                bool exon_length_distri,
                                bool exon_number_distri,
                                bool intron_length_distri)
{
  GtNodeVisitor *gv = gt_node_visitor_create(stat_visitor_class());
  StatVisitor *stat_visitor = stat_visitor_cast(gv);
  if (gene_length_distri)
    stat_visitor->gene_length_distribution = disc_distri_new();
  if (gene_score_distri)
    stat_visitor->gene_score_distribution = disc_distri_new();
  if (exon_length_distri)
    stat_visitor->exon_length_distribution = disc_distri_new();
  if (exon_number_distri)
    stat_visitor->exon_number_distribution = disc_distri_new();
  if (intron_length_distri)
    stat_visitor->intron_length_distribution = disc_distri_new();
  return gv;
}

void stat_visitor_show_stats(GtNodeVisitor *gv)
{
  StatVisitor *stat_visitor = stat_visitor_cast(gv);
  if (stat_visitor->number_of_sequence_regions) {
    printf("sequence regions: %lu (total length: %llu)\n",
           stat_visitor->number_of_sequence_regions,
           stat_visitor->total_length_of_sequence_regions);
  }
  if (stat_visitor->number_of_genes)
    printf("genes: %lu\n", stat_visitor->number_of_genes);
  if (stat_visitor->number_of_protein_coding_genes) {
    printf("protein-coding genes: %lu\n",
           stat_visitor->number_of_protein_coding_genes);
  }
  if (stat_visitor->number_of_mRNAs)
    printf("mRNAs: %lu\n", stat_visitor->number_of_mRNAs);
  if (stat_visitor->number_of_exons)
    printf("exons: %lu\n", stat_visitor->number_of_exons);
  if (stat_visitor->number_of_CDSs)
    printf("CDSs: %lu\n", stat_visitor->number_of_CDSs);
  if (stat_visitor->number_of_LTR_retrotransposons) {
    printf("LTR_retrotransposons: %lu\n",
           stat_visitor->number_of_LTR_retrotransposons);
  }
  if (stat_visitor->gene_length_distribution) {
    printf("gene length distribution:\n");
    disc_distri_show(stat_visitor->gene_length_distribution);
  }
  if (stat_visitor->gene_score_distribution) {
    printf("gene score distribution:\n");
    disc_distri_show(stat_visitor->gene_score_distribution);
  }
  if (stat_visitor->exon_length_distribution) {
    printf("exon length distribution:\n");
    disc_distri_show(stat_visitor->exon_length_distribution);
  }
  if (stat_visitor->exon_number_distribution) {
    printf("exon number distribution:\n");
    disc_distri_show(stat_visitor->exon_number_distribution);
  }
  if (stat_visitor->intron_length_distribution) {
    printf("intron length distribution:\n");
    disc_distri_show(stat_visitor->intron_length_distribution);
  }
}
