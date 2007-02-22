/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include "disc_distri.h"
#include "genome_visitor_rep.h"
#include "stat_visitor.h"

struct StatVisitor {
  const GenomeVisitor parent_instance;
  unsigned long number_of_genes,
                number_of_mRNAs,
                number_of_exons;
  DiscDistri *gene_length_distribution,
             *gene_score_distribution,
             *exon_length_distribution,
             *intron_length_distribution;
};

#define stat_visitor_cast(GV)\
        genome_visitor_cast(stat_visitor_class(), GV)

static void stat_visitor_free(GenomeVisitor *gv, Env *env)
{
  StatVisitor *stat_visitor = stat_visitor_cast(gv);
  disc_distri_delete(stat_visitor->gene_length_distribution, env);
  disc_distri_delete(stat_visitor->gene_score_distribution, env);
}

static int stat_visitor_genome_feature(GenomeVisitor *gv, GenomeFeature *gf,
                                       Env *env)
{
  StatVisitor *stat_visitor;
  env_error_check(env);
  stat_visitor = stat_visitor_cast(gv);
  switch (genome_feature_get_type(gf)) {
    case gft_gene:
      stat_visitor->number_of_genes++;
      if (stat_visitor->gene_length_distribution) {
        disc_distri_add(stat_visitor->gene_length_distribution,
                        range_length(genome_node_get_range((GenomeNode*) gf)),
                        env);
      }
      if (stat_visitor->gene_score_distribution) {
        disc_distri_add(stat_visitor->gene_score_distribution,
                        genome_feature_get_score(gf) * 100.0, env);
      }
      break;
    case gft_mRNA:
      stat_visitor->number_of_mRNAs++;
      break;
    case gft_exon:
      stat_visitor->number_of_exons++;
      if (stat_visitor->exon_length_distribution) {
        disc_distri_add(stat_visitor->exon_length_distribution,
                        range_length(genome_node_get_range((GenomeNode*) gf)),
                        env);
      }
      break;
    case gft_intron:
      if (stat_visitor->intron_length_distribution) {
        disc_distri_add(stat_visitor->intron_length_distribution,
                        range_length(genome_node_get_range((GenomeNode*) gf)),
                        env);
      }
      break;
    default: assert(1); /* nothing to do for all other types */
  }
  return 0;
}

const GenomeVisitorClass* stat_visitor_class()
{
  static const GenomeVisitorClass gvc = { sizeof (StatVisitor),
                                          stat_visitor_free,
                                          NULL,
                                          stat_visitor_genome_feature,
                                          NULL,
                                          NULL };
  return &gvc;
}

GenomeVisitor* stat_visitor_new(bool gene_length_distri,
                                bool gene_score_distri,
                                bool exon_length_distri,
                                bool intron_length_distri)
{
  GenomeVisitor *gv = genome_visitor_create(stat_visitor_class());
  StatVisitor *stat_visitor = stat_visitor_cast(gv);
  if (gene_length_distri)
    stat_visitor->gene_length_distribution = disc_distri_new();
  if (gene_score_distri)
    stat_visitor->gene_score_distribution = disc_distri_new();
  if (exon_length_distri)
    stat_visitor->exon_length_distribution = disc_distri_new();
  if (intron_length_distri)
    stat_visitor->intron_length_distribution = disc_distri_new();
  return gv;
}

void stat_visitor_show_stats(GenomeVisitor *gv)
{
  StatVisitor *stat_visitor = stat_visitor_cast(gv);
  if (stat_visitor->number_of_exons)
    printf("genes: %lu\n", stat_visitor->number_of_genes);
  if (stat_visitor->number_of_mRNAs)
    printf("mRNAs: %lu\n", stat_visitor->number_of_mRNAs);
  if (stat_visitor->number_of_exons)
    printf("exons: %lu\n", stat_visitor->number_of_exons);
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
  if (stat_visitor->intron_length_distribution) {
    printf("intron length distribution:\n");
    disc_distri_show(stat_visitor->intron_length_distribution);
  }
}
