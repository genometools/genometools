/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include "disc_distri.h"
#include "genome_visitor_rep.h"
#include "stat_visitor.h"

struct Stat_visitor {
  const GenomeVisitor parent_instance;
  unsigned long number_of_genes,
                number_of_mRNAs,
                number_of_exons;
  Disc_distri *gene_length_distribution,
              *gene_score_distribution;
};

#define stat_visitor_cast(GV)\
        genome_visitor_cast(stat_visitor_class(), GV)

static void stat_visitor_free(GenomeVisitor *gv)
{
  Stat_visitor *stat_visitor = stat_visitor_cast(gv);
  disc_distri_free(stat_visitor->gene_length_distribution);
  disc_distri_free(stat_visitor->gene_score_distribution);
}

static void stat_visitor_genome_feature(GenomeVisitor *gv, Genome_feature *gf,
                                        /*@unused@*/ Log *l)
{
  Stat_visitor *stat_visitor = stat_visitor_cast(gv);
  switch (genome_feature_get_type(gf)) {
    case gft_gene:
      stat_visitor->number_of_genes++;
     if (stat_visitor->gene_length_distribution) {
       disc_distri_add(stat_visitor->gene_length_distribution,
                       range_length(genome_node_get_range((GenomeNode*) gf)));
     }
     if (stat_visitor->gene_score_distribution) {
       disc_distri_add(stat_visitor->gene_score_distribution,
                       genome_feature_get_score(gf) * 100.0);
     }
     break;
    case gft_mRNA:
      stat_visitor->number_of_mRNAs++;
      break;
    case gft_exon:
      stat_visitor->number_of_exons++;
      break;
    default:
      assert(0);
  }
}

const GenomeVisitorClass* stat_visitor_class()
{
  static const GenomeVisitorClass gvc = { sizeof(Stat_visitor),
                                            stat_visitor_free,
                                            NULL,
                                            stat_visitor_genome_feature,
                                            NULL,
                                            NULL };
  return &gvc;
}

GenomeVisitor* stat_visitor_new(unsigned int gene_length_distri,
                                 unsigned int gene_score_distri)
{
  GenomeVisitor *gv = genome_visitor_create(stat_visitor_class());
  Stat_visitor *stat_visitor = stat_visitor_cast(gv);
  if (gene_length_distri)
    stat_visitor->gene_length_distribution = disc_distri_new();
  if (gene_score_distri)
    stat_visitor->gene_score_distribution = disc_distri_new();
  return gv;
}

void stat_visitor_show_stats(GenomeVisitor *gv)
{
  Stat_visitor *stat_visitor = stat_visitor_cast(gv);
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
}
