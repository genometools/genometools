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
#include "libgtcore/disc_distri.h"
#include "libgtcore/unused.h"
#include "libgtext/genome_visitor_rep.h"
#include "libgtext/stat_visitor.h"

struct StatVisitor {
  const GenomeVisitor parent_instance;
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
        genome_visitor_cast(stat_visitor_class(), GV)

static void stat_visitor_free(GenomeVisitor *gv)
{
  StatVisitor *stat_visitor = stat_visitor_cast(gv);
  disc_distri_delete(stat_visitor->gene_length_distribution);
  disc_distri_delete(stat_visitor->gene_score_distribution);
  disc_distri_delete(stat_visitor->exon_length_distribution);
  disc_distri_delete(stat_visitor->exon_number_distribution);
  disc_distri_delete(stat_visitor->intron_length_distribution);
}

static int add_exon_number(GenomeNode *gn, void *data, UNUSED Error *err)
{
  StatVisitor *stat_visitor = (StatVisitor*) data;
  GenomeFeature *gf = (GenomeFeature*) gn;
  error_check(err);
  assert(stat_visitor && gf);
  if (genome_feature_get_type(gf) == gft_exon)
    stat_visitor->exon_number_for_distri++;
  return 0;
}

static int compute_statistics(GenomeNode *gn, void *data, Error *e)
{
  StatVisitor *stat_visitor;
  GenomeFeature *gf;
  int rval;
  error_check(e);
  assert(data);
  stat_visitor = (StatVisitor*) data;
  gf = (GenomeFeature*) gn;
  switch (genome_feature_get_type(gf)) {
    case gft_gene:
      stat_visitor->number_of_genes++;
      if (genome_feature_has_CDS(gf))
        stat_visitor->number_of_protein_coding_genes++;
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
      if (stat_visitor->exon_length_distribution) {
        disc_distri_add(stat_visitor->exon_length_distribution,
                       range_length(genome_node_get_range((GenomeNode*) gf)));
      }
      break;
    case gft_CDS:
      stat_visitor->number_of_CDSs++;
      break;
    case gft_intron:
      if (stat_visitor->intron_length_distribution) {
        disc_distri_add(stat_visitor->intron_length_distribution,
                       range_length(genome_node_get_range((GenomeNode*) gf)));
      }
      break;
    case gft_LTR_retrotransposon:
      stat_visitor->number_of_LTR_retrotransposons++;
      break;
    default: assert(1); /* nothing to do for all other types */
  }
  if (stat_visitor->exon_number_distribution) {
    stat_visitor->exon_number_for_distri = 0;
    rval = genome_node_traverse_direct_children(gn, stat_visitor,
                                                add_exon_number, e);
    assert(!rval); /* add_exon_number() is sane */
    if (stat_visitor->exon_number_for_distri) {
      disc_distri_add(stat_visitor->exon_number_distribution,
                     stat_visitor->exon_number_for_distri);
    }
  }
  return 0;
}

static int stat_visitor_genome_feature(GenomeVisitor *gv, GenomeFeature *gf,
                                       Error *e)
{
  StatVisitor *stat_visitor;
  error_check(e);
  stat_visitor = stat_visitor_cast(gv);
  return genome_node_traverse_children((GenomeNode*) gf, stat_visitor,
                                       compute_statistics, false, e);
}

static int stat_visitor_sequence_region(GenomeVisitor *gv, SequenceRegion *sr,
                                        UNUSED Error *err)
{
  StatVisitor *stat_visitor;
  error_check(err);
  stat_visitor = stat_visitor_cast(gv);
  stat_visitor->number_of_sequence_regions++;
  stat_visitor->total_length_of_sequence_regions +=
    range_length(genome_node_get_range((GenomeNode*) sr));
  return 0;
}

const GenomeVisitorClass* stat_visitor_class()
{
  static const GenomeVisitorClass gvc = { sizeof (StatVisitor),
                                          stat_visitor_free,
                                          NULL,
                                          stat_visitor_genome_feature,
                                          stat_visitor_sequence_region };
  return &gvc;
}

GenomeVisitor* stat_visitor_new(bool gene_length_distri,
                                bool gene_score_distri,
                                bool exon_length_distri,
                                bool exon_number_distri,
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
  if (exon_number_distri)
    stat_visitor->exon_number_distribution = disc_distri_new();
  if (intron_length_distri)
    stat_visitor->intron_length_distribution = disc_distri_new();
  return gv;
}

void stat_visitor_show_stats(GenomeVisitor *gv)
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
