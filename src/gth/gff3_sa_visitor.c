/*
  Copyright (c) 2004-2011, 2013 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2004-2008       Center for Bioinformatics, University of Hamburg

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

#include "core/unused_api.h"
#include "extended/feature_type.h"
#include "extended/genome_node.h"
#include "gth/gff3_sa_visitor.h"
#include "gth/sa_visitor_rep.h"
#include "gth/region_factory.h"

#define SA_DELIMITERLINECHAR  '*'

struct GthGFF3SAVisitor {
  const GthSAVisitor parent_instance;
  GthInput *input;
  GthRegionFactory *region_factory;
  GtStr *gthsourcetag;
  GtArray *nodes;
  GtFile *outfp;
};

#define gff3_sa_visitor_cast(GV)\
        gth_sa_visitor_cast(gth_gff3_sa_visitor_class(), GV)

static void gff3_sa_visitor_free(GthSAVisitor *sa_visitor)
{
  unsigned long i;
  GthGFF3SAVisitor *visitor = gff3_sa_visitor_cast(sa_visitor);
  for (i = 0; i < gt_array_size(visitor->nodes); i++)
    gt_genome_node_delete(*(GtGenomeNode**) gt_array_get(visitor->nodes, i));
  gt_array_delete(visitor->nodes);
  gt_str_delete(visitor->gthsourcetag);
  gth_region_factory_delete(visitor->region_factory);
}

static void gff3_sa_visitor_preface(GthSAVisitor *sa_visitor)
{
  GthGFF3SAVisitor *visitor = gff3_sa_visitor_cast(sa_visitor);
  gth_region_factory_save(visitor->region_factory, visitor->nodes,
                          visitor->input);
}

static void save_sa_in_gff3(GthSA *sa, GthRegionFactory *region_factory,
                            GtArray *nodes, GtStr *gthsourcetag, bool md5ids)
{
  GtFeatureNode *gene_feature, *exon_feature;
  GtRange range;
  GtStr *seqid;
  unsigned long i;
  long offset;

  gt_assert(sa && region_factory && nodes);

  seqid = gth_region_factory_get_seqid(region_factory, gth_sa_gen_file_num(sa),
                                       gth_sa_gen_seq_num(sa));
  offset = gth_region_factory_offset(region_factory, gth_sa_gen_file_num(sa),
                                     gth_sa_gen_seq_num(sa)) - 1;
  /* create gene feature */
  range.start = gth_sa_left_genomic_exon_border(sa, 0);
  range.end   = gth_sa_right_genomic_exon_border(sa, gth_sa_num_of_exons(sa)-1);
  range = gt_range_reorder(range);
  range = gt_range_offset(&range, offset);
  gene_feature = (GtFeatureNode*)
                 gt_feature_node_new(seqid, gt_ft_gene, range.start, range.end,
                                     gth_sa_gen_strand(sa));
  gt_feature_node_set_source(gene_feature, gthsourcetag);
  gt_feature_node_set_score(gene_feature, gth_sa_score(sa));
  /* set target attribute, if possible */
  if (strlen(gth_sa_gff3_target_attribute(sa, md5ids))) {
    gt_feature_node_add_attribute(gene_feature, "Target",
                                  gth_sa_gff3_target_attribute(sa, md5ids));
  }

  for (i = 0; i < gth_sa_num_of_exons(sa); i++) {
    /* create splice sites, if necessary */
    if (i > 0) {
      GtGenomeNode *ss_feature;
      /* donor site */
      range = gth_sa_donor_site_range(sa, i-1),
      range = gt_range_offset(&range, offset);
      ss_feature = gt_feature_node_new(seqid, gt_ft_five_prime_cis_splice_site,
                                       range.start, range.end,
                                       gth_sa_gen_strand(sa));
      gt_feature_node_set_source((GtFeatureNode*) ss_feature, gthsourcetag);
      gt_feature_node_set_score((GtFeatureNode*) ss_feature,
                                gth_sa_donor_site_prob(sa, i-1));
      gt_feature_node_add_child(gene_feature, (GtFeatureNode*) ss_feature);

      /* acceptor site */
      range = gth_sa_acceptor_site_range(sa, i-1),
      range = gt_range_offset(&range, offset);
      ss_feature = gt_feature_node_new(seqid, gt_ft_three_prime_cis_splice_site,
                                       range.start, range.end,
                                       gth_sa_gen_strand(sa));
      gt_feature_node_set_source((GtFeatureNode*) ss_feature, gthsourcetag);
      gt_feature_node_set_score((GtFeatureNode*) ss_feature,
                               gth_sa_acceptor_site_prob(sa, i-1));
      gt_feature_node_add_child(gene_feature, (GtFeatureNode*) ss_feature);
    }

    /* create exon */
    range.start = gth_sa_left_genomic_exon_border(sa, i);
    range.end   = gth_sa_right_genomic_exon_border(sa, i);
    range = gt_range_reorder(range);
    range = gt_range_offset(&range, offset);
    exon_feature = (GtFeatureNode*)
                   gt_feature_node_new(seqid, gt_ft_exon, range.start,
                                       range.end, gth_sa_gen_strand(sa));
    gt_feature_node_set_source(exon_feature, gthsourcetag);
    gt_feature_node_set_score(exon_feature,
                              gth_sa_exon_score(sa, i));
    gt_feature_node_add_child(gene_feature, exon_feature);
  }
  gt_array_add(nodes, gene_feature);
}

static void gff3_sa_visitor_visit_sa(GthSAVisitor *sa_visitor, GthSA *sa)
{
  GthGFF3SAVisitor *visitor = gff3_sa_visitor_cast(sa_visitor);
  gt_assert(sa);
  save_sa_in_gff3(sa, visitor->region_factory, visitor->nodes,
                  visitor->gthsourcetag, gth_input_md5ids(visitor->input));
}

static void gff3_sa_visitor_trailer(GthSAVisitor *sa_visitor,
                                    GT_UNUSED unsigned long num_of_sas)
{
  GthGFF3SAVisitor *visitor = gff3_sa_visitor_cast(sa_visitor);
  gt_genome_nodes_sort_stable(visitor->nodes);
  gt_genome_nodes_show(visitor->nodes, visitor->outfp);
}

const GthSAVisitorClass* gth_gff3_sa_visitor_class()
{
  static const GthSAVisitorClass savc = { sizeof (GthGFF3SAVisitor),
                                          gff3_sa_visitor_free,
                                          gff3_sa_visitor_preface,
                                          gff3_sa_visitor_visit_sa,
                                          gff3_sa_visitor_trailer };
  return &savc;
}

GthSAVisitor* gth_gff3_sa_visitor_new(GthInput *input, bool use_desc_ranges,
                                      GtFile *outfp)
{
  GthSAVisitor *sa_visitor = gth_sa_visitor_create(gth_gff3_sa_visitor_class());
  GthGFF3SAVisitor *visitor = gff3_sa_visitor_cast(sa_visitor);
  visitor->input = input;
  visitor->region_factory = gth_region_factory_new(use_desc_ranges);
  visitor->gthsourcetag = gt_str_new_cstr(GTHSOURCETAG);
  visitor->nodes = gt_array_new(sizeof (GtGenomeNode*));
  visitor->outfp = outfp;
  return sa_visitor;
}
