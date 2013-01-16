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
#include "extended/cds_visitor.h"
#include "extended/feature_type.h"
#include "extended/genome_node.h"
#include "gth/ags.h"
#include "gth/gff3_pgl_visitor.h"
#include "gth/pgl_visitor_rep.h"
#include "gth/region_factory.h"

struct GthGFF3PGLVisitor {
  const GthPGLVisitor parent_instance;
  GthInput *input;
  GthRegionFactory *region_factory;
  GtStr *gthsourcetag;
  GtNodeVisitor *cds_visitor;
  GtArray *nodes;
  GtFile *outfp;
};

#define gff3_pgl_visitor_cast(GV)\
        gth_pgl_visitor_cast(gth_gff3_pgl_visitor_class(), GV)

static void gff3_pgl_visitor_free(GthPGLVisitor *pgl_visitor)
{
  unsigned long i;
  GthGFF3PGLVisitor *visitor = gff3_pgl_visitor_cast(pgl_visitor);
  for (i = 0; i < gt_array_size(visitor->nodes); i++)
    gt_genome_node_delete(*(GtGenomeNode**) gt_array_get(visitor->nodes, i));
  gt_array_delete(visitor->nodes);
  gt_node_visitor_delete(visitor->cds_visitor);
  gt_str_delete(visitor->gthsourcetag);
  gth_region_factory_delete(visitor->region_factory);
}

static void gff3_pgl_visitor_preface(GthPGLVisitor *pgl_visitor,
                                     GT_UNUSED unsigned long num_of_pgls)
{
  GthGFF3PGLVisitor *visitor = gff3_pgl_visitor_cast(pgl_visitor);
  gth_region_factory_save(visitor->region_factory, visitor->nodes,
                          visitor->input);
}

static void gff3_pgl_visitor_set_region_mapping(GthPGLVisitor *pgl_visitor,
                                                GtRegionMapping *region_mapping)
{
  GthGFF3PGLVisitor *visitor;
  gt_assert(pgl_visitor && region_mapping);
  visitor = gff3_pgl_visitor_cast(pgl_visitor);
  gt_cds_visitor_set_region_mapping(gt_cds_visitor_cast(visitor->cds_visitor),
                                    region_mapping);
}

static void add_target_attributes(GtFeatureNode *mrna_feature, GthAGS *ags,
                                  bool md5ids)
{
  GtStr *target_attribute;
  GthSA *sa;
  unsigned long i;
  gt_assert(mrna_feature && ags);
  target_attribute = gt_str_new();
  for (i = 0; i < gt_array_size(ags->alignments); i++) {
    sa = *(GthSA**) gt_array_get(ags->alignments, i);
    if (strlen(gth_sa_gff3_target_attribute(sa, md5ids))) {
      if (gt_str_length(target_attribute))
        gt_str_append_char(target_attribute, ',');
      gt_str_append_cstr(target_attribute,
                         gth_sa_gff3_target_attribute(sa, md5ids));
    }
  }
  if (gt_str_length(target_attribute)) {
    gt_feature_node_add_attribute(mrna_feature, "Target",
                                  gt_str_get(target_attribute));
  }
  gt_str_delete(target_attribute);
}

static void save_pgl_in_gff3(GthPGL *pgl, GthRegionFactory *region_factory,
                             GtNodeVisitor *cds_visitor, GtArray *nodes,
                             GtStr *gthsourcetag, bool md5ids)
{
  GthExonAGS *exon, *first_exon, *last_exon;
  GtFeatureNode *gene_feature, *mrna_feature, *exon_feature;
  unsigned long i, j;
  GtRange range;
  GtStr *seqid;
  GT_UNUSED int had_err;
  struct GthAGS *ags;
  long offset;

  gt_assert(pgl && region_factory && cds_visitor && nodes);

  seqid = gth_region_factory_get_seqid(region_factory, gth_pgl_filenum(pgl),
                                       gth_pgl_seqnum(pgl));
  offset = gth_region_factory_offset(region_factory, gth_pgl_filenum(pgl),
                                     gth_pgl_seqnum(pgl)) - 1;

  /* create gene feature */
  range = gth_pgl_genomic_range(pgl);
  range = gt_range_offset(&range, offset);
  gene_feature = (GtFeatureNode*)
                 gt_feature_node_new(seqid, gt_ft_gene, range.start, range.end,
                                     gth_pgl_genomic_strand(pgl));
  gt_feature_node_set_source(gene_feature, gthsourcetag);

  for (i = 0; i < gth_pgl_num_of_ags(pgl); i++) {
    ags = gth_pgl_get_ags(pgl, i);
    /* create mRNA feature */
    first_exon = (GthExonAGS*) gt_array_get_first(ags->exons);
    last_exon  = (GthExonAGS*) gt_array_get_last(ags->exons);
    range.start = gth_pgl_is_forward(pgl)
                  ? SHOWGENPOSAGS(first_exon->range.start)
                  : SHOWGENPOSAGS(last_exon->range.end);
    range.end   = gth_pgl_is_forward(pgl)
                  ? SHOWGENPOSAGS(last_exon->range.end)
                  : SHOWGENPOSAGS(first_exon->range.start);
    gt_assert(range.start <= range.end);
    range = gt_range_offset(&range, offset);
    mrna_feature = (GtFeatureNode*)
                   gt_feature_node_new(seqid, gt_ft_mRNA, range.start,
                                       range.end, gth_ags_genomic_strand(ags));
    gt_feature_node_set_source(mrna_feature, gthsourcetag);
    add_target_attributes(mrna_feature, ags, md5ids);
    gt_feature_node_add_child(gene_feature, mrna_feature);

    for (j = 0; j < gt_array_size(ags->exons); j++) {
      /* create splice sites, if necessary */
      if (j > 0) {
        GtGenomeNode *ss_feature;
        /* donor site */
        range = gth_ags_donor_site_range(ags, j-1),
        range = gt_range_offset(&range, offset);
        ss_feature = gt_feature_node_new(seqid,
                                         gt_ft_five_prime_cis_splice_site,
                                         range.start, range.end,
                                         gth_ags_genomic_strand(ags));
        gt_feature_node_set_source((GtFeatureNode*) ss_feature, gthsourcetag);
        gt_feature_node_set_score((GtFeatureNode*) ss_feature,
                                  gth_ags_donor_site_prob(ags, j-1));
        gt_feature_node_add_child(mrna_feature, (GtFeatureNode*) ss_feature);

        /* acceptor site */
        range = gth_ags_acceptor_site_range(ags, j-1),
        range = gt_range_offset(&range, offset);
        ss_feature = gt_feature_node_new(seqid,
                                         gt_ft_three_prime_cis_splice_site,
                                         range.start, range.end,
                                         gth_ags_genomic_strand(ags));
        gt_feature_node_set_source((GtFeatureNode*) ss_feature, gthsourcetag);
        gt_feature_node_set_score((GtFeatureNode*) ss_feature,
                                  gth_ags_acceptor_site_prob(ags, j-1));
        gt_feature_node_add_child(mrna_feature, (GtFeatureNode*) ss_feature);
      }

      /* create exon */
      exon = gt_array_get(ags->exons, j);
      range.start = gth_pgl_is_forward(pgl)
                    ? SHOWGENPOSAGS(exon->range.start)
                    : SHOWGENPOSAGS(exon->range.end);
      range.end = gth_pgl_is_forward(pgl)
                  ? SHOWGENPOSAGS(exon->range.end)
                  : SHOWGENPOSAGS(exon->range.start);
      gt_assert(range.start <= range.end);
      range = gt_range_offset(&range, offset);
      exon_feature = (GtFeatureNode*)
                     gt_feature_node_new(seqid, gt_ft_exon, range.start,
                                         range.end,
                                         gth_ags_genomic_strand(ags));
      gt_feature_node_set_source(exon_feature, gthsourcetag);
      gt_feature_node_set_score(exon_feature, exon->score);
      gt_feature_node_add_child(mrna_feature, exon_feature);
    }
  }
  had_err = gt_genome_node_accept((GtGenomeNode*) gene_feature, cds_visitor,
                                  NULL);
  gt_assert(!had_err); /* should not happen */
  gt_array_add(nodes, gene_feature);
}

static void gff3_pgl_visitor_visit_pgl(GthPGLVisitor *pgl_visitor,
                                       GthPGL *pgl,
                                       GT_UNUSED unsigned long pglnum)
{
  GthGFF3PGLVisitor *visitor = gff3_pgl_visitor_cast(pgl_visitor);
  gt_assert(pgl);
  save_pgl_in_gff3(pgl, visitor->region_factory, visitor->cds_visitor,
                   visitor->nodes, visitor->gthsourcetag,
                   gth_input_md5ids(visitor->input));
}

static void gff3_pgl_visitor_trailer(GthPGLVisitor *pgl_visitor)
{
  GthGFF3PGLVisitor *visitor = gff3_pgl_visitor_cast(pgl_visitor);
  gt_genome_nodes_sort_stable(visitor->nodes);
  gt_genome_nodes_show(visitor->nodes, visitor->outfp);
}

const GthPGLVisitorClass* gth_gff3_pgl_visitor_class()
{
  static const GthPGLVisitorClass pglvc = { sizeof (GthGFF3PGLVisitor),
                                            gff3_pgl_visitor_free,
                                            gff3_pgl_visitor_preface,
                                            gff3_pgl_visitor_set_region_mapping,
                                            gff3_pgl_visitor_visit_pgl,
                                            gff3_pgl_visitor_trailer };
  return &pglvc;
}

GthPGLVisitor* gth_gff3_pgl_visitor_new(GthInput *input, bool use_desc_ranges,
                                        unsigned long minORFlength,
                                        bool start_codon, bool final_stop_codon,
                                        GtFile *outfp)
{
  GthPGLVisitor *pgl_visitor =
    gth_pgl_visitor_create(gth_gff3_pgl_visitor_class());
  GthGFF3PGLVisitor *visitor = gff3_pgl_visitor_cast(pgl_visitor);
  visitor->input = input;
  visitor->region_factory = gth_region_factory_new(use_desc_ranges);
  visitor->gthsourcetag = gt_str_new_cstr(GTHSOURCETAG);
  visitor->cds_visitor = gt_cds_visitor_new(NULL, minORFlength,
                                            visitor->gthsourcetag, start_codon,
                                            final_stop_codon, false /* XXX */);
  visitor->nodes = gt_array_new(sizeof (GtGenomeNode*));
  visitor->outfp = outfp;
  return pgl_visitor;
}
