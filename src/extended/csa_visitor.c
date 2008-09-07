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
#include "core/log.h"
#include "core/queue.h"
#include "core/undef.h"
#include "core/unused.h"
#include "extended/csa_variable_strands.h"
#include "extended/csa_visitor.h"
#include "extended/genome_visitor_rep.h"

#define GT_CSA_SOURCE_TAG "gt csa"

struct CSAVisitor {
  const GenomeVisitor parent_instance;
  Queue *gt_genome_node_buffer;
  unsigned long join_length;
  GT_Array *cluster;
  GenomeFeature *buffered_feature;
  GT_Range first_range,
        second_range;
  Str *first_str,
      *second_str,
      *gt_csa_source_str;
};

#define csa_visitor_cast(GV)\
        genome_visitor_cast(csa_visitor_class(), GV)

static void csa_visitor_free(GenomeVisitor *gv)
{
  CSAVisitor *csa_visitor = csa_visitor_cast(gv);
  queue_delete(csa_visitor->gt_genome_node_buffer);
  gt_array_delete(csa_visitor->cluster);
  str_delete(csa_visitor->gt_csa_source_str);
}

static int csa_visitor_genome_feature(GenomeVisitor *gv, GenomeFeature *gf,
                                      UNUSED GT_Error *err)
{
  CSAVisitor *csa_visitor;
  gt_error_check(err);
  csa_visitor = csa_visitor_cast(gv);

  /* determine the first range if necessary */
  if (csa_visitor->buffered_feature) {
    csa_visitor->first_range = gt_genome_node_get_range((GT_GenomeNode*) csa_visitor
                                                     ->buffered_feature);
    csa_visitor->first_str = gt_genome_node_get_seqid((GT_GenomeNode*) csa_visitor
                                                   ->buffered_feature);
    assert(!gt_array_size(csa_visitor->cluster));
    gt_array_add(csa_visitor->cluster, csa_visitor->buffered_feature);
    csa_visitor->buffered_feature = NULL;
  }
  else if (!gt_array_size(csa_visitor->cluster)) {
    csa_visitor->first_range = gt_genome_node_get_range((GT_GenomeNode*) gf);
    csa_visitor->first_str = gt_genome_node_get_seqid((GT_GenomeNode*) gf);
    gt_array_add(csa_visitor->cluster, gf);
    return 0;
  }

  assert(!csa_visitor->buffered_feature);
  csa_visitor->second_range = gt_genome_node_get_range((GT_GenomeNode*) gf);
  csa_visitor->second_str = gt_genome_node_get_seqid((GT_GenomeNode*) gf);

  if ((str_cmp(csa_visitor->first_str, csa_visitor->second_str) == 0) &&
      (csa_visitor->first_range.end + csa_visitor->join_length >=
       csa_visitor->second_range.start)) {
      /* we are still in the cluster */
      gt_array_add(csa_visitor->cluster, gf);
      /* update first range */
      assert(csa_visitor->second_range.start >= csa_visitor->first_range.start);
      if (csa_visitor->second_range.end > csa_visitor->first_range.end)
        csa_visitor->first_range.end = csa_visitor->second_range.end;
  }
  else {
    /* end of cluster -> process it */
    log_log("process cluster");
    csa_visitor->buffered_feature = gf;
    csa_visitor_process_cluster(gv, false);
    csa_visitor->first_range = csa_visitor->second_range;
    csa_visitor->first_str = csa_visitor->second_str;
  }
  return 0;
}

static int csa_visitor_default_func(GenomeVisitor *gv, GT_GenomeNode *gn,
                                    UNUSED GT_Error *err)
{
  CSAVisitor *csa_visitor;
  gt_error_check(err);
  csa_visitor = csa_visitor_cast(gv);
  queue_add(csa_visitor->gt_genome_node_buffer, gn);
  return 0;
}

static int csa_visitor_comment(GenomeVisitor *gv, GT_Comment *c, GT_Error *err)
{
  return csa_visitor_default_func(gv, (GT_GenomeNode*) c, err);
}

static int csa_visitor_sequence_region(GenomeVisitor *gv, SequenceRegion *sr,
                                       GT_Error *err)
{
  return csa_visitor_default_func(gv, (GT_GenomeNode*) sr, err);
}

static int csa_visitor_sequence_node(GenomeVisitor *gv, SequenceNode *sn,
                                     GT_Error *err)
{
  return csa_visitor_default_func(gv, (GT_GenomeNode*) sn, err);
}

const GenomeVisitorClass* csa_visitor_class()
{
  static const GenomeVisitorClass gvc = { sizeof (CSAVisitor),
                                          csa_visitor_free,
                                          csa_visitor_comment,
                                          csa_visitor_genome_feature,
                                          csa_visitor_sequence_region,
                                          csa_visitor_sequence_node };
  return &gvc;
}

GenomeVisitor* csa_visitor_new(unsigned long join_length)
{
  GenomeVisitor *gv = genome_visitor_create(csa_visitor_class());
  CSAVisitor *csa_visitor = csa_visitor_cast(gv);
  csa_visitor->gt_genome_node_buffer = queue_new();
  csa_visitor->join_length = join_length;
  csa_visitor->cluster = gt_array_new(sizeof (GenomeFeature*));
  csa_visitor->buffered_feature = NULL;
  csa_visitor->gt_csa_source_str = str_new_cstr(GT_CSA_SOURCE_TAG);
  return gv;
}

unsigned long csa_visitor_node_buffer_size(GenomeVisitor *gv)
{
  CSAVisitor *csa_visitor = csa_visitor_cast(gv);
  return queue_size(csa_visitor->gt_genome_node_buffer);
}

GT_GenomeNode* csa_visitor_get_node(GenomeVisitor *gv)
{
  CSAVisitor *csa_visitor;
  csa_visitor = csa_visitor_cast(gv);
  return queue_get(csa_visitor->gt_genome_node_buffer);
}

static GT_Range get_genomic_range(const void *sa)
{
  GenomeFeature *gf = *(GenomeFeature**) sa;
  assert(gf && genome_feature_has_type(gf, gft_gene));
  return gt_genome_node_get_range((GT_GenomeNode*) gf);
}

static Strand get_strand(const void *sa)
{
  GenomeFeature *gf = *(GenomeFeature**) sa;
  assert(gf && genome_feature_has_type(gf, gft_gene));
  return genome_feature_get_strand(gf);
}

static int save_exon(GT_GenomeNode *gn, void *data, UNUSED GT_Error *err)
{
  GenomeFeature *gf = (GenomeFeature*) gn;
  GT_Array *exon_ranges = (GT_Array*) data;
  GT_Range range;
  gt_error_check(err);
  assert(gf && exon_ranges);
  if (genome_feature_has_type(gf, gft_exon)) {
    range = gt_genome_node_get_range(gn);
    gt_array_add(exon_ranges, range);
  }
  return 0;
}

static void get_exons(GT_Array *exon_ranges, const void *sa)
{
  GenomeFeature *gf = *(GenomeFeature**) sa;
  int had_err;
  assert(exon_ranges && gf && genome_feature_has_type(gf, gft_gene));
  had_err = gt_genome_node_traverse_children((GT_GenomeNode*) gf, exon_ranges,
                                          save_exon, false, NULL);
  /* we cannot have an error here, because save_exon() doesn't produces one. */
  assert(!had_err);
  /* we got at least one exon */
  assert(gt_array_size(exon_ranges));
  assert(ranges_are_sorted_and_do_not_overlap(exon_ranges));
}

static void add_sa_to_exon_feature_array(GT_Array *exon_nodes,
                                         GenomeFeature *sa,
                                         Str *seqid,
                                         Str *gt_csa_source_str,
                                         Strand gene_strand,
                                         GenomeFeatureType *exon_type)
{
  GT_Array *exons_from_sa;
  unsigned long i,
                exon_feature_index = 0,
                exons_from_sa_index = 0;
  GenomeFeature *exon_feature, *exons_from_sa_feature;
  GT_GenomeNode *new_feature;
  GT_Range exon_feature_range, exons_from_sa_range;

  assert(exon_nodes && sa);
  assert(gene_strand != STRAND_BOTH); /* is defined */

  exons_from_sa = gt_array_new(sizeof (GenomeFeature*));
  genome_feature_get_exons(sa, exons_from_sa);
  genome_nodes_sort(exons_from_sa);

  while (exon_feature_index < gt_array_size(exon_nodes) &&
         exons_from_sa_index < gt_array_size(exons_from_sa)) {
    exon_feature = *(GenomeFeature**)
                   gt_array_get(exon_nodes, exon_feature_index);
    exons_from_sa_feature = *(GenomeFeature**)
                            gt_array_get(exons_from_sa, exons_from_sa_index);

    exon_feature_range = gt_genome_node_get_range((GT_GenomeNode*) exon_feature);
    exons_from_sa_range =
      gt_genome_node_get_range((GT_GenomeNode*) exons_from_sa_feature);

    switch (gt_range_compare(exon_feature_range, exons_from_sa_range)) {
      case -1:
        if (gt_range_overlap(exon_feature_range, exons_from_sa_range)) {
          if (!gt_range_contains(exon_feature_range, exons_from_sa_range)) {
            assert(gt_genome_node_get_start((GT_GenomeNode*) exon_feature) <=
                   gt_genome_node_get_start((GT_GenomeNode*) exons_from_sa_feature));
            assert(gt_genome_node_get_end((GT_GenomeNode*) exon_feature) <
                   gt_genome_node_get_end((GT_GenomeNode*) exons_from_sa_feature));
            /* update right border and score */
            genome_feature_set_end(exon_feature,
                                   gt_genome_node_get_end((GT_GenomeNode*)
                                                       exons_from_sa_feature));
            if (genome_feature_score_is_defined(exons_from_sa_feature)) {
              genome_feature_set_score(exon_feature,
                               genome_feature_get_score(exons_from_sa_feature));
            }
          }
          exons_from_sa_index++;
        }
        exon_feature_index++;
        break;
      case 0:
        assert(gt_range_overlap(exon_feature_range, exons_from_sa_range));
        /* update score if necessary */
        if ((genome_feature_score_is_defined(exon_feature) &&
             genome_feature_score_is_defined(exons_from_sa_feature) &&
             genome_feature_get_score(exon_feature) <
             genome_feature_get_score(exons_from_sa_feature)) ||
            (!genome_feature_score_is_defined(exon_feature) &&
             genome_feature_score_is_defined(exons_from_sa_feature))) {
          genome_feature_set_score(exon_feature,
                               genome_feature_get_score(exons_from_sa_feature));
        }
        exon_feature_index++;
        exons_from_sa_index++;
        break;
      case 1:
        assert(gt_range_overlap(exon_feature_range, exons_from_sa_range));
        assert(gt_genome_node_get_start((GT_GenomeNode*) exon_feature) <=
               gt_genome_node_get_start((GT_GenomeNode*) exons_from_sa_feature));
        /* update right border and score, if necessary */
        if (gt_genome_node_get_end((GT_GenomeNode*) exons_from_sa_feature) >
            gt_genome_node_get_end((GT_GenomeNode*) exon_feature)) {
          genome_feature_set_end(exon_feature,
                                 gt_genome_node_get_end((GT_GenomeNode*)
                                                     exons_from_sa_feature));
          if (genome_feature_score_is_defined(exons_from_sa_feature)) {
            genome_feature_set_score(exon_feature,
                               genome_feature_get_score(exons_from_sa_feature));
          }
        }
        exon_feature_index++;
        exons_from_sa_index++;
        break;
      default: assert(0);
    }
  }

  /* add remaining exons */
  for (i = exons_from_sa_index; i < gt_array_size(exons_from_sa); i++) {
    exons_from_sa_feature = *(GenomeFeature**) gt_array_get(exons_from_sa, i);
    new_feature =
      genome_feature_new(seqid, exon_type,
                         gt_genome_node_get_range((GT_GenomeNode*)
                                               exons_from_sa_feature),
                         gene_strand);
    if (genome_feature_score_is_defined(exons_from_sa_feature)) {
      genome_feature_set_score((GenomeFeature*) new_feature,
                               genome_feature_get_score(exons_from_sa_feature));
    }
    genome_feature_set_source(new_feature, gt_csa_source_str);
    gt_array_add(exon_nodes, new_feature);
  }

  gt_array_delete(exons_from_sa);
}

#ifndef NDEBUG
static bool genome_nodes_are_sorted_and_do_not_overlap(GT_Array *exon_nodes)
{
  GT_Array *ranges = gt_array_new(sizeof (GT_Range));
  unsigned long i;
  GT_Range range;
  bool rval;
  assert(exon_nodes);
  for (i = 0; i < gt_array_size(exon_nodes); i++) {
    range = gt_genome_node_get_range(*(GT_GenomeNode**) gt_array_get(exon_nodes, i));
    gt_array_add(ranges, range);
  }
  rval = ranges_are_sorted_and_do_not_overlap(ranges);
  gt_array_delete(ranges);
  return rval;
}
#endif

static void mRNA_set_target_attribute(GenomeFeature *mRNA_feature,
                                      const CSASpliceForm *csa_splice_form)
{
  unsigned long i;
  Str *targets;
  assert(mRNA_feature && csa_splice_form);
  targets = str_new();
  for (i = 0; i < csa_splice_form_num_of_sas(csa_splice_form); i++) {
    GT_GenomeNode *sa = *(GT_GenomeNode**) csa_splice_form_get_sa(csa_splice_form, i);
    if (genome_feature_get_attribute(sa, "Target")) {
      if (str_length(targets))
        str_append_char(targets, ',');
      str_append_cstr(targets, genome_feature_get_attribute(sa, "Target"));
    }
  }
  if (str_length(targets))
    genome_feature_add_attribute(mRNA_feature, "Target", str_get(targets));
  str_delete(targets);
}

static GT_GenomeNode* create_mRNA_feature(CSASpliceForm *csa_splice_form,
                                       Str *gt_csa_source_str,
                                       GenomeFeatureType *mRNA_type,
                                       GenomeFeatureType *exon_type)
{
  GT_GenomeNode *mRNA_feature;
  GT_Array *exon_nodes;
  unsigned long i;
  Strand strand;
  Str *seqid;
  assert(csa_splice_form && gt_csa_source_str);

  /* create the mRNA feature itself */
  strand = csa_splice_form_strand(csa_splice_form),
  seqid = gt_genome_node_get_seqid(*(GT_GenomeNode**)
                               csa_splice_form_get_representative(
                                                              csa_splice_form));
  mRNA_feature = genome_feature_new(seqid, mRNA_type,
                                 csa_splice_form_genomic_range(csa_splice_form),
                                    csa_splice_form_strand(csa_splice_form));
  genome_feature_set_source(mRNA_feature, gt_csa_source_str);
  mRNA_set_target_attribute((GenomeFeature*) mRNA_feature, csa_splice_form);

  /* create exon features */
  exon_nodes = gt_array_new(sizeof (GT_GenomeNode*));
  for (i = 0; i < csa_splice_form_num_of_sas(csa_splice_form); i++) {
    add_sa_to_exon_feature_array(exon_nodes,
                                 *(GenomeFeature**)
                                 csa_splice_form_get_sa(csa_splice_form, i),
                                 seqid, gt_csa_source_str, strand, exon_type);
  }
  assert(genome_nodes_are_sorted_and_do_not_overlap(exon_nodes));

  /* add exon features to mRNA feature */
  for (i = 0; i < gt_array_size(exon_nodes); i++) {
    gt_genome_node_is_part_of_genome_node(mRNA_feature,
                                       *(GT_GenomeNode**)
                                       gt_array_get(exon_nodes, i));
  }

  gt_array_delete(exon_nodes);

  return mRNA_feature;
}

static GT_GenomeNode* create_gene_feature(CSAGene *csa_gene,
                                       Str *gt_csa_source_str,
                                       GenomeFeatureType *gene_type,
                                       GenomeFeatureType *mRNA_type,
                                       GenomeFeatureType *exon_type)
{
  GT_GenomeNode *gene_feature, *mRNA_feature;
  unsigned long i;
  assert(csa_gene && gt_csa_source_str);

  /* create top-level gene feature */
  gene_feature = genome_feature_new(gt_genome_node_get_seqid(*(GT_GenomeNode**)
                                         csa_gene_get_representative(csa_gene)),
                                    gene_type, csa_gene_genomic_range(csa_gene),
                                    csa_gene_strand(csa_gene));
  genome_feature_set_source(gene_feature, gt_csa_source_str);

  /* create mRNA features representing the splice forms */
  for (i = 0; i < csa_gene_num_of_splice_forms(csa_gene); i++) {
    CSASpliceForm *csa_splice_form = csa_gene_get_splice_form(csa_gene, i);
    mRNA_feature = create_mRNA_feature(csa_splice_form, gt_csa_source_str,
                                       mRNA_type, exon_type);
    gt_genome_node_is_part_of_genome_node(gene_feature, mRNA_feature);
  }

  return gene_feature;
}

static void process_csa_genes(Queue *gt_genome_node_buffer, GT_Array *csa_genes,
                              Str *gt_csa_source_str,
                              GenomeFeatureType *gene_type,
                              GenomeFeatureType *mRNA_type,
                              GenomeFeatureType *exon_type)
{
  unsigned long i;
  assert(csa_genes);
  for (i = 0; i < gt_array_size(csa_genes); i++) {
    GT_GenomeNode *gene_feature = create_gene_feature(*(CSAGene**)
                                                   gt_array_get(csa_genes, i),
                                                   gt_csa_source_str, gene_type,
                                                   mRNA_type, exon_type);
    queue_add(gt_genome_node_buffer, gene_feature);
  }
}

void csa_visitor_process_cluster(GenomeVisitor *gv, bool final_cluster)
{
  CSAVisitor *csa_visitor = csa_visitor_cast(gv);
  GenomeFeature *first_feature;
  GenomeFeatureType *gene_type, *mRNA_type, *exon_type;
  GT_Array *csa_genes;
  unsigned long i;

  if (final_cluster) {
    assert(!gt_array_size(csa_visitor->cluster) ||
           !csa_visitor->buffered_feature);
    if (csa_visitor->buffered_feature) {
      gt_array_add(csa_visitor->cluster, csa_visitor->buffered_feature);
      csa_visitor->buffered_feature = NULL;
    }
  }

  if (!gt_array_size(csa_visitor->cluster)) {
    assert(final_cluster);
    return;
  }

  /* compute the consensus spliced alignments */
  first_feature = *(GenomeFeature**) gt_array_get_first(csa_visitor->cluster);
  csa_genes = csa_variable_strands(gt_array_get_space(csa_visitor->cluster),
                                   gt_array_size(csa_visitor->cluster),
                                   sizeof (GenomeFeature*), get_genomic_range,
                                   get_strand, get_exons);

  gene_type = genome_feature_create_gft(first_feature, gft_gene);
  mRNA_type = genome_feature_create_gft(first_feature, gft_mRNA);
  exon_type = genome_feature_create_gft(first_feature, gft_exon);
  process_csa_genes(csa_visitor->gt_genome_node_buffer, csa_genes,
                    csa_visitor->gt_csa_source_str, gene_type, mRNA_type,
                    exon_type);

  for (i = 0; i < gt_array_size(csa_genes); i++)
    csa_gene_delete(*(CSAGene**) gt_array_get(csa_genes, i));
  gt_array_delete(csa_genes);

  /* remove the cluster genome nodes */
  for (i = 0; i < gt_array_size(csa_visitor->cluster); i++)
    gt_genome_node_rec_delete(*(GT_GenomeNode**) gt_array_get(csa_visitor->cluster, i));
  gt_array_reset(csa_visitor->cluster);
}
