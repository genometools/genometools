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
#include "libgtcore/log.h"
#include "libgtcore/queue.h"
#include "libgtcore/undef.h"
#include "libgtcore/unused.h"
#include "libgtext/csa_variable_strands.h"
#include "libgtext/csa_visitor.h"
#include "libgtext/genome_visitor_rep.h"

#define GT_CSA_SOURCE_TAG "gt csa"

struct CSAVisitor {
  const GenomeVisitor parent_instance;
  Queue *genome_node_buffer;
  unsigned long join_length;
  Array *cluster;
  GenomeFeature *buffered_feature;
  Range first_range,
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
  queue_delete(csa_visitor->genome_node_buffer);
  array_delete(csa_visitor->cluster);
  str_delete(csa_visitor->gt_csa_source_str);
}

static int csa_visitor_genome_feature(GenomeVisitor *gv, GenomeFeature *gf,
                                      UNUSED Error *err)
{
  CSAVisitor *csa_visitor;
  error_check(err);
  csa_visitor = csa_visitor_cast(gv);

  /* determine the first range if necessary */
  if (csa_visitor->buffered_feature) {
    csa_visitor->first_range = genome_node_get_range((GenomeNode*) csa_visitor
                                                     ->buffered_feature);
    csa_visitor->first_str = genome_node_get_seqid((GenomeNode*) csa_visitor
                                                   ->buffered_feature);
    assert(!array_size(csa_visitor->cluster));
    array_add(csa_visitor->cluster, csa_visitor->buffered_feature);
    csa_visitor->buffered_feature = NULL;
  }
  else if (!array_size(csa_visitor->cluster)) {
    csa_visitor->first_range = genome_node_get_range((GenomeNode*) gf);
    csa_visitor->first_str = genome_node_get_seqid((GenomeNode*) gf);
    array_add(csa_visitor->cluster, gf);
    return 0;
  }

  assert(!csa_visitor->buffered_feature);
  csa_visitor->second_range = genome_node_get_range((GenomeNode*) gf);
  csa_visitor->second_str = genome_node_get_seqid((GenomeNode*) gf);

  if ((str_cmp(csa_visitor->first_str, csa_visitor->second_str) == 0) &&
      (csa_visitor->first_range.end + csa_visitor->join_length >=
       csa_visitor->second_range.start)) {
      /* we are still in the cluster */
      array_add(csa_visitor->cluster, gf);
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

static int csa_visitor_default_func(GenomeVisitor *gv, GenomeNode *gn,
                                    UNUSED Error *err)
{
  CSAVisitor *csa_visitor;
  error_check(err);
  csa_visitor = csa_visitor_cast(gv);
  queue_add(csa_visitor->genome_node_buffer, gn);
  return 0;
}

static int csa_visitor_comment(GenomeVisitor *gv, Comment *c, Error *e)
{
  return csa_visitor_default_func(gv, (GenomeNode*) c, e);
}

static int csa_visitor_sequence_region(GenomeVisitor *gv, SequenceRegion *sr,
                                       Error *e)
{
  return csa_visitor_default_func(gv, (GenomeNode*) sr, e);
}

const GenomeVisitorClass* csa_visitor_class()
{
  static const GenomeVisitorClass gvc = { sizeof (CSAVisitor),
                                          csa_visitor_free,
                                          csa_visitor_comment,
                                          csa_visitor_genome_feature,
                                          csa_visitor_sequence_region };
  return &gvc;
}

GenomeVisitor* csa_visitor_new(unsigned long join_length)
{
  GenomeVisitor *gv = genome_visitor_create(csa_visitor_class());
  CSAVisitor *csa_visitor = csa_visitor_cast(gv);
  csa_visitor->genome_node_buffer = queue_new();
  csa_visitor->join_length = join_length;
  csa_visitor->cluster = array_new(sizeof (GenomeFeature*));
  csa_visitor->buffered_feature = NULL;
  csa_visitor->gt_csa_source_str = str_new_cstr(GT_CSA_SOURCE_TAG);
  return gv;
}

unsigned long csa_visitor_node_buffer_size(GenomeVisitor *gv)
{
  CSAVisitor *csa_visitor = csa_visitor_cast(gv);
  return queue_size(csa_visitor->genome_node_buffer);
}

GenomeNode* csa_visitor_get_node(GenomeVisitor *gv)
{
  CSAVisitor *csa_visitor;
  csa_visitor = csa_visitor_cast(gv);
  return queue_get(csa_visitor->genome_node_buffer);
}

static Range get_genomic_range(const void *sa)
{
  GenomeFeature *gf = *(GenomeFeature**) sa;
  assert(gf && genome_feature_get_type(gf) == gft_gene);
  return genome_node_get_range((GenomeNode*) gf);
}

static Strand get_strand(const void *sa)
{
  GenomeFeature *gf = *(GenomeFeature**) sa;
  assert(gf && genome_feature_get_type(gf) == gft_gene);
  return genome_feature_get_strand(gf);
}

static int save_exon(GenomeNode *gn, void *data, UNUSED Error *err)
{
  GenomeFeature *gf = (GenomeFeature*) gn;
  Array *exon_ranges = (Array*) data;
  Range range;
  error_check(err);
  assert(gf && exon_ranges);
  if (genome_feature_get_type(gf) == gft_exon) {
    range = genome_node_get_range(gn);
    array_add(exon_ranges, range);
  }
  return 0;
}

static void get_exons(Array *exon_ranges, const void *sa)
{
  GenomeFeature *gf = *(GenomeFeature**) sa;
  int had_err;
  assert(exon_ranges && gf && genome_feature_get_type(gf) == gft_gene);
  had_err = genome_node_traverse_children((GenomeNode*) gf, exon_ranges,
                                          save_exon, false, NULL);
  /* we cannot have an error here, because save_exon() doesn't produces one. */
  assert(!had_err);
  /* we got at least one exon */
  assert(array_size(exon_ranges));
  assert(ranges_are_sorted_and_do_not_overlap(exon_ranges));
}

static void add_sa_to_exon_feature_array(Array *exon_nodes,
                                         GenomeFeature *sa,
                                         Str *seqid,
                                         Str *gt_csa_source_str,
                                         Strand gene_strand)
{
  Array *exons_from_sa;
  unsigned long i,
                exon_feature_index = 0,
                exons_from_sa_index = 0;
  GenomeFeature *exon_feature, *exons_from_sa_feature;
  GenomeNode *new_feature;
  Range exon_feature_range, exons_from_sa_range;

  assert(exon_nodes && sa);
  assert(gene_strand != STRAND_BOTH); /* is defined */

  exons_from_sa = array_new(sizeof (GenomeFeature*));
  genome_feature_get_exons(sa, exons_from_sa);
  genome_nodes_sort(exons_from_sa);

  while (exon_feature_index < array_size(exon_nodes) &&
         exons_from_sa_index < array_size(exons_from_sa)) {
    exon_feature = *(GenomeFeature**)
                   array_get(exon_nodes, exon_feature_index);
    exons_from_sa_feature = *(GenomeFeature**)
                            array_get(exons_from_sa, exons_from_sa_index);

    exon_feature_range = genome_node_get_range((GenomeNode*) exon_feature);
    exons_from_sa_range =
      genome_node_get_range((GenomeNode*) exons_from_sa_feature);

    switch (range_compare(exon_feature_range, exons_from_sa_range)) {
      case -1:
        if (range_overlap(exon_feature_range, exons_from_sa_range)) {
          if (!range_contains(exon_feature_range, exons_from_sa_range)) {
            assert(genome_node_get_start((GenomeNode*) exon_feature) <=
                   genome_node_get_start((GenomeNode*) exons_from_sa_feature));
            assert(genome_node_get_end((GenomeNode*) exon_feature) <
                   genome_node_get_end((GenomeNode*) exons_from_sa_feature));
            /* update right border and score */
            genome_feature_set_end(exon_feature,
                                   genome_node_get_end((GenomeNode*)
                                                       exons_from_sa_feature));
            genome_feature_set_score(exon_feature,
                               genome_feature_get_score(exons_from_sa_feature));
          }
          exons_from_sa_index++;
        }
        exon_feature_index++;
        break;
      case 0:
        assert(range_overlap(exon_feature_range, exons_from_sa_range));
        /* update score if necessary */
        if (genome_feature_get_score(exon_feature) <
            genome_feature_get_score(exons_from_sa_feature)) {
          genome_feature_set_score(exon_feature,
                               genome_feature_get_score(exons_from_sa_feature));
        }
        exon_feature_index++;
        exons_from_sa_index++;
        break;
      case 1:
        assert(range_overlap(exon_feature_range, exons_from_sa_range));
        assert(genome_node_get_start((GenomeNode*) exon_feature) <=
               genome_node_get_start((GenomeNode*) exons_from_sa_feature));
        /* update right border and score, if necessary */
        if (genome_node_get_end((GenomeNode*) exons_from_sa_feature) >
            genome_node_get_end((GenomeNode*) exon_feature)) {
          genome_feature_set_end(exon_feature,
                                 genome_node_get_end((GenomeNode*)
                                                     exons_from_sa_feature));
          genome_feature_set_score(exon_feature,
                               genome_feature_get_score(exons_from_sa_feature));
        }
        exon_feature_index++;
        exons_from_sa_index++;
        break;
      default: assert(0);
    }
  }

  /* add remaining exons */
  for (i = exons_from_sa_index; i < array_size(exons_from_sa); i++) {
    exons_from_sa_feature = *(GenomeFeature**) array_get(exons_from_sa, i);
    new_feature =
      genome_feature_new(gft_exon,
                         genome_node_get_range((GenomeNode*)
                                               exons_from_sa_feature),
                         gene_strand, NULL, UNDEF_ULONG);
    genome_node_set_seqid(new_feature, seqid);
    genome_feature_set_score((GenomeFeature*) new_feature,
                             genome_feature_get_score(exons_from_sa_feature));
    genome_feature_set_source(new_feature, gt_csa_source_str);
    array_add(exon_nodes, new_feature);
  }

  array_delete(exons_from_sa);
}

#ifndef NDEBUG
static bool genome_nodes_are_sorted_and_do_not_overlap(Array *exon_nodes)
{
  Array *ranges = array_new(sizeof (Range));
  unsigned long i;
  Range range;
  bool rval;
  assert(exon_nodes);
  for (i = 0; i < array_size(exon_nodes); i++) {
    range = genome_node_get_range(*(GenomeNode**) array_get(exon_nodes, i));
    array_add(ranges, range);
  }
  rval = ranges_are_sorted_and_do_not_overlap(ranges);
  array_delete(ranges);
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
    GenomeNode *sa = *(GenomeNode**) csa_splice_form_get_sa(csa_splice_form, i);
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

static GenomeNode* create_mRNA_feature(CSASpliceForm *csa_splice_form,
                                       Str *gt_csa_source_str)
{
  GenomeNode *mRNA_feature;
  Array *exon_nodes;
  unsigned long i;
  Strand strand;
  Str *seqid;
  assert(csa_splice_form && gt_csa_source_str);

  /* create the mRNA feature itself */
  strand = csa_splice_form_strand(csa_splice_form),
  mRNA_feature = genome_feature_new(gft_mRNA,
                                 csa_splice_form_genomic_range(csa_splice_form),
                                    csa_splice_form_strand(csa_splice_form),
                                    NULL, UNDEF_ULONG);
  seqid = genome_node_get_seqid(*(GenomeNode**)
                               csa_splice_form_get_representative(
                                                              csa_splice_form));
  genome_node_set_seqid(mRNA_feature, seqid);
  genome_feature_set_source(mRNA_feature, gt_csa_source_str);
  mRNA_set_target_attribute((GenomeFeature*) mRNA_feature, csa_splice_form);

  /* create exon features */
  exon_nodes = array_new(sizeof (GenomeNode*));
  for (i = 0; i < csa_splice_form_num_of_sas(csa_splice_form); i++) {
    add_sa_to_exon_feature_array(exon_nodes,
                                 *(GenomeFeature**)
                                 csa_splice_form_get_sa(csa_splice_form, i),
                                 seqid, gt_csa_source_str, strand);
  }
  assert(genome_nodes_are_sorted_and_do_not_overlap(exon_nodes));

  /* add exon features to mRNA feature */
  for (i = 0; i < array_size(exon_nodes); i++) {
    genome_node_is_part_of_genome_node(mRNA_feature,
                                       *(GenomeNode**)
                                       array_get(exon_nodes, i));
  }

  array_delete(exon_nodes);

  return mRNA_feature;
}

static GenomeNode* create_gene_feature(CSAGene *csa_gene,
                                       Str *gt_csa_source_str)
{
  GenomeNode *gene_feature, *mRNA_feature;
  unsigned long i;
  assert(csa_gene && gt_csa_source_str);

  /* create top-level gene feature */
  gene_feature = genome_feature_new(gft_gene, csa_gene_genomic_range(csa_gene),
                                    csa_gene_strand(csa_gene), NULL,
                                    UNDEF_ULONG);
  genome_node_set_seqid(gene_feature,
                        genome_node_get_seqid(*(GenomeNode**)
                                              csa_gene_get_representative(
                                                                    csa_gene)));
  genome_feature_set_source(gene_feature, gt_csa_source_str);

  /* create mRNA features representing the splice forms */
  for (i = 0; i < csa_gene_num_of_splice_forms(csa_gene); i++) {
    CSASpliceForm *csa_splice_form = csa_gene_get_splice_form(csa_gene, i);
    mRNA_feature = create_mRNA_feature(csa_splice_form, gt_csa_source_str);
    genome_node_is_part_of_genome_node(gene_feature, mRNA_feature);
  }

  return gene_feature;
}

static void process_csa_genes(Queue *genome_node_buffer, Array *csa_genes,
                              Str *gt_csa_source_str)
{
  unsigned long i;
  assert(csa_genes);
  for (i = 0; i < array_size(csa_genes); i++) {
    GenomeNode *gene_feature = create_gene_feature(*(CSAGene**)
                                                   array_get(csa_genes, i),
                                                   gt_csa_source_str);
    queue_add(genome_node_buffer, gene_feature);
  }
}

void csa_visitor_process_cluster(GenomeVisitor *gv, bool final_cluster)
{
  CSAVisitor *csa_visitor = csa_visitor_cast(gv);
  Array *csa_genes;
  unsigned long i;

  if (final_cluster) {
    assert(!array_size(csa_visitor->cluster) || !csa_visitor->buffered_feature);
    if (csa_visitor->buffered_feature) {
      array_add(csa_visitor->cluster, csa_visitor->buffered_feature);
      csa_visitor->buffered_feature = NULL;
    }
  }

  if (!array_size(csa_visitor->cluster)) {
    assert(final_cluster);
    return;
  }

  /* compute the consensus spliced alignments */
  csa_genes = csa_variable_strands(array_get_space(csa_visitor->cluster),
                                   array_size(csa_visitor->cluster),
                                   sizeof (GenomeFeature*), get_genomic_range,
                                   get_strand, get_exons);

  process_csa_genes(csa_visitor->genome_node_buffer, csa_genes,
                    csa_visitor->gt_csa_source_str);

  for (i = 0; i < array_size(csa_genes); i++)
    csa_gene_delete(*(CSAGene**) array_get(csa_genes, i));
  array_delete(csa_genes);

  /* remove the cluster genome nodes */
  for (i = 0; i < array_size(csa_visitor->cluster); i++)
    genome_node_rec_delete(*(GenomeNode**) array_get(csa_visitor->cluster, i));
  array_reset(csa_visitor->cluster);
}
