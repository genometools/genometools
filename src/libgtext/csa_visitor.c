/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg

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
#include "libgtcore/queue.h"
#include "libgtcore/undef.h"
#include "libgtext/consensus_sa.h"
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
      *gth_csa_source_str;
};

typedef struct {
  bool is_first_splice_form;
  GenomeNode* gene_feature;
  Str *seqid,
      *gth_csa_source_str;
  Strand gene_strand;
} Process_splice_form_info;

#define csa_visitor_cast(GV)\
        genome_visitor_cast(csa_visitor_class(), GV)

static void csa_visitor_free(GenomeVisitor *gv, Env *env)
{
  CSAVisitor *csa_visitor = csa_visitor_cast(gv);
  queue_delete(csa_visitor->genome_node_buffer, env);
  array_delete(csa_visitor->cluster, env);
  str_delete(csa_visitor->gth_csa_source_str, env);
}

static int csa_visitor_genome_feature(GenomeVisitor *gv, GenomeFeature *gf,
                                      Env *env)
{
  CSAVisitor *csa_visitor;
  env_error_check(env);
  csa_visitor = csa_visitor_cast(gv);

  /* determine the first range if necessary */
  if (csa_visitor->buffered_feature) {
    csa_visitor->first_range = genome_node_get_range((GenomeNode*) csa_visitor
                                                     ->buffered_feature);
    csa_visitor->first_str = genome_node_get_seqid((GenomeNode*) csa_visitor
                                                   ->buffered_feature);
    assert(!array_size(csa_visitor->cluster));
    array_add(csa_visitor->cluster, csa_visitor->buffered_feature, env);
    csa_visitor->buffered_feature = NULL;
  }
  else if (!array_size(csa_visitor->cluster)) {
    csa_visitor->first_range = genome_node_get_range((GenomeNode*) gf);
    csa_visitor->first_str = genome_node_get_seqid((GenomeNode*) gf);
    array_add(csa_visitor->cluster, gf, env);
    return 0;
  }

  assert(!csa_visitor->buffered_feature);
  csa_visitor->second_range = genome_node_get_range((GenomeNode*) gf);
  csa_visitor->second_str = genome_node_get_seqid((GenomeNode*) gf);

  if ((str_cmp(csa_visitor->first_str, csa_visitor->second_str) == 0) &&
      (csa_visitor->first_range.end + csa_visitor->join_length >=
       csa_visitor->second_range.start)) {
      /* we are still in the cluster */
      array_add(csa_visitor->cluster, gf, env);
      /* update first range */
      assert(csa_visitor->second_range.start >= csa_visitor->first_range.start);
      if (csa_visitor->second_range.end > csa_visitor->first_range.end)
        csa_visitor->first_range.end = csa_visitor->second_range.end;
  }
  else {
    /* end of cluster -> process it */
    env_log_log(env, "process cluster");
    csa_visitor->buffered_feature = gf;
    csa_visitor_process_cluster(gv, false, env);
    csa_visitor->first_range = csa_visitor->second_range;
    csa_visitor->first_str = csa_visitor->second_str;
  }
  return 0;
}

static int csa_visitor_default_func(GenomeVisitor *gv, GenomeNode *gn, Env *env)
{
  CSAVisitor *csa_visitor;
  env_error_check(env);
  csa_visitor = csa_visitor_cast(gv);
  queue_add(csa_visitor->genome_node_buffer, gn, env);
  return 0;
}

static int csa_visitor_comment(GenomeVisitor *gv, Comment *c, Env *env)
{
  return csa_visitor_default_func(gv, (GenomeNode*) c, env);
}

static int csa_visitor_sequence_region(GenomeVisitor *gv, SequenceRegion *sr,
                                       Env *env)
{
  return csa_visitor_default_func(gv, (GenomeNode*) sr, env);
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

GenomeVisitor* csa_visitor_new(unsigned long join_length, Env *env)
{
  GenomeVisitor *gv = genome_visitor_create(csa_visitor_class(), env);
  CSAVisitor *csa_visitor = csa_visitor_cast(gv);
  csa_visitor->genome_node_buffer = queue_new(sizeof (GenomeNode*), env);
  csa_visitor->join_length = join_length;
  csa_visitor->cluster = array_new(sizeof (GenomeFeature*), env);
  csa_visitor->buffered_feature = NULL;
  csa_visitor->gth_csa_source_str = str_new_cstr(GT_CSA_SOURCE_TAG, env);
  return gv;
}

unsigned long csa_visitor_node_buffer_size(GenomeVisitor *gv)
{
  CSAVisitor *csa_visitor = csa_visitor_cast(gv);
  return queue_size(csa_visitor->genome_node_buffer);
}

GenomeNode* csa_visitor_get_node(GenomeVisitor *gv)
{
  CSAVisitor *csa_visitor = csa_visitor_cast(gv);
  return *(GenomeNode**) queue_get(csa_visitor->genome_node_buffer);
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

static int save_exon(GenomeNode *gn, void *data, Env *env)
{
  GenomeFeature *gf = (GenomeFeature*) gn;
  Array *exon_ranges = (Array*) data;
  Range range;
  env_error_check(env);
  assert(gf && exon_ranges);
  if (genome_feature_get_type(gf) == gft_exon) {
    range = genome_node_get_range(gn);
    array_add(exon_ranges, range, env);
  }
  return 0;
}

static void get_exons(Array *exon_ranges, const void *sa, Env *env)
{
  GenomeFeature *gf = *(GenomeFeature**) sa;
  int had_err;
  assert(exon_ranges && gf && genome_feature_get_type(gf) == gft_gene);
  had_err = genome_node_traverse_children((GenomeNode*) gf, exon_ranges,
                                          save_exon, false, env);
  /* we cannot have an error here, because save_exon() doesn't produces one. */
  assert(!had_err);
  /* we got at least one exon */
  assert(array_size(exon_ranges));
  assert(ranges_are_sorted_and_do_not_overlap(exon_ranges));
}

static void add_sa_to_exon_feature_array(Array *exon_nodes,
                                         GenomeFeature *sa,
                                         Str *seqid,
                                         Str *gth_csa_source_str,
                                         Strand gene_strand, Env *env)
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

  exons_from_sa = array_new(sizeof (GenomeFeature*), env);
  genome_feature_get_exons(sa, exons_from_sa, env);
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
                         gene_strand, NULL, UNDEF_ULONG, env);
    genome_node_set_seqid(new_feature, seqid, env);
    genome_feature_set_score((GenomeFeature*) new_feature,
                             genome_feature_get_score(exons_from_sa_feature));
    genome_node_set_source(new_feature, gth_csa_source_str);
    array_add(exon_nodes, new_feature, env);
  }

  array_delete(exons_from_sa, env);
}

#ifndef NDEBUG
static bool genome_nodes_are_sorted_and_do_not_overlap(Array *exon_nodes,
                                                       Env *env)
{
  Array *ranges = array_new(sizeof (Range), env);
  unsigned long i;
  Range range;
  bool rval;
  assert(exon_nodes);
  for (i = 0; i < array_size(exon_nodes); i++) {
    range = genome_node_get_range(*(GenomeNode**) array_get(exon_nodes, i));
    array_add(ranges, range, env);
  }
  rval = ranges_are_sorted_and_do_not_overlap(ranges);
  array_delete(ranges, env);
  return rval;
}
#endif

static void process_splice_form(Array *spliced_alignments_in_form,
                                const void *set_of_sas,
                                unsigned long number_of_sas,
                                size_t size_of_sa,
                                void *userdata, Env *env)
{
  Process_splice_form_info *info = (Process_splice_form_info*) userdata;
  GenomeNode *mRNA_feature, *sa_node;
  Array *exon_nodes;
  Range gene_range, mRNA_range, tmp_range;
  unsigned long i, sa_num;

  if (info->is_first_splice_form) {
    assert(!info->gene_feature);
    /* compute gene range */
    gene_range.start = ~0UL;
    gene_range.end   = 0UL;
    /* for the gene range we have to iterate over all spliced alignments */
    for (i = 0; i < number_of_sas; i++) {
      sa_node = *(GenomeNode**) (set_of_sas + i * size_of_sa);

      /* gene range */
      tmp_range = genome_node_get_range(sa_node);
      if (tmp_range.start < gene_range.start)
        gene_range.start = tmp_range.start;
      if (tmp_range.end > gene_range.end)
        gene_range.end = tmp_range.end;

      /* determine seqid */
      if (!i)
        info->seqid = genome_node_get_seqid(*(GenomeNode**) (set_of_sas));
      else {
        assert(!str_cmp(info->seqid,
                        genome_node_get_seqid(*(GenomeNode**)
                                              (set_of_sas + i * size_of_sa))));
      }
    }

    /* determine strand (we have to iterate over all spliced alignments in this
       splice form */
    for (i = 0; i < array_size(spliced_alignments_in_form); i++) {
      sa_node = *(GenomeNode**)
                (set_of_sas +
                 *(unsigned long*) array_get(spliced_alignments_in_form, i) *
                 size_of_sa);
      info->gene_strand =
        strand_join(info->gene_strand,
                    genome_feature_get_strand((GenomeFeature*) sa_node));
    }

    assert(info->gene_strand != STRAND_BOTH);
    info->gene_feature = genome_feature_new(gft_gene, gene_range,
                                            info->gene_strand, NULL,
                                            UNDEF_ULONG, env);
    genome_node_set_seqid(info->gene_feature, info->seqid, env);
    genome_node_set_source(info->gene_feature, info->gth_csa_source_str);
    info->is_first_splice_form = false;
  }

  exon_nodes = array_new(sizeof (GenomeNode*), env);

  for (i = 0; i < array_size(spliced_alignments_in_form); i++) {
    sa_num = *(unsigned long*) array_get(spliced_alignments_in_form, i);
    assert(sa_num < number_of_sas);
    add_sa_to_exon_feature_array(exon_nodes,
                                 *(GenomeFeature**)
                                 (set_of_sas + sa_num * size_of_sa),
                                 info->seqid, info->gth_csa_source_str,
                                 info->gene_strand, env);
  }
  assert(genome_nodes_are_sorted_and_do_not_overlap(exon_nodes, env));

  mRNA_range.start = genome_node_get_start(*(GenomeNode**)
                                           array_get(exon_nodes, 0));
  mRNA_range.end   = genome_node_get_end(*(GenomeNode**)
                                         array_get(exon_nodes,
                                         array_size(exon_nodes) - 1));
  assert(info->gene_strand != STRAND_BOTH);
  mRNA_feature = genome_feature_new(gft_mRNA, mRNA_range, info->gene_strand,
                                    NULL, UNDEF_ULONG, env);
  genome_node_set_seqid(mRNA_feature, info->seqid, env);
  genome_node_set_source(mRNA_feature, info->gth_csa_source_str);
  genome_node_is_part_of_genome_node(info->gene_feature, mRNA_feature, env);
  for (i = 0; i < array_size(exon_nodes); i++) {
    genome_node_is_part_of_genome_node(mRNA_feature,
                                       *(GenomeNode**)
                                       array_get(exon_nodes, i), env);
  }

  array_delete(exon_nodes, env);
}

void csa_visitor_process_cluster(GenomeVisitor *gv, bool final_cluster,
                                 Env *env)
{
  CSAVisitor *csa_visitor = csa_visitor_cast(gv);
  Process_splice_form_info info;
  unsigned long i;

  if (final_cluster) {
    assert(!array_size(csa_visitor->cluster) || !csa_visitor->buffered_feature);
    if (csa_visitor->buffered_feature) {
      array_add(csa_visitor->cluster, csa_visitor->buffered_feature, env);
      csa_visitor->buffered_feature = NULL;
    }
  }

  /* compute the consenus spliced alignments */
  info.is_first_splice_form = true;
  info.gene_feature = NULL;
  info.seqid = NULL;
  info.gth_csa_source_str = csa_visitor->gth_csa_source_str;
  info.gene_strand = STRAND_BOTH; /* undefined */
  if (!array_size(csa_visitor->cluster)) {
    assert(final_cluster);
    return;
  }
  consensus_sa(array_get_space(csa_visitor->cluster),
               array_size(csa_visitor->cluster),
               sizeof (GenomeFeature*),
               get_genomic_range,
               get_strand,
               get_exons,
               process_splice_form,
               &info,
               env);
  assert(info.gene_feature);
  queue_add(csa_visitor->genome_node_buffer, info.gene_feature, env);

  /* remove the cluster genome nodes */
  for (i = 0; i < array_size(csa_visitor->cluster); i++) {
    genome_node_rec_delete(*(GenomeNode**) array_get(csa_visitor->cluster, i),
                           env);
  }
  array_reset(csa_visitor->cluster);
}
