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
#include "core/orf.h"
#include "core/translate.h"
#include "core/undef.h"
#include "extended/cds_visitor.h"
#include "extended/genome_visitor_rep.h"
#include "extended/splicedseq.h"

struct CDSVisitor {
  const GenomeVisitor parent_instance;
  Str *source;
  Splicedseq *splicedseq; /* the (spliced) sequence of the currently considered
                             gene */
  RegionMapping *region_mapping;
};

#define cds_visitor_cast(GV)\
        genome_visitor_cast(cds_visitor_class(), GV)

static void cds_visitor_free(GenomeVisitor *gv)
{
  CDSVisitor *cds_visitor = cds_visitor_cast(gv);
  assert(cds_visitor);
  str_delete(cds_visitor->source);
  splicedseq_delete(cds_visitor->splicedseq);
  region_mapping_delete(cds_visitor->region_mapping);
}

static int extract_cds_if_necessary(GT_GenomeNode *gn, void *data, GT_Error *err)
{
  CDSVisitor *v = (CDSVisitor*) data;
  GT_GenomeFeature *gf;
  GT_Range range;
  const char *raw_sequence;
  unsigned long raw_sequence_length;
  int had_err = 0;

  gt_error_check(err);
  gf = gt_genome_node_cast(gt_genome_feature_class(), gn);
  assert(gf);

  if (gt_genome_feature_has_type(gf, gft_exon) &&
      (gt_genome_feature_get_strand(gf) == GT_STRAND_FORWARD ||
       gt_genome_feature_get_strand(gf) == GT_STRAND_REVERSE)) {
    had_err = region_mapping_get_raw_sequence(v->region_mapping, &raw_sequence,
                                              gt_genome_node_get_seqid(gn), err);
    if (!had_err) {
      range = gt_genome_node_get_range(gn);
      assert(range.start && range.end); /* 1-based coordinates */
      had_err = region_mapping_get_raw_sequence_length(v->region_mapping,
                                                       &raw_sequence_length,
                                                      gt_genome_node_get_seqid(gn),
                                                       err);
    }
    if (!had_err) {
      assert(range.end <= raw_sequence_length);
      splicedseq_add(v->splicedseq, range.start - 1, range.end - 1,
                     raw_sequence);
    }
  }
  return had_err;
}

static int extract_spliced_seq(GT_GenomeNode *gn, CDSVisitor *visitor, GT_Error *err)
{
  gt_error_check(err);
  assert(gn && visitor);
  /* traverse the direct children */
  splicedseq_reset(visitor->splicedseq);
  return gt_genome_node_traverse_direct_children(gn, visitor,
                                              extract_cds_if_necessary, err);
}

static GT_Array* determine_ORFs_for_all_three_frames(Splicedseq *ss)
{
  Str *pr_0, *pr_1, *pr_2;
  GT_Array *orfs;
  assert(ss);

  pr_0 = str_new();
  pr_1 = str_new();
  pr_2 = str_new();
  orfs = gt_array_new(sizeof (GT_Range));

  translate_dna(pr_0, splicedseq_get(ss), splicedseq_length(ss), 0);
  translate_dna(pr_1, splicedseq_get(ss), splicedseq_length(ss), 1);
  translate_dna(pr_2, splicedseq_get(ss), splicedseq_length(ss), 2);
  determine_ORFs(orfs, 0, str_get(pr_0), str_length(pr_0));
  determine_ORFs(orfs, 1, str_get(pr_1), str_length(pr_1));
  determine_ORFs(orfs, 2, str_get(pr_2), str_length(pr_2));

  str_delete(pr_2);
  str_delete(pr_1);
  str_delete(pr_0);

  return orfs;
}

static void create_CDS_features_for_ORF(GT_Range orf, CDSVisitor *v,
                                        GT_GenomeNode *gn)
{
  GT_GenomeFeatureType *cds_type;
  GT_GenomeNode *cds_feature;
  unsigned long i;
  GT_Range cds;
  GT_Strand strand = gt_genome_feature_get_strand((GT_GenomeFeature*) gn);

  assert(gt_range_length(orf) >= 3);
  /* the first CDS feature */
  cds.start = splicedseq_map(v->splicedseq, strand == GT_STRAND_FORWARD
                             ? orf.start : orf.end) + 1;
  cds.end = splicedseq_map(v->splicedseq, strand == GT_STRAND_FORWARD
                           ? orf.end : orf.start) + 1;
  cds_type = gt_genome_feature_create_gft((GT_GenomeFeature*) gn, gft_CDS);
  assert(cds_type);
  cds_feature = gt_genome_feature_new(gt_genome_node_get_seqid(gn), cds_type, cds,
                                   gt_genome_feature_get_strand((GT_GenomeFeature*)
                                                             gn));
  gt_genome_feature_set_source(cds_feature, v->source);
  gt_genome_feature_set_phase(cds_feature, PHASE_ZERO);
  /* all CDS features in between */
  for (i = strand == GT_STRAND_FORWARD ? orf.start : orf.end;
       strand == GT_STRAND_FORWARD ? i < orf.end : i > orf.start;
       strand == GT_STRAND_FORWARD ? i++ : i--) {
    if (splicedseq_pos_is_border(v->splicedseq, i)) {
      gt_genome_feature_set_end((GT_GenomeFeature*) cds_feature,
                             splicedseq_map(v->splicedseq, i) + 1);
      gt_genome_node_is_part_of_genome_node(gn, cds_feature);
      if (strand == GT_STRAND_FORWARD)
        orf.start = i + 1;
      else
        orf.end = i - 1;
      cds.start = splicedseq_map(v->splicedseq, strand == GT_STRAND_FORWARD
                                 ? orf.start : orf.end) + 1;
      cds.end = splicedseq_map(v->splicedseq, strand == GT_STRAND_FORWARD
                               ? orf.end : orf.start) + 1;
      cds_feature = gt_genome_feature_new(gt_genome_node_get_seqid(gn), cds_type, cds,
                                gt_genome_feature_get_strand((GT_GenomeFeature*) gn));
      gt_genome_feature_set_source(cds_feature, v->source);
      /* XXX correct this */
      gt_genome_feature_set_phase(cds_feature, (Phase)
                               splicedseq_map(v->splicedseq, orf.start) % 3);
    }
  }
  /* set the end of the last CDS feature and store it */
  gt_genome_feature_set_end((GT_GenomeFeature*) cds_feature,
                         splicedseq_map(v->splicedseq,
                                        strand == GT_STRAND_FORWARD
                                        ? orf.end : orf.start) + 1);
  gt_genome_node_is_part_of_genome_node(gn, cds_feature);
}

static void create_CDS_features_for_longest_ORF(GT_Array *orfs, CDSVisitor *v,
                                                GT_GenomeNode *gn)
{
  if (gt_array_size(orfs)) {
    /* sort ORFs according to length */
    ranges_sort_by_length_stable(orfs);

    /* create CDS features from the longest ORF */
    create_CDS_features_for_ORF(*(GT_Range*) gt_array_get_first(orfs), v, gn);
  }
}

static int add_cds_if_necessary(GT_GenomeNode *gn, void *data, GT_Error *err)
{
  CDSVisitor *v = (CDSVisitor*) data;
  GT_GenomeFeature *gf;
  int had_err;

  gt_error_check(err);
  gf = gt_genome_node_cast(gt_genome_feature_class(), gn);
  assert(gf);

  had_err = extract_spliced_seq(gn, v, err);
  if (!had_err && splicedseq_length(v->splicedseq) > 2) {
    GT_Array *orfs;

    if (gt_genome_feature_get_strand(gf) == GT_STRAND_REVERSE) {
      if (splicedseq_reverse(v->splicedseq, err))
        return -1;
    }

    orfs = determine_ORFs_for_all_three_frames(v->splicedseq);
    create_CDS_features_for_longest_ORF(orfs, v, gn);

    gt_array_delete(orfs);
  }
  return had_err;
}

static int cds_visitor_genome_feature(GenomeVisitor *gv, GT_GenomeFeature *gf,
                                      GT_Error *err)
{
  CDSVisitor *v = cds_visitor_cast(gv);
  gt_error_check(err);
  return gt_genome_node_traverse_children((GT_GenomeNode*) gf, v,
                                       add_cds_if_necessary, false, err);

}

const GenomeVisitorClass* cds_visitor_class()
{
  static const GenomeVisitorClass gvc = { sizeof (CDSVisitor),
                                          cds_visitor_free,
                                          NULL,
                                          cds_visitor_genome_feature,
                                          NULL,
                                          NULL };
  return &gvc;
}

GenomeVisitor* cds_visitor_new(RegionMapping *region_mapping, Str *source)
{
  GenomeVisitor *gv;
  CDSVisitor *cds_visitor;
  assert(region_mapping);
  gv = genome_visitor_create(cds_visitor_class());
  cds_visitor = cds_visitor_cast(gv);
  cds_visitor->source = str_ref(source);
  cds_visitor->splicedseq = splicedseq_new();
  cds_visitor->region_mapping = region_mapping;
  return gv;
}
