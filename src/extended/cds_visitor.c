/*
  Copyright (c) 2006-2010 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#include "core/assert_api.h"
#include "core/codon.h"
#include "core/codon_iterator_simple.h"
#include "core/orf.h"
#include "core/translator.h"
#include "core/undef.h"
#include "extended/cds_visitor.h"
#include "extended/node_visitor_rep.h"
#include "extended/splicedseq.h"

struct GtCDSVisitor {
  const GtNodeVisitor parent_instance;
  unsigned minorflen;
  GtStr *source;
  Splicedseq *splicedseq; /* the (spliced) sequence of the currently considered
                             gene */
  GtRegionMapping *region_mapping;
  unsigned long offset;
  bool start_codon,
       final_stop_codon;
};

static void cds_visitor_free(GtNodeVisitor *nv)
{
  GtCDSVisitor *cds_visitor = gt_cds_visitor_cast(nv);
  gt_assert(cds_visitor);
  gt_str_delete(cds_visitor->source);
  gt_splicedseq_delete(cds_visitor->splicedseq);
  gt_region_mapping_delete(cds_visitor->region_mapping);
}

static int extract_cds_if_necessary(GtGenomeNode *gn, void *data,
                                    GtError *err)
{
  GtCDSVisitor *v = (GtCDSVisitor*) data;
  GtFeatureNode *fn;
  GtRange range;
  const char *raw_sequence;
  unsigned long raw_sequence_length;
  int had_err = 0;

  gt_error_check(err);
  fn = gt_genome_node_cast(gt_feature_node_class(), gn);
  gt_assert(fn);

  if (gt_feature_node_has_type(fn, gt_ft_exon) &&
      (gt_feature_node_get_strand(fn) == GT_STRAND_FORWARD ||
       gt_feature_node_get_strand(fn) == GT_STRAND_REVERSE)) {
    range = gt_genome_node_get_range(gn);
    gt_assert(v->region_mapping);
    had_err = gt_region_mapping_get_raw_sequence(v->region_mapping,
                                                 &raw_sequence,
                                                 &raw_sequence_length,
                                                 &v->offset,
                                                 gt_genome_node_get_seqid(gn),
                                                 &range, err);
    if (!had_err) {
      gt_assert(range.start && range.end); /* 1-based coordinates */
      gt_assert(range.end - v->offset < raw_sequence_length);
      gt_splicedseq_add(v->splicedseq, range.start - v->offset,
                        range.end - v->offset, raw_sequence);
    }
  }
  return had_err;
}

static int extract_spliced_seq(GtGenomeNode *gn, GtCDSVisitor *visitor,
                               GtError *err)
{
  gt_error_check(err);
  gt_assert(gn && visitor);
  /* traverse the direct children */
  gt_splicedseq_reset(visitor->splicedseq);
  return gt_genome_node_traverse_direct_children(gn, visitor,
                                                 extract_cds_if_necessary, err);
}

static GtArray* determine_ORFs_for_all_three_frames(Splicedseq *ss,
                                                    bool start_codon,
                                                    bool final_stop_codon)
{
  GtStr *pr[3];
  GtArray *orfs;
  GtTranslator *tr;
  char translated;
  int rval;
  unsigned int frame;
  GtCodonIterator *ci;
  gt_assert(ss);
  ci = gt_codon_iterator_simple_new(gt_splicedseq_get(ss),
                                    gt_splicedseq_length(ss),
                                    NULL);

  pr[0] = gt_str_new();
  pr[1] = gt_str_new();
  pr[2] = gt_str_new();
  orfs = gt_array_new(sizeof (GtRange));

  gt_assert(ci);
  tr = gt_translator_new(ci);
  rval = gt_translator_next(tr, &translated, &frame, NULL);
  while (!rval && translated) {
    gt_str_append_char(pr[frame], translated);
    rval = gt_translator_next(tr, &translated, &frame, NULL);
  }
  gt_determine_ORFs(orfs, 0, gt_str_get(pr[0]), gt_str_length(pr[0]),
                    start_codon, final_stop_codon);
  gt_determine_ORFs(orfs, 1, gt_str_get(pr[1]), gt_str_length(pr[1]),
                    start_codon, final_stop_codon);
  gt_determine_ORFs(orfs, 2, gt_str_get(pr[2]), gt_str_length(pr[2]),
                    start_codon, final_stop_codon);

  gt_str_delete(pr[2]);
  gt_str_delete(pr[1]);
  gt_str_delete(pr[0]);
  gt_translator_delete(tr);
  gt_codon_iterator_delete(ci);

  return orfs;
}

static void create_CDS_features_for_ORF(GtRange orf, GtCDSVisitor *v,
                                        GtGenomeNode *gn)
{
  GtFeatureNode *cds_feature;
  unsigned long i;
  GtRange cds;
  GtStrand strand = gt_feature_node_get_strand((GtFeatureNode*) gn);

  gt_assert(gt_range_length(&orf) >= 3);
  /* the first CDS feature */
  cds.start = gt_splicedseq_map(v->splicedseq,strand == GT_STRAND_FORWARD
                                ? orf.start : orf.end) + v->offset;
  cds.end = gt_splicedseq_map(v->splicedseq, strand == GT_STRAND_FORWARD
                              ? orf.end : orf.start) + v->offset;
  cds_feature = (GtFeatureNode*)
                gt_feature_node_new(gt_genome_node_get_seqid(gn), gt_ft_CDS,
                                    cds.start, cds.end,
                                    gt_feature_node_get_strand((GtFeatureNode*)
                                                               gn));
  gt_feature_node_set_source(cds_feature, v->source);
  gt_feature_node_set_phase(cds_feature, GT_PHASE_ZERO);
  /* all CDS features in between */
  for (i = strand == GT_STRAND_FORWARD ? orf.start : orf.end;
       strand == GT_STRAND_FORWARD ? i < orf.end : i > orf.start;
       strand == GT_STRAND_FORWARD ? i++ : i--) {
    if (gt_splicedseq_pos_is_border(v->splicedseq, i)) {
      gt_feature_node_set_end((GtFeatureNode*) cds_feature,
                              gt_splicedseq_map(v->splicedseq, i) + v->offset);
      gt_feature_node_add_child(gt_feature_node_cast(gn), cds_feature);
      if (strand == GT_STRAND_FORWARD)
        orf.start = i + 1;
      else
        orf.end = i - 1;
      cds.start = gt_splicedseq_map(v->splicedseq, strand == GT_STRAND_FORWARD
                                    ? orf.start : orf.end) + v->offset;
      cds.end = gt_splicedseq_map(v->splicedseq, strand == GT_STRAND_FORWARD
                                  ? orf.end : orf.start) + v->offset;
      cds_feature = (GtFeatureNode*)
                    gt_feature_node_new(gt_genome_node_get_seqid(gn), gt_ft_CDS,
                                        cds.start, cds.end,
                               gt_feature_node_get_strand((GtFeatureNode*) gn));
      gt_feature_node_set_source(cds_feature, v->source);
      /* XXX correct this */
      gt_feature_node_set_phase(cds_feature, (GtPhase)
                               gt_splicedseq_map(v->splicedseq, orf.start) % 3);
    }
  }
  /* set the end of the last CDS feature and store it */
  gt_feature_node_set_end((GtFeatureNode*) cds_feature,
                          gt_splicedseq_map(v->splicedseq,
                                            strand == GT_STRAND_FORWARD
                                            ? orf.end : orf.start) + v->offset);
  gt_feature_node_add_child(gt_feature_node_cast(gn), cds_feature);
}

static void create_CDS_features_for_longest_ORF(GtArray *orfs, GtCDSVisitor *v,
                                                GtGenomeNode *gn)
{
  if (gt_array_size(orfs)) {
    /* sort ORFs according to length */
    gt_ranges_sort_by_length_stable(orfs);

    /* create CDS features from the longest ORF */
    if (gt_range_length(gt_array_get_first(orfs)) >=
        v->minorflen * GT_CODON_LENGTH) {
      create_CDS_features_for_ORF(*(GtRange*) gt_array_get_first(orfs), v, gn);
    }
  }
}

static int add_cds_if_necessary(GtGenomeNode *gn, void *data, GtError *err)
{
  GtCDSVisitor *v = (GtCDSVisitor*) data;
  GtFeatureNode *fn;
  int had_err;

  gt_error_check(err);
  fn = gt_genome_node_cast(gt_feature_node_class(), gn);
  gt_assert(fn);

  had_err = extract_spliced_seq(gn, v, err);
  if (!had_err && gt_splicedseq_length(v->splicedseq) > 2) {
    GtArray *orfs;

    if (gt_feature_node_get_strand(fn) == GT_STRAND_REVERSE) {
      if (gt_splicedseq_reverse(v->splicedseq, err))
        return -1;
    }

    orfs = determine_ORFs_for_all_three_frames(v->splicedseq, v->start_codon,
                                               v->final_stop_codon);
    create_CDS_features_for_longest_ORF(orfs, v, gn);

    gt_array_delete(orfs);
  }
  return had_err;
}

static int cds_visitor_feature_node(GtNodeVisitor *nv, GtFeatureNode *fn,
                                    GtError *err)
{
  GtCDSVisitor *v = gt_cds_visitor_cast(nv);
  gt_error_check(err);
  return gt_genome_node_traverse_children((GtGenomeNode*) fn, v,
                                          add_cds_if_necessary, false, err);
}

const GtNodeVisitorClass* gt_cds_visitor_class()
{
  static const GtNodeVisitorClass *nvc = NULL;
  if (!nvc) {
    nvc = gt_node_visitor_class_new(sizeof (GtCDSVisitor),
                                    cds_visitor_free,
                                    NULL,
                                    cds_visitor_feature_node,
                                    NULL,
                                    NULL);
  }
  return nvc;
}

GtNodeVisitor* gt_cds_visitor_new(GtRegionMapping *region_mapping,
                                  unsigned int minorflen, GtStr *source,
                                  bool start_codon, bool final_stop_codon)
{
  GtNodeVisitor *nv;
  GtCDSVisitor *cds_visitor;
  nv = gt_node_visitor_create(gt_cds_visitor_class());
  cds_visitor = gt_cds_visitor_cast(nv);
  cds_visitor->minorflen = minorflen;
  cds_visitor->source = gt_str_ref(source);
  cds_visitor->splicedseq = gt_splicedseq_new();
  cds_visitor->region_mapping = region_mapping;
  cds_visitor->start_codon = start_codon;
  cds_visitor->final_stop_codon = final_stop_codon;
  return nv;
}

void gt_cds_visitor_set_region_mapping(GtCDSVisitor *cds_visitor,
                                       GtRegionMapping *region_mapping)
{
  gt_assert(cds_visitor && region_mapping);
  gt_region_mapping_delete(cds_visitor->region_mapping);
  cds_visitor->region_mapping = gt_region_mapping_ref(region_mapping);
}
