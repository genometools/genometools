/*
  Copyright (c) 2006-2011 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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
#include "core/class_alloc_lock.h"
#include "core/codon_api.h"
#include "core/codon_iterator_simple_api.h"
#include "core/ma.h"
#include "core/orf.h"
#include "core/translator.h"
#include "core/trans_table.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "extended/cds_visitor.h"
#include "extended/feature_node.h"
#include "extended/feature_type.h"
#include "extended/node_visitor_api.h"
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
       final_stop_codon,
       generic_start_codons;
};

static void cds_visitor_free(GtNodeVisitor *nv)
{
  GtCDSVisitor *cds_visitor = gt_cds_visitor_cast(nv);
  gt_assert(cds_visitor);
  gt_str_delete(cds_visitor->source);
  gt_splicedseq_delete(cds_visitor->splicedseq);
  gt_region_mapping_delete(cds_visitor->region_mapping);
}

static int extract_cds_if_necessary(GtFeatureNode *fn, void *data,
                                    GtError *err)
{
  GtCDSVisitor *v = (GtCDSVisitor*) data;
  GtRange range;
  char *outsequence;
  int had_err = 0;

  gt_error_check(err);
  gt_assert(fn);

  if (gt_feature_node_has_type(fn, gt_ft_exon) &&
      (gt_feature_node_get_strand(fn) == GT_STRAND_FORWARD ||
       gt_feature_node_get_strand(fn) == GT_STRAND_REVERSE)) {
    range = gt_genome_node_get_range((GtGenomeNode*) fn);
    gt_assert(v->region_mapping);
    had_err = gt_region_mapping_get_sequence(v->region_mapping, &outsequence,
                                   gt_genome_node_get_seqid((GtGenomeNode*) fn),
                                   range.start, range.end, err);
    if (!had_err) {
      gt_assert(range.start && range.end); /* 1-based coordinates */
      gt_splicedseq_add(v->splicedseq, range.start, range.end, outsequence);
      gt_free(outsequence);
    }
  }
  return had_err;
}

static int extract_spliced_seq(GtFeatureNode *fn, GtCDSVisitor *visitor,
                               GtError *err)
{
  gt_error_check(err);
  gt_assert(fn && visitor);
  /* traverse the direct children */
  gt_splicedseq_reset(visitor->splicedseq);
  return gt_feature_node_traverse_direct_children(fn, visitor,
                                                  extract_cds_if_necessary,
                                                  err);
}

static void save_orf(void *data, GtRange *orf, GT_UNUSED unsigned long framenum,
                     GT_UNUSED const char *frame,
                     GT_UNUSED bool ends_with_stop_codon)
{
  GtArray *ranges = data;
  gt_assert(ranges && orf);
  gt_array_add(ranges, *orf);
}

static GtArray* determine_ORFs_for_all_three_frames(Splicedseq *ss,
                                                    bool start_codon,
                                                    bool final_stop_codon,
                                                    bool generic_start_codons)
{
  GtTranslatorStatus status;
  GtStr *pr[3], *start_codons[3];
  GtArray *orfs = NULL;
  GtTranslator *tr;
  char translated;
  int had_err = 0;
  unsigned int frame;
  GtCodonIterator *ci;
  bool start;
  gt_assert(ss);
  ci = gt_codon_iterator_simple_new(gt_splicedseq_get(ss),
                                    gt_splicedseq_length(ss),
                                    NULL);

  pr[0] = gt_str_new();
  pr[1] = gt_str_new();
  pr[2] = gt_str_new();
  start_codons[0] = generic_start_codons ? gt_str_new() : NULL;
  start_codons[1] = generic_start_codons ? gt_str_new() : NULL;
  start_codons[2] = generic_start_codons ? gt_str_new() : NULL;

  gt_assert(ci);
  tr = gt_translator_new(ci);
  status = gt_translator_next_with_start(tr, &translated, &frame, &start, NULL);
  while (status == GT_TRANSLATOR_OK) {
    gt_str_append_char(pr[frame], translated);
    if (generic_start_codons)
      gt_str_append_char(start_codons[frame], start ? GT_START_AMINO : '-');
    status = gt_translator_next_with_start(tr, &translated, &frame, &start,
                                           NULL);
  }
  if (status == GT_TRANSLATOR_ERROR)
    had_err = -1;

  if (!had_err) {
    orfs = gt_array_new(sizeof (GtRange));
    gt_determine_ORFs(save_orf, orfs, 0, gt_str_get(pr[0]),
                      gt_str_length(pr[0]), start_codon, final_stop_codon,
                      false, generic_start_codons ? gt_str_get(start_codons[0])
                                                  : NULL);
    gt_determine_ORFs(save_orf, orfs, 1, gt_str_get(pr[1]),
                      gt_str_length(pr[1]), start_codon, final_stop_codon,
                      false, generic_start_codons ? gt_str_get(start_codons[1])
                                                  : NULL);
    gt_determine_ORFs(save_orf, orfs, 2, gt_str_get(pr[2]),
                      gt_str_length(pr[2]), start_codon, final_stop_codon,
                      false, generic_start_codons ? gt_str_get(start_codons[2])
                                                  : NULL);
  }

  gt_str_delete(start_codons[2]);
  gt_str_delete(start_codons[1]);
  gt_str_delete(start_codons[0]);
  gt_str_delete(pr[2]);
  gt_str_delete(pr[1]);
  gt_str_delete(pr[0]);
  gt_translator_delete(tr);
  gt_codon_iterator_delete(ci);

  return orfs;
}

static void set_phases(GtArray *cds_features)
{
  GtPhase phase = GT_PHASE_ZERO;
  GtFeatureNode *cds_feature;
  unsigned long i, length;
  gt_assert(cds_features && gt_array_size(cds_features));
  for (i = 0; i < gt_array_size(cds_features); i++) {
    cds_feature = *(GtFeatureNode**) gt_array_get(cds_features, i);
    gt_feature_node_set_phase(cds_feature, phase);
    length = gt_genome_node_get_length((GtGenomeNode*) cds_feature);
    phase = (3 - (length - phase) % 3) % 3;
    gt_assert(phase <= GT_PHASE_TWO);
  }
}

static void create_CDS_features_for_ORF(GtRange orf, GtCDSVisitor *v,
                                        GtFeatureNode *fn)
{
  GtFeatureNode *representative, *cds_feature;
  GtArray *cds_features;
  unsigned long i;
  GtRange cds;
  GtStrand strand = gt_feature_node_get_strand(fn);

  gt_assert(gt_range_length(&orf) >= 3);
  /* the first CDS feature */
  cds.start = gt_splicedseq_map(v->splicedseq,strand == GT_STRAND_FORWARD
                                ? orf.start : orf.end) + v->offset;
  cds.end = gt_splicedseq_map(v->splicedseq, strand == GT_STRAND_FORWARD
                              ? orf.end : orf.start) + v->offset;
  cds_features = gt_array_new(sizeof (GtFeatureNode*));
  cds_feature = (GtFeatureNode*)
                gt_feature_node_new(gt_genome_node_get_seqid((GtGenomeNode*)
                                                             fn),
                                    gt_ft_CDS, cds.start, cds.end,
                                    gt_feature_node_get_strand(fn));
  gt_feature_node_set_source(cds_feature, v->source);
  gt_feature_node_set_phase(cds_feature, GT_PHASE_ZERO);
  gt_feature_node_make_multi_representative(cds_feature);
  representative = cds_feature;
  /* all CDS features in between */
  for (i = strand == GT_STRAND_FORWARD ? orf.start : orf.end;
       strand == GT_STRAND_FORWARD ? i < orf.end : i > orf.start;
       strand == GT_STRAND_FORWARD ? i++ : i--) {
    GT_UNUSED unsigned long prev_length;
    if (gt_splicedseq_pos_is_border(v->splicedseq, i)) {
      gt_feature_node_set_end(cds_feature,
                              gt_splicedseq_map(v->splicedseq, i) + v->offset);
      prev_length = gt_genome_node_get_length((GtGenomeNode*) cds_feature);
      gt_feature_node_add_child(fn, cds_feature);
      gt_array_add(cds_features, cds_feature);
      if (strand == GT_STRAND_FORWARD)
        orf.start = i + 1;
      else
        orf.end = i - 1;
      cds.start = gt_splicedseq_map(v->splicedseq, strand == GT_STRAND_FORWARD
                                    ? orf.start : orf.end) + v->offset;
      cds.end = gt_splicedseq_map(v->splicedseq, strand == GT_STRAND_FORWARD
                                  ? orf.end : orf.start) + v->offset;
      cds_feature = (GtFeatureNode*)
                    gt_feature_node_new(gt_genome_node_get_seqid((GtGenomeNode*)
                                                                 fn),
                                        gt_ft_CDS, cds.start, cds.end,
                                        gt_feature_node_get_strand(fn));
      gt_feature_node_set_source(cds_feature, v->source);
      gt_feature_node_set_multi_representative(cds_feature, representative);
    }
  }
  /* set the end of the last CDS feature and store it */
  gt_feature_node_set_end((GtFeatureNode*) cds_feature,
                          gt_splicedseq_map(v->splicedseq,
                                            strand == GT_STRAND_FORWARD
                                            ? orf.end : orf.start) + v->offset);
  gt_feature_node_add_child(fn, cds_feature);
  gt_array_add(cds_features, cds_feature);
  if (strand == GT_STRAND_REVERSE)
    gt_array_reverse(cds_features);
  set_phases(cds_features);
  gt_array_delete(cds_features);
}

static void create_CDS_features_for_longest_ORF(GtArray *orfs, GtCDSVisitor *v,
                                                GtFeatureNode *fn)
{
  if (gt_array_size(orfs)) {
    /* sort ORFs according to length */
    gt_ranges_sort_by_length_stable(orfs);

    /* create CDS features from the longest ORF */
    if (gt_range_length(gt_array_get_first(orfs)) >=
        v->minorflen * GT_CODON_LENGTH) {
      create_CDS_features_for_ORF(*(GtRange*) gt_array_get_first(orfs), v, fn);
    }
  }
}

static int add_cds_if_necessary(GtFeatureNode *fn, void *data, GtError *err)
{
  GtCDSVisitor *v = (GtCDSVisitor*) data;
  int had_err;

  gt_error_check(err);
  gt_assert(fn);

  had_err = extract_spliced_seq(fn, v, err);
  if (!had_err && gt_splicedseq_length(v->splicedseq) > 2) {
    GtArray *orfs;

    if (gt_feature_node_get_strand(fn) == GT_STRAND_REVERSE) {
      if (gt_splicedseq_reverse(v->splicedseq, err))
        return -1;
    }

    orfs = determine_ORFs_for_all_three_frames(v->splicedseq, v->start_codon,
                                               v->final_stop_codon,
                                               v->generic_start_codons);
    if (!orfs)
      had_err = -1;
    if (!had_err)
      create_CDS_features_for_longest_ORF(orfs, v, fn);
    gt_array_delete(orfs);
  }
  return had_err;
}

static int cds_visitor_feature_node(GtNodeVisitor *nv, GtFeatureNode *fn,
                                    GtError *err)
{
  GtCDSVisitor *v = gt_cds_visitor_cast(nv);
  gt_error_check(err);
  return gt_feature_node_traverse_children(fn, v, add_cds_if_necessary, false,
                                           err);
}

const GtNodeVisitorClass* gt_cds_visitor_class()
{
  static const GtNodeVisitorClass *nvc = NULL;
  gt_class_alloc_lock_enter();
  if (!nvc) {
    nvc = gt_node_visitor_class_new(sizeof (GtCDSVisitor),
                                    cds_visitor_free,
                                    NULL,
                                    cds_visitor_feature_node,
                                    NULL,
                                    NULL,
                                    NULL);
  }
  gt_class_alloc_lock_leave();
  return nvc;
}

GtNodeVisitor* gt_cds_visitor_new(GtRegionMapping *region_mapping,
                                  unsigned int minorflen, GtStr *source,
                                  bool start_codon, bool final_stop_codon,
                                  bool generic_start_codons)
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
  cds_visitor->generic_start_codons = generic_start_codons;
  return nv;
}

void gt_cds_visitor_set_region_mapping(GtCDSVisitor *cds_visitor,
                                       GtRegionMapping *region_mapping)
{
  gt_assert(cds_visitor && region_mapping);
  gt_region_mapping_delete(cds_visitor->region_mapping);
  cds_visitor->region_mapping = gt_region_mapping_ref(region_mapping);
}
