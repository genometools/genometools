/*
  Copyright (c) 2006-2011 Gordon Gremme <gordon@gremme.org>
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

#include "core/codon_iterator_api.h"
#include "core/codon_iterator_simple_api.h"
#include "core/ma.h"
#include "core/trans_table_api.h"
#include "core/translator_api.h"
#include "extended/feature_node_iterator_api.h"
#include "extended/feature_node.h"
#include "extended/genome_node.h"
#include "extended/gff3_defines.h"
#include "extended/gff3_parser.h"
#include "extended/region_mapping_api.h"
#include "extended/reverse_api.h"

static int extract_join_feature(GtGenomeNode *gn, const char *type,
                                GtRegionMapping *region_mapping,
                                GtStr *sequence, bool *reverse_strand,
                                GtError *err)
{
  char *outsequence;
  GtFeatureNode *fn;
  GtRange range;
  int had_err = 0;

  gt_error_check(err);
  fn = gt_feature_node_cast(gn);
  gt_assert(fn);

  if (gt_feature_node_has_type(fn, type)) {
    range = gt_genome_node_get_range(gn);
    had_err = gt_region_mapping_get_sequence(region_mapping, &outsequence,
                                             gt_genome_node_get_seqid(gn),
                                             range.start, range.end, err);
    if (!had_err) {
      gt_str_append_cstr_nt(sequence, outsequence, gt_range_length(&range));
      gt_free(outsequence);
      if (gt_feature_node_get_strand(fn) == GT_STRAND_REVERSE)
        *reverse_strand = true;
    }
  }
  return had_err;
}

static int gt_extract_feature_sequence_generic(GtStr *sequence,
                                GtGenomeNode *gn,
                                const char *type, bool join, GtStr *seqid,
                                GtStrArray *target_ids,
                                unsigned int *out_phase_offset,
                                GtRegionMapping *region_mapping, GtError *err)
{
  GtFeatureNode *fn;
  GtRange range;
  unsigned int phase_offset = 0;
  char *outsequence;
  const char *target;
  int had_err = 0;

  gt_error_check(err);
  fn = gt_genome_node_cast(gt_feature_node_class(), gn);
  gt_assert(fn);

  if (seqid)
    gt_str_append_str(seqid, gt_genome_node_get_seqid(gn));
  if (target_ids &&
      (target = gt_feature_node_get_attribute(fn, GT_GFF_TARGET))) {
    had_err = gt_gff3_parser_parse_all_target_attributes(target, false,
                                                         target_ids, NULL,
                                                         NULL, "", 0, err);
  }
  if (!had_err) {
    if (join) {
      GtFeatureNodeIterator *fni;
      GtFeatureNode *child;
      bool reverse_strand = false,
           first_child = true;
      /* in this case we have to traverse the children */
      fni = gt_feature_node_iterator_new_direct(gt_feature_node_cast(gn));
      while (!had_err && (child = gt_feature_node_iterator_next(fni))) {
        GtPhase phase = gt_feature_node_get_phase(child);
        if (gt_feature_node_get_strand(gt_feature_node_cast(gn))
              == GT_STRAND_REVERSE  || first_child) {
          phase_offset = (int) phase;
        }
        if (first_child) {
          if (target_ids &&
               (target = gt_feature_node_get_attribute(child, GT_GFF_TARGET))) {
            gt_str_array_reset(target_ids);
            had_err = gt_gff3_parser_parse_all_target_attributes(target, false,
                                                                 target_ids,
                                                                 NULL,
                                                                 NULL, "", 0,
                                                                 err);
          }
          first_child = false;
        }
        if (!had_err) {
          if (extract_join_feature((GtGenomeNode*) child, type, region_mapping,
                                   sequence, &reverse_strand, err)) {
            had_err = -1;
          }
        }
      }
      gt_feature_node_iterator_delete(fni);
      gt_assert(phase_offset <= (unsigned int) GT_PHASE_UNDEFINED);
      if (!had_err && gt_str_length(sequence)) {
        if (reverse_strand) {
          had_err = gt_reverse_complement(gt_str_get(sequence),
                                          gt_str_length(sequence), err);
        }
      }
    }
    else if (gt_feature_node_get_type(fn) == type) {
      GtPhase phase = gt_feature_node_get_phase(fn);
      gt_assert(!had_err);
      if (phase != GT_PHASE_UNDEFINED)
        phase_offset = (unsigned int) phase;
      /* otherwise we only have to look this feature */
      range = gt_genome_node_get_range(gn);
      gt_assert(range.start); /* 1-based coordinates */
      had_err = gt_region_mapping_get_sequence(region_mapping, &outsequence,
                                               gt_genome_node_get_seqid(gn),
                                               range.start, range.end, err);
      if (!had_err) {
        gt_str_append_cstr_nt(sequence, outsequence, gt_range_length(&range));
        gt_free(outsequence);
        if (gt_feature_node_get_strand(fn) == GT_STRAND_REVERSE) {
          had_err = gt_reverse_complement(gt_str_get(sequence),
                                          gt_str_length(sequence), err);
        }
      }
    }
  }
  if (out_phase_offset && phase_offset != GT_PHASE_UNDEFINED) {
    *out_phase_offset = phase_offset;
  }
  return had_err;
}

int gt_extract_feature_sequence(GtStr *sequence, GtGenomeNode *gn,
                                const char *type, bool join, GtStr *seqid,
                                GtStrArray *target_ids,
                                GtRegionMapping *region_mapping, GtError *err)
{
  return gt_extract_feature_sequence_generic(sequence, gn, type, join, seqid,
                                             target_ids, NULL, region_mapping,
                                             err);
}

int gt_extract_and_translate_feature_sequence(GtFeatureNode *feature_node,
                                              const char *type,
                                              bool join,
                                              GtRegionMapping *rm,
                                              GtTransTable *ttable,
                                              GtStr *translation_fr1,
                                              GtStr *translation_fr2,
                                              GtStr *translation_fr3,
                                              GtError *err)
{
  GtTranslator *tr = NULL;
  GtTranslatorStatus status;
  GtCodonIterator *ci = NULL;
  unsigned int frame, phase_offset = 0;
  char translated;
  int had_err = 0;
  GtStr *sequence = gt_str_new();
  gt_assert(feature_node && type);

  had_err = gt_extract_feature_sequence_generic(sequence,
                                                (GtGenomeNode*) feature_node,
                                                type, join, NULL, NULL,
                                                &phase_offset, rm, err);

  if (!had_err && gt_str_length(sequence) > phase_offset) {
    ci = gt_codon_iterator_simple_new(gt_str_get(sequence) + phase_offset,
                                      gt_str_length(sequence) - phase_offset,
                                      NULL);
    tr = gt_translator_new(ci);
    if (ttable)
      gt_translator_set_translation_table(tr, ttable);
    status = gt_translator_next(tr, &translated, &frame, NULL);
    while (status == GT_TRANSLATOR_OK) {
      if (frame == 0 && translation_fr1)
        gt_str_append_char(translation_fr1, translated);
      else if (frame == 1 && translation_fr2)
        gt_str_append_char(translation_fr2, translated);
      else if (frame == 2 && translation_fr3)
        gt_str_append_char(translation_fr3, translated);
      status = gt_translator_next(tr, &translated, &frame, NULL);
    }
    if (status == GT_TRANSLATOR_ERROR)
      had_err = -1;
  }
  gt_translator_delete(tr);
  gt_codon_iterator_delete(ci);
  gt_str_delete(sequence);

  return had_err;
}
