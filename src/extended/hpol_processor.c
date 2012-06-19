/*
  Copyright (c) 2012 Giorgio Gonnella <gonnella@zbh.uni-hamburg.de>
  Copyright (c) 2012 Center for Bioinformatics, University of Hamburg

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

#include "core/disc_distri.h"
#include "core/fastq.h"
#include "core/log_api.h"
#include "core/ma.h"
#include "core/minmax.h"
#include "extended/reverse.h"
#include "core/undef_api.h"
#include "extended/feature_type.h"
#include "extended/seqpos_classifier.h"
#include "extended/aligned_segments_pile.h"
#include "extended/hpol_processor.h"

struct GtHpolProcessor
{
  GtEncseq *encseq;
  unsigned long hmin, max_hlen_diff;
  GtDiscDistri *hdist;
  GtSeqposClassifier *cds_oracle;
  GtAlignedSegmentsPile *asp;
  unsigned long nof_complete_edited, nof_complete_not_edited,
                nof_skipped, nof_unmapped;
  bool adjust_s_hlen;
  double min_alt_consensus;
  GtFile *outfp;
  GtAlphabet *alpha;
};

GtHpolProcessor *gt_hpol_processor_new(GtEncseq *encseq, unsigned long hmin)
{
  GtHpolProcessor *hpp;
  hpp = gt_malloc(sizeof (GtHpolProcessor));
  gt_assert(encseq != NULL);
  gt_assert(hmin > 0);
  hpp->encseq = encseq;
  hpp->hmin = hmin;
  hpp->hdist = gt_disc_distri_new();
  hpp->cds_oracle = NULL;
  hpp->asp = NULL;
  hpp->nof_complete_edited = 0;
  hpp->nof_complete_not_edited = 0;
  hpp->nof_skipped = 0;
  hpp->nof_unmapped = 0;
  hpp->max_hlen_diff = GT_UNDEF_ULONG;
#ifndef S_SPLINT_S
  hpp->min_alt_consensus = 2.0D;
#endif
  hpp->outfp = NULL;
  hpp->alpha = gt_alphabet_new_dna();
  hpp->adjust_s_hlen = false;
  return hpp;
}

void gt_hpol_processor_restrict_to_feature_type(GtHpolProcessor *hpp,
    GtSeqposClassifier *spc)
{
  gt_assert(hpp != NULL);
  gt_assert(spc != NULL);
  hpp->cds_oracle = spc;
}

static void gt_hpol_processor_output_segment(GtAlignedSegment *as,
    bool may_be_gapped, GtFile *outfp)
{
  unsigned long slen;
  if (may_be_gapped)
    gt_aligned_segment_ungap_seq_and_qual(as);
  slen = (unsigned long)strlen(gt_aligned_segment_seq(as));
  gt_assert(slen == (unsigned long)strlen(gt_aligned_segment_qual(as)));
  if (gt_aligned_segment_is_reverse(as))
  {
    GtError *err = gt_error_new();
    char *q = gt_aligned_segment_qual(as), tmp;
    unsigned long i;
    for (i = 0; i < (slen + 1UL) >> 1; i++)
    {
      tmp = q[i];
      q[i] = q[slen - i - 1UL];
      q[slen - i - 1UL] = tmp;
    }
    gt_assert((unsigned long)strlen(gt_aligned_segment_qual(as)) == slen);
    if (gt_reverse_complement(gt_aligned_segment_seq(as), slen, err) != 0)
    {
      fprintf(stderr, "error: %s", gt_error_get(err));
      exit(EXIT_FAILURE);
    }
    gt_error_delete(err);

  }
  gt_fastq_show_entry(gt_aligned_segment_description(as),
      gt_aligned_segment_seq(as), gt_aligned_segment_qual(as),
      slen, 0, false, outfp);
}

static void gt_hpol_processor_process_complete_segment(
    GtAlignedSegment *as, void *data)
{
  GtHpolProcessor *hpp = data;
  gt_assert(hpp != NULL);
  if (gt_aligned_segment_seq_edited(as))
  {
    hpp->nof_complete_edited++;
  }
  else
  {
    hpp->nof_complete_not_edited++;
  }
  if (hpp->outfp != NULL)
    gt_hpol_processor_output_segment(as, gt_aligned_segment_has_indels(as),
        hpp->outfp);
}

static void gt_hpol_processor_process_skipped_segment(
    GtAlignedSegment *as, void *data)
{
  GtHpolProcessor *hpp = data;
  gt_assert(hpp != NULL);
  hpp->nof_skipped++;
  if (hpp->outfp != NULL)
    gt_hpol_processor_output_segment(as, gt_aligned_segment_has_indels(as),
        hpp->outfp);
}

static void gt_hpol_processor_process_unmapped_segment(
    GtAlignedSegment *as, void *data)
{
  GtHpolProcessor *hpp = data;
  gt_assert(hpp != NULL);
  hpp->nof_unmapped++;
  if (hpp->outfp != NULL)
    gt_hpol_processor_output_segment(as, false, hpp->outfp);
}

static void gt_hpol_processor_refregioncheck(
    GtAlignedSegment *as, void *data)
{
  GtHpolProcessor *hpp = data;
  unsigned long rlen, i, startpos;
  char *r;
  bool haserr = false;
  gt_assert(hpp != NULL);
  gt_aligned_segment_ungap_refregion(as);
  gt_aligned_segment_assign_refregion_chars(as, hpp->encseq);
  r = gt_aligned_segment_refregion(as);
  rlen = (unsigned long)strlen(r);
  startpos = gt_aligned_segment_refregion_startpos(as);
  for (i = 0; i < rlen; i++)
  {
    if (r[i] != gt_encseq_get_decoded_char(hpp->encseq, startpos + i,
          GT_READMODE_FORWARD))
    {
      haserr = true;
    }
  }
  if (haserr)
  {
    char *s;
    gt_aligned_segment_ungap_seq_and_qual(as);
    s = gt_aligned_segment_seq(as);
    fprintf(stderr, "segment   = %s\n", s);
    fprintf(stderr, "refregion = %s\n", r);
    fprintf(stderr, "            ");
    for (i = 0; i < rlen; i++)
    {
      if (r[i] != gt_encseq_get_decoded_char(hpp->encseq, startpos + i,
            GT_READMODE_FORWARD))
        fprintf(stderr, "X");
      else
        fprintf(stderr, " ");
    }
    fprintf(stderr, "\n");
    fprintf(stderr, "reference = ");
    for (i = 0; i < rlen; i++)
    {
      fprintf(stderr, "%c", gt_encseq_get_decoded_char(hpp->encseq,
            startpos + i, GT_READMODE_FORWARD));
    }
    fprintf(stderr, "\n");
    fprintf(stderr, "[startpos = %lu]\n", startpos);
    exit(EXIT_FAILURE);
  }
}

void gt_hpol_processor_enable_aligned_segments_refregionscheck(
    GtHpolProcessor *hpp, GtAlignedSegmentsPile *asp)
{
  gt_assert(hpp != NULL);
  gt_assert(asp != NULL);
  gt_assert(hpp->asp == NULL);
  hpp->asp = asp;
  gt_aligned_segments_pile_register_process_complete(hpp->asp,
      gt_hpol_processor_refregioncheck, hpp);
  gt_aligned_segments_pile_register_process_skipped(hpp->asp,
      gt_hpol_processor_refregioncheck, hpp);
}

void gt_hpol_processor_enable_segments_hlen_adjustment(GtHpolProcessor *hpp,
    GtAlignedSegmentsPile *asp, unsigned long max_hlen_diff,
    double min_alt_consensus)
{
  gt_assert(hpp != NULL);
  gt_assert(asp != NULL);
  gt_assert(hpp->asp == NULL);
  hpp->asp = asp;
  hpp->max_hlen_diff = max_hlen_diff;
  hpp->min_alt_consensus = min_alt_consensus;
  hpp->adjust_s_hlen = true;
  gt_aligned_segments_pile_register_process_complete(hpp->asp,
      gt_hpol_processor_process_complete_segment, hpp);
  gt_aligned_segments_pile_register_process_skipped(hpp->asp,
      gt_hpol_processor_process_skipped_segment, hpp);
  gt_aligned_segments_pile_register_process_unmapped(hpp->asp,
      gt_hpol_processor_process_unmapped_segment, hpp);
}

void gt_hpol_processor_enable_segments_output(
    GtHpolProcessor *hpp, GtFile *outfile)
{
  gt_assert(hpp != NULL);
  gt_assert(outfile != NULL);
  hpp->outfp = outfile;
}

GT_UNUSED
static unsigned long gt_hpol_processor_determine_hlen_backwards(char *s,
    char *q, unsigned long pos, char c, unsigned long *q_sum,
    unsigned long *gaps)
{
  unsigned long s_hlen = 0;
  pos++;
  while (pos > 0)
  {
    pos--;
    if (s[pos] == c)
    {
      s_hlen++;
      if (q_sum != NULL)
        *q_sum += q[pos];
    }
    else if (s[pos] == '-')
    {
      if (gaps != NULL)
        (*gaps)++;
    }
    else
      break;
  }
  return s_hlen;
}

static unsigned long gt_hpol_processor_determine_hlen_forwards(char *s, char *q,
    unsigned long pos, unsigned long maxpos, char c, unsigned long *q_sum,
    unsigned long *gaps)
{
  unsigned long s_hlen = 0;
  while (pos <= maxpos)
  {
    if (s[pos] == c)
    {
      s_hlen++;
      if (q_sum != NULL)
        *q_sum += q[pos];
    }
    else if (s[pos] == '-')
    {
      if (gaps != NULL)
        (*gaps)++;
    }
    else
      break;
    pos++;
  }
  return s_hlen;
}

static void gt_hpol_processor_subst_in_range(char *s, char *q,
    unsigned long left, unsigned long right, unsigned long length, char old,
    char new, char qual)
{
  unsigned long pos;
  for (pos = right + 1UL; pos > left && length > 0; /**/)
  {
    pos--;
    if (s[pos] == old)
    {
      s[pos] = new;
      q[pos] = qual;
      length--;
    }
  }
  gt_assert(length == 0);
}

static void gt_hpol_processor_shrink_hpol(char *s, char *q, unsigned long left,
    unsigned long right, unsigned long length, char c)
{
  gt_hpol_processor_subst_in_range(s, q, left, right, length, c, '-',
      GT_UNDEF_CHAR);
}

static void gt_hpol_processor_enlarge_hpol(char *s, char *q, unsigned long left,
    unsigned long right, unsigned long length, char c, char qual)
{
  gt_hpol_processor_subst_in_range(s, q, left, right, length, '-', c, qual);
}

static void gt_hpol_processor_adjust_hlen_of_a_segment(GtAlignedSegment *as,
    char c, unsigned long r_hstart, unsigned long r_hlen,
    unsigned long max_hlen_diff)
{
  unsigned long left, right, s_hlen, q_sum = 0, s_free = 0;
  char *s, *q;
  left = gt_aligned_segment_offset_for_refpos(as, r_hstart);
  right = gt_aligned_segment_offset_for_refpos(as, r_hstart + r_hlen);
  if (left == GT_UNDEF_ULONG || right == GT_UNDEF_ULONG)
    return;
  right--;
  s = gt_aligned_segment_seq(as);
  q = gt_aligned_segment_qual(as);
#ifdef GG_DEBUG
  if (left > 0)
    s_hlen = gt_hpol_processor_determine_hlen_backwards(s, q, left - 1UL,
          c, &q_sum, &s_free);
  gt_assert(s_hlen == 0);
#endif
  s_hlen = gt_hpol_processor_determine_hlen_forwards(s, q, left, right, c,
      &q_sum, &s_free);
  if (s_hlen == 0)
    return;
#ifdef GG_DEBUG
  if (s_hlen != r_hlen)
  {
    unsigned long i;
    gt_aligned_segment_show(as, NULL);
    printf("   ");
    for (i = 0; i < left; i++)
      printf(" ");
    printf("|");
    for (/**/; i < right - 1UL; i++)
      printf(".");
    printf("|\n");
    printf("[s_hlen=%lu, r_hlen=%lu]\n", s_hlen, r_hlen);
    printf("\n");
  }
#endif
  /* do the editing if necessary and possible */
  if (s_hlen < r_hlen)
  {
    if (s_free > 0)
    {
      unsigned long hlen_diff = r_hlen - s_hlen;
      if (hlen_diff < max_hlen_diff)
      {
        gt_hpol_processor_enlarge_hpol(s, q, left, right,
            MIN(s_free, hlen_diff), c, (char)(q_sum / s_hlen));
      }
      gt_aligned_segment_seq_set_edited(as);
    }
  }
  else if (s_hlen > r_hlen)
  {
    unsigned long hlen_diff = s_hlen - r_hlen;
    if (hlen_diff < max_hlen_diff)
    {
      gt_hpol_processor_shrink_hpol(s, q, left, right, hlen_diff, c);
    }
    gt_aligned_segment_seq_set_edited(as);
  }
}

static void gt_hpol_processor_adjust_hlen_of_all_segments(
    GtAlignedSegmentsPile *asp, char c, unsigned long r_hstart,
    unsigned long r_hlen, unsigned long max_hlen_diff)
{
  GtDlistelem *dlistelem;
  gt_assert(asp != NULL);
  for (dlistelem = gt_dlist_first(gt_aligned_segments_pile_get(asp));
        dlistelem != NULL; dlistelem = gt_dlistelem_next(dlistelem))
  {
    GtAlignedSegment *as;
    as = gt_dlistelem_get_data(dlistelem);
    if (gt_aligned_segment_has_indels(as))
      gt_hpol_processor_adjust_hlen_of_a_segment(as, c, r_hstart, r_hlen,
          max_hlen_diff);
  }
}

static void gt_hpol_processor_determine_alternative_consensus(
    GtHpolProcessor *hpp, char c, unsigned long r_hstart, unsigned long r_hlen,
    unsigned long *c_s_hlen, unsigned long *c_support, unsigned long *piled,
    unsigned long *r_hlen_support)
{
  unsigned long *occ;
  GtDlistelem *dlistelem;
  unsigned long s_hlen_max = r_hlen << 1, i;
  *c_support = 0;
  *piled = 0;
  occ = gt_calloc((size_t)(s_hlen_max + 1UL), sizeof (*occ));
  for (dlistelem = gt_dlist_first(gt_aligned_segments_pile_get(hpp->asp));
      dlistelem != NULL; dlistelem = gt_dlistelem_next(dlistelem))
  {
    GtAlignedSegment *as;
    unsigned long left, right, s_hlen;
    char *s;
    as = gt_dlistelem_get_data(dlistelem);
#ifdef GG_DEBUG
    if (gt_aligned_segment_has_indels(as))
      gt_aligned_segment_assign_refregion_chars(as, hpp->encseq);
#endif
    left = gt_aligned_segment_offset_for_refpos(as, r_hstart);
    right = gt_aligned_segment_offset_for_refpos(as, r_hstart + r_hlen);
    if (left == GT_UNDEF_ULONG || right == GT_UNDEF_ULONG)
      continue;
    (*piled)++;
    right--;
    s = gt_aligned_segment_seq(as);
    s_hlen = gt_hpol_processor_determine_hlen_forwards(s, NULL, left, right, c,
        NULL, NULL);
    occ[MIN(s_hlen, s_hlen_max)]++;
  }
  *r_hlen_support = occ[r_hlen];
  *c_s_hlen = 0;
  *c_support = occ[0];
  for (i = 0; i <= s_hlen_max; i++)
  {
    if (i != r_hlen && occ[i] > *c_support)
    {
      *c_support = occ[i];
      *c_s_hlen = i;
    }
  }
  gt_free(occ);
}

static void gt_hpol_processor_process_hpol_end(GtHpolProcessor *hpp,
    GtUchar c, unsigned long endpos, unsigned long hlen,
    unsigned long max_hlen_diff, double min_alt_consensus)
{
  gt_disc_distri_add(hpp->hdist, hlen);
  if (hpp->adjust_s_hlen)
  {
    char ch = gt_alphabet_decode(hpp->alpha, c);
    unsigned long min_alt_support = ULONG_MAX, alt_support = 0, piled,
                  r_hlen_support = 0;
    gt_aligned_segments_pile_move_over_position(hpp->asp, endpos + 1UL);
    piled = gt_aligned_segments_pile_size(hpp->asp);
    if (piled > 0)
    {
#if !defined(GG_DEBUG) && !defined(S_SPLINT_S)
      if (min_alt_consensus <= 1.0D)
#endif
      {
        unsigned long alt_s_hlen;
        gt_hpol_processor_determine_alternative_consensus(hpp, ch,
            endpos + 1UL - hlen, hlen, &alt_s_hlen, &alt_support, &piled,
            &r_hlen_support);
        if (piled > 0)
        {
          min_alt_support = (unsigned long)(min_alt_consensus * (double)piled);
        }
      }
      if (piled > 0 && r_hlen_support != piled && alt_support < min_alt_support)
        gt_hpol_processor_adjust_hlen_of_all_segments(hpp->asp, ch,
            endpos + 1UL - hlen, hlen, max_hlen_diff);
    }
  }
}

static void gt_hpol_processor_show_hdist_elem(unsigned long key,
    unsigned long long value, void *data)
{
  GtLogger *logger = data;
  gt_logger_log(logger, "%lu: %llu", key, value);
}

static void gt_hpol_processor_show_hdist(GtHpolProcessor *hpp, GtLogger *logger)
{
  gt_assert(hpp != NULL);
  gt_assert(hpp->hdist != NULL);
  gt_logger_log(logger, "Distribution of homopolymers of length >= %lu %s",
      hpp->hmin, (hpp->cds_oracle != NULL ?  "in coding sequences" :
        "in whole reference sequence"));
  gt_disc_distri_foreach(hpp->hdist, gt_hpol_processor_show_hdist_elem,
      logger);
  if (hpp->cds_oracle != NULL)
  {
    gt_logger_log(logger, "coding sequences: %lu",
        gt_seqpos_classifier_nof_features_found(hpp->cds_oracle));
  }
  if (hpp->adjust_s_hlen)
  {
    gt_logger_log(logger, "segments in SAM file: %lu",
        hpp->nof_complete_edited + hpp->nof_complete_not_edited +
        hpp->nof_skipped + hpp->nof_unmapped);
    gt_logger_log(logger, "- processed and not edited: %lu",
        hpp->nof_complete_not_edited);
    gt_logger_log(logger, "- processed and edited: %lu",
        hpp->nof_complete_edited);
    gt_logger_log(logger, "- not processed: %lu",
        hpp->nof_skipped);
    gt_logger_log(logger, "- not mapping: %lu",
        hpp->nof_unmapped);
  }
}

int gt_hpol_processor_run(GtHpolProcessor *hpp, GtLogger *logger, GtError *err)
{
  int had_err = 0;
  unsigned long i, hlen, tlen;
  GtUchar prev, c;
  bool coding = false;
  bool end_of_annotation = true;
  GtEncseqReader *esr;
  gt_assert(hpp != NULL);
  gt_assert(hpp->encseq != NULL);
  esr = gt_encseq_create_reader_with_readmode(hpp->encseq,
      GT_READMODE_FORWARD, 0);
  tlen = gt_encseq_total_length(hpp->encseq);
  prev = gt_encseq_reader_next_encoded_char(esr);
  hlen = 1UL;
  if (hpp->cds_oracle != NULL)
    had_err = gt_seqpos_classifier_position_is_inside_feature(
        hpp->cds_oracle, 0, &coding, &end_of_annotation, err);
  for (i = 1UL; i < tlen && !had_err; i++)
  {
    if (hpp->cds_oracle != NULL)
    {
      GT_UNUSED bool was_not_coding = !coding;
      had_err = gt_seqpos_classifier_position_is_inside_feature(
          hpp->cds_oracle, i, &coding, &end_of_annotation, err);
      /* break at boundary:
      if (was_not_coding && coding)
        hlen = 0; */
    }
    if (!had_err)
    {
      c = gt_encseq_reader_next_encoded_char(esr);
      if (prev == c)
      {
        hlen++;
      }
      else
      {
        if (hlen >= hpp->hmin && (hpp->cds_oracle == NULL || coding))
          gt_hpol_processor_process_hpol_end(hpp, prev, i - 1UL, hlen,
              hpp->max_hlen_diff, hpp->min_alt_consensus);
        hlen = 1UL;
      }
      prev = c;
    }
  }
  if (!had_err)
  {
    if (hlen >= hpp->hmin && (hpp->cds_oracle == NULL || coding))
      gt_hpol_processor_process_hpol_end(hpp, prev, i - 1UL, hlen,
          hpp->max_hlen_diff, hpp->min_alt_consensus);
  }
  gt_encseq_reader_delete(esr);
  if (logger != NULL && !had_err)
    gt_hpol_processor_show_hdist(hpp, logger);
  return had_err;
}

void gt_hpol_processor_delete(GtHpolProcessor *hpp)
{
  if (hpp != NULL)
  {
    gt_disc_distri_delete(hpp->hdist);
    gt_alphabet_delete(hpp->alpha);
    gt_free(hpp);
  }
}
