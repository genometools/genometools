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

#include "core/complement.h"
#include "core/disc_distri.h"
#include "core/fastq.h"
#include "core/hashmap_api.h"
#include "core/log_api.h"
#include "core/ma.h"
#include "core/minmax.h"
#include "core/splitter.h"
#include "core/undef_api.h"
#include "core/warning_api.h"
#include "extended/aligned_segments_pile.h"
#include "extended/feature_type.h"
#include "extended/hpol_processor.h"
#include "extended/reverse_api.h"
#include "extended/seqpos_classifier.h"

struct GtHpolProcessor
{
  GtEncseq *encseq;
  unsigned long hmin, max_hlen_diff;
  GtDiscDistri *hdist, *hdist_e;
  GtSeqposClassifier *cds_oracle;
  GtAlignedSegmentsPile *asp;
  unsigned long nof_complete, nof_complete_edited, nof_complete_not_edited,
                nof_skipped, nof_unmapped, nof_h, nof_h_e, hlen_max;
  bool adjust_s_hlen;
  double min_alt_consensus;
  GtFile *outfp_segments, *outfp_stats;
  GtSeqIterator *reads_iter;
  GtHashmap *processed_segments;
  GtAlphabet *alpha;
  bool output_segments, output_stats;
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
  hpp->nof_h = 0;
  hpp->hdist_e = gt_disc_distri_new();
  hpp->nof_h_e = 0;
  hpp->hlen_max = 0;
  hpp->cds_oracle = NULL;
  hpp->asp = NULL;
  hpp->nof_complete = 0;
  hpp->nof_complete_edited = 0;
  hpp->nof_complete_not_edited = 0;
  hpp->nof_skipped = 0;
  hpp->nof_unmapped = 0;
  hpp->max_hlen_diff = GT_UNDEF_ULONG;
  hpp->min_alt_consensus = (double) 2.0;
  hpp->outfp = NULL;
  hpp->alpha = gt_alphabet_new_dna();
  hpp->adjust_s_hlen = false;
  hpp->output_segments = false;
  hpp->outfp_segments = NULL;
  hpp->output_stats = false;
  hpp->outfp_stats = NULL;
  hpp->processed_segments = NULL;
  hpp->reads_iter = NULL;
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
    bool may_be_gapped, GtFile *outfp, const char *desc)
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
  gt_fastq_show_entry((desc != NULL) ? desc :
      gt_aligned_segment_description(as), gt_aligned_segment_seq(as),
      gt_aligned_segment_qual(as), slen, 0, false, outfp);
}

static void gt_hpol_processor_process_complete_segment(
    GtAlignedSegment *as, void *data)
{
  GtHpolProcessor *hpp = data;
  gt_assert(hpp != NULL);
  (hpp->nof_complete)++;
  if (gt_aligned_segment_seq_edited(as))
  {
    (hpp->nof_complete_edited)++;
  }
  else
  {
    (hpp->nof_complete_not_edited)++;
  }
  if (hpp->output_segments)
    gt_hpol_processor_output_segment(as, gt_aligned_segment_has_indels(as),
        hpp->outfp_segments, NULL);
  if (hpp->processed_segments != NULL)
    gt_hashmap_add(hpp->processed_segments,
        (void*)gt_aligned_segment_description(as), as);
}

static void gt_hpol_processor_process_skipped_segment(
    GtAlignedSegment *as, void *data)
{
  GtHpolProcessor *hpp = data;
  gt_assert(hpp != NULL);
  hpp->nof_skipped++;
  if (hpp->output_segments)
    gt_hpol_processor_output_segment(as, gt_aligned_segment_has_indels(as),
        hpp->outfp_segments, NULL);
  if (hpp->processed_segments != NULL)
    gt_hashmap_add(hpp->processed_segments,
        (void*)gt_aligned_segment_description(as), as);
}

static void gt_hpol_processor_process_unmapped_segment(
    GtAlignedSegment *as, void *data)
{
  GtHpolProcessor *hpp = data;
  gt_assert(hpp != NULL);
  hpp->nof_unmapped++;
  if (hpp->output_segments)
    gt_hpol_processor_output_segment(as, false, hpp->outfp_segments, NULL);
  if (hpp->processed_segments != NULL)
    gt_hashmap_add(hpp->processed_segments,
        (void*)gt_aligned_segment_description(as), as);
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

void gt_hpol_processor_enable_segments_output(GtHpolProcessor *hpp,
    GtFile *outfile)
{
  gt_assert(hpp != NULL);
  hpp->output_segments = true;
  hpp->outfp_segments = outfile;
}

void gt_hpol_processor_sort_segments_output(GtHpolProcessor *hpp,
    GtSeqIterator *reads_iter)
{
  gt_assert(hpp != NULL);
  gt_assert(hpp->output_segments);
  hpp->output_segments = false;
  gt_aligned_segments_pile_disable_segment_deletion(hpp->asp);
  hpp->processed_segments = gt_hashmap_new(GT_HASH_STRING, NULL,
      (GtFree)gt_aligned_segment_delete);
  hpp->reads_iter = reads_iter;
}

static void gt_hpol_processor_output_stats_header(GtFile *outfp)
{
  gt_file_xprintf(outfp, "# correction statistics\n");
  gt_file_xprintf(outfp, "# edit =     edit operation (I or D)\n");
  gt_file_xprintf(outfp, "# r_hpos =   start pos of hpol on reference\n");
  gt_file_xprintf(outfp, "# r_hlen =   length of hpol on reference\n");
  gt_file_xprintf(outfp, "# s_hpos =   start pos of hpol on read\n");
  gt_file_xprintf(outfp, "# s_hend =   end pos of hpol on read\n");
  gt_file_xprintf(outfp, "# s_hlen =   length of hpol in read\n");
  gt_file_xprintf(outfp, "# c_len =    correction length\n");
  gt_file_xprintf(outfp, "# s_char =   hpol character in read\n");
  gt_file_xprintf(outfp, "# s_or =     orientation of read "
      "(+ or -; + = same as reference)\n");
  gt_file_xprintf(outfp, "# s_q_ave =  average quality of read in the hpol "
      "positions\n");
  gt_file_xprintf(outfp, "# s_q_bef =  quality of base before the hpol\n");
  gt_file_xprintf(outfp, "# s_q_aft =  quality of base after the hpol\n");
  gt_file_xprintf(outfp, "# s_qual =   quality string in read for the hpol "
      "positions\n");
  gt_file_xprintf(outfp, "# s_id =     read identifier\n");
  gt_file_xprintf(outfp, "# coordinates are 1-based\n");
  gt_file_xprintf(outfp, "#\n");
  gt_file_xprintf(outfp, "# edit\tr_hpos\tr_hlen\ts_hpos\ts_hend\ts_hlen"
      "\tc_len\ts_char\ts_or\ts_q_ave\ts_q_bef\ts_q_aft\ts_qual\ts_id\n");
}

#define GT_HPOL_PROCESSOR_PHREDOFFSET 33UL

static void gt_hpol_processor_output_stats(GtAlignedSegment *as,
    unsigned long r_pos, unsigned long r_hlen, unsigned long s_hlen,
    char c, unsigned long s_q_sum, unsigned long c_len, GtFile *outfp)
{
  unsigned long i, pos, s_pos = 0, s_offset, s_q_bef, s_q_aft;
  double s_q_ave;
  char *s_qual, *q;
  const char *r_desc;
  r_desc = gt_aligned_segment_description(as);
  q = gt_aligned_segment_qual(as);
  if (gt_aligned_segment_is_reverse(as))
  {
    GtError *err = gt_error_new();
    (void)gt_complement(&c, c, err);
    gt_error_delete(err);
  }
  gt_assert(s_hlen > 0);
  gt_assert(s_q_sum >= GT_HPOL_PROCESSOR_PHREDOFFSET * s_hlen);
  s_q_ave = (double)(s_q_sum - GT_HPOL_PROCESSOR_PHREDOFFSET * s_hlen)
    / (double)s_hlen;
  s_pos = gt_aligned_segment_orig_seqpos_for_refpos(as, r_pos);
  s_offset = gt_aligned_segment_offset_for_refpos(as, r_pos);
  s_qual = gt_malloc(sizeof (*s_qual) * (s_hlen + 1UL));
  s_q_bef = GT_UNDEF_ULONG;
  for (i = s_offset; i > 0; /**/)
  {
    i--;
    if (q[i] != GT_UNDEF_CHAR)
    {
      s_q_bef = (unsigned long)q[i] - GT_HPOL_PROCESSOR_PHREDOFFSET;
      break;
    }
  }
  gt_assert(s_q_bef != GT_UNDEF_ULONG);
  if (!gt_aligned_segment_is_reverse(as))
  {
    for (i = s_offset, pos = 0; pos < s_hlen; i++)
    {
      if (q[i] != GT_UNDEF_CHAR)
      {
        s_qual[pos] = q[i];
        pos++;
      }
    }
  }
  else
  {
    for (i = s_offset, pos = s_hlen; pos > 0; i++)
    {
      if (q[i] != GT_UNDEF_CHAR)
      {
        s_qual[pos - 1UL] = q[i];
        pos--;
      }
    }
  }
  s_q_aft = GT_UNDEF_ULONG;
  for (/**/; i < gt_aligned_segment_length(as); i++)
  {
    if (q[i] != GT_UNDEF_CHAR)
    {
      s_q_aft = (unsigned long)q[i] - GT_HPOL_PROCESSOR_PHREDOFFSET;
      break;
    }
  }
  gt_assert(s_q_aft != GT_UNDEF_ULONG);
  s_qual[s_hlen] = '\0';
  gt_assert(r_hlen != s_hlen);
  s_pos++;
  r_pos++;
  gt_file_xprintf(outfp,
      "%c\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%c\t%c\t%.2f\t%lu\t%lu\t%s\t%s\n",
      r_hlen > s_hlen ? 'I' : 'D', r_pos, r_hlen,
      gt_aligned_segment_is_reverse(as) ? s_pos - s_hlen + 1UL : s_pos,
      gt_aligned_segment_is_reverse(as) ? s_pos : s_pos + s_hlen - 1UL,
      s_hlen, c_len, c, gt_aligned_segment_is_reverse(as) ? '-' : '+',
      s_q_ave, gt_aligned_segment_is_reverse(as) ? s_q_aft : s_q_bef,
      gt_aligned_segment_is_reverse(as) ? s_q_bef : s_q_aft, s_qual, r_desc);
  gt_free(s_qual);
}

void gt_hpol_processor_enable_statistics_output(GtHpolProcessor *hpp,
    GtFile *outfile)
{
  gt_assert(hpp != NULL);
  hpp->output_stats = true;
  hpp->outfp_stats = outfile;
  gt_hpol_processor_output_stats_header(hpp->outfp_stats);
  gt_aligned_segments_pile_enable_edit_tracking(hpp->asp);
}

#ifdef GG_DEBUG
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
#endif

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

static bool gt_hpol_processor_adjust_hlen_of_a_segment(GtAlignedSegment *as,
    char c, unsigned long r_hstart, unsigned long r_hlen,
    unsigned long max_hlen_diff, bool output_stats, GtFile *outfp_stats)
{
  unsigned long left, right, s_hlen = 0, q_sum = 0, s_free = 0;
  char *s, *q;
  left = gt_aligned_segment_offset_for_refpos(as, r_hstart);
  right = gt_aligned_segment_offset_for_refpos(as, r_hstart + r_hlen);
  if (left == GT_UNDEF_ULONG || left == 0 ||
      right == GT_UNDEF_ULONG || right == gt_aligned_segment_length(as))
    return false;
  right--;
  s = gt_aligned_segment_seq(as);
  q = gt_aligned_segment_qual(as);
#ifdef GG_DEBUG
  if (left > 0)
    s_hlen = gt_hpol_processor_determine_hlen_backwards(s, q, left - 1UL,
          c, &q_sum, &s_free);
  if (s_hlen > 0)
    printf("backwards_hlen = %lu\n", s_hlen);
  s_free = 0;
  q_sum = 0;
#endif
  s_hlen = gt_hpol_processor_determine_hlen_forwards(s, q, left, right, c,
      &q_sum, &s_free);
  if (s_hlen == 0)
    return false;
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
  else
  {
    printf("s_hlen == r_hlen = %lu\n\n", s_hlen);
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
        if (output_stats)
        {
          gt_hpol_processor_output_stats(as, r_hstart, r_hlen, s_hlen,
              c, q_sum, MIN(s_free, hlen_diff), outfp_stats);
        }
        gt_aligned_segment_seq_set_edited(as);
        gt_log_log("edited: %s", gt_aligned_segment_description(as));
        gt_hpol_processor_enlarge_hpol(s, q, left, right,
            MIN(s_free, hlen_diff), c, (char)(q_sum / s_hlen));
      }
    }
    return true;
  }
  else if (s_hlen > r_hlen)
  {
    unsigned long hlen_diff = s_hlen - r_hlen;
    if (hlen_diff < max_hlen_diff)
    {
      if (output_stats)
      {
        gt_hpol_processor_output_stats(as, r_hstart, r_hlen, s_hlen,
            c, q_sum, hlen_diff, outfp_stats);
      }
      gt_aligned_segment_seq_set_edited(as);
      gt_log_log("edited: %s", gt_aligned_segment_description(as));
      gt_hpol_processor_shrink_hpol(s, q, left, right, hlen_diff, c);
    }
    return true;
  }
  return false;
}

static bool gt_hpol_processor_adjust_hlen_of_all_segments(
    GtAlignedSegmentsPile *asp, char c, unsigned long r_hstart,
    unsigned long r_hlen, unsigned long max_hlen_diff, bool output_stats,
    GtFile *outfp_stats)
{
  GtDlistelem *dlistelem;
  bool any_edited = false, edited;
  gt_assert(asp != NULL);
  for (dlistelem = gt_dlist_first(gt_aligned_segments_pile_get(asp));
        dlistelem != NULL; dlistelem = gt_dlistelem_next(dlistelem))
  {
    GtAlignedSegment *as;
    as = gt_dlistelem_get_data(dlistelem);
    if (gt_aligned_segment_has_indels(as))
    {
      edited = gt_hpol_processor_adjust_hlen_of_a_segment(as, c, r_hstart,
          r_hlen, max_hlen_diff, output_stats, outfp_stats);
      if (edited)
        any_edited = true;
    }
  }
  return any_edited;
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
  bool edited = false;
  gt_disc_distri_add(hpp->hdist, hlen);
  hpp->nof_h++;
  if (hlen > hpp->hlen_max)
    hpp->hlen_max = hlen;
  if (hpp->adjust_s_hlen)
  {
    char ch = gt_alphabet_decode(hpp->alpha, c);
    unsigned long min_alt_support = ULONG_MAX, alt_support = 0, piled,
                  r_hlen_support = 0;
    gt_aligned_segments_pile_move_over_position(hpp->asp, endpos + 1UL);
    piled = gt_aligned_segments_pile_size(hpp->asp);
    if (piled > 0)
    {
#ifndef GG_DEBUG
      if (min_alt_consensus <= (double) 1.0)
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
      {
        edited = gt_hpol_processor_adjust_hlen_of_all_segments(hpp->asp, ch,
            endpos + 1UL - hlen, hlen, max_hlen_diff, hpp->output_stats,
            hpp->outfp_stats);
      }
    }
  }
  if (edited)
  {
    hpp->nof_h_e++;
    gt_disc_distri_add(hpp->hdist_e, hlen);
  }
}

static void gt_hpol_processor_show_hdist(GtHpolProcessor *hpp, GtLogger *logger)
{
  unsigned long i;
  gt_assert(hpp != NULL);
  gt_assert(hpp->hdist != NULL);
  gt_logger_log(logger, "Distribution of homopolymers of length >= %lu %s",
      hpp->hmin, (hpp->cds_oracle != NULL ?  "in coding sequences" :
        "in whole reference sequence"));
  gt_logger_log(logger, "length\toccurrences\tedited");
  for (i = hpp->hmin; i < hpp->hlen_max; i++)
  {
    gt_logger_log(logger, "%-6lu\t%-11lu\t%-6lu\t(%.2f%%)", i,
        (unsigned long)gt_disc_distri_get(hpp->hdist, i),
        (unsigned long)gt_disc_distri_get(hpp->hdist_e, i),
        (double)gt_disc_distri_get(hpp->hdist_e, i) * 100 /
        (double)gt_disc_distri_get(hpp->hdist, i));
  }
  gt_logger_log(logger, "total \t%-11lu\t%-6lu\t(%.2f%%)",
      hpp->nof_h, hpp->nof_h_e, (double)hpp->nof_h_e * 100 /
      (double)hpp->nof_h);
  if (hpp->cds_oracle != NULL)
  {
    gt_logger_log(logger, "coding sequences: %lu",
        gt_seqpos_classifier_nof_features_found(hpp->cds_oracle));
  }
  if (hpp->adjust_s_hlen)
  {
    unsigned long total_segments =
        hpp->nof_complete_edited + hpp->nof_complete_not_edited +
        hpp->nof_skipped + hpp->nof_unmapped;
    gt_logger_log(logger, "segments in SAM file:       %lu",
        total_segments);
    gt_logger_log(logger, "- processed:                %-7lu (%.2f%%)",
        hpp->nof_complete, (double)hpp->nof_complete * 100 / total_segments);
    gt_logger_log(logger, "  ... and not edited:       %-7lu (%.2f%%)",
        hpp->nof_complete_not_edited,
        (double)hpp->nof_complete_not_edited * 100 / total_segments);
    gt_logger_log(logger, "  ... and edited:           %-7lu (%.2f%%)",
        hpp->nof_complete_edited,
        (double)hpp->nof_complete_edited * 100 / total_segments);
    gt_logger_log(logger, "- not processed:            %-7lu (%.2f%%)",
        hpp->nof_skipped, (double)hpp->nof_skipped * 100 / total_segments);
    gt_logger_log(logger, "- not mapping:              %-7lu (%.2f%%)",
        hpp->nof_unmapped, (double)hpp->nof_unmapped * 100 / total_segments);
  }
}

static int gt_hpol_processor_output_sorted_segments(GtHpolProcessor *hpp,
    GtError *err)
{
  const GtUchar *s;
  char *d;
  int next_rval;
  unsigned long len;
  GtStr *d_str = NULL;
  GtAlignedSegment *as;
  GtSplitter *spl = gt_splitter_new();
  gt_assert(hpp != NULL);
  gt_assert(hpp->processed_segments != NULL);
  gt_assert(hpp->reads_iter != NULL);
  d_str = gt_str_new();
  while ((next_rval = gt_seqiterator_next(hpp->reads_iter, &s, &len, &d, err))
      > 0)
  {
    gt_str_set(d_str, d);
    gt_splitter_split(spl, d, (unsigned long)strlen(d), ' ');
    if ((as = gt_hashmap_get(hpp->processed_segments,
            gt_splitter_get_token(spl, 0))) != NULL)
    {
      gt_hpol_processor_output_segment(as, true, hpp->outfp_segments,
          gt_str_get(d_str));
    }
    else
    {
      gt_warning("ID not found: %s", gt_splitter_get_token(spl, 0));
    }
  }
  gt_str_delete(d_str);
  gt_splitter_delete(spl);
  return next_rval;
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
      had_err = gt_seqpos_classifier_position_is_inside_feature(
          hpp->cds_oracle, i, &coding, &end_of_annotation, err);
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
  gt_hpol_processor_process_hpol_end(hpp, prev, i - 1UL, hlen,
              hpp->max_hlen_diff, hpp->min_alt_consensus);
  gt_aligned_segments_pile_flush(hpp->asp, true);
  if (!had_err && hpp->processed_segments != NULL)
    had_err = gt_hpol_processor_output_sorted_segments(hpp, err);
  if (logger != NULL && !had_err)
    gt_hpol_processor_show_hdist(hpp, logger);
  return had_err;
}

void gt_hpol_processor_delete(GtHpolProcessor *hpp)
{
  if (hpp != NULL)
  {
    gt_disc_distri_delete(hpp->hdist);
    gt_disc_distri_delete(hpp->hdist_e);
    gt_hashmap_delete(hpp->processed_segments);
    gt_alphabet_delete(hpp->alpha);
    gt_free(hpp);
  }
}
