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
#include "core/disc_distri_api.h"
#include "core/fastq.h"
#include "core/hashmap_api.h"
#include "core/log_api.h"
#include "core/ma.h"
#include "core/minmax.h"
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
  unsigned long hmin, clenmax, read_hmin, mapqmin, covmin;
  bool allow_partial, allow_multiple;
  double altmax, refmin;
  GtDiscDistri *hdist, *hdist_e;
  GtSeqposClassifier *cds_oracle;
  GtAlignedSegmentsPile *asp;
  unsigned long nof_complete_edited, nof_complete_not_edited,
                nof_skipped, nof_unmapped, nof_h, nof_h_e, hlen_max,
                nof_multihits, nof_replaced;
  bool adjust_s_hlen;
  GtFile *outfp_segments, *outfp_stats, **outfiles;
  GtSeqIterator **reads_iters;
  unsigned long nfiles;
  GtHashmap *processed_segments;
  GtAlphabet *alpha;
  bool output_segments, output_stats, output_multihit_stats;
};

GtHpolProcessor *gt_hpol_processor_new(GtEncseq *encseq, unsigned long hmin)
{
  GtHpolProcessor *hpp;
  hpp = gt_malloc(sizeof (GtHpolProcessor));
  gt_assert(encseq != NULL);
  gt_assert(hmin > 0);
  hpp->encseq = encseq;
  hpp->hmin = hmin;
  hpp->read_hmin = 0;
  hpp->mapqmin = 0;
  hpp->covmin = 0;
  hpp->allow_partial = false;
  hpp->allow_multiple = false;
  hpp->hdist = gt_disc_distri_new();
  hpp->nof_h = 0;
  hpp->hdist_e = gt_disc_distri_new();
  hpp->nof_h_e = 0;
  hpp->hlen_max = 0;
  hpp->cds_oracle = NULL;
  hpp->asp = NULL;
  hpp->nof_complete_edited = 0;
  hpp->nof_complete_not_edited = 0;
  hpp->nof_skipped = 0;
  hpp->nof_unmapped = 0;
  hpp->nof_multihits = 0;
  hpp->nof_replaced = 0;
  hpp->clenmax = GT_UNDEF_ULONG;
  hpp->altmax = (double) 1.0;
  hpp->refmin = (double) 0.0;
  hpp->alpha = gt_alphabet_new_dna();
  hpp->adjust_s_hlen = false;
  hpp->output_segments = false;
  hpp->outfp_segments = NULL;
  hpp->output_stats = false;
  hpp->output_multihit_stats = false;
  hpp->outfp_stats = NULL;
  hpp->processed_segments = NULL;
  hpp->reads_iters = NULL;
  hpp->outfiles = NULL;
  hpp->nfiles = 0;
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

typedef enum {
  GT_HPOL_PROCESSOR_NEW_RECORD,
  GT_HPOL_PROCESSOR_REPLACED,
  GT_HPOL_PROCESSOR_NOT_REPLACED,
} GtHpolProcessorAddToHashResult;

static GtHpolProcessorAddToHashResult gt_hpol_processor_add_segment_to_hashmap(
    GtHpolProcessor *hpp, GtAlignedSegment *as)
{
  GtAlignedSegment *stored_as;
  if ((stored_as = gt_hashmap_get(hpp->processed_segments,
          gt_aligned_segment_description(as))) != NULL)
  {
    hpp->nof_multihits++;
    if (!gt_aligned_segment_seq_edited(stored_as) &&
        gt_aligned_segment_seq_edited(as))
    {
      hpp->nof_replaced++;
      /* change with newly edited as */
      gt_hashmap_remove(hpp->processed_segments,
          (void*)gt_aligned_segment_description(as));
      gt_hashmap_add(hpp->processed_segments,
          (void*)gt_aligned_segment_description(as), as);
      return GT_HPOL_PROCESSOR_REPLACED;
    }
    /* otherwise discard (todo: implement combination of edits) */
    else
    {
      return GT_HPOL_PROCESSOR_NOT_REPLACED;
    }
  }
  else
  {
    gt_hashmap_add(hpp->processed_segments,
        (void*)gt_aligned_segment_description(as), as);
    return GT_HPOL_PROCESSOR_NEW_RECORD;
  }
}

static void gt_hpol_processor_process_complete_segment(
    GtAlignedSegment *as, void *data)
{
  GtHpolProcessorAddToHashResult multihit = GT_HPOL_PROCESSOR_NEW_RECORD;
  GtHpolProcessor *hpp = data;
  gt_assert(hpp != NULL);
  if (hpp->output_segments)
    gt_hpol_processor_output_segment(as, gt_aligned_segment_has_indels(as),
        hpp->outfp_segments, NULL);
  if (hpp->processed_segments != NULL)
    multihit = gt_hpol_processor_add_segment_to_hashmap(hpp, as);
  if (multihit == GT_HPOL_PROCESSOR_NEW_RECORD)
  {
    if (gt_aligned_segment_seq_edited(as))
      (hpp->nof_complete_edited)++;
    else
      (hpp->nof_complete_not_edited)++;
  }
  else if (multihit == GT_HPOL_PROCESSOR_REPLACED)
  {
    gt_assert(gt_aligned_segment_seq_edited(as));
    (hpp->nof_complete_edited)++;
    gt_assert(hpp->nof_complete_not_edited > 0);
    (hpp->nof_complete_not_edited)--;
  }
  else if (multihit == GT_HPOL_PROCESSOR_NOT_REPLACED)
  {
    gt_aligned_segment_delete(as);
  }
}

static void gt_hpol_processor_process_skipped_segment(
    GtAlignedSegment *as, void *data)
{
  GtHpolProcessorAddToHashResult multihit = GT_HPOL_PROCESSOR_NEW_RECORD;
  GtHpolProcessor *hpp = data;
  gt_assert(hpp != NULL);
  if (hpp->output_segments)
    gt_hpol_processor_output_segment(as, gt_aligned_segment_has_indels(as),
        hpp->outfp_segments, NULL);
  if (hpp->processed_segments != NULL)
    multihit = gt_hpol_processor_add_segment_to_hashmap(hpp, as);
  gt_assert(multihit != GT_HPOL_PROCESSOR_REPLACED);
  if (multihit == GT_HPOL_PROCESSOR_NEW_RECORD)
    hpp->nof_skipped++;
  else if (multihit == GT_HPOL_PROCESSOR_NOT_REPLACED)
    gt_aligned_segment_delete(as);
}

static void gt_hpol_processor_process_unmapped_segment(
    GtAlignedSegment *as, void *data)
{
  GT_UNUSED GtHpolProcessorAddToHashResult multihit =
    GT_HPOL_PROCESSOR_NEW_RECORD;
  GtHpolProcessor *hpp = data;
  gt_assert(hpp != NULL);
  if (hpp->output_segments)
    gt_hpol_processor_output_segment(as, false, hpp->outfp_segments, NULL);
  if (hpp->processed_segments != NULL)
    multihit = gt_hpol_processor_add_segment_to_hashmap(hpp, as);
  gt_assert(multihit == GT_HPOL_PROCESSOR_NEW_RECORD);
  hpp->nof_unmapped++;
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
    GtAlignedSegmentsPile *asp, unsigned long read_hmin, double altmax,
    double refmin, unsigned long mapqmin, unsigned long covmin,
    bool allow_partial, bool allow_multiple, unsigned long clenmax)
{
  gt_assert(hpp != NULL);
  gt_assert(asp != NULL);
  gt_assert(hpp->asp == NULL);
  hpp->adjust_s_hlen = true;
  hpp->asp = asp;
  hpp->read_hmin = read_hmin;
  hpp->altmax = altmax;
  hpp->refmin = refmin;
  hpp->mapqmin = mapqmin;
  hpp->covmin = covmin;
  hpp->allow_partial = allow_partial;
  hpp->allow_multiple = allow_multiple;
  hpp->clenmax = clenmax;
  gt_aligned_segments_pile_register_process_complete(hpp->asp,
      gt_hpol_processor_process_complete_segment, hpp);
  gt_aligned_segments_pile_register_process_skipped(hpp->asp,
      gt_hpol_processor_process_skipped_segment, hpp);
  gt_aligned_segments_pile_register_process_unmapped(hpp->asp,
      gt_hpol_processor_process_unmapped_segment, hpp);
}

void gt_hpol_processor_enable_direct_segments_output(GtHpolProcessor *hpp,
    GtFile *outfile)
{
  gt_assert(hpp != NULL);
  hpp->output_segments = true;
  hpp->outfp_segments = outfile;
}

void gt_hpol_processor_enable_sorted_segments_output(GtHpolProcessor *hpp,
    unsigned long nfiles, GtSeqIterator **reads_iters, GtFile **outfiles)
{
  gt_assert(hpp != NULL);
  gt_aligned_segments_pile_disable_segment_deletion(hpp->asp);
  hpp->processed_segments = gt_hashmap_new(GT_HASH_STRING, NULL,
      (GtFree)gt_aligned_segment_delete);
  hpp->reads_iters = reads_iters;
  hpp->outfiles = outfiles;
  hpp->nfiles = nfiles;
}

static void gt_hpol_processor_output_stats_header(GtFile *outfp)
{
  gt_file_xprintf(outfp, "# correction statistics\n");
  gt_file_xprintf(outfp, "# r_hpos =    start pos of hpol on reference\n");
  gt_file_xprintf(outfp, "# edit =      edit operation on the read (I or D)\n");
  gt_file_xprintf(outfp, "# s_hpos =    first pos of hpol on read\n");
  gt_file_xprintf(outfp, "# s_hend =    last pos of hpol on read\n");
  gt_file_xprintf(outfp, "# s_char =    hpol character in read\n");
  gt_file_xprintf(outfp, "# s_or =      orientation of read "
      "(+ or -; + = same as reference)\n");
  gt_file_xprintf(outfp, "# c_len =     correction length\n");
  gt_file_xprintf(outfp, "# coverage =  number of reads over entire hpol\n");
  gt_file_xprintf(outfp, "# r_hlen =    length of hpol on reference\n");
  gt_file_xprintf(outfp, "# r_supp =    %% reads with ref hpol length\n");
  gt_file_xprintf(outfp, "# s_hlen =    length of hpol in read\n");
  gt_file_xprintf(outfp, "# a_hlen =    alt consensus hpol length in reads\n");
  gt_file_xprintf(outfp, "# a_supp =    %% reads with alt hpol length\n");
  gt_file_xprintf(outfp, "# s_mapq =    mapping quality of read\n");
  gt_file_xprintf(outfp, "# s_q_bef =   quality of base before the hpol\n");
  gt_file_xprintf(outfp, "# s_q_first = quality of first hpol base\n");
  gt_file_xprintf(outfp, "# s_q_min =   min quality among hpol bases\n");
  gt_file_xprintf(outfp, "# s_q_ave =   average quality of read in the hpol "
      "positions\n");
  gt_file_xprintf(outfp, "# s_q_max =   max quality among hpol bases\n");
  gt_file_xprintf(outfp, "# s_q_range = s_q_max - s_q_min + 1\n");
  gt_file_xprintf(outfp, "# s_q_last =  quality of last hpol base\n");
  gt_file_xprintf(outfp, "# s_q_aft =   quality of base after the hpol\n");
  gt_file_xprintf(outfp, "# s_qual =    quality string in read for the hpol "
      "positions\n");
  gt_file_xprintf(outfp, "# s_id =      read identifier\n");
  gt_file_xprintf(outfp, "# coordinates are 1-based\n");
  gt_file_xprintf(outfp, "#\n");
  gt_file_xprintf(outfp, "# r_hpos\tedit\ts_hpos\ts_hend\ts_char\ts_or\tc_len\t"
      "coverage\tr_hlen\tr_supp\ts_hlen\ta_hlen\ta_supp\ts_mapq\ts_q_bef\t"
      "s_q_first\ts_q_min\ts_q_ave\ts_q_max\ts_q_range\ts_q_last\ts_q_aft\t"
      "s_qual\ts_id\n");
}

#define GT_HPOL_PROCESSOR_PHREDOFFSET 33UL
#define GT_HPOL_PROCESSOR_QUAL(QCHAR) \
  (unsigned long)(QCHAR) - GT_HPOL_PROCESSOR_PHREDOFFSET;

static void gt_hpol_processor_output_stats(GtAlignedSegment *as,
    unsigned long r_hpos, unsigned long coverage, unsigned long r_hlen,
    unsigned long r_supp, unsigned long s_hlen, unsigned long a_hlen,
    unsigned long a_supp, char s_char, unsigned long s_q_sum,
    unsigned long c_len, GtFile *outfp)
{
  unsigned long i, pos, s_hpos = 0, s_offset, s_q_bef, s_q_aft, s_q_value,
                s_q_min, s_q_max, s_q_range, s_q_first, s_q_last = 0, s_hend,
                s_mapq;
  double s_q_ave;
  char *s_qual, *q, edit, s_or;
  const char *s_id;
  gt_assert(r_hlen != s_hlen);
  edit = r_hlen > s_hlen ? 'I' : 'D';
  gt_assert(coverage > 0);
  r_supp = r_supp * 100 / coverage;
  a_supp = a_supp * 100 / coverage;
  s_id = gt_aligned_segment_description(as);
  s_mapq = gt_aligned_segment_mapping_quality(as);
  q = gt_aligned_segment_qual(as);
  gt_assert(s_hlen > 0);
  gt_assert(s_q_sum >= GT_HPOL_PROCESSOR_PHREDOFFSET * s_hlen);
  s_q_ave = (double)(s_q_sum - GT_HPOL_PROCESSOR_PHREDOFFSET * s_hlen)
    / (double)s_hlen;
  s_hpos = gt_aligned_segment_orig_seqpos_for_refpos(as, r_hpos);
  s_offset = gt_aligned_segment_offset_for_refpos(as, r_hpos);
  s_qual = gt_malloc(sizeof (*s_qual) * (s_hlen + 1UL));
  s_q_bef = GT_UNDEF_ULONG;
  for (i = s_offset; i > 0; /**/)
  {
    i--;
    if (q[i] != GT_UNDEF_CHAR)
    {
      s_q_bef = GT_HPOL_PROCESSOR_QUAL(q[i]);
      break;
    }
  }
  gt_assert(s_q_bef != GT_UNDEF_ULONG);
  s_q_min = ULONG_MAX;
  s_q_max = 0;
  s_q_first = GT_UNDEF_ULONG;
  if (!gt_aligned_segment_is_reverse(as))
  {
    for (i = s_offset, pos = 0; pos < s_hlen; i++)
    {
      if (q[i] != GT_UNDEF_CHAR)
      {
        s_qual[pos] = q[i];
        pos++;
        s_q_value = GT_HPOL_PROCESSOR_QUAL(q[i]);
        if (s_q_value < s_q_min)
          s_q_min = s_q_value;
        if (s_q_value > s_q_max)
          s_q_max = s_q_value;
        if (s_q_first == GT_UNDEF_ULONG)
          s_q_first = s_q_value;
        s_q_last = s_q_value;
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
        s_q_value = GT_HPOL_PROCESSOR_QUAL(q[i]);
        if (s_q_value < s_q_min)
          s_q_min = s_q_value;
        if (s_q_value > s_q_max)
          s_q_max = s_q_value;
        if (s_q_first == GT_UNDEF_ULONG)
          s_q_first = s_q_value;
        s_q_last = s_q_value;
      }
    }
  }
  s_qual[s_hlen] = '\0';
  s_q_aft = GT_UNDEF_ULONG;
  for (/**/; i < gt_aligned_segment_length(as); i++)
  {
    if (q[i] != GT_UNDEF_CHAR)
    {
      s_q_aft = GT_HPOL_PROCESSOR_QUAL(q[i]);
      break;
    }
  }
  gt_assert(s_q_aft != GT_UNDEF_ULONG);
  gt_assert(s_q_min < ULONG_MAX);
  gt_assert(s_q_max >= s_q_min);
  s_q_range = s_q_max - s_q_min + 1UL;
  /* convert to 1-based coordinates */
  r_hpos++;
  s_hpos++;
  /* handle reverse alignments */
  if (gt_aligned_segment_is_reverse(as))
  {
    /* complement char */
    GtError *err = gt_error_new();
    (void)gt_complement(&s_char, s_char, err);
    gt_error_delete(err);
    /* correct coords on s */
    s_hend = s_hpos;
    s_hpos = s_hpos - s_hlen + 1UL;
    s_or = '-';
    /* swap q values */
    s_q_value = s_q_aft;
    s_q_aft = s_q_bef;
    s_q_bef = s_q_value;
    s_q_value = s_q_last;
    s_q_last = s_q_first;
    s_q_first = s_q_value;
  }
  else
  {
    s_hend = s_hpos + s_hlen - 1UL;
    s_or = '+';
  }
  gt_file_xprintf(outfp,
      "%lu\t%c\t%lu\t%lu\t%c\t%c\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t"
      "%lu\t%lu\t%lu\t%.2f\t%lu\t%lu\t%lu\t%lu\t%s\t%s\n",
      r_hpos, edit, s_hpos, s_hend, s_char, s_or, c_len, coverage, r_hlen,
      r_supp, s_hlen, a_hlen, a_supp, s_mapq, s_q_bef, s_q_first, s_q_min,
      s_q_ave, s_q_max, s_q_range, s_q_last, s_q_aft, s_qual, s_id);
  gt_free(s_qual);
}

void gt_hpol_processor_enable_statistics_output(GtHpolProcessor *hpp,
    bool output_multihit_stats, GtFile *outfile)
{
  gt_assert(hpp != NULL);
  hpp->output_stats = true;
  hpp->outfp_stats = outfile;
  hpp->output_multihit_stats = output_multihit_stats;
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
    char c, unsigned long r_hstart, unsigned long coverage,
    unsigned long r_hlen, unsigned long r_supp, unsigned long a_hlen,
    unsigned long a_supp, unsigned long clenmax, bool allow_partial,
    bool allow_multiple, unsigned long s_hmin, bool output_stats,
    GtFile *outfp_stats)
{
  unsigned long left, right, s_hlen = 0, q_sum = 0, s_free = 0;
  char *s, *q;
  bool edited = false;
  left = gt_aligned_segment_offset_for_refpos(as, r_hstart);
  right = gt_aligned_segment_offset_for_refpos(as, r_hstart + r_hlen);
  if (left == GT_UNDEF_ULONG || left == 0 ||
      right == GT_UNDEF_ULONG || right == gt_aligned_segment_length(as))
    return edited;
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
  if (s_hlen == 0 || s_hlen < s_hmin)
    return edited;
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
      if (hlen_diff < clenmax &&
          (s_free >= hlen_diff || allow_partial) &&
          (!gt_aligned_segment_seq_edited(as) || allow_multiple))
      {
        if (output_stats)
        {
          gt_hpol_processor_output_stats(as, r_hstart, coverage, r_hlen, r_supp,
              s_hlen, a_hlen, a_supp, c, q_sum, MIN(s_free, hlen_diff),
              outfp_stats);
        }
        gt_aligned_segment_seq_set_edited(as);
        gt_hpol_processor_enlarge_hpol(s, q, left, right,
            MIN(s_free, hlen_diff), c, (char)(q_sum / s_hlen));
        edited = true;
      }
    }
  }
  else if (s_hlen > r_hlen)
  {
    unsigned long hlen_diff = s_hlen - r_hlen;
    if (hlen_diff < clenmax &&
       (!gt_aligned_segment_seq_edited(as) || allow_multiple))
    {
      if (output_stats)
      {
        gt_hpol_processor_output_stats(as, r_hstart, coverage, r_hlen, r_supp,
            s_hlen, a_hlen, a_supp, c, q_sum, hlen_diff, outfp_stats);
      }
      gt_aligned_segment_seq_set_edited(as);
      gt_hpol_processor_shrink_hpol(s, q, left, right, hlen_diff, c);
      edited = true;
    }
  }
  return edited;
}

static bool gt_hpol_processor_adjust_hlen_of_all_segments(
    GtAlignedSegmentsPile *asp, char c, unsigned long r_hstart,
    unsigned long coverage, unsigned long r_hlen, unsigned long r_supp,
    unsigned long a_hlen, unsigned long a_supp, unsigned long clenmax,
    bool allow_partial, bool allow_multiple, unsigned long s_hmin,
    unsigned long mapqmin, bool output_stats, GtFile *outfp_stats,
    GtHashmap *processed_segments, bool output_multihit_stats)
{
  GtDlistelem *dlistelem;
  bool any_edited = false, edited;
  gt_assert(asp != NULL);
  for (dlistelem = gt_dlist_first(gt_aligned_segments_pile_get(asp));
        dlistelem != NULL; dlistelem = gt_dlistelem_next(dlistelem))
  {
    GtAlignedSegment *as;
    as = gt_dlistelem_get_data(dlistelem);
    if (gt_aligned_segment_has_indels(as) &&
        gt_aligned_segment_mapping_quality(as) >= mapqmin)
    {
      bool output_stats_for_as = output_stats;
      if (output_stats_for_as && !output_multihit_stats && processed_segments)
      {
        GtAlignedSegment *stored_as;
        /* if processed_segment exists, disable output_stats if stored_as
         * exists and was edited */
        if ((stored_as = gt_hashmap_get(processed_segments,
                gt_aligned_segment_description(as))) != NULL &&
            gt_aligned_segment_seq_edited(stored_as))
          output_stats_for_as = false;
      }
      edited = gt_hpol_processor_adjust_hlen_of_a_segment(as, c, r_hstart,
          coverage, r_hlen, r_supp, a_hlen, a_supp, clenmax, allow_partial,
          allow_multiple, s_hmin, output_stats_for_as, outfp_stats);
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
    GtUchar c, unsigned long endpos, unsigned long hlen)
{
  bool edited = false;
  gt_disc_distri_add(hpp->hdist, hlen);
  hpp->nof_h++;
  if (hlen > hpp->hlen_max)
    hpp->hlen_max = hlen;
  if (hpp->adjust_s_hlen)
  {
    char ch = gt_alphabet_decode(hpp->alpha, c);
    unsigned long a_supp_max = ULONG_MAX, a_supp = 0, piled, r_supp = 0, a_hlen,
                  r_supp_min = 0;
    gt_aligned_segments_pile_move_over_position(hpp->asp, endpos + 1UL);
    piled = gt_aligned_segments_pile_size(hpp->asp);
    if (piled >= hpp->covmin)
    {
      gt_hpol_processor_determine_alternative_consensus(hpp, ch,
          endpos + 1UL - hlen, hlen, &a_hlen, &a_supp, &piled,
          &r_supp);
      a_supp_max = (unsigned long)(hpp->altmax * (double)piled);
      r_supp_min = (unsigned long)(hpp->refmin * (double)piled);
      if (r_supp < piled && r_supp >= r_supp_min && a_supp <= a_supp_max)
      {
        edited = gt_hpol_processor_adjust_hlen_of_all_segments(hpp->asp, ch,
            endpos + 1UL - hlen, piled, hlen, r_supp, a_hlen, a_supp,
            hpp->clenmax, hpp->allow_partial, hpp->allow_multiple,
            hpp->read_hmin, hpp->mapqmin, hpp->output_stats, hpp->outfp_stats,
            hpp->processed_segments, hpp->output_multihit_stats);
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
        "in reference sequence"));
  gt_logger_log(logger, "length\toccurrences\tedited");
  for (i = hpp->hmin; i <= hpp->hlen_max; i++)
  {
    unsigned long n = (unsigned long)gt_disc_distri_get(hpp->hdist, i);
    if (n > 0)
    {
      gt_logger_log(logger, "%-6lu\t%-11lu\t%-6lu\t(%.2f%%)", i, n,
          (unsigned long)gt_disc_distri_get(hpp->hdist_e, i),
          (double)gt_disc_distri_get(hpp->hdist_e, i) * 100 / n);
    }
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
    unsigned long nof_complete = hpp->nof_complete_edited +
      hpp->nof_complete_not_edited;
    unsigned long total_segments =
      nof_complete + hpp->nof_skipped + hpp->nof_unmapped;
    gt_logger_log(logger, "segments in SAM file:       %lu",
        total_segments);
    gt_logger_log(logger, "- processed:                %-7lu (%.2f%%)",
        nof_complete, (double)nof_complete * 100 / total_segments);
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
    if (hpp->processed_segments != NULL)
    {
      gt_logger_log(logger, "- multiple hits:            %-7lu",
          hpp->nof_multihits);
      gt_log_log("replacements in hashmap: %lu", hpp->nof_replaced);
    }
  }
}

static int gt_hpol_processor_output_sorted_segments(GtHpolProcessor *hpp,
    GtSeqIterator *reads_iter, GtFile *outfile, GtError *err)
{
  const GtUchar *s;
  char *d;
  int next_rval;
  unsigned long len, i;
  GtStr *d_str = NULL;
  GtAlignedSegment *as;
  gt_assert(hpp != NULL);
  gt_assert(hpp->processed_segments != NULL);
  gt_assert(reads_iter != NULL);
  d_str = gt_str_new();
  while ((next_rval = gt_seq_iterator_next(reads_iter, &s, &len, &d, err))
      > 0)
  {
    gt_str_set(d_str, d);
    for (i = 0; i < (unsigned long)strlen(d); i++)
    {
      if (d[i] == ' ')
        d[i] = '\0';
    }
    if ((as = gt_hashmap_get(hpp->processed_segments, d)) != NULL)
    {
      gt_hpol_processor_output_segment(as, true, outfile, gt_str_get(d_str));
    }
    else
    {
      gt_warning("ID not found: %s", d);
    }
  }
  gt_str_delete(d_str);
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
          gt_hpol_processor_process_hpol_end(hpp, prev, i - 1UL, hlen);
        hlen = 1UL;
      }
      prev = c;
    }
  }
  if (!had_err)
  {
    if (hlen >= hpp->hmin && (hpp->cds_oracle == NULL || coding))
      gt_hpol_processor_process_hpol_end(hpp, prev, i - 1UL, hlen);
  }
  gt_encseq_reader_delete(esr);
  gt_aligned_segments_pile_flush(hpp->asp, true);
  if (!had_err && hpp->processed_segments != NULL)
  {
    for (i = 0; i < hpp->nfiles; i++)
      had_err = gt_hpol_processor_output_sorted_segments(hpp,
          hpp->reads_iters[i], hpp->outfiles[i], err);
  }
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
