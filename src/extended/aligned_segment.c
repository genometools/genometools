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

#include <string.h>
#include "core/log_api.h"
#include "core/ma.h"
#include "core/undef_api.h"
#include "extended/aligned_segment.h"

struct GtAlignedSegment
{
  char *s, *q, *r, *d; /* segment sequence, segment qualities, reference region,
                          segment description */
  unsigned long alen;
  /* leftmost and rightmost coordinates of ref region on ref sequence */
  unsigned long r_left, r_right;
  /* true if ref region is on the strand opposite to ref sequence
   * (left and right still apply to the ref sequence strand) */
  bool r_reverse;
  bool has_indels;
  bool s_edited, r_edited;
  char *s_orig;
  unsigned long orig_seqlen;
  unsigned long mapq;
};

static unsigned long gt_aligned_segment_cigar2alen(GtSamAlignment *sa)
{
  unsigned long alen;
  uint16_t clen, opnum;
  unsigned char opcode;
  alen = 0;
  clen = gt_sam_alignment_cigar_length(sa);
  for (opnum = 0; opnum < clen; opnum++)
  {
    opcode = gt_sam_alignment_cigar_i_operation(sa, opnum);
    /* ignore hard clipping and padding operations */
    if (opcode != (unsigned char)'H' || opcode != (unsigned char)'P')
      alen += (unsigned long)gt_sam_alignment_cigar_i_length(sa, opnum);
  }
  gt_assert(alen >= gt_sam_alignment_read_length(sa));
  return alen;
}

/* s and q are stored at the end of the buffers */
static void gt_aligned_segment_fetch_s_and_q_from_sa(GtAlignedSegment *as,
    GtSamAlignment *sa)
{
  unsigned long alen = as->alen, slen = gt_sam_alignment_read_length(sa);
  GtUchar *s = (GtUchar*) (as->s + (alen - slen)),
          *q = (GtUchar*) (as->q + (alen - slen));
  GtAlphabet *alpha;
  gt_sam_alignment_sequence_external_buffer(sa, &s, &alen);
  gt_assert(as->alen == alen);
  gt_assert(s == (GtUchar*) (as->s + (alen - slen)));
  gt_sam_alignment_qualitystring_external_buffer(sa, &q,
      &alen);
  gt_assert(as->alen == alen);
  gt_assert(q == (GtUchar*) (as->q + (alen - slen)));
  alpha = gt_alphabet_new_dna();
  gt_alphabet_decode_seq_to_cstr(alpha, (char*) s, s, slen);
  gt_alphabet_delete(alpha);
}

static void gt_aligned_segment_align_using_cigar(GtAlignedSegment *as,
    GtSamAlignment *sa)
{
  unsigned long srcpos, pos = 0;
  uint16_t clen, opnum, oplen, i;
  unsigned char opcode;
  gt_assert(as->alen >= gt_sam_alignment_read_length(sa));
  srcpos = as->alen - gt_sam_alignment_read_length(sa);
  clen = gt_sam_alignment_cigar_length(sa);
  for (opnum = 0; opnum < clen; opnum++)
  {
    opcode = gt_sam_alignment_cigar_i_operation(sa, opnum);
    oplen = gt_sam_alignment_cigar_i_length(sa, opnum);
    switch (opcode)
    {
      case 'S':
        if (opnum == 0)
          as->r_left -= oplen;
        else if (opnum == clen - 1UL)
          as->r_right += oplen;
        else
          gt_assert(false); /* not allowed by SAM format specification */
        /*@fallthrough@*/
      case 'X': /*@fallthrough@*/
      case '=': /*@fallthrough@*/
      case 'M':
        for (i = 0; i < oplen; i++)
        {
          if (pos != srcpos)
          {
            gt_assert(pos < srcpos);
            gt_assert(as->s[srcpos] != '-');
            gt_assert(as->q[srcpos] != GT_UNDEF_CHAR);
            as->s[pos] = as->s[srcpos];
            as->q[pos] = as->q[srcpos];
          }
          as->r[pos] = (opcode == '=') ? as->s[srcpos] : '?';
          pos++;
          srcpos++;
        }
        break;
      case 'I':
        for (i = 0; i < oplen; i++)
        {
          if (pos != srcpos)
          {
            gt_assert(pos < srcpos);
            gt_assert(as->s[srcpos] != '-');
            gt_assert(as->q[srcpos] != GT_UNDEF_CHAR);
            as->s[pos] = as->s[srcpos];
            as->q[pos] = as->q[srcpos];
          }
          as->r[pos] = '-';
          pos++;
          srcpos++;
        }
        as->has_indels = true;
        break;
      case 'N': /*@fallthrough@*/
      case 'D':
        for (i = 0; i < oplen; i++)
        {
          as->s[pos] = '-';
          as->q[pos] = GT_UNDEF_CHAR;
          as->r[pos] = '?';
          pos++;
        }
        as->has_indels = true;
        break;
      case 'P': /*@fallthrough@*/
      case 'H':
        /* ignore padding and hard clipping */
        break;
      default:
        exit(GT_EXIT_PROGRAMMING_ERROR);
    }
  }
}

static void gt_aligned_segment_init_from_unmapped_sa(GtAlignedSegment *as,
    GtSamAlignment *sa)
{
  as->r_left = GT_UNDEF_ULONG;
  as->r_right = GT_UNDEF_ULONG;
  as->alen = gt_sam_alignment_read_length(sa);
  as->s = gt_malloc(sizeof (*as->s) * (as->alen + 1UL));
  as->q = gt_malloc(sizeof (*as->q) * (as->alen + 1UL));
  as->r = NULL;
  as->s[as->alen] = 0;
  as->q[as->alen] = 0;
  gt_aligned_segment_fetch_s_and_q_from_sa(as, sa);
}

static void gt_aligned_segment_init_from_mapped_sa(GtAlignedSegment *as,
    GtSamAlignment *sa, GtSamfileEncseqMapping *sem)
{
  as->r_left = gt_samfile_encseq_mapping_seqpos(sem,
      gt_sam_alignment_ref_num(sa), gt_sam_alignment_pos(sa));
  as->r_right = gt_samfile_encseq_mapping_seqpos(sem,
      gt_sam_alignment_ref_num(sa), gt_sam_alignment_rightmost_pos(sa));
  as->alen = gt_aligned_segment_cigar2alen(sa);
  as->s = gt_malloc(sizeof (*as->s) * (as->alen + 1UL));
  as->q = gt_malloc(sizeof (*as->q) * (as->alen + 1UL));
  as->r = gt_malloc(sizeof (*as->r) * (as->alen + 1UL));
  as->s[as->alen] = 0;
  as->q[as->alen] = 0;
  as->r[as->alen] = 0;
  gt_aligned_segment_fetch_s_and_q_from_sa(as, sa);
  gt_aligned_segment_align_using_cigar(as, sa);
}

GtAlignedSegment *gt_aligned_segment_new_from_sa(GtSamAlignment *sa,
    GtSamfileEncseqMapping *sem)
{
  GtAlignedSegment *as;
  size_t dlen;
  as = gt_malloc(sizeof (GtAlignedSegment));
  gt_assert(sa != NULL);
  as->r_reverse = gt_sam_alignment_is_reverse(sa);
  as->has_indels = false;
  if (gt_sam_alignment_is_unmapped(sa))
    gt_aligned_segment_init_from_unmapped_sa(as, sa);
  else
    gt_aligned_segment_init_from_mapped_sa(as, sa, sem);
  dlen = strlen(gt_sam_alignment_identifier(sa)) + 1UL;
  as->d = gt_malloc(sizeof (*as->d) * dlen);
  (void)memcpy(as->d, gt_sam_alignment_identifier(sa), sizeof (*as->d) * dlen);
  as->s_edited = false;
  as->r_edited = false;
  as->s_orig = NULL;
  as->mapq = gt_sam_alignment_mapping_quality(sa);
  as->orig_seqlen = gt_sam_alignment_read_length(sa);
  return as;
}

void gt_aligned_segment_enable_edit_tracking(GtAlignedSegment *as)
{
  gt_assert(as != NULL);
  gt_assert(as->s_orig == NULL);
  as->s_orig = gt_malloc(sizeof (*as->s_orig) * (as->alen + 1UL));
  memcpy(as->s_orig, as->s, (size_t)(as->alen + 1UL));
}

unsigned long gt_aligned_segment_mapping_quality(GtAlignedSegment *as)
{
  gt_assert(as != NULL);
  return as->mapq;
}

const char *gt_aligned_segment_orig_seq(GtAlignedSegment *as)
{
  gt_assert(as != NULL);
  return as->s_orig;
}

char *gt_aligned_segment_seq(GtAlignedSegment *as)
{
  gt_assert(as != NULL);
  return as->s;
}

char *gt_aligned_segment_qual(GtAlignedSegment *as)
{
  gt_assert(as != NULL);
  return as->q;
}

char *gt_aligned_segment_refregion(GtAlignedSegment *as)
{
  gt_assert(as != NULL);
  return as->r;
}

unsigned long gt_aligned_segment_length(const GtAlignedSegment *as)
{
  gt_assert(as != NULL);
  return as->alen;
}

unsigned long gt_aligned_segment_refregion_startpos(const GtAlignedSegment *as)
{
  gt_assert(as != NULL);
  return as->r_left;
}

unsigned long gt_aligned_segment_refregion_endpos(const GtAlignedSegment *as)
{
  gt_assert(as != NULL);
  return as->r_right;
}

bool gt_aligned_segment_is_reverse(const GtAlignedSegment *as)
{
  gt_assert(as != NULL);
  return as->r_reverse;
}

unsigned long gt_aligned_segment_offset_for_refpos(const GtAlignedSegment *as,
    unsigned long refpos)
{
  unsigned long r_offset, pos, ungapped_pos;
  gt_assert(as != NULL);
  if (refpos < as->r_left || refpos > as->r_right)
    return GT_UNDEF_ULONG;
  r_offset = refpos - as->r_left;
  pos = 0;
  ungapped_pos = 0;
  while (ungapped_pos < r_offset)
  {
    if (as->r[pos] != '-')
      ungapped_pos++;
    pos++;
  }
  gt_assert(pos <= as->alen);
  return pos;
}

unsigned long gt_aligned_segment_orig_seqlen(const GtAlignedSegment *as)
{
  gt_assert(as != NULL);
  return as->orig_seqlen;
}

unsigned long gt_aligned_segment_orig_seqpos_for_refpos(
    const GtAlignedSegment *as, unsigned long refpos)
{
  unsigned long r_offset, gapped_pos, ungapped_pos_on_r, ungapped_pos_on_s_orig;
  gt_assert(as != NULL);
  gt_assert(as->s_orig != NULL);
  if (refpos < as->r_left || refpos > as->r_right)
    return GT_UNDEF_ULONG;
  r_offset = refpos - as->r_left;
  gapped_pos = 0;
  ungapped_pos_on_r = 0;
  ungapped_pos_on_s_orig = 0;
  while (ungapped_pos_on_r < r_offset)
  {
    if (as->r[gapped_pos] != '-')
      ungapped_pos_on_r++;
    if (as->s_orig[gapped_pos] != '-')
      ungapped_pos_on_s_orig++;
    gapped_pos++;
  }
  gt_assert(gapped_pos <= as->alen);
  if (as->r_reverse)
    return as->orig_seqlen - 1UL - ungapped_pos_on_s_orig;
  else
    return ungapped_pos_on_s_orig;
}

void gt_aligned_segment_ungap_refregion(GtAlignedSegment *as)
{
  unsigned long pos, srcpos;
  gt_assert(as != NULL);
  gt_assert(as->r != NULL);
  for (srcpos = 0, pos = 0; srcpos < as->alen; srcpos++)
  {
    if (as->r[srcpos] != '-')
    {
      if (pos != srcpos)
      {
        gt_assert(pos < srcpos);
        as->r[pos] = as->r[srcpos];
      }
      pos++;
    }
  }
  gt_assert(pos <= as->alen + 1UL);
  if (pos <= as->alen)
    as->r[pos] = '\0';
}

void gt_aligned_segment_ungap_seq_and_qual(GtAlignedSegment *as)
{
  unsigned long pos, srcpos;
  gt_assert(as != NULL);
  for (srcpos = 0, pos = 0; srcpos < as->alen; srcpos++)
  {
    if (as->s[srcpos] != '-')
    {
      if (pos != srcpos)
      {
        gt_assert(pos < srcpos);
        as->s[pos] = as->s[srcpos];
        gt_assert(as->q[srcpos] != GT_UNDEF_CHAR);
        as->q[pos] = as->q[srcpos];
      }
      pos++;
    }
  }
  gt_assert(pos <= as->alen + 1UL);
  if (pos <= as->alen)
  {
    as->s[pos] = '\0';
    as->q[pos] = '\0';
  }
}

const char *gt_aligned_segment_description(const GtAlignedSegment *as)
{
  gt_assert(as != NULL);
  return as->d;
}

void gt_aligned_segment_seq_set_edited(GtAlignedSegment *as)
{
  gt_assert(as != NULL);
  as->s_edited = true;
}

bool gt_aligned_segment_seq_edited(const GtAlignedSegment *as)
{
  gt_assert(as != NULL);
  return as->s_edited;
}

void gt_aligned_segment_refregion_set_edited(GtAlignedSegment *as)
{
  gt_assert(as != NULL);
  as->r_edited = true;
}

bool gt_aligned_segment_refregion_edited(const GtAlignedSegment *as)
{
  gt_assert(as != NULL);
  return as->r_edited;
}

bool gt_aligned_segment_has_indels(const GtAlignedSegment *as)
{
  gt_assert(as != NULL);
  return as->has_indels;
}

void gt_aligned_segment_delete(GtAlignedSegment *as)
{
  if (as != NULL)
  {
    gt_free(as->s);
    gt_free(as->q);
    gt_free(as->r);
    gt_free(as->d);
    gt_free(as->s_orig);
    gt_free(as);
  }
}

void gt_aligned_segment_assign_refregion_chars(GtAlignedSegment *as,
    GtEncseq *encseq)
{
  unsigned long i, pos;
  gt_assert(as != NULL);
  gt_assert(as->r != NULL);
  for (pos = as->r_left, i = 0; i < as->alen; i++)
  {
    if (as->r[i] == '?')
    {
      as->r[i] = gt_encseq_get_decoded_char(encseq, pos, GT_READMODE_FORWARD);
    }
    if (as->r[i] != '-')
      pos++;
  }
}

void gt_aligned_segment_show(GtAlignedSegment *as, GtFile *outfp)
{
  gt_assert(as != NULL);
  if (as->d != NULL)
    gt_file_xprintf(outfp, "D: %s\n", as->d);
  if (as->r != NULL)
    gt_file_xprintf(outfp, "R: %s\n", as->r);
  if (as->s_orig != NULL)
    gt_file_xprintf(outfp, "O: %s\n", as->s_orig);
  gt_file_xprintf(outfp, "S: %s\n", as->s);
  gt_file_xprintf(outfp, "Q: %s\n", as->q);
}
