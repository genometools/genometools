/*
  Copyright (c) 2011      Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
  Copyright (c)      2012 Giorgio Gonnella <gonnella@zbh.uni-hamburg.de>
  Copyright (c) 2011-2012 Center for Bioinformatics, University of Hamburg

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

#include <sam.h>
#include "core/alphabet_api.h"
#include "core/ensure.h"
#include "core/ma_api.h"
#include "core/types_api.h"
#include "core/undef_api.h"
#include "extended/sam_alignment_rep.h"
#include "extended/sam_alignment.h"

#define BAMBASEA 1
#define BAMBASEC 2
#define BAMBASEG 4
#define BAMBASET 8
#define BAMBASEN 15

#define PHREDOFFSET 33

GtSamAlignment *gt_sam_alignment_new(GtAlphabet *alphabet)
{
  GtSamAlignment *sam_alignment;
  gt_assert(alphabet != NULL);

  sam_alignment = gt_malloc(sizeof (GtSamAlignment));
  sam_alignment->s_alignment = bam_init1();
  sam_alignment->alphabet = gt_alphabet_ref(alphabet);
  sam_alignment->seq_buffer = NULL;
  sam_alignment->qual_buffer = NULL;
  sam_alignment->s_bufsize = 0;
  sam_alignment->q_bufsize = 0;
  sam_alignment->rightmost = GT_UNDEF_ULONG;
  return sam_alignment;
}

GtSamAlignment *gt_sam_alignment_clone(GtSamAlignment *sam_alignment)
{
  GtSamAlignment *clone;
  gt_assert(sam_alignment != NULL);
  clone = gt_sam_alignment_new(sam_alignment->alphabet);
  (void)bam_copy1(clone->s_alignment, sam_alignment->s_alignment);
  clone->rightmost = sam_alignment->rightmost;
  return clone;
}

void gt_sam_alignment_delete(GtSamAlignment *sam_alignment)
{
  if (sam_alignment != NULL) {
    bam_destroy1(sam_alignment->s_alignment);
    gt_alphabet_delete(sam_alignment->alphabet);
    gt_free(sam_alignment->qual_buffer);
    gt_free(sam_alignment->seq_buffer);
    gt_free(sam_alignment);
  }
}

static GtUchar bambase2gtbase(uint8_t base, GtAlphabet *alphabet)
{
  switch (base)
  {
    case BAMBASEA:
      return gt_alphabet_encode(alphabet, 'A');
    case BAMBASEC:
      return gt_alphabet_encode(alphabet, 'C');
    case BAMBASEG:
      return gt_alphabet_encode(alphabet, 'G');
    case BAMBASET:
      return gt_alphabet_encode(alphabet, 'T');
    default:
      return gt_alphabet_encode(alphabet,
                                (char) gt_alphabet_wildcard_show(alphabet));
  }
}

unsigned long gt_sam_alignment_pos(GtSamAlignment *sam_alignment)
{
  gt_assert(sam_alignment != NULL);
  return (unsigned long) sam_alignment->s_alignment->core.pos;
}

unsigned long gt_sam_alignment_read_length(GtSamAlignment *sam_alignment)
{
  gt_assert(sam_alignment != NULL);
  return (unsigned long) sam_alignment->s_alignment->core.l_qseq;
}

unsigned long gt_sam_alignment_mapping_quality(GtSamAlignment *sam_alignment)
{
  gt_assert(sam_alignment != NULL);
  return (unsigned long) sam_alignment->s_alignment->core.qual;
}

unsigned long gt_sam_alignment_rightmost_pos(GtSamAlignment *sam_alignment)
{
  gt_assert(sam_alignment != NULL);
  if (sam_alignment->rightmost == GT_UNDEF_ULONG)
  {
    sam_alignment->rightmost = (unsigned long)bam_calend(
        &sam_alignment->s_alignment->core,
        bam1_cigar(sam_alignment->s_alignment));
  }
  return sam_alignment->rightmost;
}

const char *gt_sam_alignment_identifier(GtSamAlignment *sam_alignment)
{
  gt_assert(sam_alignment != NULL);
  return bam1_qname(sam_alignment->s_alignment);
}

int32_t gt_sam_alignment_ref_num(GtSamAlignment *sam_alignment)
{
  gt_assert(sam_alignment != NULL);
  return sam_alignment->s_alignment->core.tid;
}

void gt_sam_alignment_sequence_external_buffer(GtSamAlignment *sam_alignment,
    GtUchar **seq_buffer, unsigned long *bufsize)
{
  unsigned long query_len, idx;
  uint8_t *bam_seq;

  query_len = (unsigned long) sam_alignment->s_alignment->core.l_qseq;
  if ((*bufsize) < query_len) {
    (*seq_buffer) = gt_realloc(*seq_buffer, sizeof (**seq_buffer) *
                                   (query_len + 1UL));
    (*bufsize) = query_len;
  }
  gt_assert((*seq_buffer) != NULL);

  bam_seq = bam1_seq(sam_alignment->s_alignment);

  for (idx = 0UL; idx < query_len; idx++) {
    (*seq_buffer)[idx] = bambase2gtbase((uint8_t) bam1_seqi(bam_seq, idx),
                                      sam_alignment->alphabet);
  }
  (*seq_buffer)[query_len] = (GtUchar)'\0';
}

const GtUchar *gt_sam_alignment_sequence(GtSamAlignment *sam_alignment)
{
  gt_sam_alignment_sequence_external_buffer(sam_alignment,
      &sam_alignment->seq_buffer, &sam_alignment->s_bufsize);
  return sam_alignment->seq_buffer;
}

void gt_sam_alignment_qualitystring_external_buffer(
    GtSamAlignment *sam_alignment, GtUchar **qual_buffer,
    unsigned long *bufsize)
{
  unsigned long query_len, idx;
  uint8_t *qual;

  gt_assert(qual_buffer != NULL);
  gt_assert(bufsize != NULL);
  query_len = (unsigned long) sam_alignment->s_alignment->core.l_qseq;
  if ((*bufsize) < query_len) {
    (*qual_buffer) = gt_realloc((*qual_buffer), sizeof (**qual_buffer) *
        (query_len + 1UL));
    (*bufsize) = query_len;
  }
  gt_assert((*qual_buffer) != NULL);
  qual = bam1_qual(sam_alignment->s_alignment);
  for (idx = 0UL; idx < query_len; idx++) {
    (*qual_buffer)[idx] = (GtUchar)(qual[idx] + PHREDOFFSET);
  }
  (*qual_buffer)[query_len] = (GtUchar)'\0';
}

const GtUchar *gt_sam_alignment_qualitystring(GtSamAlignment *sam_alignment)
{
  gt_sam_alignment_qualitystring_external_buffer(sam_alignment,
      &sam_alignment->qual_buffer, &sam_alignment->q_bufsize);
  return sam_alignment->qual_buffer;
}

uint16_t gt_sam_alignment_cigar_length(GtSamAlignment *sam_alignment)
{
  gt_assert(sam_alignment != NULL);
  return sam_alignment->s_alignment->core.n_cigar;
}

uint32_t gt_sam_alignment_cigar_i_length(GtSamAlignment *sam_alignment,
                                         uint16_t i)
{
  gt_assert(sam_alignment != NULL);
  return bam1_cigar(sam_alignment->s_alignment)[i]>>BAM_CIGAR_SHIFT;
}

unsigned char gt_sam_alignment_cigar_i_operation(GtSamAlignment *sam_alignment,
                                                 uint16_t i)
{
  gt_assert(sam_alignment != NULL);
  switch ((unsigned char) bam1_cigar(sam_alignment->s_alignment)[i]
                          & BAM_CIGAR_MASK) {
    case BAM_CMATCH:
      return 'M';
    case BAM_CINS:
      return 'I';
    case BAM_CDEL:
      return 'D';
    case BAM_CREF_SKIP:
      return 'N';
    case BAM_CSOFT_CLIP:
      return 'S';
    case BAM_CHARD_CLIP:
      return 'H';
    case BAM_CPAD:
      return 'P';
    case BAM_CEQUAL:
      return '=';
    case BAM_CDIFF:
      return 'X';
    default:
      exit(GT_EXIT_PROGRAMMING_ERROR);
  }
}

uint32_t gt_sam_alignment_flag(GtSamAlignment *sam_alignment)
{
  gt_assert(sam_alignment != NULL);
  return sam_alignment->s_alignment->core.flag;
}

bool gt_sam_alignment_is_paired(GtSamAlignment *sam_alignment)
{
  gt_assert(sam_alignment != NULL);
  if (sam_alignment->s_alignment->core.flag & BAM_FPAIRED)
    return true;
  else
    return false;
}

bool gt_sam_alignment_is_proper_paired(GtSamAlignment *sam_alignment)
{
  gt_assert(sam_alignment != NULL);
  if (sam_alignment->s_alignment->core.flag & BAM_FPROPER_PAIR)
    return true;
  else
    return false;
}

bool gt_sam_alignment_is_unmapped(GtSamAlignment *sam_alignment)
{
  gt_assert(sam_alignment != NULL);
  if (sam_alignment->s_alignment->core.flag & BAM_FUNMAP)
    return true;
  else
    return false;
}

bool gt_sam_alignment_mate_is_unmapped(GtSamAlignment *sam_alignment)
{
  gt_assert(sam_alignment != NULL);
  if (sam_alignment->s_alignment->core.flag & BAM_FMUNMAP)
    return true;
  else
    return false;
}

bool gt_sam_alignment_is_reverse(GtSamAlignment *sam_alignment)
{
  gt_assert(sam_alignment != NULL);
  if (sam_alignment->s_alignment->core.flag & BAM_FREVERSE)
    return true;
  else
    return false;
}

bool gt_sam_alignment_mate_is_reverse(GtSamAlignment *sam_alignment)
{
  gt_assert(sam_alignment != NULL);
  if (sam_alignment->s_alignment->core.flag & BAM_FMREVERSE)
    return true;
  else
    return false;
}

bool gt_sam_alignment_is_read1(GtSamAlignment *sam_alignment)
{
  gt_assert(sam_alignment != NULL);
  if (sam_alignment->s_alignment->core.flag & BAM_FREAD1)
    return true;
  else
    return false;
}

bool gt_sam_alignment_is_read2(GtSamAlignment *sam_alignment)
{
  gt_assert(sam_alignment != NULL);
  if (sam_alignment->s_alignment->core.flag & BAM_FREAD2)
    return true;
  else
    return false;
}

bool gt_sam_alignment_is_secondary(GtSamAlignment *sam_alignment)
{
  gt_assert(sam_alignment != NULL);
  if (sam_alignment->s_alignment->core.flag & BAM_FSECONDARY)
    return true;
  else
    return false;
}

bool gt_sam_alignment_has_qc_failure(GtSamAlignment *sam_alignment)
{
  gt_assert(sam_alignment != NULL);
  if (sam_alignment->s_alignment->core.flag & BAM_FQCFAIL)
    return true;
  else
    return false;
}

bool
gt_sam_alignment_is_optical_pcr_duplicate(GtSamAlignment *sam_alignment)
{
  gt_assert(sam_alignment != NULL);
  if (sam_alignment->s_alignment->core.flag & BAM_FDUP)
    return true;
  else
    return false;
}
