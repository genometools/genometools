/*
  Copyright (c) 2011 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
  Copyright (c) 2011 Center for Bioinformatics, University of Hamburg

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

#include "core/alphabet_api.h"
#include "core/ensure.h"
#include "core/ma_api.h"
#include "core/types_api.h"
#include "core/undef_api.h"
#include "extended/sam_alignment_rep.h"
#include "extended/sam_alignment.h"
#include <sam.h>

#define BAMBASEA 1
#define BAMBASEC 2
#define BAMBASEG 4
#define BAMBASET 8
#define BAMBASEN 15

#define PHREDOFFSET 33

GtSamAlignment *gt_sam_alignment_new(GtAlphabet *alphabet)
{
  GtSamAlignment *gt_s_alignment;
  gt_s_alignment = gt_malloc(sizeof (GtSamAlignment));
  gt_s_alignment->s_alignment = bam_init1();
  gt_s_alignment->alphabet = gt_alphabet_ref(alphabet);
  gt_s_alignment->seq_buffer = NULL;
  gt_s_alignment->qual_buffer = NULL;
  gt_s_alignment->buffsize = 0;
  return gt_s_alignment;
}

void gt_sam_alignment_delete(GtSamAlignment *gt_s_alignment)
{
  bam_destroy1(gt_s_alignment->s_alignment);
  gt_alphabet_delete(gt_s_alignment->alphabet);
  if (gt_s_alignment->buffsize > 0) {
    gt_free(gt_s_alignment->qual_buffer);
    gt_free(gt_s_alignment->seq_buffer);
  }
  gt_free(gt_s_alignment);
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
      return gt_alphabet_encode(alphabet, gt_alphabet_wildcard_show(alphabet));
  }
  return GT_UNDEF_UCHAR;
}

unsigned long gt_sam_alignment_pos(GtSamAlignment *gt_s_alignment)
{
  return (unsigned long) gt_s_alignment->s_alignment->core.pos;
}

unsigned long gt_sam_alignment_read_length(GtSamAlignment *gt_s_alignment)
{
  return (unsigned long) gt_s_alignment->s_alignment->core.l_qseq;
}

const char *gt_sam_alignment_identifier(GtSamAlignment *gt_s_alignment)
{
  return bam1_qname(gt_s_alignment->s_alignment);
}

int32_t gt_sam_alignment_ref_num(GtSamAlignment *gt_s_alignment)
{
  return gt_s_alignment->s_alignment->core.tid;
}

const GtUchar *gt_sam_alignment_sequence(GtSamAlignment *gt_s_alignment)
{
  unsigned long query_len, idx;
  uint8_t *bam_seq;

  query_len = (unsigned long) gt_s_alignment->s_alignment->core.l_qseq;
  if (gt_s_alignment->buffsize < query_len + 1) {
    gt_s_alignment->seq_buffer = gt_realloc(gt_s_alignment->seq_buffer,
                                   sizeof (gt_s_alignment->seq_buffer) *
                                   (query_len + 1));
    gt_s_alignment->qual_buffer = gt_realloc(gt_s_alignment->qual_buffer,
                                    sizeof (gt_s_alignment->qual_buffer) *
                                    (query_len + 1));
    gt_s_alignment->buffsize = query_len + 1;
    gt_assert(gt_s_alignment->qual_buffer);
    gt_assert(gt_s_alignment->seq_buffer);
  }

  bam_seq = bam1_seq(gt_s_alignment->s_alignment);

  for (idx = 0UL; idx < query_len; idx++) {
    gt_s_alignment->seq_buffer[idx] = bambase2gtbase(bam1_seqi(bam_seq, idx),
                                      gt_s_alignment->alphabet);
  }
  gt_s_alignment->seq_buffer[query_len] = '\0';

  return gt_s_alignment->seq_buffer;
}

const GtUchar *gt_sam_alignment_qualitystring(GtSamAlignment *gt_s_alignment)
{
  unsigned long query_len, idx;
  uint8_t *qual;

  query_len = (unsigned long) gt_s_alignment->s_alignment->core.l_qseq;
  if (gt_s_alignment->buffsize < query_len) {
    gt_s_alignment->seq_buffer = gt_realloc(gt_s_alignment->seq_buffer,
                                   sizeof (gt_s_alignment->seq_buffer) *
                                   (query_len + 1));
    gt_s_alignment->qual_buffer = gt_realloc(gt_s_alignment->qual_buffer,
                                    sizeof (gt_s_alignment->qual_buffer) *
                                    (query_len + 1));
    gt_s_alignment->buffsize = query_len;
    gt_assert(gt_s_alignment->qual_buffer);
    gt_assert(gt_s_alignment->seq_buffer);
  }
  qual = bam1_qual(gt_s_alignment->s_alignment);
  for (idx = 0UL; idx < query_len; idx++) {
    gt_s_alignment->qual_buffer[idx] = (char) qual[idx] + PHREDOFFSET;
  }
  gt_s_alignment->qual_buffer[query_len] = '\0';

  return gt_s_alignment->qual_buffer;
}

uint16_t gt_sam_alignment_cigar_length(GtSamAlignment *gt_s_alignment)
{
  return gt_s_alignment->s_alignment->core.n_cigar;
}

uint32_t gt_sam_alignment_cigar_i_length(GtSamAlignment *gt_s_alignment,
                                         uint16_t i)
{
  return bam1_cigar(gt_s_alignment->s_alignment)[i]>>BAM_CIGAR_SHIFT;
}

unsigned char gt_sam_alignment_cigar_i_operation(GtSamAlignment *gt_s_alignment,
                                                 uint16_t i)
{
  switch ((unsigned char) bam1_cigar(gt_s_alignment->s_alignment)[i]
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
      return 'M';
    case BAM_CDIFF:
      return 'X';
    default:
      exit(GT_EXIT_PROGRAMMING_ERROR);
  }
}

uint32_t gt_sam_alignment_flag(GtSamAlignment *gt_s_alignment)
{
  return gt_s_alignment->s_alignment->core.flag;
}

bool gt_sam_alignment_is_paired(GtSamAlignment *gt_s_alignment)
{
  if (gt_s_alignment->s_alignment->core.flag & BAM_FPAIRED)
    return true;
  else
    return false;
}

bool gt_sam_alignment_is_proper_paired(GtSamAlignment *gt_s_alignment)
{
  if (gt_s_alignment->s_alignment->core.flag & BAM_FPROPER_PAIR)
    return true;
  else
    return false;
}

bool gt_sam_alignment_is_unmapped(GtSamAlignment *gt_s_alignment)
{
  if (gt_s_alignment->s_alignment->core.flag & BAM_FUNMAP)
    return true;
  else
    return false;
}

bool gt_sam_alignment_mate_is_unmapped(GtSamAlignment *gt_s_alignment)
{
  if (gt_s_alignment->s_alignment->core.flag & BAM_FMUNMAP)
    return true;
  else
    return false;
}

bool gt_sam_alignment_is_reverse(GtSamAlignment *gt_s_alignment)
{
  if (gt_s_alignment->s_alignment->core.flag & BAM_FREVERSE)
    return true;
  else
    return false;
}

bool gt_sam_alignment_mate_is_reverse(GtSamAlignment *gt_s_alignment)
{
  if (gt_s_alignment->s_alignment->core.flag & BAM_FMREVERSE)
    return true;
  else
    return false;
}

bool gt_sam_alignment_is_read1(GtSamAlignment *gt_s_alignment)
{
  if (gt_s_alignment->s_alignment->core.flag & BAM_FREAD1)
    return true;
  else
    return false;
}

bool gt_sam_alignment_is_read2(GtSamAlignment *gt_s_alignment)
{
  if (gt_s_alignment->s_alignment->core.flag & BAM_FREAD2)
    return true;
  else
    return false;
}

bool gt_sam_alignment_is_secondary(GtSamAlignment *gt_s_alignment)
{
  if (gt_s_alignment->s_alignment->core.flag & BAM_FSECONDARY)
    return true;
  else
    return false;
}

bool gt_sam_alignment_has_qc_failure(GtSamAlignment *gt_s_alignment)
{
  if (gt_s_alignment->s_alignment->core.flag & BAM_FQCFAIL)
    return true;
  else
    return false;
}

bool
gt_sam_alignment_is_optical_pcr_duplicate(GtSamAlignment *gt_s_alignment)
{
  if (gt_s_alignment->s_alignment->core.flag & BAM_FDUP)
    return true;
  else
    return false;
}
