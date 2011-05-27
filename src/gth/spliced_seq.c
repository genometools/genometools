/*
  Copyright (c) 2004-2009 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2004-2008 Center for Bioinformatics, University of Hamburg

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

#include "core/ma.h"
#include "core/range.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "gth/gthoutput.h"
#include "gth/spliced_seq.h"

static void fillsplicedseq(unsigned char *splicedseq,
                           const unsigned char *origseq, GtArray *ranges)
{
  const unsigned char *genomicptr;
  unsigned long i;
  gt_assert(ranges);
  for (i = 0; i < gt_array_size(ranges); i++) {
    for (genomicptr = origseq + ((GtRange*) gt_array_get(ranges, i))->start;
         genomicptr <= origseq +
                       ((GtRange*) gt_array_get(ranges, i))->end;
         *splicedseq++ = *genomicptr++);
  }
}

static void computepositionmapping(unsigned long *positionmapping,
                                   GtArray *ranges,
                                   GT_UNUSED unsigned long splicedseqlen)
{
  unsigned long i, rangecounter, templatepos = 0;
  for (rangecounter = 0; rangecounter < gt_array_size(ranges); rangecounter++) {
    for (i = ((GtRange*) gt_array_get(ranges, rangecounter))->start;
         i <= ((GtRange*) gt_array_get(ranges, rangecounter))->end;
         i++) {
      positionmapping[templatepos++] = i;
    }
  }
  gt_assert(templatepos == splicedseqlen);
}

GthSplicedSeq* gth_spliced_seq_new(const unsigned char *sequence,
                                   GtArray *ranges)
{
  GthSplicedSeq *spliced_seq = gt_malloc(sizeof *spliced_seq);

  gt_assert(sequence && ranges);

  spliced_seq->origseq = sequence;
  spliced_seq->ranges  = ranges;

  gt_assert(gt_ranges_are_consecutive(ranges));

  /* save total length of ranges */
  spliced_seq->splicedseqlen = gt_ranges_total_length(ranges);

  /* allocate space */
  spliced_seq->splicedseq = gt_malloc(sizeof (unsigned char) *
                                      spliced_seq->splicedseqlen);
  spliced_seq->positionmapping = gt_malloc(sizeof (unsigned long) *
                                           spliced_seq->splicedseqlen);

  /* processing */
  fillsplicedseq(spliced_seq->splicedseq, spliced_seq->origseq,
                 spliced_seq->ranges);

  /* compute position mapping */
  computepositionmapping(spliced_seq->positionmapping, spliced_seq->ranges,
                         spliced_seq->splicedseqlen);

  return spliced_seq;
}

GthSplicedSeq* gth_spliced_seq_new_with_comments(const unsigned char *sequence,
                                                 GtArray *ranges, bool comments,
                                                 GtFile *outfp)
{
  GthSplicedSeq *spliced_seq;
  gt_assert(sequence && ranges);
  spliced_seq = gth_spliced_seq_new(sequence, ranges);
  gt_assert(spliced_seq);
  if (comments) {
    double fraction = ((double) spliced_seq->splicedseqlen /
                       (double) gt_ranges_spanned_length(ranges)) * 100.0;
    /* fraction is valid percent value */
    gt_assert(fraction >= 0.0 && fraction <= 100.0);
    gt_file_xprintf(outfp,
                       "%c spliced sequence is %.2f%% of original sequence\n",
                       COMMENTCHAR, fraction);
  }
  return spliced_seq;
}

void gth_spliced_seq_delete(GthSplicedSeq *spliced_seq)
{
  if (!spliced_seq) return;
  gt_free(spliced_seq->splicedseq);
  gt_free(spliced_seq->positionmapping);
  gt_free(spliced_seq);
}

bool gth_spliced_seq_pos_is_border(const GthSplicedSeq *spliced_seq,
                                   unsigned long position)
{
  gt_assert(spliced_seq);
  /* position is legal */
  gt_assert(position < spliced_seq->splicedseqlen);
  if ((position + 1 < spliced_seq->splicedseqlen) &&
      (spliced_seq->positionmapping[position] + 1 !=
       spliced_seq->positionmapping[position+1])) {
    return true;
  }
  return false;
}

unsigned long gth_spliced_seq_border_length(const GthSplicedSeq *spliced_seq,
                                            unsigned long position)
{
  gt_assert(gth_spliced_seq_pos_is_border(spliced_seq, position));
  return spliced_seq->positionmapping[position+1] -
         spliced_seq->positionmapping[position] - 1;
}

unsigned long gth_spliced_seq_num_of_borders(const GthSplicedSeq *spliced_seq)
{
  unsigned long i, borders = 0;
  gt_assert(spliced_seq);
  for (i = 0; i < spliced_seq->splicedseqlen; i++) {
    if (gth_spliced_seq_pos_is_border(spliced_seq, i))
      borders++;
  }
  return borders;
}

static int cmpulong(const void *u1, const void *u2)
{
  return *(unsigned long*) u1 - *(unsigned long*) u2;
}

unsigned long gth_spliced_seq_orig_to_spliced_pos(const GthSplicedSeq
                                                  *spliced_seq,
                                                  unsigned long orig_pos)
{
  unsigned long *splicedposptr ;
  gt_assert(spliced_seq);
  splicedposptr = bsearch(&orig_pos, spliced_seq->positionmapping,
                          spliced_seq->splicedseqlen,
                          sizeof (unsigned long), cmpulong);
  if (splicedposptr)
    return splicedposptr - spliced_seq->positionmapping;
  return GT_UNDEF_ULONG;
}
