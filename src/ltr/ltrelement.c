/*
  Copyright (c) 2008 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2008 Center for Bioinformatics, University of Hamburg

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
#include "core/cstr.h"
#include "core/ensure.h"
#include "core/ma.h"
#include "core/range.h"
#include "core/unused_api.h"
#include "extended/reverse.h"
#include "ltr/ltrelement.h"
#include "match/alphadef.h"

unsigned long gt_ltrelement_length(GtLTRElement *e)
{
  gt_assert(e && (e->leftLTR_3 >= e->leftLTR_5));
  return e->rightLTR_3 - e->leftLTR_5 + 1;
}

unsigned long gt_ltrelement_leftltrlen(GtLTRElement *e)
{
  gt_assert(e && (e->leftLTR_3 >= e->leftLTR_5));
  return e->leftLTR_3-e->leftLTR_5 + 1;
}

char* gt_ltrelement_get_sequence(unsigned long start, unsigned long end,
                                         GtStrand strand, Encodedsequence *seq,
                                         Seqinfo *seqinfo, GtError *err)
{
  char *out;
  Uchar *symbolstring;
  const SfxAlphabet *alpha;
  Encodedsequencescanstate *ess;
  unsigned long len, i;

  gt_assert(seq && end >= start);
  gt_error_check(err);

  ess = newEncodedsequencescanstate();
  alpha = getencseqAlphabet(seq);
  len = end - start + 1;

  out          = gt_malloc((len + 1) * sizeof (char));
  symbolstring = gt_malloc((len + 1) * sizeof (Uchar));

  initEncodedsequencescanstate(ess, seq, Forwardmode,
                               seqinfo->seqstartpos + start);
  for (i=0;i<len;i++)
  {
    symbolstring[i] = sequentialgetencodedchar(seq, ess,
                                               seqinfo->seqstartpos + start + i,
                                               Forwardmode);
  }
  sprintfsymbolstring(out, alpha, symbolstring, len);
  gt_free(symbolstring);
  freeEncodedsequencescanstate(&ess);

  if (strand == GT_STRAND_REVERSE)
    (void) gt_reverse_complement(out, len, err);
  out[len]='\0';
  return out;
}

unsigned long gt_ltrelement_rightltrlen(GtLTRElement *e)
{
  gt_assert(e && (e->rightLTR_3 >= e->rightLTR_5));
  return e->rightLTR_3 - e->rightLTR_5 + 1;
}

void gt_ltrelement_offset2pos(GtLTRElement *e, GtRange *rng,
                              unsigned long radius,
                              enum GtOffset o,
                              GtStrand strand)
{
  unsigned long len = gt_range_length(rng);
  switch (strand)
  {
    case GT_STRAND_FORWARD:
      switch (o)
      {
        case GT_OFFSET_BEGIN_LEFT_LTR:
          rng->start = e->leftLTR_5 - radius + rng->start;
          break;
        case GT_OFFSET_END_RIGHT_LTR:
          rng->start = e->rightLTR_3 - radius + rng->start;
          break;
        case GT_OFFSET_END_LEFT_LTR:
          rng->start = e->leftLTR_3 - radius + rng->start;
          break;
        case GT_OFFSET_BEGIN_RIGHT_LTR:
          rng->start = e->rightLTR_5 - radius + rng->start;
          break;
      }
      break;
    case GT_STRAND_REVERSE:
      switch (o)
      {
        case GT_OFFSET_END_RIGHT_LTR:
          rng->start = e->leftLTR_5 + radius - rng->end - 1;
          break;
        case GT_OFFSET_BEGIN_LEFT_LTR:
          rng->start = e->rightLTR_3 + radius - rng->end;
          break;
        case GT_OFFSET_END_LEFT_LTR:
          rng->start = e->rightLTR_5 + radius - rng->end  - 1;
          break;
        case GT_OFFSET_BEGIN_RIGHT_LTR:
          rng->start = e->leftLTR_3 + radius - rng->end;
          break;
      }
      break;
    default:
      break;
  }
  rng->end = rng->start + len -1 ;
}

int gt_ltrelement_format_description(GtLTRElement *e, unsigned int seqnamelen,
                                     char *buf, size_t buflen)
{
  int ret;
  char *tmpstr;
  gt_assert(buf && e);
  tmpstr = gt_calloc(seqnamelen+1, sizeof (char));
  memset(tmpstr,0,sizeof (char) * (seqnamelen + 1));
  (void) snprintf(tmpstr, seqnamelen, "%s", e->seqid);
  gt_cstr_rep(tmpstr, ' ', '_');
  ret = snprintf(buf, buflen, "%s_%lu_%lu",
                 tmpstr, e->leftLTR_5, e->rightLTR_3);
  gt_free(tmpstr);
  return ret;
}

/*
int gt_ltrelement_unit_test(GtError *err)
{
  int had_err = 0;
  GtLTRElement element;
  GtRange rng1;
  GtSeq *seq;
  char *testseq = "ATCGAGGGGTCGAAT", *cseq;
  GtAlpha *alpha;
  unsigned long radius = 30;

  gt_error_check(err);

  alpha = gt_alpha_new_dna();
  seq = gt_seq_new(testseq, 15, alpha);

  element.leftLTR_5 = 100;
  element.leftLTR_3 = 150;
  element.rightLTR_5 = 450;
  element.rightLTR_3 = 600;

  rng1.start = 2;
  rng1.end = 28;

  gt_ltrelement_offset2pos(&element, &rng1, radius, GT_OFFSET_END_LEFT_LTR,
                           GT_STRAND_FORWARD);
  ensure(had_err, 122 == rng1.start);
  ensure(had_err, 148 == rng1.end);

  rng1.start = 2;
  rng1.end = 28;

  gt_ltrelement_offset2pos(&element, &rng1, radius, GT_OFFSET_BEGIN_RIGHT_LTR,
                           GT_STRAND_FORWARD);
  ensure(had_err, 422 == rng1.start);
  ensure(had_err, 448 == rng1.end);

  rng1.start = 2;
  rng1.end = 28;

  gt_ltrelement_offset2pos(&element, &rng1, radius, GT_OFFSET_END_LEFT_LTR,
                           GT_STRAND_REVERSE);
  ensure(had_err, 451 == rng1.start);
  ensure(had_err, 477 == rng1.end);

  cseq = gt_ltrelement_get_sequence(2, 6, GT_STRAND_FORWARD,
                                    gt_seq_get_orig(seq),
                                    gt_seq_length(seq), err);
  ensure(had_err,
         !strcmp("CGAGG\0", cseq));
  gt_free(cseq);

  cseq = gt_ltrelement_get_sequence(2, 6, GT_STRAND_REVERSE,
                                    gt_seq_get_orig(seq),
                                    gt_seq_length(seq), err);
  ensure(had_err,
         !strcmp("CCTCG\0", cseq));
  gt_free(cseq);

  gt_alpha_delete(alpha);
  gt_seq_delete(seq);

  return had_err;
}
*/
