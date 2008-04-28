/*
  Copyright (c) 2008 Sascha Steinbiss <ssteinbiss@zbh.uni-hamburg.de>
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

#include <assert.h>
#include <string.h>
#include "libgtcore/cstr.h"
#include "libgtcore/ensure.h"
#include "libgtcore/ma.h"
#include "libgtcore/range.h"
#include "libgtext/reverse.h"
#include "ltrelement.h"

unsigned long ltrelement_length(LTRElement *e)
{
  assert(e && (e->leftLTR_3 >= e->leftLTR_5));
  return e->rightLTR_3 - e->leftLTR_5;
}

unsigned long ltrelement_leftltrlen(LTRElement *e)
{
  assert(e && (e->leftLTR_3 >= e->leftLTR_5));
  return e->leftLTR_3-e->leftLTR_5+1;
}

char* ltrelement_get_sequence(unsigned long start, unsigned long end,
                              Strand strand, Seq *seq, Error *err)
{
  char *out;
  unsigned long len;

  assert(seq && end >= start && end <= seq_length(seq));

  len = end - start + 1;
  out = ma_malloc(sizeof (char) * (len + 1));
  memcpy(out, seq_get_orig(seq)+start, sizeof (char) * len);
  if (strand == STRAND_REVERSE)
    reverse_complement(out, len, err);
  out[len]='\0';
  return out;
}

unsigned long ltrelement_rightltrlen(LTRElement *e)
{
  assert(e && (e->rightLTR_3 >= e->rightLTR_5));
  return e->rightLTR_3 - e->rightLTR_5 + 1;
}

void ltrelement_offset2pos(LTRElement *e, Range *rng,
                               unsigned long radius,
                               enum Offset o,
                               Strand strand)
{
  unsigned long len = range_length(*rng);
  switch (strand)
  {
    case STRAND_FORWARD:
      switch (o)
      {
        case OFFSET_BEGIN_LEFT_LTR:
          rng->start = e->leftLTR_5 - radius + rng->start;
          break;
        case OFFSET_END_RIGHT_LTR:
          rng->start = e->rightLTR_3 - radius + rng->start;
          break;
        case OFFSET_END_LEFT_LTR:
          rng->start = e->leftLTR_3 - radius + rng->start;
          break;
        case OFFSET_BEGIN_RIGHT_LTR:
          rng->start = e->rightLTR_5 - radius + rng->start;
          break;
      }
      break;
    case STRAND_REVERSE:
      switch (o)
      {
        case OFFSET_END_RIGHT_LTR:
          rng->start = e->leftLTR_5 + radius - rng->end - 1;
          break;
        case OFFSET_BEGIN_LEFT_LTR:
          rng->start = e->rightLTR_3 + radius - rng->end;
          break;
        case OFFSET_END_LEFT_LTR:
          rng->start = e->rightLTR_5 + radius - rng->end  - 1;
          break;
        case OFFSET_BEGIN_RIGHT_LTR:
          rng->start = e->leftLTR_3 + radius - rng->end;
          break;
      }
      break;
    default:
      break;
  }
  rng->end = rng->start + len -1 ;
}

int ltrelement_format_description(LTRElement *e, unsigned int seqnamelen,
                                  char *buf, size_t buflen)
{
  int ret;
  char *tmpstr;
  assert(buf && e);
  tmpstr = ma_malloc(sizeof (char) * (seqnamelen + 1));
  memset(tmpstr,0,sizeof (char) * (seqnamelen + 1));
  snprintf(tmpstr, seqnamelen, "%s", e->seqid);
  cstr_rep(tmpstr, ' ', '_');
  ret = snprintf(buf, buflen, "%s_%lu_%lu",
                  tmpstr, e->leftLTR_5, e->rightLTR_3);
  ma_free(tmpstr);
  return ret;
}

int ltrelement_unit_test(Error *err)
{
  int had_err = 0;
  error_check(err);

  LTRElement element;
  Range rng1;
  Seq *seq;
  char *testseq = "ATCGAGGGGTCGAAT", *cseq;
  Alpha *alpha;
  unsigned long radius = 30;

  alpha = alpha_new_dna();
  seq = seq_new(testseq, 15, alpha);

  element.leftLTR_5 = 100;
  element.leftLTR_3 = 150;
  element.rightLTR_5 = 450;
  element.rightLTR_3 = 600;

  rng1.start = 2;
  rng1.end = 28;

  ltrelement_offset2pos(&element, &rng1, radius, OFFSET_END_LEFT_LTR,
                        STRAND_FORWARD);
  ensure(had_err, 122 == rng1.start);
  ensure(had_err, 148 == rng1.end);

  rng1.start = 2;
  rng1.end = 28;

  ltrelement_offset2pos(&element, &rng1, radius, OFFSET_BEGIN_RIGHT_LTR,
                        STRAND_FORWARD);
  ensure(had_err, 422 == rng1.start);
  ensure(had_err, 448 == rng1.end);

  rng1.start = 2;
  rng1.end = 28;

  ltrelement_offset2pos(&element, &rng1, radius, OFFSET_END_LEFT_LTR,
                        STRAND_REVERSE);
  ensure(had_err, 451 == rng1.start);
  ensure(had_err, 477 == rng1.end);

  cseq = ltrelement_get_sequence(2, 6, STRAND_FORWARD, seq, err);
  ensure(had_err,
         !strcmp("CGAGG\0", cseq));
  ma_free(cseq);

  cseq = ltrelement_get_sequence(2, 6, STRAND_REVERSE, seq, err);
  ensure(had_err,
         !strcmp("CCTCG\0", cseq));
  ma_free(cseq);

  alpha_delete(alpha);
  seq_delete(seq);

  return had_err;
}
