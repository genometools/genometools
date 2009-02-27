/*
  Copyright (c) 2008-2009 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2008-2009 Center for Bioinformatics, University of Hamburg

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
          rng->start = e->leftLTR_5 - (radius - 1) + rng->start;
          break;
        case GT_OFFSET_END_LEFT_LTR:
          rng->start = e->leftLTR_5 + (gt_ltrelement_leftltrlen(e) - 1)
                                    - (radius - 1) + rng->start;
          break;
        case GT_OFFSET_BEGIN_RIGHT_LTR:
          rng->start = e->leftLTR_5 + (gt_ltrelement_length(e) - 1)
                                    - (gt_ltrelement_rightltrlen(e) - 1)
                                    - (radius - 1) + rng->start;
          break;
        case GT_OFFSET_END_RIGHT_LTR:
          rng->start = e->leftLTR_5 + (gt_ltrelement_length(e) - 1)
                                    - (radius - 1) + rng->start;
          break;
      }
      break;
    case GT_STRAND_REVERSE:
      switch (o)
      {
        case GT_OFFSET_BEGIN_LEFT_LTR:
          rng->start = e->rightLTR_3 + (radius - 1) - rng->end;
          break;
        case GT_OFFSET_END_LEFT_LTR:
          rng->start = e->rightLTR_3 - (gt_ltrelement_rightltrlen(e) - 1)
                                    + (radius - 1) - rng->end;
          break;
        case GT_OFFSET_BEGIN_RIGHT_LTR:
          rng->start = e->rightLTR_3 - (gt_ltrelement_length(e) - 1)
                                     + (gt_ltrelement_leftltrlen(e) - 1)
                                     + (radius - 1) - rng->end;
          break;
        case GT_OFFSET_END_RIGHT_LTR:
          rng->start = e->rightLTR_3 - (gt_ltrelement_length(e) - 1)
                                     + (radius - 1) - rng->end;
          break;

      }
      break;
    default:
      break;
  }
  rng->end = rng->start + len - 1;
}

int gt_ltrelement_format_description(GtLTRElement *e, unsigned int seqnamelen,
                                     char *buf, size_t buflen)
{
  int ret;
  char *tmpstr;
  gt_assert(buf && e);
  tmpstr = gt_calloc(seqnamelen+1, sizeof (char));
  (void) snprintf(tmpstr, seqnamelen, "%s", e->seqid);
  gt_cstr_rep(tmpstr, ' ', '_');
  ret = snprintf(buf, buflen, "%s_%lu_%lu", tmpstr, e->leftLTR_5+1,
                 e->rightLTR_3+1);
  gt_free(tmpstr);
  return ret;
}

int gt_ltrelement_unit_test(GtError *err)
{
  int had_err = 0;
  GtLTRElement element;
  char tmp[BUFSIZ];
  const char *fullseq =                           "aaaaaaaaaaaaaaaaaaaa"
                    "tatagcactgcatttcgaatatagtttcgaatatagcactgcatttcgaa"
                    "tatagcactgcatttcgaatatagtttcgaatatagcactgcatttcgaa"
                    "acatactaggatgctagaatatagtttcgaatatagcactgcatttcgaa"
                    "tatagcactgcatttcgaatatagtttcgaatatagcactgcatttcgaa"
                    "tatagcactgcatttcgaatatagtttcgaatatagcactgcatttcgaa"
                    "tatagcactgcatttcgaatatagtttcgaatatagcactgcatttcgaa"
                    "tatagcactgcatttcgaatatagtttcgaatatagcactgcatttcgaa"
                    "tatagcactgcatttcgaatatagtttcgaatatagcactgcatttcgaa"
                    "tatagcactgcatttcgaatatagtttcgaatatagcactgcatttcgaa"
                    "tatagcactgcatttcgaatatagtttcgaatataggatcctaaggctac"
                    "tatagcactgcatttcgaatatagtttcgaatatagcactgcatttcgaa"
                    "tatagcactgcatttcgaatatagtttcgaatatagcactgcatttcgaa"
                    "aaaaaaaaaaaaaaaaaaaa";

  gt_error_check(err);

  element.leftLTR_5 = 20;
  element.leftLTR_3 = 119;
  ensure(had_err, gt_ltrelement_leftltrlen(&element) == 100);
  memset(tmp, 0, BUFSIZ);
  memcpy(tmp, fullseq + (element.leftLTR_5 * sizeof (char)),
         (element.leftLTR_3 - element.leftLTR_5+ 1) * sizeof (char));
  ensure(had_err, strcmp(tmp, "tatagcactgcatttcgaatatagtttcgaatatagcactgcatttcg"
                              "aatatagcactgcatttcgaatatagtttcgaatatagcactgcattt"
                              "cgaa" ) == 0);

  /* check right LTR */
  element.rightLTR_5 = 520;
  element.rightLTR_3 = 619;
  ensure(had_err, gt_ltrelement_rightltrlen(&element) == 100);
  memset(tmp, 0, BUFSIZ);
  memcpy(tmp, fullseq + (element.rightLTR_5 * sizeof (char)),
         (element.rightLTR_3 - element.rightLTR_5+ 1) * sizeof (char));
  ensure(had_err, strcmp(tmp, "tatagcactgcatttcgaatatagtttcgaatatagcactgcatttcg"
                              "aatatagcactgcatttcgaatatagtttcgaatatagcactgcattt"
                              "cgaa" ) == 0);

  return had_err;
}
