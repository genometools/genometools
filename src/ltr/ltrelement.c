/*
  Copyright (c) 2008-2010 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2008-2010 Center for Bioinformatics, University of Hamburg

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
#include "core/alphabet.h"
#include "core/cstr_api.h"
#include "core/ensure.h"
#include "core/ma.h"
#include "core/minmax.h"
#include "extended/reverse_api.h"
#include "ltr/ltrelement.h"

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
                                 GtStrand strand, GtEncseq *seq,
                                 unsigned long seqstartpos, GtError *err)
{
  char *out;
  int had_err = 0;
  GtUchar *symbolstring;
  const GtAlphabet *alpha;
  GtEncseqReader *esr;
  unsigned long len, i;

  gt_assert(seq && end >= start);
  gt_error_check(err);

  alpha = gt_encseq_alphabet(seq);
  len = end - start + 1;

  out          = gt_malloc((size_t) (len + 1) * sizeof (char));
  symbolstring = gt_malloc((size_t) (len + 1) * sizeof (GtUchar));

  esr = gt_encseq_create_reader_with_readmode(seq, GT_READMODE_FORWARD,
                                              seqstartpos + start);
  for (i=0;i<len;i++)
  {
    symbolstring[i] = gt_encseq_reader_next_encoded_char(esr);
  }
  gt_alphabet_decode_seq_to_cstr(alpha, out, symbolstring, len);
  gt_free(symbolstring);
  gt_encseq_reader_delete(esr);

  if (strand == GT_STRAND_REVERSE)
    had_err = gt_reverse_complement(out, len, err);
  out[len]='\0';
  if (had_err)
  {
    gt_free(out);
    out = NULL;
  }
  return out;
}

unsigned long gt_ltrelement_rightltrlen(GtLTRElement *e)
{
  gt_assert(e && (e->rightLTR_3 >= e->rightLTR_5));
  return e->rightLTR_3 - e->rightLTR_5 + 1;
}

int gt_ltrelement_format_description(GtLTRElement *e, unsigned int seqnamelen,
                                     char *buf, size_t buflen)
{
  int ret;
  char tmpstr[BUFSIZ];
  gt_assert(buf && e);

  (void) snprintf(tmpstr, MIN(BUFSIZ, (size_t) seqnamelen+1), "%s", e->seqid);
  tmpstr[seqnamelen+1] = '\0';
  gt_cstr_rep(tmpstr, ' ', '_');
  ret = snprintf(buf, buflen, "%s_%lu_%lu", tmpstr, e->leftLTR_5+1,
                 e->rightLTR_3+1);
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

  /* check left LTR */
  element.leftLTR_5 = 20UL;
  element.leftLTR_3 = 119UL;
  gt_ensure(had_err, gt_ltrelement_leftltrlen(&element) == 100UL);
  memset(tmp, 0, BUFSIZ);
  memcpy(tmp, fullseq + (element.leftLTR_5 * sizeof (char)),
         (size_t) (element.leftLTR_3 - element.leftLTR_5+ 1) * sizeof (char));
  gt_ensure(had_err,
            strcmp(tmp, "tatagcactgcatttcgaatatagtttcgaatatagcactgcatttcg"
                        "aatatagcactgcatttcgaatatagtttcgaatatagcactgcattt"
                        "cgaa" ) == 0);

  /* check right LTR */
  element.rightLTR_5 = 520UL;
  element.rightLTR_3 = 619UL;
  gt_ensure(had_err, gt_ltrelement_rightltrlen(&element) == 100UL);
  memset(tmp, 0, BUFSIZ);
  memcpy(tmp, fullseq + (element.rightLTR_5 * sizeof (char)),
         (size_t) (element.rightLTR_3 - element.rightLTR_5+ 1) * sizeof (char));
  gt_ensure(had_err,
            strcmp(tmp, "tatagcactgcatttcgaatatagtttcgaatatagcactgcatttcg"
                        "aatatagcactgcatttcgaatatagtttcgaatatagcactgcattt"
                        "cgaa" ) == 0);

  return had_err;
}
