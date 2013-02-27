/*
  Copyright (c) 2010-2011 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2010-2011 Center for Bioinformatics, University of Hamburg

  Permission to use, copy, modify, and distribute this software for any
  purpose with or without fee is hereby granted, provided that the above
  copyright notice and this permission notice appear in all copies.

  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
  WHAciOEVER RESULTING FROM LOSS OF USE, DATA OR PROFIci, WHETHER IN AN
  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/

#include "core/class_alloc_lock.h"
#include "core/encseq.h"
#include "core/ensure.h"
#include "core/ma.h"
#include "core/codon_api.h"
#include "core/codon_iterator_encseq_api.h"
#include "core/codon_iterator_rep.h"

struct GtCodonIteratorEncseq {
  const GtCodonIterator parent_instance;
  GtReadmode readmode;
  GtEncseq *encseq;
};

#define gt_codon_iterator_encseq_cast(CI)\
        gt_codon_iterator_cast(gt_codon_iterator_encseq_class(), CI)

static GtCodonIteratorStatus gt_codon_iterator_encseq_next(GtCodonIterator *ci,
                                                           char *n1,
                                                           char *n2,
                                                           char *n3,
                                                           unsigned int *frame,
                                                           GT_UNUSED GtError
                                                                           *err)
{
  GtCodonIteratorEncseq *cie;
  gt_assert(n1 && n2 && n3 && frame);
  gt_error_check(err);
  if (ci->pvt->curpos+2 >= ci->pvt->length)
    return GT_CODON_ITERATOR_END;
  cie = gt_codon_iterator_encseq_cast(ci);
  *n1 = gt_encseq_get_decoded_char(cie->encseq,
                                   ci->pvt->startpos+ci->pvt->curpos,
                                   cie->readmode);
  *n2 = gt_encseq_get_decoded_char(cie->encseq,
                                   ci->pvt->startpos+ci->pvt->curpos+1,
                                   cie->readmode);
  *n3 = gt_encseq_get_decoded_char(cie->encseq,
                                   ci->pvt->startpos+ci->pvt->curpos+2,
                                   cie->readmode);
  *frame = ci->pvt->curpos % GT_CODON_LENGTH;
  ci->pvt->curpos++;
  return GT_CODON_ITERATOR_OK;
}

static unsigned long gt_codon_iterator_encseq_current_position(GtCodonIterator
                                                                            *ci)
{
  gt_assert(ci);
  return ci->pvt->curpos;
}

static unsigned long gt_codon_iterator_encseq_length(GtCodonIterator *ci)
{
  gt_assert(ci);
  return ci->pvt->length;
}

static void gt_codon_iterator_encseq_rewind(GtCodonIterator *ci)
{
  ci->pvt->curpos = 0;
}

static void gt_codon_iterator_encseq_delete(GtCodonIterator *ci)
{
  GtCodonIteratorEncseq *cie;
  if (!ci) return;
  cie = gt_codon_iterator_encseq_cast(ci);
  gt_encseq_delete(cie->encseq);
}

const GtCodonIteratorClass* gt_codon_iterator_encseq_class(void)
{
  static const GtCodonIteratorClass *cic = NULL;
  gt_class_alloc_lock_enter();
  if (!cic)
  {
    cic = gt_codon_iterator_class_new(sizeof (GtCodonIteratorEncseq),
                                      gt_codon_iterator_encseq_delete,
                                      gt_codon_iterator_encseq_current_position,
                                      gt_codon_iterator_encseq_length,
                                      gt_codon_iterator_encseq_rewind,
                                      gt_codon_iterator_encseq_next);
  }
  gt_class_alloc_lock_leave();
  return cic;
}

GtCodonIterator* gt_codon_iterator_encseq_new_with_readmode(GtEncseq *encseq,
                                                         unsigned long startpos,
                                                         unsigned long length,
                                                         GtReadmode readmode,
                                                         GT_UNUSED GtError *err)
{
  GtCodonIteratorEncseq *cie;
  GtCodonIterator *ci;
  gt_assert(encseq && startpos + length - 1 < gt_encseq_total_length(encseq));
  gt_error_check(err);
  ci = gt_codon_iterator_create(gt_codon_iterator_encseq_class());
  cie = gt_codon_iterator_encseq_cast(ci);
  cie->encseq = gt_encseq_ref(encseq);
  cie->readmode = readmode;
  ci->pvt->length = length;
  ci->pvt->curpos = 0;
  ci->pvt->startpos = startpos;
  return ci;
}

GtCodonIterator* gt_codon_iterator_encseq_new(GtEncseq *encseq,
                                              unsigned long startpos,
                                              unsigned long length,
                                              GT_UNUSED GtError *err)
{
  return gt_codon_iterator_encseq_new_with_readmode(encseq, startpos, length,
                                                    GT_READMODE_FORWARD, err);
}

static int gt_codon_iterator_encseq_single_test(GtEncseq *encseq,
                                                const char *testseq,
                                                const char *testseq_cmp,
                                                GtReadmode readmode,
                                                GT_UNUSED GtError *err)
{
  unsigned long j, k, i;
  GtCodonIterator *ci;
  char n1, n2, n3;
  int had_err = 0;
  unsigned int frame;
  gt_error_check(err);

  for (j = 0; !had_err && j < strlen(testseq); j++) {
    for (k = j; !had_err && k < strlen(testseq); k++) {
      GtCodonIteratorStatus s;
      ci = gt_codon_iterator_encseq_new_with_readmode(encseq, j,
                                                      strlen(testseq) - k,
                                                      readmode,
                                                      NULL);
      i = j;
      while (!had_err && !(s = gt_codon_iterator_next(ci, &n1, &n2, &n3,
                                                      &frame, NULL))) {
        gt_ensure(had_err, n1 == testseq_cmp[i]);
        gt_ensure(had_err, n2 == testseq_cmp[i+1]);
        gt_ensure(had_err, n3 == testseq_cmp[i+2]);
        i++;
      }
      gt_codon_iterator_delete(ci);
    }
  }
  return had_err;
}

int gt_codon_iterator_encseq_unit_test(GtError *err)
{
  int had_err = 0,
      i, j;
  const char *testseq    = "gctgatcgactgaacatagctagcacggccgcgcgatcgtacgatg",
             *testseq_rc = "catcgtacgatcgcgcggccgtgctagctatgttcagtcgatcagc",
             *testseq_rv = "gtagcatgctagcgcgccggcacgatcgatacaagtcagctagtcg",
             *testseq_cm = "cgactagctgacttgtatcgatcgtgccggcgcgctagcatgctac";
  GtEncseq *encseq;
  GtEncseqBuilder *eb;
  GtCodonIterator *ci;
  GtAlphabet *alpha;
  char n1, n2, n3;
  unsigned int frame;
  gt_error_check(err);

  alpha = gt_alphabet_new_dna();
  eb = gt_encseq_builder_new(alpha);
  gt_encseq_builder_add_cstr(eb, testseq, strlen(testseq), "foo");
  encseq = gt_encseq_builder_build(eb, NULL);

  /* forward tests */
  had_err = gt_codon_iterator_encseq_single_test(encseq, testseq, testseq,
                                                 GT_READMODE_FORWARD, err);

  /* complement tests */
  had_err = gt_codon_iterator_encseq_single_test(encseq, testseq, testseq_cm,
                                                 GT_READMODE_COMPL, err);

  /* revcompl tests */
  had_err = gt_codon_iterator_encseq_single_test(encseq, testseq, testseq_rc,
                                                 GT_READMODE_REVCOMPL, err);

  /* reverse tests */
  had_err = gt_codon_iterator_encseq_single_test(encseq, testseq, testseq_rv,
                                                 GT_READMODE_REVERSE, err);

  /* lengths < 3 */
  for (j = 0; !had_err && j < 3; j++) {
    ci = gt_codon_iterator_encseq_new_with_readmode(encseq, 10, j,
                                                    GT_READMODE_REVCOMPL, NULL);
    i = 10;
    while (!(gt_codon_iterator_next(ci, &n1, &n2, &n3, &frame, NULL))) {
      gt_ensure(had_err, false);
    }
    gt_ensure(had_err, i == 10);
    gt_codon_iterator_delete(ci);
  }

  gt_encseq_delete(encseq);
  gt_encseq_builder_delete(eb);
  gt_alphabet_delete(alpha);
  return had_err;
}
