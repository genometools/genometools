/*
  Copyright (c) 2001-2003 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2009-2010 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2009-2010 Center for Bioinformatics, University of Hamburg

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
#include "core/assert_api.h"
#include "core/class_alloc.h"
#include "core/codon.h"
#include "core/codon_iterator_simple.h"
#include "core/ensure.h"
#include "core/ma.h"
#include "core/log.h"
#include "core/translator.h"
#include "core/unused_api.h"

struct GtTranslator {
  GtTransTable *transtable;
  bool owntable;
  GtCodonIterator *ci;
};

GtTranslator* gt_translator_new_with_table(GtTransTable *tt,
                                           GtCodonIterator *ci)
{
  GtTranslator *tr;
  gt_assert(tt && ci);
  tr = gt_calloc(1, sizeof (GtTranslator));
  tr->transtable = tt;
  tr->ci = ci;
  tr->owntable = false;
  return tr;
}

GtTranslator* gt_translator_new(GtCodonIterator *ci)
{
  GtTranslator *tr;
  gt_assert(ci);
  tr = gt_translator_new_with_table(gt_trans_table_new_standard(NULL), ci);
  tr->owntable = true;
  return tr;
}

void gt_translator_reset(GtTranslator *tr, GtCodonIterator *ci)
{
  gt_assert(tr && ci);
  tr->ci = ci;
}

void gt_translator_set_codon_iterator(GtTranslator *tr,
                                      GtCodonIterator *ci)
{
  gt_assert(tr && ci);
  tr->ci = ci;
}

void gt_translator_set_translation_table(GtTranslator *tr,
                                         GtTransTable *tt)
{
  gt_assert(tr && tt);
  if (tr->owntable) {
    gt_trans_table_delete(tr->transtable);
    tr->owntable = false;
  }
  tr->transtable = tt;
}

GtTranslatorStatus gt_translator_find_startcodon(GtTranslator *tr,
                                                 unsigned long *pos,
                                                 GtError *err)
{
  char n1, n2, n3;
  unsigned int frame;
  GtCodonIteratorStatus retval;
  gt_assert(tr && pos);
  gt_error_check(err);

  while (!(retval =
             gt_codon_iterator_next(tr->ci, &n1, &n2, &n3, &frame, err))) {
    if (gt_trans_table_is_start_codon(tr->transtable, n1, n2, n3)) {
      *pos = gt_codon_iterator_current_position(tr->ci)-1;
      return GT_TRANSLATOR_OK;
    }
  }
  if (retval == GT_CODON_ITERATOR_END)
    return GT_TRANSLATOR_END;
  else
    return GT_TRANSLATOR_ERROR;
}

GtTranslatorStatus gt_translator_find_stopcodon(GtTranslator *tr,
                                                unsigned long *pos,
                                                GtError *err)
{
  char n1, n2, n3;
  unsigned int frame;
  GtCodonIteratorStatus retval;
  gt_assert(tr && pos);
  gt_error_check(err);

  while (!(retval =
             gt_codon_iterator_next(tr->ci, &n1, &n2, &n3, &frame, err))) {
    if (gt_trans_table_is_stop_codon(tr->transtable, n1, n2, n3)) {
      *pos = gt_codon_iterator_current_position(tr->ci)-1;
      return GT_TRANSLATOR_OK;
    }
  }
  if (retval == GT_CODON_ITERATOR_END)
    return GT_TRANSLATOR_END;
  else
    return GT_TRANSLATOR_ERROR;
}

GtTranslatorStatus gt_translator_find_codon(GtTranslator *tr,
                                            GtStrArray *codons,
                                            unsigned long *pos,
                                            GtError *err)
{
  char n1, n2, n3;
  unsigned int frame;
  unsigned long i;
  GtCodonIteratorStatus retval;
  gt_assert(tr && codons && pos);
  gt_error_check(err);

  for (i = 0; i<gt_str_array_size(codons); i++) {
    int len;
    if ((len = strlen(gt_str_array_get(codons, i))) != GT_CODON_LENGTH) {
      gt_error_set(err, "invalid codon length for codon %s: %d",
                   gt_str_array_get(codons, i), len);
      return GT_TRANSLATOR_ERROR;
    }
  }

  while (!(retval =
             gt_codon_iterator_next(tr->ci, &n1, &n2, &n3, &frame, err))) {
    for (i = 0; i<gt_str_array_size(codons); i++) {
      const char *codon;
      codon = gt_str_array_get(codons, i);
      if (n1 == codon[0] && n2 == codon[1] && n3 == codon[2]) {
        *pos = gt_codon_iterator_current_position(tr->ci)-1;
        return GT_TRANSLATOR_OK;
      }
    }
  }
  if (retval == GT_CODON_ITERATOR_END)
    return GT_TRANSLATOR_END;
  else
    return GT_TRANSLATOR_ERROR;
}

GtTranslatorStatus gt_translator_next(GtTranslator *tr,
                                      char *translated,
                                      unsigned int *frame,
                                      GtError *err)
{
  char n1, n2, n3;
  int retval;
  gt_assert(tr && translated && frame);
  gt_error_check(err);

  retval = gt_codon_iterator_next(tr->ci, &n1, &n2, &n3, frame, err);
  if (retval == GT_CODON_ITERATOR_END)
    return GT_TRANSLATOR_END;
  if (retval == GT_CODON_ITERATOR_ERROR)
    return GT_TRANSLATOR_ERROR;
  retval = gt_trans_table_translate_codon(tr->transtable, n1, n2, n3,
                                          translated, err);
  if (retval)
    return GT_TRANSLATOR_ERROR;
  gt_assert(*frame < GT_CODON_LENGTH);
  return GT_TRANSLATOR_OK;
}

void gt_translator_delete(GtTranslator *tr)
{
  if (!tr) return;
  if (tr->owntable)
    gt_trans_table_delete(tr->transtable);
  gt_free(tr);
}

int gt_translator_unit_test(GtError *err)
{
  int had_err = 0, test_errnum = 0;
  GtTranslator *tr;
  GtCodonIterator *ci;
  GtError *test_err;
  GtStrArray *codons, *invalidcodons;
  const char *seq = "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGT"
                    "GGATTAAAAAAAGAGTGTCTGATAGCAGCTTCTGAACTGGT"
                    "TACCTGCCGTGAGTAAATTAAAATTTTATTGACTTAGG";
  const char *no_startcodon = "AAAAAAAAAATCATCTCCCCATTTTTTT";
  const char *invalidseq  = "ZAGCTTTTCATTCTGACTGCAAATATGTCTCTGTGT";
  const char *invalidseq2 = "AGCTTTTCATTCTGACZTGCAAATATGTCTCTGTGT";

  char translated;
  unsigned int frame;
  unsigned long pos = 0;
  GtStr *protein[3];
  gt_error_check(err);

  test_err = gt_error_new();
  ci = gt_codon_iterator_simple_new(seq, strlen(seq), test_err);
  tr = gt_translator_new(ci);
  protein[0] = gt_str_new();
  protein[1] = gt_str_new();
  protein[2] = gt_str_new();
  codons = gt_str_array_new();
  gt_str_array_add_cstr(codons, "ACG");
  gt_str_array_add_cstr(codons, "ACT");
  invalidcodons = gt_str_array_new();
  gt_str_array_add_cstr(invalidcodons, "ACG");
  gt_str_array_add_cstr(invalidcodons, "AC");

  /* do 3-frame translation */
  gt_error_unset(test_err);
  test_errnum = gt_translator_next(tr, &translated, &frame, test_err);
  while (!test_errnum && translated) {
    gt_str_append_char(protein[frame], translated);
    test_errnum = gt_translator_next(tr, &translated, &frame, test_err);
    ensure(had_err,
           test_errnum != GT_TRANSLATOR_ERROR && !gt_error_is_set(test_err));
  }
  ensure(had_err,
         test_errnum == GT_TRANSLATOR_END && !gt_error_is_set(test_err));

  /* check 3-frame translation */
  ensure(had_err, strcmp(gt_str_get(protein[0]),
                         "SFSF*LQRAICLCVD*KKSV**QLLNWLPAVSKLKFY*LR") == 0);
  ensure(had_err, strcmp(gt_str_get(protein[1]),
                         "AFHSDCNGQYVSVWIKKRVSDSSF*TGYLP*VN*NFIDL") == 0);
  ensure(had_err, strcmp(gt_str_get(protein[2]),
                         "LFILTATGNMSLCGLKKECLIAASELVTCRE*IKILLT*") == 0);

  /* find start codon -- positive */
  gt_error_unset(test_err);
  gt_codon_iterator_rewind(ci);
  test_errnum = gt_translator_find_startcodon(tr, &pos, test_err);
  ensure(had_err, !test_errnum && !gt_error_is_set(test_err));
  ensure(had_err, pos == 11);

  /* find stop codon -- positive */
  gt_error_unset(test_err);
  gt_codon_iterator_rewind(ci);
  test_errnum = gt_translator_find_stopcodon(tr, &pos, test_err);
  ensure(had_err, !test_errnum && !gt_error_is_set(test_err));
  ensure(had_err, pos == 12);

  /* find arbitrary codons -- positive */
  gt_error_unset(test_err);
  gt_codon_iterator_rewind(ci);
  test_errnum = gt_translator_find_codon(tr, codons, &pos, test_err);
  ensure(had_err, !test_errnum && !gt_error_is_set(test_err));
  ensure(had_err, pos == 14);

  /* find arbitrary codons -- negative (invalid codons) */
  gt_error_unset(test_err);
  gt_codon_iterator_rewind(ci);
  test_errnum = gt_translator_find_codon(tr, invalidcodons, &pos, test_err);
  ensure(had_err,
         test_errnum == GT_TRANSLATOR_ERROR && gt_error_is_set(test_err));

  gt_error_unset(test_err);
  gt_codon_iterator_delete(ci);
  ci = gt_codon_iterator_simple_new(invalidseq, strlen(invalidseq), test_err);
  ensure(had_err, ci && !gt_error_is_set(test_err));
  gt_translator_reset(tr, ci);
  /* check translation of sequence with invalid beginning */
  test_errnum = gt_translator_next(tr, &translated, &frame, test_err);
  ensure(had_err, test_errnum && gt_error_is_set(test_err));

  /* check translation of sequence with invalid character within */
  gt_error_unset(test_err);
  gt_codon_iterator_delete(ci);
  ci = gt_codon_iterator_simple_new(invalidseq2, strlen(invalidseq2), test_err);
  ensure(had_err, ci && !gt_error_is_set(test_err));
  gt_translator_reset(tr, ci);
  test_errnum = gt_translator_next(tr, &translated, &frame, test_err);
  while (!test_errnum && translated) {
    gt_str_append_char(protein[frame], translated);
    test_errnum = gt_translator_next(tr, &translated, &frame, test_err);
  }
  ensure(had_err,
         test_errnum == GT_TRANSLATOR_ERROR && gt_error_is_set(test_err));

  /* find start codon -- fail */
  gt_error_unset(test_err);
  gt_codon_iterator_delete(ci);
  ci = gt_codon_iterator_simple_new(no_startcodon, strlen(no_startcodon),
                                    test_err);
  ensure(had_err, ci && !gt_error_is_set(test_err));
  gt_translator_reset(tr, ci);
  test_errnum = gt_translator_find_startcodon(tr, &pos, test_err);
  ensure(had_err,
         test_errnum == GT_TRANSLATOR_END && !gt_error_is_set(test_err));

  /* find stop codon -- fail */
  gt_error_unset(test_err);
  gt_codon_iterator_rewind(ci);
  test_errnum = gt_translator_find_stopcodon(tr, &pos, test_err);
  ensure(had_err,
         test_errnum == GT_TRANSLATOR_END && !gt_error_is_set(test_err));

  /* find arbitrary codons -- negative (none there) */
  gt_error_unset(test_err);
  gt_codon_iterator_rewind(ci);
  test_errnum = gt_translator_find_codon(tr, codons, &pos, test_err);
  ensure(had_err,
         test_errnum == GT_TRANSLATOR_END && !gt_error_is_set(test_err));

  gt_codon_iterator_delete(ci);
  gt_translator_delete(tr);
  gt_str_delete(protein[0]);
  gt_str_delete(protein[1]);
  gt_str_delete(protein[2]);
  gt_str_array_delete(codons);
  gt_str_array_delete(invalidcodons);
  gt_error_delete(test_err);

  return had_err;
}
