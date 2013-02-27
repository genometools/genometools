/*
  Copyright (c) 2001-2003 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2009-2010 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c)      2010 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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
#include "core/codon_api.h"
#include "core/codon_iterator_simple_api.h"
#include "core/ensure.h"
#include "core/ma.h"
#include "core/translator.h"

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
  tr = gt_calloc((size_t) 1, sizeof (GtTranslator));
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

void gt_translator_reset(GtTranslator *translator, GtCodonIterator *ci)
{
  gt_assert(translator && ci);
  translator->ci = ci;
}

void gt_translator_set_codon_iterator(GtTranslator *translator,
                                      GtCodonIterator *ci)
{
  gt_assert(translator && ci);
  translator->ci = ci;
}

void gt_translator_set_translation_table(GtTranslator *translator,
                                         GtTransTable *tt)
{
  gt_assert(translator && tt);
  if (translator->owntable) {
    gt_trans_table_delete(translator->transtable);
    translator->owntable = false;
  }
  translator->transtable = tt;
}

GtTranslatorStatus gt_translator_find_startcodon(GtTranslator *translator,
                                                 unsigned long *pos,
                                                 GtError *err)
{
  char n1, n2, n3;
  unsigned int frame;
  GtCodonIteratorStatus retval;
  gt_assert(translator && pos);
  gt_error_check(err);

  while (!(retval =
             gt_codon_iterator_next(translator->ci,
                                    &n1, &n2, &n3, &frame, err))) {
    if (gt_trans_table_is_start_codon(translator->transtable, n1, n2, n3)) {
      *pos = gt_codon_iterator_current_position(translator->ci)-1;
      return GT_TRANSLATOR_OK;
    }
  }
  if (retval == GT_CODON_ITERATOR_END)
    return GT_TRANSLATOR_END;
  else
    return GT_TRANSLATOR_ERROR;
}

GtTranslatorStatus gt_translator_find_stopcodon(GtTranslator *translator,
                                                unsigned long *pos,
                                                GtError *err)
{
  char n1, n2, n3;
  unsigned int frame;
  GtCodonIteratorStatus retval;
  gt_assert(translator && pos);
  gt_error_check(err);

  while (!(retval =
             gt_codon_iterator_next(translator->ci,
                                    &n1, &n2, &n3, &frame, err))) {
    if (gt_trans_table_is_stop_codon(translator->transtable, n1, n2, n3)) {
      *pos = gt_codon_iterator_current_position(translator->ci)-1;
      return GT_TRANSLATOR_OK;
    }
  }
  if (retval == GT_CODON_ITERATOR_END)
    return GT_TRANSLATOR_END;
  else
    return GT_TRANSLATOR_ERROR;
}

GtTranslatorStatus gt_translator_find_codon(GtTranslator *translator,
                                            GtStrArray *codons,
                                            unsigned long *pos,
                                            GtError *err)
{
  char n1, n2, n3;
  unsigned int frame;
  unsigned long i;
  GtCodonIteratorStatus retval;
  gt_assert(translator && codons && pos);
  gt_error_check(err);

  for (i = 0; i<gt_str_array_size(codons); i++) {
    int len;
    if ((len = (int) strlen(gt_str_array_get(codons, i))) != GT_CODON_LENGTH) {
      gt_error_set(err, "invalid codon length for codon %s: %d",
                   gt_str_array_get(codons, i), len);
      return GT_TRANSLATOR_ERROR;
    }
  }

  while (!(retval =
             gt_codon_iterator_next(translator->ci,
                                    &n1, &n2, &n3, &frame, err))) {
    for (i = 0; i<gt_str_array_size(codons); i++) {
      const char *codon;
      codon = gt_str_array_get(codons, i);
      if (n1 == codon[0] && n2 == codon[1] && n3 == codon[2]) {
        *pos = gt_codon_iterator_current_position(translator->ci)-1;
        return GT_TRANSLATOR_OK;
      }
    }
  }
  if (retval == GT_CODON_ITERATOR_END)
    return GT_TRANSLATOR_END;
  else
    return GT_TRANSLATOR_ERROR;
}

GtTranslatorStatus gt_translator_next(GtTranslator *translator,
                                      char *translated,
                                      unsigned int *frame,
                                      GtError *err)
{
  return gt_translator_next_with_start(translator, translated, frame, NULL,
                                       err);
}

GtTranslatorStatus gt_translator_next_with_start(GtTranslator *translator,
                                                 char *translated,
                                                 unsigned int *frame,
                                                 bool *start,
                                                 GtError *err)
{
  char n1, n2, n3;
  GtCodonIteratorStatus status;
  int retval;
  gt_assert(translator && translated && frame);
  gt_error_check(err);

  status = gt_codon_iterator_next(translator->ci, &n1, &n2, &n3, frame, err);
  if (status == GT_CODON_ITERATOR_END)
    return GT_TRANSLATOR_END;
  if (status == GT_CODON_ITERATOR_ERROR)
    return GT_TRANSLATOR_ERROR;
  retval = gt_trans_table_translate_codon(translator->transtable, n1, n2, n3,
                                          translated, err);
  if (retval)
    return GT_TRANSLATOR_ERROR;
  if (start != NULL)
    *start = gt_trans_table_is_start_codon(translator->transtable, n1, n2, n3);
  gt_assert(*frame < (unsigned int) GT_CODON_LENGTH);
  return GT_TRANSLATOR_OK;
}

void gt_translator_delete(GtTranslator *translator)
{
  if (!translator) return;
  if (translator->owntable)
    gt_trans_table_delete(translator->transtable);
  gt_free(translator);
}

int gt_translator_unit_test(GtError *err)
{
  int had_err = 0;
  GtTranslatorStatus test_errnum;
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
  ci = gt_codon_iterator_simple_new(seq, (unsigned long) strlen(seq), test_err);
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
    gt_ensure(had_err,
           test_errnum != GT_TRANSLATOR_ERROR && !gt_error_is_set(test_err));
  }
  gt_ensure(had_err,
         test_errnum == GT_TRANSLATOR_END && !gt_error_is_set(test_err));

  /* check 3-frame translation */
  gt_ensure(had_err, strcmp(gt_str_get(protein[0]),
                         "SFSF*LQRAICLCVD*KKSV**QLLNWLPAVSKLKFY*LR") == 0);
  gt_ensure(had_err, strcmp(gt_str_get(protein[1]),
                         "AFHSDCNGQYVSVWIKKRVSDSSF*TGYLP*VN*NFIDL") == 0);
  gt_ensure(had_err, strcmp(gt_str_get(protein[2]),
                         "LFILTATGNMSLCGLKKECLIAASELVTCRE*IKILLT*") == 0);

  /* find start codon -- positive */
  gt_error_unset(test_err);
  gt_codon_iterator_rewind(ci);
  test_errnum = gt_translator_find_startcodon(tr, &pos, test_err);
  gt_ensure(had_err, !test_errnum && !gt_error_is_set(test_err));
  gt_ensure(had_err, pos == 11UL);

  /* find stop codon -- positive */
  gt_error_unset(test_err);
  gt_codon_iterator_rewind(ci);
  test_errnum = gt_translator_find_stopcodon(tr, &pos, test_err);
  gt_ensure(had_err, !test_errnum && !gt_error_is_set(test_err));
  gt_ensure(had_err, pos == 12UL);

  /* find arbitrary codons -- positive */
  gt_error_unset(test_err);
  gt_codon_iterator_rewind(ci);
  test_errnum = gt_translator_find_codon(tr, codons, &pos, test_err);
  gt_ensure(had_err, !test_errnum && !gt_error_is_set(test_err));
  gt_ensure(had_err, pos == 14UL);

  /* find arbitrary codons -- negative (invalid codons) */
  gt_error_unset(test_err);
  gt_codon_iterator_rewind(ci);
  test_errnum = gt_translator_find_codon(tr, invalidcodons, &pos, test_err);
  gt_ensure(had_err,
         test_errnum == GT_TRANSLATOR_ERROR && gt_error_is_set(test_err));

  gt_error_unset(test_err);
  gt_codon_iterator_delete(ci);
  ci = gt_codon_iterator_simple_new(invalidseq,
                                    (unsigned long) strlen(invalidseq),
                                    test_err);
  gt_ensure(had_err, ci && !gt_error_is_set(test_err));
  gt_translator_reset(tr, ci);
  /* check translation of sequence with invalid beginning */
  test_errnum = gt_translator_next(tr, &translated, &frame, test_err);
  gt_ensure(had_err, test_errnum && gt_error_is_set(test_err));

  /* check translation of sequence with invalid character within */
  gt_error_unset(test_err);
  gt_codon_iterator_delete(ci);
  ci = gt_codon_iterator_simple_new(invalidseq2,
                                    (unsigned long) strlen(invalidseq2),
                                    test_err);
  gt_ensure(had_err, ci && !gt_error_is_set(test_err));
  gt_translator_reset(tr, ci);
  test_errnum = gt_translator_next(tr, &translated, &frame, test_err);
  while (!test_errnum && translated) {
    gt_str_append_char(protein[frame], translated);
    test_errnum = gt_translator_next(tr, &translated, &frame, test_err);
  }
  gt_ensure(had_err,
         test_errnum == GT_TRANSLATOR_ERROR && gt_error_is_set(test_err));

  /* find start codon -- fail */
  gt_error_unset(test_err);
  gt_codon_iterator_delete(ci);
  ci = gt_codon_iterator_simple_new(no_startcodon,
                                    (unsigned long) strlen(no_startcodon),
                                    test_err);
  gt_ensure(had_err, ci && !gt_error_is_set(test_err));
  gt_translator_reset(tr, ci);
  test_errnum = gt_translator_find_startcodon(tr, &pos, test_err);
  gt_ensure(had_err,
         test_errnum == GT_TRANSLATOR_END && !gt_error_is_set(test_err));

  /* find stop codon -- fail */
  gt_error_unset(test_err);
  gt_codon_iterator_rewind(ci);
  test_errnum = gt_translator_find_stopcodon(tr, &pos, test_err);
  gt_ensure(had_err,
         test_errnum == GT_TRANSLATOR_END && !gt_error_is_set(test_err));

  /* find arbitrary codons -- negative (none there) */
  gt_error_unset(test_err);
  gt_codon_iterator_rewind(ci);
  test_errnum = gt_translator_find_codon(tr, codons, &pos, test_err);
  gt_ensure(had_err,
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
