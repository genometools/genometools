/*
  Copyright (c) 2010 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2010 Center for Bioinformatics, University of Hamburg

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
#include "core/ensure.h"
#include "core/ma.h"
#include "core/codon_api.h"
#include "core/codon_iterator_simple_api.h"
#include "core/codon_iterator_rep.h"
#include "core/unused_api.h"

struct GtCodonIteratorSimple {
  const GtCodonIterator parent_instance;
  const char *dnaseq;
};

#define gt_codon_iterator_simple_cast(CI)\
        gt_codon_iterator_cast(gt_codon_iterator_simple_class(), CI)

static GtCodonIteratorStatus gt_codon_iterator_simple_next(GtCodonIterator *ci,
                                                           char *n1,
                                                           char *n2,
                                                           char *n3,
                                                           unsigned int *frame,
                                                           GT_UNUSED GtError
                                                                           *err)
{
  GtCodonIteratorSimple *cis;
  gt_assert(n1 && n2 && n3 && frame);
  gt_error_check(err);
  if (ci->pvt->curpos+2 >= ci->pvt->length)
    return GT_CODON_ITERATOR_END;
  cis = gt_codon_iterator_simple_cast(ci);
  *n1 = *(cis->dnaseq+ci->pvt->curpos);
  *n2 = *(cis->dnaseq+ci->pvt->curpos+1);
  *n3 = *(cis->dnaseq+ci->pvt->curpos+2);
  *frame = ci->pvt->curpos % GT_CODON_LENGTH;
  ci->pvt->curpos++;
  return GT_CODON_ITERATOR_OK;
}

static unsigned long gt_codon_iterator_simple_current_position(GtCodonIterator
                                                                            *ci)
{
  gt_assert(ci);
  return ci->pvt->curpos;
}

static unsigned long gt_codon_iterator_simple_length(GtCodonIterator *ci)
{
  gt_assert(ci);
  return ci->pvt->length;
}

static void gt_codon_iterator_simple_rewind(GtCodonIterator *ci)
{
  ci->pvt->curpos = ci->pvt->startpos;
}

const GtCodonIteratorClass* gt_codon_iterator_simple_class(void)
{
  static const GtCodonIteratorClass *cic = NULL;
  gt_class_alloc_lock_enter();
  if (!cic)
  {
    cic = gt_codon_iterator_class_new(sizeof (GtCodonIteratorSimple),
                                      NULL,
                                      gt_codon_iterator_simple_current_position,
                                      gt_codon_iterator_simple_length,
                                      gt_codon_iterator_simple_rewind,
                                      gt_codon_iterator_simple_next);
  }
  gt_class_alloc_lock_leave();
  return cic;
}

GtCodonIterator* gt_codon_iterator_simple_new(const char *seq,
                                              unsigned long length,
                                              GT_UNUSED GtError *err)
{
  GtCodonIteratorSimple *cis;
  GtCodonIterator *ci;
  gt_assert(seq && length >= GT_CODON_LENGTH);
  gt_error_check(err);
  ci = gt_codon_iterator_create(gt_codon_iterator_simple_class());
  cis = gt_codon_iterator_simple_cast(ci);
  cis->dnaseq = seq;
  ci->pvt->length = length;
  ci->pvt->curpos = 0;
  ci->pvt->startpos = 0;
  return ci;
}

int gt_codon_iterator_simple_unit_test(GtError *err)
{
  int had_err = 0,
      i;
  const char *testseq = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
  GtCodonIterator *ci;
  char n1, n2, n3;
  unsigned int frame;
  gt_error_check(err);

  ci = gt_codon_iterator_simple_new(testseq, 26, NULL);
  i = 0;
  while (!gt_codon_iterator_next(ci, &n1, &n2, &n3, &frame, NULL)) {
    gt_ensure(had_err, n1 == testseq[i]);
    gt_ensure(had_err, n2 == testseq[i+1]);
    gt_ensure(had_err, n3 == testseq[i+2]);
    i++;
  }
  gt_ensure(had_err, i == 24);
  gt_codon_iterator_delete(ci);
  return had_err;
}
