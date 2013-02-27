/*
  Copyright (c) 2010 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
  Copyright (c) 2010 Center for Bioinformatics, University of Hamburg

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

#include "core/assert_api.h"
#include "core/class_alloc_lock.h"
#include "core/error_api.h"
#include "core/seq_iterator_api.h"
#include "core/seq_iterator_sequence_buffer_api.h"
#include "core/unused_api.h"
#include "extended/fasta_header_iterator.h"
#include "extended/cstr_iterator.h"
#include "extended/cstr_iterator_rep.h"

struct GtFastaHeaderIterator {
  const GtCstrIterator parent_instance;
  GtSeqIterator *seq_iter;
  GtStrArray *filenametab;
};

static int gt_fasta_header_iterator_next(GtCstrIterator *cstr_iterator,
                                   const char **desc_str,
                                   GtError *err)
{
  GtFastaHeaderIterator *fhi = gt_fasta_header_iterator_cast(cstr_iterator);
  GT_UNUSED unsigned long len;
  return gt_seq_iterator_next(fhi->seq_iter, NULL, &len,
                             (char **) desc_str, err);
}

static int gt_fasta_header_iterator_reset(GtCstrIterator *cstr_iterator,
                                    GtError *err)
{
  GtFastaHeaderIterator *fhi = gt_fasta_header_iterator_cast(cstr_iterator);
  gt_error_check(err);
  gt_seq_iterator_delete(fhi->seq_iter);
  fhi->seq_iter = gt_seq_iterator_sequence_buffer_new(fhi->filenametab, err);
  if (fhi->seq_iter == NULL)
    return -1;
  gt_seq_iterator_set_sequence_output(fhi->seq_iter, false);
  return 0;
}

static void gt_fasta_header_iterator_delete(GtCstrIterator *cstr_iterator)
{
  GtFastaHeaderIterator *fhi = gt_fasta_header_iterator_cast(cstr_iterator);
  gt_seq_iterator_delete(fhi->seq_iter);
  gt_str_array_delete(fhi->filenametab);
}

/* map static local method to interface */
const GtCstrIteratorClass* gt_fasta_header_iterator_class(void)
{
  static const GtCstrIteratorClass *sic = NULL;
  gt_class_alloc_lock_enter();
  if (sic == NULL) {
    sic = gt_cstr_iterator_class_new(sizeof (GtFastaHeaderIterator),
                                     gt_fasta_header_iterator_next,
                                     gt_fasta_header_iterator_reset,
                                     gt_fasta_header_iterator_delete);
  }
  gt_class_alloc_lock_leave();
  return sic;
}

GtCstrIterator* gt_fasta_header_iterator_new(GtStrArray *filenametab,
                                               GtError *err)
{
  GtCstrIterator *cstr_iterator =
    gt_cstr_iterator_create(gt_fasta_header_iterator_class());
  GtFastaHeaderIterator *fhi = gt_fasta_header_iterator_cast(cstr_iterator);
  fhi->filenametab = gt_str_array_ref(filenametab);
  if (gt_fasta_header_iterator_reset(cstr_iterator, err) != 0)
    return NULL;
  return cstr_iterator;
}
