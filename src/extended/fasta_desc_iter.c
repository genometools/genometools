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
#include "core/error_api.h"
#include "core/seqiterator.h"
#include "core/seqiterator_sequence_buffer.h"
#include "core/unused_api.h"
#include "extended/fasta_desc_iter.h"
#include "extended/string_iter.h"
#include "extended/string_iter_rep.h"

struct GtFastaDescIter {
  const GtStringIter parent_instance;
  GtSeqIterator *seq_iter;
  GtStrArray *filenametab;
};

static int gt_fasta_desc_iter_next(GtStringIter *str_iter,
                                   const char **desc_str,
                                   GtError *err)
{
  GtFastaDescIter *fdi = gt_fasta_desc_iter_cast(str_iter);
  GT_UNUSED unsigned long len;
  return gt_seqiterator_next(fdi->seq_iter, NULL, &len, (char **) desc_str, err);
}

static int gt_fasta_desc_iter_reset(GtStringIter *str_iter,
                                    GtError *err)
{
  GtFastaDescIter *fdi = gt_fasta_desc_iter_cast(str_iter);
  gt_error_check(err);
  gt_seqiterator_delete(fdi->seq_iter);
  fdi->seq_iter = gt_seqiterator_sequence_buffer_new(fdi->filenametab, err);
  if (fdi->seq_iter == NULL)
    return -1;
  gt_seqiterator_set_sequence_output(fdi->seq_iter, false);
  return 0;
}

static void gt_fasta_desc_iter_delete(GtStringIter *str_iter)
{
  GtFastaDescIter *fdi = gt_fasta_desc_iter_cast(str_iter);
  gt_seqiterator_delete(fdi->seq_iter);
  gt_str_array_delete(fdi->filenametab);
}

/* map static local method to interface */
const GtStringIterClass* gt_fasta_desc_iter_class(void)
{
  static const GtStringIterClass *sic = NULL;
  if (sic == NULL) {
    sic = gt_string_iter_class_new(sizeof (GtFastaDescIter),
                                   gt_fasta_desc_iter_next,
                                   gt_fasta_desc_iter_reset,
                                   gt_fasta_desc_iter_delete);
  }
  return sic;
}

GtStringIter* gt_fasta_desc_iter_new(GtStrArray *filenametab,
                                     GtError *err)
{
  GtStringIter *str_iter = gt_string_iter_create(gt_fasta_desc_iter_class());
  GtFastaDescIter *fdi = gt_fasta_desc_iter_cast(str_iter);
  fdi->filenametab = gt_str_array_ref(filenametab);
  if (gt_fasta_desc_iter_reset(str_iter, err) != 0)
    return NULL;
  return str_iter;
}
