/*
  Copyright (c) 2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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
#include "core/fasta_reader_seqit.h"
#include "core/fasta_reader_rep.h"
#include "core/ma.h"
#include "core/seqiterator.h"
#include "core/strarray.h"

struct GT_FastaReaderSeqIt {
  const GT_FastaReader parent_instance;
  GT_StrArray *filenametab;
  SeqIterator *seqit;
};

#define gt_fasta_reader_seqit_cast(FR)\
        gt_fasta_reader_cast(gt_fasta_reader_seqit_class(), FR)

static int gt_fasta_reader_seqit_run(GT_FastaReader *fasta_reader,
                                     GT_FastaReaderProcDescription
                                     proc_description,
                                     GT_FastaReaderProcSequencePart
                                     proc_sequence_part,
                                     GT_FastaReaderProcSequenceLength
                                     proc_sequence_length,
                                     void *data, GT_Error *err)
{
  GT_FastaReaderSeqIt *gt_fasta_reader_seqit =
    gt_fasta_reader_seqit_cast(fasta_reader);
  const Uchar *sequence;
  unsigned long len;
  char *desc;
  int rval, had_err = 0;
  gt_error_check(err);

  /* at least one function has to be defined */
  assert(proc_description || proc_sequence_part || proc_sequence_length);

  while ((rval = seqiterator_next(gt_fasta_reader_seqit->seqit, &sequence, &len,
                                  &desc, err))) {
    if (rval < 0) {
      had_err = -1;
      break;
    }

    if (proc_description)
      had_err = proc_description(desc, strlen(desc), data, err);
    if (!had_err && proc_sequence_part)
      had_err = proc_sequence_part((char*) sequence, len, data, err);
    if (!had_err && proc_sequence_length)
      had_err = proc_sequence_length(len, data, err);

    gt_free(desc);
    if (had_err)
      break;
  }

  return had_err;
}

static void gt_fasta_reader_seqit_free(GT_FastaReader *fr)
{
  GT_FastaReaderSeqIt *gt_fasta_reader_seqit = gt_fasta_reader_seqit_cast(fr);
  gt_strarray_delete(gt_fasta_reader_seqit->filenametab);
  seqiterator_delete(gt_fasta_reader_seqit->seqit);
}

const GT_FastaReaderClass* gt_fasta_reader_seqit_class(void)
{
  static const GT_FastaReaderClass frc = { sizeof (GT_FastaReaderSeqIt),
                                        gt_fasta_reader_seqit_run,
                                        gt_fasta_reader_seqit_free };
  return &frc;
}

GT_FastaReader* gt_fasta_reader_seqit_new(GT_Str *sequence_filename)
{
  GT_FastaReader *fr;
  GT_FastaReaderSeqIt *gt_fasta_reader_seqit;
  assert(sequence_filename);
  fr = gt_fasta_reader_create(gt_fasta_reader_seqit_class());
  gt_fasta_reader_seqit = gt_fasta_reader_seqit_cast(fr);
  gt_fasta_reader_seqit->filenametab = gt_strarray_new();
  gt_strarray_add_cstr(gt_fasta_reader_seqit->filenametab,
                       gt_str_get(sequence_filename));
  gt_fasta_reader_seqit->seqit = seqiterator_new(gt_fasta_reader_seqit
                                                 ->filenametab, NULL, true);
  return fr;
}
