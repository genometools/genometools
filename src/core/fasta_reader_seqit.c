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
#include "core/seq_iterator_sequence_buffer_api.h"
#include "core/str_array.h"

struct GtFastaReaderSeqIt {
  const GtFastaReader parent_instance;
  GtStrArray *filenametab;
  GtSeqIterator *seqit;
};

#define gt_fasta_reader_seqit_cast(FR)\
        gt_fasta_reader_cast(gt_fasta_reader_seqit_class(), FR)

static int gt_fasta_reader_seqit_run(GtFastaReader *fasta_reader,
                                     GtFastaReaderProcDescription
                                     proc_description,
                                     GtFastaReaderProcSequencePart
                                     proc_sequence_part,
                                     GtFastaReaderProcSequenceLength
                                     proc_sequence_length,
                                     void *data, GtError *err)
{
  GtFastaReaderSeqIt *gt_fasta_reader_seqit =
    gt_fasta_reader_seqit_cast(fasta_reader);
  const GtUchar *sequence;
  unsigned long len;
  char *desc;
  int rval, had_err = 0;
  gt_error_check(err);

  /* at least one function has to be defined */
  gt_assert(proc_description || proc_sequence_part || proc_sequence_length);

  while ((rval = gt_seq_iterator_next(gt_fasta_reader_seqit->seqit, &sequence,
                                     &len, &desc, err))) {
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

static void gt_fasta_reader_seqit_free(GtFastaReader *fr)
{
  GtFastaReaderSeqIt *gt_fasta_reader_seqit = gt_fasta_reader_seqit_cast(fr);
  gt_str_array_delete(gt_fasta_reader_seqit->filenametab);
  gt_seq_iterator_delete(gt_fasta_reader_seqit->seqit);
}

const GtFastaReaderClass* gt_fasta_reader_seqit_class(void)
{
  static const GtFastaReaderClass frc = { sizeof (GtFastaReaderSeqIt),
                                        gt_fasta_reader_seqit_run,
                                        gt_fasta_reader_seqit_free };
  return &frc;
}

GtFastaReader* gt_fasta_reader_seqit_new(GtStr *sequence_filename)
{
  GtFastaReader *fr;
  GtFastaReaderSeqIt *gt_fasta_reader_seqit;
  gt_assert(sequence_filename);
  fr = gt_fasta_reader_create(gt_fasta_reader_seqit_class());
  gt_fasta_reader_seqit = gt_fasta_reader_seqit_cast(fr);
  gt_fasta_reader_seqit->filenametab = gt_str_array_new();
  gt_str_array_add_cstr(gt_fasta_reader_seqit->filenametab,
                       gt_str_get(sequence_filename));
  gt_fasta_reader_seqit->seqit = gt_seq_iterator_sequence_buffer_new(
                                             gt_fasta_reader_seqit->filenametab,
                                             NULL);
  return fr;
}
