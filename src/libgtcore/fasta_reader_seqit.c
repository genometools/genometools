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
#include "libgtcore/fasta_reader_seqit.h"
#include "libgtcore/fasta_reader_rep.h"
#include "libgtcore/ma.h"
#include "libgtcore/seqiterator.h"
#include "libgtcore/strarray.h"

struct FastaReaderSeqIt {
  const FastaReader parent_instance;
  StrArray *filenametab;
  SeqIterator *seqit;
};

#define fasta_reader_seqit_cast(FR)\
        fasta_reader_cast(fasta_reader_seqit_class(), FR)

static int fasta_reader_seqit_run(FastaReader *fasta_reader,
                                  FastaReaderProcDescription proc_description,
                                  FastaReaderProcSequencePart
                                  proc_sequence_part,
                                  FastaReaderProcSequenceLength
                                  proc_sequence_length,
                                  void *data, Error *err)
{
  FastaReaderSeqIt *fasta_reader_seqit = fasta_reader_seqit_cast(fasta_reader);
  const Uchar *sequence;
  unsigned long len;
  char *desc;
  int rval, had_err = 0;
  error_check(err);

  /* at least one function has to be defined */
  assert(proc_description || proc_sequence_part || proc_sequence_length);

  while ((rval = seqiterator_next(fasta_reader_seqit->seqit, &sequence, &len,
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

    ma_free(desc);
    if (had_err)
      break;
  }

  return had_err;
}

static void fasta_reader_seqit_free(FastaReader *fr)
{
  FastaReaderSeqIt *fasta_reader_seqit = fasta_reader_seqit_cast(fr);
  strarray_delete(fasta_reader_seqit->filenametab);
  seqiterator_delete(fasta_reader_seqit->seqit);
}

const FastaReaderClass* fasta_reader_seqit_class(void)
{
  static const FastaReaderClass frc = { sizeof (FastaReaderSeqIt),
                                        fasta_reader_seqit_run,
                                        fasta_reader_seqit_free };
  return &frc;
}

FastaReader* fasta_reader_seqit_new(Str *sequence_filename)
{
  FastaReader *fr;
  FastaReaderSeqIt *fasta_reader_seqit;
  assert(sequence_filename);
  fr = fasta_reader_create(fasta_reader_seqit_class());
  fasta_reader_seqit = fasta_reader_seqit_cast(fr);
  fasta_reader_seqit->filenametab = strarray_new();
  strarray_add_cstr(fasta_reader_seqit->filenametab,
                    str_get(sequence_filename));
  fasta_reader_seqit->seqit = seqiterator_new(fasta_reader_seqit->filenametab,
                                              NULL, true);
  return fr;
}
