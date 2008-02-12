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

#include "libgtcore/fasta_reader_rep.h"
#include "libgtcore/ma.h"
#include "libgtcore/unused.h"

FastaReader* fasta_reader_create(const FastaReaderClass *frc)
{
  FastaReader *fr;
  assert(frc && frc->size);
  fr = ma_calloc(1, frc->size);
  fr->c_class = frc;
  return fr;
}

void fasta_reader_delete(FastaReader *fr)
{
  if (!fr) return;
  assert(fr->c_class && fr->c_class->free);
  fr->c_class->free(fr);
  ma_free(fr);
}

int fasta_reader_run(FastaReader *fr,
                     FastaReaderProcDescription proc_description,
                     FastaReaderProcSequencePart proc_sequence_part,
                     FastaReaderProcSequenceLength proc_sequence_length,
                     void *data, Error *err)
{
  error_check(err);
  assert(fr && fr->c_class && fr->c_class->run);
  return fr->c_class->run(fr, proc_description, proc_sequence_part,
                          proc_sequence_length, data, err);
}

void* fasta_reader_cast(UNUSED const FastaReaderClass *frc, FastaReader *fr)
{
  assert(frc && fr && fr->c_class == frc);
  return fr;
}
