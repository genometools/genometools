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

#include "core/fasta_reader_rep.h"
#include "core/ma.h"
#include "core/unused_api.h"

GtFastaReader* gt_fasta_reader_create(const GtFastaReaderClass *frc)
{
  GtFastaReader *fr;
  gt_assert(frc && frc->size);
  fr = gt_calloc(1, frc->size);
  fr->c_class = frc;
  return fr;
}

void gt_fasta_reader_delete(GtFastaReader *fr)
{
  if (!fr) return;
  gt_assert(fr->c_class && fr->c_class->free);
  fr->c_class->free(fr);
  gt_free(fr);
}

int gt_fasta_reader_run(GtFastaReader *fr,
                     GtFastaReaderProcDescription proc_description,
                     GtFastaReaderProcSequencePart proc_sequence_part,
                     GtFastaReaderProcSequenceLength proc_sequence_length,
                     void *data, GtError *err)
{
  gt_error_check(err);
  gt_assert(fr && fr->c_class && fr->c_class->run);
  return fr->c_class->run(fr, proc_description, proc_sequence_part,
                          proc_sequence_length, data, err);
}

void* gt_fasta_reader_cast(GT_UNUSED const GtFastaReaderClass *frc,
                           GtFastaReader *fr)
{
  gt_assert(frc && fr && fr->c_class == frc);
  return fr;
}
