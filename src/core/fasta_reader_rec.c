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
#include "core/fasta_reader_rec.h"
#include "core/fasta_reader_rep.h"
#include "core/fasta_separator.h"
#include "core/io.h"

struct FastaReaderRec {
  const FastaReader parent_instance;
  IO *seqio;
};

#define fasta_reader_rec_cast(FR)\
        fasta_reader_cast(fasta_reader_rec_class(), FR)

static int parse_fasta_description(GT_Str *description, IO *seqio, GT_Error *err)
{
  int rval;
  char cc;
  gt_error_check(err);
  assert(description && seqio);
  rval = io_get_char(seqio, &cc);
  assert(!rval); /* was checked earlier */
  /* make sure we got a proper fasta description */
  if (cc != FASTA_SEPARATOR) {
    gt_error_set(err, "the first character of fasta file \"%s\" has to be '%c'",
              io_get_filename(seqio), FASTA_SEPARATOR);

    return -1;
  }
  /* read description */
  while (!io_get_char(seqio, &cc) && cc != '\n')
    gt_str_append_char(description, cc);
  return 0;
}

static int parse_fasta_sequence(GT_Str *sequence, IO *seqio, GT_Error *err)
{
  char cc;
  gt_error_check(err);
  assert(sequence && seqio);
  assert(!gt_str_length(sequence));
  /* read sequence */
  while (!io_get_char(seqio, &cc) && cc != FASTA_SEPARATOR) {
    if (cc != '\n' && cc != ' ')
      gt_str_append_char(sequence, cc);
  }
  if (!gt_str_length(sequence)) {
    gt_error_set(err, "empty sequence given in line %lu",
              io_get_line_number(seqio));
    return -1;
  }
  if (cc == FASTA_SEPARATOR)
    io_unget_char(seqio, FASTA_SEPARATOR);
  return 0;
}

static int parse_fasta_entry(GT_Str *description, GT_Str *sequence, IO *seqio,
                             GT_Error *err)
{
  int had_err;
  gt_error_check(err);
  assert(description && sequence && seqio);
  had_err = parse_fasta_description(description, seqio, err);
  if (!had_err)
    had_err = parse_fasta_sequence(sequence, seqio, err);
  return had_err;
}

static int fasta_reader_rec_run(FastaReader *fasta_reader,
                                FastaReaderProcDescription proc_description,
                                FastaReaderProcSequencePart proc_sequence_part,
                                FastaReaderProcSequenceLength
                                proc_sequence_length, void *data, GT_Error *err)
{
  FastaReaderRec *fr = fasta_reader_rec_cast(fasta_reader);
  GT_Str *description, *sequence;
  int had_err = 0;
  gt_error_check(err);

  /* at least one function has to be defined */
  assert(proc_description || proc_sequence_part || proc_sequence_length);

  /* init */
  description = gt_str_new();
  sequence    = gt_str_new();

  /* make sure file is not empty */
  if (!io_has_char(fr->seqio)) {
    gt_error_set(err, "sequence file \"%s\" is empty", io_get_filename(fr->seqio));
    had_err = -1;
  }

  /* parse file */
  while (!had_err && io_has_char(fr->seqio)) {
    /* reset */
    gt_str_reset(description);
    gt_str_reset(sequence);

    /* parse entry */
    had_err = parse_fasta_entry(description, sequence, fr->seqio, err);

    /* process entry */
    if (!had_err && proc_description) {
      had_err = proc_description(gt_str_get(description), gt_str_length(description),
                                 data, err);
    }
    if (!had_err && proc_sequence_part) {
      had_err = proc_sequence_part(gt_str_get(sequence), gt_str_length(sequence),
                                   data, err);
    }
    if (!had_err && proc_sequence_length)
      had_err = proc_sequence_length(gt_str_length(sequence), data, err);
  }

  /* free */
  gt_str_delete(description);
  gt_str_delete(sequence);

  return had_err;
}

static void fasta_reader_rec_free(FastaReader *fr)
{
  FastaReaderRec *fasta_reader_rec = fasta_reader_rec_cast(fr);
  io_delete(fasta_reader_rec->seqio);
}

const FastaReaderClass* fasta_reader_rec_class(void)
{
  static const FastaReaderClass frc = { sizeof (FastaReaderRec),
                                        fasta_reader_rec_run,
                                        fasta_reader_rec_free };
  return &frc;
}

FastaReader* fasta_reader_rec_new(GT_Str *sequence_filename)
{
  FastaReader *fr = fasta_reader_create(fasta_reader_rec_class());
  FastaReaderRec *fasta_reader_rec = fasta_reader_rec_cast(fr);
  fasta_reader_rec->seqio = io_new(sequence_filename
                                   ? gt_str_get(sequence_filename) : NULL,
                                   "r");
  return fr;
}
