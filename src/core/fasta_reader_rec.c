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
#include "core/unused_api.h"

struct GtFastaReaderRec {
  const GtFastaReader parent_instance;
  GtIO *seqio;
};

#define gt_fasta_reader_rec_cast(FR)\
        gt_fasta_reader_cast(gt_fasta_reader_rec_class(), FR)

static int parse_fasta_description(GtStr *description, GtIO *seqio,
                                   GtError *err)
{
  GT_UNUSED int rval;
  char cc;
  gt_error_check(err);
  gt_assert(description && seqio);
  rval = gt_io_get_char(seqio, &cc);
  gt_assert(!rval); /* was checked earlier */
  /* make sure we got a proper fasta description */
  if (cc != GT_FASTA_SEPARATOR) {
    gt_error_set(err, "the first character of fasta file \"%s\" has to be '%c'",
                 gt_io_get_filename(seqio), GT_FASTA_SEPARATOR);
    return -1;
  }
  /* read description */
  while (!gt_io_get_char(seqio, &cc) && cc != '\n') {
    if (cc != '\r')
      gt_str_append_char(description, cc);
  }
  return 0;
}

static int parse_fasta_sequence(GtStr *sequence, GtIO *seqio, GtError *err)
{
  char cc;
  gt_error_check(err);
  gt_assert(sequence && seqio);
  gt_assert(!gt_str_length(sequence));
  /* read sequence */
  while (!gt_io_get_char(seqio, &cc) && cc != GT_FASTA_SEPARATOR) {
    if (cc != '\n' && cc != ' ' && cc != '\r')
      gt_str_append_char(sequence, cc);
  }
  if (!gt_str_length(sequence)) {
    gt_error_set(err, "empty sequence given in line %lu",
              gt_io_get_line_number(seqio));
    return -1;
  }
  if (cc == GT_FASTA_SEPARATOR)
    gt_io_unget_char(seqio, GT_FASTA_SEPARATOR);
  return 0;
}

static int parse_fasta_entry(GtStr *description, GtStr *sequence,
                             GtIO *seqio, GtError *err)
{
  int had_err;
  gt_error_check(err);
  gt_assert(description && sequence && seqio);
  had_err = parse_fasta_description(description, seqio, err);
  if (!had_err)
    had_err = parse_fasta_sequence(sequence, seqio, err);
  return had_err;
}

static int gt_fasta_reader_rec_run(GtFastaReader *fasta_reader,
                                   GtFastaReaderProcDescription
                                   proc_description,
                                   GtFastaReaderProcSequencePart
                                   proc_sequence_part,
                                   GtFastaReaderProcSequenceLength
                                   proc_sequence_length, void *data,
                                   GtError *err)
{
  GtFastaReaderRec *fr = gt_fasta_reader_rec_cast(fasta_reader);
  GtStr *description, *sequence;
  int had_err = 0;
  gt_error_check(err);

  /* at least one function has to be defined */
  gt_assert(proc_description || proc_sequence_part || proc_sequence_length);

  /* init */
  description = gt_str_new();
  sequence    = gt_str_new();

  /* make sure file is not empty */
  if (!gt_io_has_char(fr->seqio)) {
    gt_error_set(err, "sequence file \"%s\" is empty",
                 gt_io_get_filename(fr->seqio));
    had_err = -1;
  }

  /* parse file */
  while (!had_err && gt_io_has_char(fr->seqio)) {
    /* reset */
    gt_str_reset(description);
    gt_str_reset(sequence);

    /* parse entry */
    had_err = parse_fasta_entry(description, sequence, fr->seqio, err);

    /* process entry */
    if (!had_err && proc_description) {
      had_err = proc_description(gt_str_get(description),
                                 gt_str_length(description), data, err);
    }
    if (!had_err && proc_sequence_part) {
      had_err = proc_sequence_part(gt_str_get(sequence),
                                   gt_str_length(sequence), data, err);
    }
    if (!had_err && proc_sequence_length)
      had_err = proc_sequence_length(gt_str_length(sequence), data, err);
  }

  /* free */
  gt_str_delete(description);
  gt_str_delete(sequence);

  return had_err;
}

static void gt_fasta_reader_rec_free(GtFastaReader *fr)
{
  GtFastaReaderRec *gt_fasta_reader_rec = gt_fasta_reader_rec_cast(fr);
  gt_io_delete(gt_fasta_reader_rec->seqio);
}

const GtFastaReaderClass* gt_fasta_reader_rec_class(void)
{
  static const GtFastaReaderClass frc = { sizeof (GtFastaReaderRec),
                                        gt_fasta_reader_rec_run,
                                        gt_fasta_reader_rec_free };
  return &frc;
}

GtFastaReader* gt_fasta_reader_rec_new(GtStr *sequence_filename)
{
  GtFastaReader *fr = gt_fasta_reader_create(gt_fasta_reader_rec_class());
  GtFastaReaderRec *gt_fasta_reader_rec = gt_fasta_reader_rec_cast(fr);
  gt_fasta_reader_rec->seqio = gt_io_new(sequence_filename
                                   ? gt_str_get(sequence_filename) : NULL,
                                   "r");
  return fr;
}
