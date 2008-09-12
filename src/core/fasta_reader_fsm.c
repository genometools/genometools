/*
  Copyright (c) 2006-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2008 Center for Bioinformatics, University of Hamburg

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

#include "core/fasta_reader_fsm.h"
#include "core/fasta_reader_rep.h"
#include "core/fasta_separator.h"

struct GT_FastaReaderFSM {
  const GT_FastaReader parent_instance;
  GtStr *sequence_filename;
  GT_GenFile *sequence_file;
};

typedef enum {
  EXPECTING_SEPARATOR,
  READING_DESCRIPTION,
  READING_SEQUENCE_AFTER_NEWLINE,
  READING_SEQUENCE
} GT_FastaReaderState;

#define gt_fasta_reader_fsm_cast(FR)\
        gt_fasta_reader_cast(gt_fasta_reader_fsm_class(), FR)

static int gt_fasta_reader_fsm_run(GT_FastaReader *fasta_reader,
                                   GT_FastaReaderProcDescription
                                   proc_description,
                                   GT_FastaReaderProcSequencePart
                                   proc_sequence_part,
                                   GT_FastaReaderProcSequenceLength
                                   proc_sequence_length,
                                   void *data, GT_Error *err)
{
  GT_FastaReaderFSM *fr = gt_fasta_reader_fsm_cast(fasta_reader);
  unsigned char cc;
  GT_FastaReaderState state = EXPECTING_SEPARATOR;
  unsigned long sequence_length = 0, line_counter = 1;
  GtStr *description, *sequence;
  int had_err = 0;

  gt_error_check(err);
  assert(fr);

  /* init */
  description = gt_str_new();
  sequence    = gt_str_new();

  /* at least one function has to be defined */
  assert(proc_description || proc_sequence_part || proc_sequence_length);

  /* rewind sequence file (to allow multiple calls) */
  if (fr->sequence_file)
    gt_genfile_xrewind(fr->sequence_file);

  /* reading */
  while (!had_err && gt_genfile_xread(fr->sequence_file, &cc, 1) != 0) {
    switch (state) {
      case EXPECTING_SEPARATOR:
        if (cc != FASTA_SEPARATOR) {
          gt_error_set(err,
                    "the first character of fasta file \"%s\" has to be '%c'",
                    gt_str_get(fr->sequence_filename), FASTA_SEPARATOR);
          had_err = -1;
        }
        else
          state = READING_DESCRIPTION;
        break;
      case READING_DESCRIPTION:
        if (cc == '\n') {
          if (proc_description) {
            had_err = proc_description(gt_str_get(description),
                                       gt_str_length(description), data, err);
            if (!had_err)
              gt_str_reset(description);
          }
          if (!had_err) {
            sequence_length = 0;
            line_counter++;
            state = READING_SEQUENCE_AFTER_NEWLINE;
          }
        }
        else if (proc_description)
          gt_str_append_char(description, cc);
        break;
      case READING_SEQUENCE_AFTER_NEWLINE:
        if (cc == FASTA_SEPARATOR) {
          if (!sequence_length) {
            assert(line_counter);
            gt_error_set(err, "empty sequence after description given in line "
                              "%lu", line_counter - 1);
            had_err = -1;
            break;
          }
          else {
            if (proc_sequence_part) {
              assert(gt_str_length(sequence));
              had_err = proc_sequence_part(gt_str_get(sequence),
                                           gt_str_length(sequence), data, err);
            }
            if (had_err)
              break;
            gt_str_reset(sequence);
            if (proc_sequence_length)
              had_err = proc_sequence_length(sequence_length, data, err);
            if (had_err)
              break;
            state = READING_DESCRIPTION;
            continue;
          }
        }
        /*@fallthrough@*/
      case READING_SEQUENCE:
        if (cc == '\n') {
          line_counter++;
          state = READING_SEQUENCE_AFTER_NEWLINE;
        }
        else {
          sequence_length++;
          if (proc_sequence_part) {
            if (gt_str_length(sequence) == BUFSIZ) {
              had_err = proc_sequence_part(gt_str_get(sequence),
                                           gt_str_length(sequence), data, err);
              if (had_err)
                break;
              gt_str_reset(sequence);
            }
            if (cc != ' ')
              gt_str_append_char(sequence, cc);
          }
        }
        break;
    }
  }

  if (!had_err) {
    /* checks after reading */
    switch (state) {
      case EXPECTING_SEPARATOR:
        gt_error_set(err, "sequence file \"%s\" is empty",
                  gt_str_get(fr->sequence_filename));
        had_err = -1;
        break;
      case READING_DESCRIPTION:
        gt_error_set(err, "unfinished fasta entry in line %lu of sequence file "
                  "\"%s\"", line_counter, gt_str_get(fr->sequence_filename));
        had_err = -1;
        break;
      case READING_SEQUENCE_AFTER_NEWLINE:
      case READING_SEQUENCE:
        if (!sequence_length) {
          assert(line_counter);
          gt_error_set(err, "empty sequence after description given in line "
                            "%lu", line_counter - 1);
          had_err = -1;
        }
        else {
          if (proc_sequence_part) {
            assert(gt_str_length(sequence));
            had_err = proc_sequence_part(gt_str_get(sequence),
                                         gt_str_length(sequence), data, err);
          }
          if (!had_err && proc_sequence_length)
            had_err = proc_sequence_length(sequence_length, data, err);
        }
    }
  }

  /* free */
  gt_str_delete(sequence);
  gt_str_delete(description);

  return had_err;
}

static void gt_fasta_reader_fsm_free(GT_FastaReader *fr)
{
  GT_FastaReaderFSM *gt_fasta_reader_fsm = gt_fasta_reader_fsm_cast(fr);
  gt_str_delete(gt_fasta_reader_fsm->sequence_filename);
  gt_genfile_close(gt_fasta_reader_fsm->sequence_file);
}

const GT_FastaReaderClass* gt_fasta_reader_fsm_class(void)
{
  static const GT_FastaReaderClass frc = { sizeof (GT_FastaReaderFSM),
                                        gt_fasta_reader_fsm_run,
                                        gt_fasta_reader_fsm_free };
  return &frc;
}

GT_FastaReader* gt_fasta_reader_fsm_new(GtStr *sequence_filename)
{
  GT_FastaReader *fr = gt_fasta_reader_create(gt_fasta_reader_fsm_class());
  GT_FastaReaderFSM *gt_fasta_reader_fsm = gt_fasta_reader_fsm_cast(fr);
  gt_fasta_reader_fsm->sequence_filename = gt_str_ref(sequence_filename);
  if (sequence_filename) {
    gt_fasta_reader_fsm->sequence_file =
      gt_genfile_xopen(gt_str_get(sequence_filename), "r");
  }
  else
    gt_fasta_reader_fsm->sequence_filename = gt_str_new_cstr("stdin");
  return fr;
}
