/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include <stdio.h>
#include "error.h"
#include "fasta.h"
#include "fasta_reader.h"
#include "xansi.h"

struct FastaReader {
  Str *sequence_filename;
  FILE *sequence_file;
};

typedef enum {
  EXPECTING_SEPARATOR,
  READING_DESCRIPTION,
  READING_SEQUENCE_AFTER_NEWLINE,
  READING_SEQUENCE
} FastaReader_state;

FastaReader* fasta_reader_new(Str *sequence_filename)
{
  FastaReader *fs = xmalloc(sizeof(FastaReader));
  fs->sequence_filename = str_ref(sequence_filename);
  fs->sequence_file = xfopen(str_get(sequence_filename), "r");
  return fs;
}

void fasta_reader_run(FastaReader *fr,
                      FastaReader_proc_description proc_description,
                      FastaReader_proc_character proc_character,
                      FastaReader_proc_sequence_length proc_sequence_length,
                      void *data)
{
  unsigned char cc;
  FastaReader_state state = EXPECTING_SEPARATOR;
  unsigned long sequence_length = 0, line_counter = 1;
  Str *description;

  assert(fr);

  /* init */
  description = str_new();

  /* at least one function has to be defined */
  assert(proc_description || proc_character || proc_sequence_length);

  /* rewind sequence file (to allow multiple calls) */
  rewind(fr->sequence_file);

  /* reading */
  while (xfread(&cc, sizeof(unsigned char), 1, fr->sequence_file) != 0) {
    switch (state) {
      case EXPECTING_SEPARATOR:
        if (cc != FASTA_SEPARATOR)
          error("the first character of fasta file \"%s\" has to be '%c'",
                str_get(fr->sequence_filename), FASTA_SEPARATOR);
        state = READING_DESCRIPTION;
        break;
      case READING_DESCRIPTION:
        if (cc == '\n') {
          if (proc_description) {
            proc_description(description, data);
            str_reset(description);
          }
          sequence_length = 0;
          line_counter++;
          state = READING_SEQUENCE_AFTER_NEWLINE;
        }
        else if (proc_description)
          str_append_char(description, cc);
        break;
      case READING_SEQUENCE_AFTER_NEWLINE:
        if (cc == FASTA_SEPARATOR) {
          if (!sequence_length) {
            assert(line_counter);
            error("empty sequence after description given in line %lu",
                  line_counter - 1);
          }
          if (proc_sequence_length)
            proc_sequence_length(sequence_length, data);
          state = READING_DESCRIPTION;
          continue;
        }
        /*@fallthrough@*/
      case READING_SEQUENCE:
        if (cc == '\n') {
          line_counter++;
          state = READING_SEQUENCE_AFTER_NEWLINE;
        }
        else {
          sequence_length++;
          if (proc_character)
            proc_character(cc, data);
        }
        break;
    }
  }

  /* checks after reading */
  switch (state) {
    case EXPECTING_SEPARATOR:
      error("sequence file \"%s\" is empty", str_get(fr->sequence_filename));
      /* error() does not return, no break statement to get full coverage */
      /*@fallthrough@*/
    case READING_DESCRIPTION:
      error("unfinished fasta entry in line %lu of sequence file \"%s\"",
            line_counter, str_get(fr->sequence_filename));
      /* error() does not return, no break statement to get full coverage */
      /*@fallthrough@*/
    case READING_SEQUENCE_AFTER_NEWLINE:
    case READING_SEQUENCE:
      if (!sequence_length) {
        assert(line_counter);
        error("empty sequence after description given in line %lu",
              line_counter - 1);
      }
      if (proc_sequence_length)
        proc_sequence_length(sequence_length, data);
  }

  /* free */
  str_free(description);
}

void fasta_reader_free(FastaReader *fr)
{
  if (!fr) return;
  str_free(fr->sequence_filename);
  xfclose(fr->sequence_file);
  free(fr);
}
