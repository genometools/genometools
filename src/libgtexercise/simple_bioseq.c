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

#include <stdbool.h>
#include "libgtcore/array.h"
#include "libgtcore/fasta_separator.h"
#include "libgtcore/str.h"
#include "libgtcore/xansi.h"
#include "libgtexercise/simple_bioseq.h"

typedef struct {
  Str *description,
      *sequence;
} FastaEntry;

struct SimpleBioseq {
  Array *entries;
};

static bool has_char(FILE *fp)
{
  int cc = xfgetc(fp);
  if (cc == EOF)
    return false;
  xungetc(cc, fp);
  return true;
}

static void parse_fasta_description(Str *description, FILE *fp,
                                    const char *fasta_file)
{
  int cc;
  assert(description && fp && fasta_file);
  cc = xfgetc(fp);
  if (cc != FASTA_SEPARATOR) {
    fprintf(stderr, "the first character of fasta file \"%s\" has to be '%c'\n",
            fasta_file, FASTA_SEPARATOR);
    exit(EXIT_FAILURE);
  }
  /* read description */
  while ((cc = xfgetc(fp)) != EOF && cc != '\n')
    str_append_char(description, cc);
}

static void parse_fasta_sequence(Str *sequence, FILE *fp,
                                 const char *fasta_file)
{
  int cc;
  assert(sequence && fp && fasta_file);
  /* read sequence */
  while ((cc = xfgetc(fp)) != EOF && cc != FASTA_SEPARATOR) {
    if (cc != '\n' && cc != ' ') /* skip newlines and blanks */
      str_append_char(sequence, cc);
  }
  if (!str_length(sequence)) {
    fprintf(stderr, "empty sequence entry in fasta file \"%s\"\n", fasta_file);
    exit(EXIT_FAILURE);
  }
  if (cc == FASTA_SEPARATOR) /* lookahead */
    xungetc(FASTA_SEPARATOR, fp);
}

static void parse_fasta_entry(Array *entries, FILE *fp, const char *fasta_file)
{
  FastaEntry new_entry;
  assert(entries && fp && fasta_file);
  new_entry.description = str_new();
  new_entry.sequence = str_new();
  parse_fasta_description(new_entry.description, fp, fasta_file);
  parse_fasta_sequence(new_entry.sequence, fp, fasta_file);
  array_add(entries, new_entry);
}

static void parse_fasta_file(Array *entries, const char *fasta_file)
{
  FILE *fp;
  assert(entries && fasta_file);
  fp = xfopen(fasta_file, "r");
  if (!has_char(fp)) {
    fprintf(stderr, "sequence file \"%s\" is empty\n", fasta_file);
    exit(EXIT_FAILURE);
  }
  while (has_char(fp))
    parse_fasta_entry(entries, fp, fasta_file);
  xfclose(fp);
}

SimpleBioseq* simple_bioseq_new(const char *fasta_file)
{
  SimpleBioseq *sbs;
  assert(fasta_file);
  sbs = xmalloc(sizeof *sbs);
  sbs->entries = array_new(sizeof (FastaEntry));
  parse_fasta_file(sbs->entries, fasta_file); /* top-down / recursive descent */
  return sbs;
}

void simple_bioseq_delete(SimpleBioseq *sbs)
{
  unsigned long i;
  if (!sbs) return;
  for (i = 0; i < array_size(sbs->entries); i++) {
    FastaEntry *fasta_entry = array_get(sbs->entries, i);
    str_delete(fasta_entry->description);
    str_delete(fasta_entry->sequence);
  }
  array_delete(sbs->entries);
  free(sbs);
}

const char* simple_bioseq_get_description(SimpleBioseq *sbs,
                                          unsigned long index)
{
  assert(sbs);
  return str_get(((FastaEntry*) array_get(sbs->entries, index))->description);
}

const char* simple_bioseq_get_sequence(SimpleBioseq *sbs, unsigned long index)
{
  assert(sbs);
  return str_get(((FastaEntry*) array_get(sbs->entries, index))->sequence);
}

unsigned long simple_bioseq_get_sequence_length(SimpleBioseq *sbs,
                                                unsigned long index)
{
  assert(sbs);
  return str_length(((FastaEntry*) array_get(sbs->entries, index))->sequence);
}

unsigned long simple_bioseq_number_of_sequences(SimpleBioseq *sbs)
{
  assert(sbs);
  return array_size(sbs->entries);
}
