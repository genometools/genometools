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
#include "core/array.h"
#include "core/fasta_separator.h"
#include "core/str.h"
#include "core/xansi.h"
#include "exercise/simple_bioseq.h"

typedef struct {
  GtStr *description,
         *sequence;
} FastaEntry;

struct GtSimpleBioseq {
  GtArray *entries;
};

static bool has_char(FILE *fp)
{
  int cc = gt_xfgetc(fp);
  if (cc == EOF)
    return false;
  gt_xungetc(cc, fp);
  return true;
}

static void parse_fasta_description(GtStr *description, FILE *fp,
                                    const char *fasta_file)
{
  int cc;
  gt_assert(description && fp && fasta_file);
  cc = gt_xfgetc(fp);
  if (cc != FASTA_SEPARATOR) {
    fprintf(stderr, "the first character of fasta file \"%s\" has to be '%c'\n",
            fasta_file, FASTA_SEPARATOR);
    exit(EXIT_FAILURE);
  }
  /* read description */
  while ((cc = gt_xfgetc(fp)) != EOF && cc != '\n')
    gt_str_append_char(description, cc);
}

static void parse_fasta_sequence(GtStr *sequence, FILE *fp,
                                 const char *fasta_file)
{
  int cc;
  gt_assert(sequence && fp && fasta_file);
  /* read sequence */
  while ((cc = gt_xfgetc(fp)) != EOF && cc != FASTA_SEPARATOR) {
    if (cc != '\n' && cc != ' ') /* skip newlines and blanks */
      gt_str_append_char(sequence, cc);
  }
  if (!gt_str_length(sequence)) {
    fprintf(stderr, "empty sequence entry in fasta file \"%s\"\n", fasta_file);
    exit(EXIT_FAILURE);
  }
  if (cc == FASTA_SEPARATOR) /* lookahead */
    gt_xungetc(FASTA_SEPARATOR, fp);
}

static void parse_fasta_entry(GtArray *entries, FILE *fp,
                              const char *fasta_file)
{
  FastaEntry new_entry;
  gt_assert(entries && fp && fasta_file);
  new_entry.description = gt_str_new();
  new_entry.sequence = gt_str_new();
  parse_fasta_description(new_entry.description, fp, fasta_file);
  parse_fasta_sequence(new_entry.sequence, fp, fasta_file);
  gt_array_add(entries, new_entry);
}

static void parse_fasta_file(GtArray *entries, const char *fasta_file)
{
  FILE *fp;
  gt_assert(entries && fasta_file);
  fp = gt_xfopen(fasta_file, "r");
  if (!has_char(fp)) {
    fprintf(stderr, "sequence file \"%s\" is empty\n", fasta_file);
    exit(EXIT_FAILURE);
  }
  while (has_char(fp))
    parse_fasta_entry(entries, fp, fasta_file);
  gt_xfclose(fp);
}

GtSimpleBioseq* gt_simple_bioseq_new(const char *fasta_file)
{
  GtSimpleBioseq *sbs;
  gt_assert(fasta_file);
  sbs = gt_xmalloc(sizeof *sbs);
  sbs->entries = gt_array_new(sizeof (FastaEntry));
  parse_fasta_file(sbs->entries, fasta_file); /* top-down / recursive descent */
  return sbs;
}

void gt_simple_bioseq_delete(GtSimpleBioseq *sbs)
{
  unsigned long i;
  if (!sbs) return;
  for (i = 0; i < gt_array_size(sbs->entries); i++) {
    FastaEntry *fasta_entry = gt_array_get(sbs->entries, i);
    gt_str_delete(fasta_entry->description);
    gt_str_delete(fasta_entry->sequence);
  }
  gt_array_delete(sbs->entries);
  free(sbs);
}

const char* gt_simple_bioseq_get_description(GtSimpleBioseq *sbs,
                                          unsigned long index)
{
  gt_assert(sbs);
  return gt_str_get(((FastaEntry*) gt_array_get(sbs->entries,
                    index))->description);
}

const char* gt_simple_bioseq_get_sequence(GtSimpleBioseq *sbs,
                                          unsigned long index)
{
  gt_assert(sbs);
  return gt_str_get(((FastaEntry*) gt_array_get(sbs->entries,
                                                index))->sequence);
}

unsigned long gt_simple_bioseq_get_sequence_length(GtSimpleBioseq *sbs,
                                                   unsigned long index)
{
  gt_assert(sbs);
  return gt_str_length(((FastaEntry*) gt_array_get(sbs->entries,
                                                   index))->sequence);
}

unsigned long gt_simple_bioseq_number_of_sequences(GtSimpleBioseq *sbs)
{
  gt_assert(sbs);
  return gt_array_size(sbs->entries);
}
