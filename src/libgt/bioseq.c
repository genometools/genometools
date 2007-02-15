/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include "array.h"
#include "bioseq.h"
#include "error.h"
#include "fasta.h"
#include "fasta_reader.h"
#include "fileutils.h"
#include "grep.h"
#include "range.h"
#include "sig.h"
#include "str.h"
#include "undef.h"
#include "xansi.h"
#include "xposix.h"

struct Bioseq {
  Str *sequence_file;
  Seq **seqs;
  Array *descriptions,
        *sequence_ranges;
  char *raw_sequence;
  size_t raw_sequence_length;
  Alpha *alpha;
  bool bioseq_is_filled;
};

typedef struct {
  FILE *bioseq_index,
       *bioseq_raw;
  unsigned long offset;
} Construct_bioseq_files_info;

Bioseq* bioseq_new(const char *sequence_file)
{
  Bioseq *bs = xcalloc(1, sizeof(Bioseq));
  if (!file_exists(sequence_file)) {
    error("sequence file \"%s\" does not exist or is not readable",
          sequence_file);
  }
  bs->sequence_file = str_new_cstr(sequence_file);
  bs->descriptions = array_new(sizeof(char*));
  bs->sequence_ranges = array_new(sizeof(Range));
  return bs;
}

Bioseq* bioseq_new_str(Str *sequence_file)
{
  Bioseq *bs = xcalloc(1, sizeof(Bioseq));
  if (!file_exists(str_get(sequence_file))) {
    error("sequence file \"%s\" does not exist or is not readable",
          str_get(sequence_file));
  }
  bs->sequence_file = str_ref(sequence_file);
  bs->descriptions = array_new(sizeof(char*));
  bs->sequence_ranges = array_new(sizeof(Range));
  return bs;
}

static void proc_description(Str *description, void *data)
{
  Construct_bioseq_files_info *info = (Construct_bioseq_files_info*) data;
  if (str_length(description))
    xfputs(str_get(description), info->bioseq_index);
  xfputc('\n', info->bioseq_index);
}

static void proc_character(char character, void *data)
{
  Construct_bioseq_files_info *info = (Construct_bioseq_files_info*) data;
  xfputc(character, info->bioseq_raw);
}

static void proc_sequence_length(unsigned long sequence_length, void *data)
{
  Construct_bioseq_files_info *info = (Construct_bioseq_files_info*) data;
  fprintf(info->bioseq_index, "%lu\n", info->offset);
  assert(sequence_length);
  fprintf(info->bioseq_index, "%lu\n", info->offset + sequence_length - 1);
  info->offset += sequence_length;
}

/* this global variables are necessary for the signal handler below */
static Construct_bioseq_files_info bioseq_files_info;
static const char *bioseq_index_filename,
                  *bioseq_raw_filename;

/* removes the incomplete bioseq files */
static void remove_bioseq_files(int sigraised)
{
  /* we don't care if fclose() succeeds, xunlink() will take care of it */
  (void) fclose(bioseq_files_info.bioseq_index);
  (void) fclose(bioseq_files_info.bioseq_raw);
  xunlink(bioseq_index_filename);
  xunlink(bioseq_raw_filename);
  (void) xsignal(sigraised, SIG_DFL);
  xraise(sigraised);
}

static void construct_bioseq_files(Str *bioseq_index_file, Str *bioseq_raw_file,
                                   Str *sequence_file)
{
  FastaReader *fasta_reader;

  /* open files & init */
  bioseq_files_info.bioseq_index = xfopen(str_get(bioseq_index_file), "w");
  bioseq_files_info.bioseq_raw = xfopen(str_get(bioseq_raw_file), "w");
  bioseq_files_info.offset = 0;

  /* register the signal handler to remove incomplete files upon termination */
  bioseq_index_filename = str_get(bioseq_index_file);
  bioseq_raw_filename = str_get(bioseq_raw_file);
  sig_register_all(remove_bioseq_files);

  /* read fasta file */
  fasta_reader = fasta_reader_new(sequence_file);
  fasta_reader_run(fasta_reader, proc_description, proc_character,
                   proc_sequence_length, &bioseq_files_info);
  fasta_reader_free(fasta_reader);

  /* unregister the signal handler */
  sig_unregister_all();

  /* close files */
  xfclose(bioseq_files_info.bioseq_index);
  xfclose(bioseq_files_info.bioseq_raw);
}

static void fill_bioseq(Bioseq *bs, const char *index_filename,
                        const char *raw_filename)
{
  FILE *index_file;
  Str *index_line;
  unsigned long line_number = 1;
  char *description;
  Range range;

  /* parse the index file and fill the sequence to index mapping */
  index_line = str_new();
  index_file = xfopen(index_filename, "r");

  while (str_read_next_line(index_line, index_file) != EOF) {
    switch (line_number % 3) {
      case 1:
        /* process description */
        description = xstrdup(str_get(index_line));
        array_add(bs->descriptions, description);
        break;
      case 2:
        /* process sequence start */
        if (sscanf(str_get(index_line), "%lu", &range.start) != 1) {
          error("could not parse bioseq start in line %lu of file \"%s\"",
                line_number, index_filename);
        }
        break;
      case 0:
        /* process sequence end */
        if (sscanf(str_get(index_line), "%lu", &range.end) != 1) {
          error("could not parse bioseq end in line %lu of file \"%s\"",
                line_number, index_filename);
        }
        assert(range.start <= range.end); /* XXX */
        array_add(bs->sequence_ranges, range);
        break;
    }
    line_number++;
    str_reset(index_line);
  }
  xfclose(index_file);
  str_free(index_line);

  /* the number of descriptions equals the number of sequence ranges */
  assert(array_size(bs->descriptions) == array_size(bs->sequence_ranges));

  /* map the raw file */
  bs->raw_sequence = xmap_read(raw_filename, &bs->raw_sequence_length);
}

void bioseq_fill(Bioseq *bs, unsigned int recreate)
{
  Str *bioseq_index_file,
      *bioseq_raw_file;

  if (bs->bioseq_is_filled)
    return;

  assert(!bs->raw_sequence);

  /* construct file names */
  bioseq_index_file = str_clone(bs->sequence_file);
  str_append_cstr(bioseq_index_file, GT_BIOSEQ_INDEX);
  bioseq_raw_file = str_clone(bs->sequence_file);
  str_append_cstr(bioseq_raw_file, GT_BIOSEQ_RAW);

  /* construct the bioseq files if necessary */
  if (recreate ||
      !file_exists(str_get(bioseq_index_file)) ||
      !file_exists(str_get(bioseq_raw_file)) ||
      file_is_newer(str_get(bs->sequence_file), str_get(bioseq_index_file)) ||
      file_is_newer(str_get(bs->sequence_file), str_get(bioseq_raw_file))) {
    construct_bioseq_files(bioseq_index_file, bioseq_raw_file,
                           bs->sequence_file);
  }

  /* fill the bioseq */
  fill_bioseq(bs, str_get(bioseq_index_file), str_get(bioseq_raw_file));

  /* free */
  str_free(bioseq_index_file);
  str_free(bioseq_raw_file);

  bs->bioseq_is_filled = true;
}

Seq* bioseq_get_seq(Bioseq *bs, unsigned long idx)
{
  assert(bs);
  if (!bs->bioseq_is_filled)
    bioseq_fill(bs, 0);
  assert(idx < array_size(bs->descriptions));
  if (!bs->seqs)
    bs->seqs = xcalloc(array_size(bs->descriptions), sizeof(Seq*));
  if (!bs->alpha) {
    bs->alpha = alpha_guess(bioseq_get_raw_sequence(bs),
                            bioseq_get_raw_sequence_length(bs));
  }
  if (!bs->seqs[idx]) {
    bs->seqs[idx] = seq_new(bioseq_get_sequence(bs, idx),
                            bioseq_get_sequence_length(bs, idx), bs->alpha);
    seq_set_description(bs->seqs[idx], bioseq_get_description(bs, idx));
  }
  return bs->seqs[idx];
}

const char* bioseq_get_description(Bioseq *bs, unsigned long idx)
{
  assert(bs);
  if (!bs->bioseq_is_filled)
    bioseq_fill(bs, 0);
  return *(char**) array_get(bs->descriptions, idx);
}

const char* bioseq_get_sequence(Bioseq *bs, unsigned long idx)
{
  Range sequence_range;
  assert(bs);
  if (!bs->bioseq_is_filled)
    bioseq_fill(bs, 0);
  sequence_range = *(Range*) array_get(bs->sequence_ranges, idx);
  return bs->raw_sequence + sequence_range.start;
}

const char* bioseq_get_raw_sequence(Bioseq *bs)
{
  assert(bs);
  if (!bs->bioseq_is_filled)
    bioseq_fill(bs, 0);
  return bs->raw_sequence;
}

static unsigned long get_seqnum_with_desc(Bioseq *bs, const char *description)
{
  unsigned long i, seqnum = UNDEFULONG;
  Str *pattern;
  assert(bs && description);
  pattern = str_new();
  str_append_char(pattern, '^');
  str_append_cstr(pattern, description);
  for (i = 0; i < array_size(bs->descriptions); i++) {
    if (grep(str_get(pattern), *(char**) array_get(bs->descriptions, i))) {
      seqnum = i;
      break;
    }
  }
  str_free(pattern);
  return seqnum;
}

const char* bioseq_get_sequence_with_desc(Bioseq *bs, const char *description,
                                          unsigned long *seqnum)
{
  unsigned long rseqnum;
  assert(bs && description);
  rseqnum = get_seqnum_with_desc(bs, description);
  if (seqnum)
    *seqnum = rseqnum;
  if (rseqnum != UNDEFULONG)
    return bioseq_get_sequence(bs, rseqnum);
  return NULL;
}

unsigned long bioseq_get_sequence_length(Bioseq *bs, unsigned long idx)
{
  Range sequence_range;
  assert(bs);
  if (!bs->bioseq_is_filled)
    bioseq_fill(bs, 0);
  sequence_range = *(Range*) array_get(bs->sequence_ranges, idx);
  return range_length(sequence_range);
}

unsigned long bioseq_get_raw_sequence_length(Bioseq *bs)
{
  assert(bs);
  if (!bs->bioseq_is_filled)
    bioseq_fill(bs, 0);
  return bs->raw_sequence_length;
}

unsigned long bioseq_number_of_sequences(Bioseq *bs)
{
  assert(bs);
  if (!bs->bioseq_is_filled)
    bioseq_fill(bs, 0);
  return array_size(bs->descriptions);
}

bool bioseq_contains_sequence(Bioseq *bs, /*@unused@*/ const char *sequence)
{
  assert(bs);
  if (!bs->bioseq_is_filled)
    bioseq_fill(bs, 0);
  if (get_seqnum_with_desc(bs, sequence) != UNDEFULONG)
    return true;
  return false;
}

void bioseq_free(Bioseq *bs)
{
  unsigned long i;

  if (!bs) return;
  str_free(bs->sequence_file);
  if (bs->seqs) {
    for (i = 0; i < array_size(bs->descriptions); i++)
      seq_free(bs->seqs[i]);
    free(bs->seqs);
  }
  for (i = 0; i < array_size(bs->descriptions); i++)
    free(*(char**) array_get(bs->descriptions, i));
  array_free(bs->descriptions);
  array_free(bs->sequence_ranges);
  xmunmap(bs->raw_sequence, bs->raw_sequence_length);
  alpha_free(bs->alpha);
  free(bs);
}

void bioseq_show_as_fasta(Bioseq *bs, unsigned long width)
{
  unsigned long i;

  assert(bs);

  for (i = 0; i < bioseq_number_of_sequences(bs); i++) {
    fasta_show_entry(bioseq_get_description(bs, i), bioseq_get_sequence(bs, i),
                     bioseq_get_sequence_length(bs, i), width);
  }
}

void bioseq_show_sequence_as_fasta(Bioseq *bs, unsigned long seqnum,
                                   unsigned long width)
{
  assert(bs);
  assert(seqnum < bioseq_number_of_sequences(bs));

  fasta_show_entry(bioseq_get_description(bs, seqnum),
                   bioseq_get_sequence(bs, seqnum),
                   bioseq_get_sequence_length(bs, seqnum), width);

}

void bioseq_show_stat(Bioseq *bs)
{
  unsigned long i, num_of_seqs;
  assert(bs);
  num_of_seqs = bioseq_number_of_sequences(bs);
  printf("showing statistics for sequence file \"%s\"\n",
         str_get(bs->sequence_file));
  printf("number of sequences: %lu\n", num_of_seqs);
  printf("total length: %lu\n", bioseq_get_raw_sequence_length(bs));
  for (i = 0; i < num_of_seqs; i++) {
    printf("sequence #%lu length: %lu\n", i+1,
           bioseq_get_sequence_length(bs, i));
  }
}
