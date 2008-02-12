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

#include <assert.h>
#include "libgtcore/array.h"
#include "libgtcore/bioseq.h"
#include "libgtcore/cstr.h"
#include "libgtcore/discdistri.h"
#include "libgtcore/dynalloc.h"
#include "libgtcore/error.h"
#include "libgtcore/fa.h"
#include "libgtcore/fasta.h"
#include "libgtcore/fasta_reader.h"
#include "libgtcore/fasta_reader_fsm.h"
#include "libgtcore/fasta_reader_rec.h"
#include "libgtcore/fasta_reader_seqit.h"
#include "libgtcore/fileutils.h"
#include "libgtcore/gc_content.h"
#include "libgtcore/grep.h"
#include "libgtcore/ma.h"
#include "libgtcore/parseutils.h"
#include "libgtcore/range.h"
#include "libgtcore/sig.h"
#include "libgtcore/str.h"
#include "libgtcore/undef.h"
#include "libgtcore/unused.h"
#include "libgtcore/xansi.h"
#include "libgtcore/xposix.h"

struct Bioseq {
  bool use_stdin;
  Str *sequence_file;
  Seq **seqs;
  Array *descriptions,
        *sequence_ranges;
  char *raw_sequence;
  size_t raw_sequence_length,
         allocated;
  Alpha *alpha;
};

typedef struct {
  FILE *bioseq_index,
       *bioseq_raw;
  unsigned long offset;
  Bioseq *bs;
} Construct_bioseq_files_info;

static int proc_description(const char *description, unsigned long length,
                            void *data, UNUSED Error *err)
{
  Construct_bioseq_files_info *info = (Construct_bioseq_files_info*) data;
  char *description_cstr;
  error_check(err);
  if (info->bs->use_stdin) {
    description_cstr = cstr_dup(description);
    array_add(info->bs->descriptions, description_cstr);
  }
  else {
    if (length)
      xfputs(description, info->bioseq_index);
    xfputc('\n', info->bioseq_index);
  }
  return 0;
}

static int proc_sequence_part(const char *seqpart, unsigned long length,
                              void *data, UNUSED Error *err)
{
  Construct_bioseq_files_info *info = (Construct_bioseq_files_info*) data;
  error_check(err);
  assert(seqpart);
  if (info->bs->use_stdin) {
    info->bs->raw_sequence = dynalloc(info->bs->raw_sequence,
                                      &info->bs->allocated,
                                      info->bs->raw_sequence_length + length);
    memcpy(info->bs->raw_sequence + info->bs->raw_sequence_length, seqpart,
           length);
    info->bs->raw_sequence_length += length;
  }
  else
    xfputs(seqpart, info->bioseq_raw);
  return 0;
}

static int proc_sequence_length(unsigned long sequence_length, void *data,
                                UNUSED Error *err)
{
  Construct_bioseq_files_info *info = (Construct_bioseq_files_info*) data;
  Range range;
  error_check(err);
  if (info->bs->use_stdin) {
    range.start = info->offset;
    range.end = info->offset + sequence_length - 1;
    array_add(info->bs->sequence_ranges, range);
  }
  else {
    fprintf(info->bioseq_index, "%lu\n", info->offset);
    assert(sequence_length);
    fprintf(info->bioseq_index, "%lu\n", info->offset + sequence_length - 1);
  }
  info->offset += sequence_length;
  return 0;
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

static int fill_bioseq(Bioseq *bs, const char *index_filename,
                       const char *raw_filename, Error *err)
{
  FILE *index_file;
  Str *index_line;
  unsigned long line_number = 1;
  char *description;
  Range range;
  int had_err = 0;

  error_check(err);

  /* parse the index file and fill the sequence to index mapping */
  index_line = str_new();
  index_file = fa_xfopen(index_filename, "r");

  while (!had_err && str_read_next_line(index_line, index_file) != EOF) {
    switch (line_number % 3) {
      case 1:
        /* process description */
        description = cstr_dup(str_get(index_line));
        array_add(bs->descriptions, description);
        break;
      case 2:
        /* process sequence start */
        if (parse_ulong(&range.start, str_get(index_line))) {
          error_set(err, "could not parse bioseq start in line %lu of file "
                         "\"%s\"", line_number, index_filename);
          had_err = -1;
        }
        break;
      case 0:
        /* process sequence end */
        if (parse_ulong(&range.end, str_get(index_line))) {
          error_set(err,
                    "could not parse bioseq end in line %lu of file \"%s\"",
                    line_number, index_filename);
          had_err = -1;
        }
        else {
          assert(range.start <= range.end); /* XXX */
          array_add(bs->sequence_ranges, range);
        }
        break;
    }
    line_number++;
    str_reset(index_line);
  }

  if (!had_err) {
    /* the number of descriptions equals the number of sequence ranges */
    assert(array_size(bs->descriptions) == array_size(bs->sequence_ranges));
    /* map the raw file */
    bs->raw_sequence = fa_xmmap_read(raw_filename, &bs->raw_sequence_length);
  }

  fa_xfclose(index_file);
  str_delete(index_line);

  return had_err;
}

static int construct_bioseq_files(Bioseq *bs, Str *bioseq_index_file,
                                  Str *bioseq_raw_file,
                                  FastaReaderType fasta_reader_type, Error *err)
{
  FastaReader *fasta_reader = NULL;
  Str *sequence_filename;
  int had_err;

  error_check(err);

  /* open files & init */
  if (!bs->use_stdin) {
    bioseq_files_info.bioseq_index = fa_xfopen((const char *)
                                               str_get(bioseq_index_file), "w");
    bioseq_files_info.bioseq_raw = fa_xfopen(str_get(bioseq_raw_file), "w");
  }
  bioseq_files_info.offset = 0;
  bioseq_files_info.bs = bs;

  /* register the signal handler to remove incomplete files upon termination */
  if (!bs->use_stdin) {
    bioseq_index_filename = str_get(bioseq_index_file);
    bioseq_raw_filename = str_get(bioseq_raw_file);
    sig_register_all(remove_bioseq_files);
  }

  /* read fasta file */
  sequence_filename = bs->use_stdin ? NULL : bs->sequence_file;
  switch (fasta_reader_type) {
    case FASTA_READER_REC:
      fasta_reader = fasta_reader_rec_new(sequence_filename);
      break;
    case FASTA_READER_FSM:
      fasta_reader = fasta_reader_fsm_new(sequence_filename);
      break;
    case FASTA_READER_SEQIT:
      fasta_reader = fasta_reader_seqit_new(sequence_filename);
      break;
    default: assert(0);
  }
  had_err = fasta_reader_run(fasta_reader, proc_description, proc_sequence_part,
                             proc_sequence_length, &bioseq_files_info, err);
  fasta_reader_delete(fasta_reader);

  /* unregister the signal handler */
  if (!bs->use_stdin)
    sig_unregister_all();

  /* close files */
  if (!bs->use_stdin) {
    fa_xfclose(bioseq_files_info.bioseq_index);
    fa_xfclose(bioseq_files_info.bioseq_raw);
    if (had_err) {
      xunlink(bioseq_index_filename);
      xunlink(bioseq_raw_filename);
    }
  }

  return had_err;
}

static int bioseq_fill(Bioseq *bs, bool recreate,
                       FastaReaderType fasta_reader_type, Error *err)
{
  Str *bioseq_index_file = NULL,
      *bioseq_raw_file = NULL;
  int had_err = 0;

  assert(!bs->raw_sequence);

  /* construct file names */
  if (!bs->use_stdin) {
    bioseq_index_file = str_clone(bs->sequence_file);
    str_append_cstr(bioseq_index_file, GT_BIOSEQ_INDEX);
    bioseq_raw_file = str_clone(bs->sequence_file);
    str_append_cstr(bioseq_raw_file, GT_BIOSEQ_RAW);
  }

  /* construct the bioseq files if necessary */
  if (recreate || bs->use_stdin ||
      !file_exists(str_get(bioseq_index_file)) ||
      !file_exists(str_get(bioseq_raw_file)) ||
      file_is_newer(str_get(bs->sequence_file), str_get(bioseq_index_file)) ||
      file_is_newer(str_get(bs->sequence_file), str_get(bioseq_raw_file))) {
    had_err = construct_bioseq_files(bs, bioseq_index_file, bioseq_raw_file,
                                     fasta_reader_type, err);
  }

  if (!had_err && !bs->use_stdin) {
    /* fill the bioseq */
    had_err = fill_bioseq(bs, str_get(bioseq_index_file),
                          str_get(bioseq_raw_file), err);
  }

  /* free */
  str_delete(bioseq_index_file);
  str_delete(bioseq_raw_file);

  return had_err;
}

static Bioseq* bioseq_new_with_recreate_and_type(Str *sequence_file,
                                                 bool recreate,
                                                 FastaReaderType
                                                 fasta_reader_type,
                                                 Error *err)
{
  Bioseq *bs;
  int had_err = 0;
  error_check(err);
  bs = ma_calloc(1, sizeof (Bioseq));
  if (!strcmp(str_get(sequence_file), "-"))
    bs->use_stdin = true;
  if (!bs->use_stdin && !file_exists(str_get(sequence_file))) {
    error_set(err, "sequence file \"%s\" does not exist or is not readable",
              str_get(sequence_file));
    had_err = -1;
  }
  if (!had_err) {
    bs->sequence_file = str_ref(sequence_file);
    bs->descriptions = array_new(sizeof (char*));
    bs->sequence_ranges = array_new(sizeof (Range));
    had_err = bioseq_fill(bs, recreate, fasta_reader_type, err);
  }
  if (had_err) {
    bioseq_delete(bs);
    return NULL;
  }
  return bs;
}

Bioseq* bioseq_new(const char *sequence_file, Error *err)
{
  Bioseq *bs;
  Str *seqfile;
  error_check(err);
  seqfile = str_new_cstr(sequence_file);
  bs = bioseq_new_with_recreate_and_type(seqfile, false, FASTA_READER_REC, err);
  str_delete(seqfile);
  return bs;
}

Bioseq* bioseq_new_recreate(const char *sequence_file, Error *err)
{
  Bioseq *bs;
  Str *seqfile;
  error_check(err);
  seqfile = str_new_cstr(sequence_file);
  bs = bioseq_new_with_recreate_and_type(seqfile, true, FASTA_READER_REC, err);
  str_delete(seqfile);
  return bs;
}

Bioseq* bioseq_new_str(Str *sequence_file, Error *err)
{
  return bioseq_new_with_recreate_and_type(sequence_file, false,
                                           FASTA_READER_REC, err);
}

Bioseq* bioseq_new_with_fasta_reader(const char *sequence_file,
                                     FastaReaderType fasta_reader, Error *err)
{
  Bioseq *bs;
  Str *seqfile;
  error_check(err);
  seqfile = str_new_cstr(sequence_file);
  bs = bioseq_new_with_recreate_and_type(seqfile, true, fasta_reader, err);
  str_delete(seqfile);
  return bs;
}

static void determine_alpha_if_necessary(Bioseq *bs)
{
  assert(bs);
  if (!bs->alpha) {
    bs->alpha = alpha_guess(bioseq_get_raw_sequence(bs),
                            bioseq_get_raw_sequence_length(bs));
  }
}

Alpha* bioseq_get_alpha(Bioseq *bs)
{
  assert(bs);
  determine_alpha_if_necessary(bs);
  assert(bs->alpha);
  return bs->alpha;
}

Seq* bioseq_get_seq(Bioseq *bs, unsigned long idx)
{
  assert(bs);
  assert(idx < array_size(bs->descriptions));
  if (!bs->seqs)
    bs->seqs = ma_calloc(array_size(bs->descriptions), sizeof (Seq*));
  determine_alpha_if_necessary(bs);
  if (!bs->seqs[idx]) {
    bs->seqs[idx] = seq_new(bioseq_get_sequence(bs, idx),
                            bioseq_get_sequence_length(bs, idx),
                            bs->alpha);
    seq_set_description(bs->seqs[idx], bioseq_get_description(bs, idx));
  }
  return bs->seqs[idx];
}

const char* bioseq_get_description(Bioseq *bs, unsigned long idx)
{
  assert(bs);
  return *(char**) array_get(bs->descriptions, idx);
}

const char* bioseq_get_sequence(Bioseq *bs, unsigned long idx)
{
  Range sequence_range;
  assert(bs);
  sequence_range = *(Range*) array_get(bs->sequence_ranges, idx);
  return bs->raw_sequence + sequence_range.start;
}

const char* bioseq_get_raw_sequence(Bioseq *bs)
{
  assert(bs);
  return bs->raw_sequence;
}

unsigned long bioseq_get_sequence_length(Bioseq *bs, unsigned long idx)
{
  Range sequence_range;
  assert(bs);
  sequence_range = *(Range*) array_get(bs->sequence_ranges, idx);
  return range_length(sequence_range);
}

unsigned long bioseq_get_raw_sequence_length(Bioseq *bs)
{
  assert(bs);
  return bs->raw_sequence_length;
}

unsigned long bioseq_number_of_sequences(Bioseq *bs)
{
  assert(bs);
  return array_size(bs->descriptions);
}

void bioseq_delete(Bioseq *bs)
{
  unsigned long i;
  if (!bs) return;
  str_delete(bs->sequence_file);
  if (bs->seqs) {
    for (i = 0; i < array_size(bs->descriptions); i++)
      seq_delete(bs->seqs[i]);
    ma_free(bs->seqs);
  }
  for (i = 0; i < array_size(bs->descriptions); i++)
    ma_free(*(char**) array_get(bs->descriptions, i));
  array_delete(bs->descriptions);
  array_delete(bs->sequence_ranges);
  if (bs->use_stdin)
    ma_free(bs->raw_sequence);
  else
    fa_xmunmap(bs->raw_sequence);
  alpha_delete(bs->alpha);
  ma_free(bs);
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

void bioseq_show_gc_content(Bioseq *bs)
{
  Alpha *dna_alpha;
  assert(bs);
  determine_alpha_if_necessary(bs);
  dna_alpha = alpha_new_dna();
  if (alpha_is_compatible_with_alpha(bs->alpha, dna_alpha)) {
    printf("showing GC-content for sequence file \"%s\"\n",
           str_get(bs->sequence_file));
    gc_content_show(bioseq_get_raw_sequence(bs),
                    bioseq_get_raw_sequence_length(bs), bs->alpha);
  }
  alpha_delete(dna_alpha);
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

void bioseq_show_seqlengthdistri(Bioseq *bs)
{
  DiscDistri *d;
  unsigned long i;
  assert(bs);
  d = discdistri_new();
  for (i = 0; i < bioseq_number_of_sequences(bs); i++)
    discdistri_add(d, bioseq_get_sequence_length(bs, i));
  printf("sequence length distribution:\n");
  discdistri_show(d);
  discdistri_delete(d);
}
