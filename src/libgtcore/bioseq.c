/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include <libgtcore/array.h>
#include <libgtcore/bioseq.h>
#include <libgtcore/cstr.h>
#include <libgtcore/dynalloc.h>
#include <libgtcore/error.h>
#include <libgtcore/fasta.h>
#include <libgtcore/fasta_reader.h>
#include <libgtcore/fileutils.h>
#include <libgtcore/gc_content.h>
#include <libgtcore/grep.h>
#include <libgtcore/range.h>
#include <libgtcore/sig.h>
#include <libgtcore/str.h>
#include <libgtcore/undef.h>
#include <libgtcore/xansi.h>
#include <libgtcore/xposix.h>

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

static int proc_description(Str *description, void *data, Env *env)
{
  Construct_bioseq_files_info *info = (Construct_bioseq_files_info*) data;
  char *description_cstr;
  env_error_check(env);
  if (info->bs->use_stdin) {
    description_cstr = cstr_dup(str_get(description), env);
    array_add(info->bs->descriptions, description_cstr, env);
  }
  else {
    if (str_length(description))
      xfputs(str_get(description), info->bioseq_index);
    xfputc('\n', info->bioseq_index);
  }
  return 0;
}

static int proc_character(char character, void *data, Env *env)
{
  Construct_bioseq_files_info *info = (Construct_bioseq_files_info*) data;
  env_error_check(env);
  if (info->bs->use_stdin) {
    info->bs->raw_sequence = dynalloc(info->bs->raw_sequence,
                                      &info->bs->allocated,
                                      info->bs->raw_sequence_length + 1, env);
    info->bs->raw_sequence[info->bs->raw_sequence_length++] = character;
  }
  else
    xfputc(character, info->bioseq_raw);
  return 0;
}

static int proc_sequence_length(unsigned long sequence_length, void *data,
                                Env *env)
{
  Construct_bioseq_files_info *info = (Construct_bioseq_files_info*) data;
  Range range;
  env_error_check(env);
  if (info->bs->use_stdin) {
    range.start = info->offset;
    range.end = info->offset + sequence_length - 1;
    array_add(info->bs->sequence_ranges, range, env);
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
                       const char *raw_filename, Env *env)
{
  FILE *index_file;
  Str *index_line;
  unsigned long line_number = 1;
  char *description;
  Range range;
  int had_err = 0;

  env_error_check(env);

  /* parse the index file and fill the sequence to index mapping */
  index_line = str_new(env);
  index_file = env_fa_xfopen(env, index_filename, "r");

  while (!had_err && str_read_next_line(index_line, index_file, env) != EOF) {
    switch (line_number % 3) {
      case 1:
        /* process description */
        description = cstr_dup(str_get(index_line), env);
        array_add(bs->descriptions, description, env);
        break;
      case 2:
        /* process sequence start */
        if (sscanf(str_get(index_line), "%lu", &range.start) != 1) {
          env_error_set(env, "could not parse bioseq start in line %lu of file "
                    "\"%s\"", line_number, index_filename);
          had_err = -1;
        }
        break;
      case 0:
        /* process sequence end */
        if (sscanf(str_get(index_line), "%lu", &range.end) != 1) {
          env_error_set(env, "could not parse bioseq end in line %lu of file "
                    "\"%s\"", line_number, index_filename);
          had_err = -1;
        }
        else {
          assert(range.start <= range.end); /* XXX */
          array_add(bs->sequence_ranges, range, env);
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
    bs->raw_sequence = env_fa_xmmap_read(env, raw_filename,
                                         &bs->raw_sequence_length);
  }

  env_fa_xfclose(index_file, env);
  str_delete(index_line, env);

  return had_err;
}

static int construct_bioseq_files(Bioseq *bs, Str *bioseq_index_file,
                                  Str *bioseq_raw_file, Env *env)
{
  FastaReader *fasta_reader;
  int had_err;

  env_error_check(env);

  /* open files & init */
  if (!bs->use_stdin) {
    bioseq_files_info.bioseq_index = env_fa_xfopen(env,
                                                   (const char *) str_get(bioseq_index_file),
                                                   "w");
    bioseq_files_info.bioseq_raw = env_fa_xfopen(env, str_get(bioseq_raw_file),
                                                 "w");
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
  fasta_reader = fasta_reader_new(bs->use_stdin ? NULL : bs->sequence_file,
                                  env);
  had_err = fasta_reader_run(fasta_reader, proc_description, proc_character,
                             proc_sequence_length, &bioseq_files_info, env);
  fasta_reader_delete(fasta_reader, env);

  /* unregister the signal handler */
  if (!bs->use_stdin)
    sig_unregister_all();

  /* close files */
  if (!bs->use_stdin) {
    env_fa_xfclose(bioseq_files_info.bioseq_index, env);
    env_fa_xfclose(bioseq_files_info.bioseq_raw, env);
  }

  return had_err;
}

static int bioseq_fill(Bioseq *bs, bool recreate, Env *env)
{
  Str *bioseq_index_file = NULL,
      *bioseq_raw_file = NULL;
  int had_err = 0;

  assert(!bs->raw_sequence);

  /* construct file names */
  if (!bs->use_stdin) {
    bioseq_index_file = str_clone(bs->sequence_file, env);
    str_append_cstr(bioseq_index_file, GT_BIOSEQ_INDEX, env);
    bioseq_raw_file = str_clone(bs->sequence_file, env);
    str_append_cstr(bioseq_raw_file, GT_BIOSEQ_RAW, env);
  }

  /* construct the bioseq files if necessary */
  if (recreate || bs->use_stdin ||
      !file_exists(str_get(bioseq_index_file)) ||
      !file_exists(str_get(bioseq_raw_file)) ||
      file_is_newer(str_get(bs->sequence_file), str_get(bioseq_index_file)) ||
      file_is_newer(str_get(bs->sequence_file), str_get(bioseq_raw_file))) {
    had_err = construct_bioseq_files(bs, bioseq_index_file, bioseq_raw_file,
                                     env);
  }

  if (!had_err && !bs->use_stdin) {
    /* fill the bioseq */
    had_err = fill_bioseq(bs, str_get(bioseq_index_file),
                          str_get(bioseq_raw_file), env);
  }

  /* free */
  str_delete(bioseq_index_file, env);
  str_delete(bioseq_raw_file, env);

  return had_err;
}

static Bioseq* bioseq_new_with_recreate(Str *sequence_file, bool recreate,
                                        Env *env)
{
  Bioseq *bs;
  int had_err = 0;
  env_error_check(env);
  bs = env_ma_calloc(env, 1, sizeof (Bioseq));
  if (!strcmp(str_get(sequence_file), "-"))
    bs->use_stdin = true;
  if (!bs->use_stdin && !file_exists(str_get(sequence_file))) {
    env_error_set(env, "sequence file \"%s\" does not exist or is not readable",
              str_get(sequence_file));
    had_err = -1;
  }
  if (!had_err) {
    bs->sequence_file = str_ref(sequence_file);
    bs->descriptions = array_new(sizeof (char*), env);
    bs->sequence_ranges = array_new(sizeof (Range), env);
    had_err = bioseq_fill(bs, recreate, env);
  }
  if (had_err) {
    bioseq_delete(bs, env);
    return NULL;
  }
  return bs;
}

Bioseq* bioseq_new(const char *sequence_file, Env *env)
{
  Bioseq *bs;
  Str *seqfile;
  env_error_check(env);
  seqfile = str_new_cstr(sequence_file, env);
  bs = bioseq_new_with_recreate(seqfile, false, env);
  str_delete(seqfile, env);
  return bs;
}

Bioseq* bioseq_new_recreate(const char *sequence_file, Env *env)
{
  Bioseq *bs;
  Str *seqfile;
  env_error_check(env);
  seqfile = str_new_cstr(sequence_file, env);
  bs = bioseq_new_with_recreate(seqfile, true, env);
  str_delete(seqfile, env);
  return bs;
}

Bioseq* bioseq_new_str(Str *sequence_file, Env *env)
{
  return bioseq_new_with_recreate(sequence_file, false, env);
}

static void determine_alpha_if_necessary(Bioseq *bs, Env *env)
{
  assert(bs);
  if (!bs->alpha) {
    bs->alpha = alpha_guess(bioseq_get_raw_sequence(bs),
                            bioseq_get_raw_sequence_length(bs), env);
  }
}

Alpha* bioseq_get_alpha(Bioseq *bs, Env *env)
{
  env_error_check(env);
  assert(bs);
  determine_alpha_if_necessary(bs, env);
  assert(bs->alpha);
  return bs->alpha;
}

Seq* bioseq_get_seq(Bioseq *bs, unsigned long idx, Env *env)
{
  assert(bs);
  assert(idx < array_size(bs->descriptions));
  if (!bs->seqs)
    bs->seqs = env_ma_calloc(env, array_size(bs->descriptions), sizeof (Seq*));
  determine_alpha_if_necessary(bs, env);
  if (!bs->seqs[idx]) {
    bs->seqs[idx] = seq_new(bioseq_get_sequence(bs, idx),
                            bioseq_get_sequence_length(bs, idx),
                            bs->alpha, env);
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

void bioseq_delete(Bioseq *bs, Env *env)
{
  unsigned long i;
  if (!bs) return;
  str_delete(bs->sequence_file, env);
  if (bs->seqs) {
    for (i = 0; i < array_size(bs->descriptions); i++)
      seq_delete(bs->seqs[i], env);
    env_ma_free(bs->seqs, env);
  }
  for (i = 0; i < array_size(bs->descriptions); i++)
    env_ma_free(*(char**) array_get(bs->descriptions, i), env);
  array_delete(bs->descriptions, env);
  array_delete(bs->sequence_ranges, env);
  if (bs->use_stdin)
    env_ma_free(bs->raw_sequence, env);
  else
    env_fa_xmunmap(bs->raw_sequence, env);
  alpha_delete(bs->alpha, env);
  env_ma_free(bs, env);
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

void bioseq_show_gc_content(Bioseq *bs, Env *env)
{
  Alpha *dna_alpha;
  env_error_check(env);
  assert(bs);
  determine_alpha_if_necessary(bs, env);
  dna_alpha = alpha_new_dna(env);
  if (alpha_is_compatible_with_alpha(bs->alpha, dna_alpha)) {
    printf("showing GC-content for sequence file \"%s\"\n",
           str_get(bs->sequence_file));
    gc_content_show(bioseq_get_raw_sequence(bs),
                    bioseq_get_raw_sequence_length(bs), bs->alpha, env);
  }
  alpha_delete(dna_alpha, env);
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
