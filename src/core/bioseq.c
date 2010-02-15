/*
  Copyright (c) 2006-2009 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#include "core/assert_api.h"
#include "core/array.h"
#include "core/bioseq.h"
#include "core/cstr_api.h"
#include "core/disc_distri.h"
#include "core/dynalloc.h"
#include "core/error.h"
#include "core/fa.h"
#include "core/fasta.h"
#include "core/fasta_reader.h"
#include "core/fasta_reader_fsm.h"
#include "core/fasta_reader_rec.h"
#include "core/fasta_reader_seqit.h"
#include "core/fileutils_api.h"
#include "core/gc_content.h"
#include "core/grep_api.h"
#include "core/ma.h"
#include "core/md5_fingerprint.h"
#include "core/parseutils.h"
#include "core/range.h"
#include "core/sig.h"
#include "core/str.h"
#include "core/undef.h"
#include "core/unused_api.h"
#include "core/xansi.h"
#include "core/xposix.h"

typedef struct {
  GtStrArray *md5_fingerprints;
} GtBioseqFingerprints;

struct GtBioseq {
  bool use_stdin;
  GtStr *sequence_file;
  GtSeq **seqs;
  GtArray *descriptions,
        *sequence_ranges;
  char *raw_sequence;
  size_t raw_sequence_length,
         allocated;
  GtAlphabet *alphabet;
  GtBioseqFingerprints *fingerprints;
};

static bool read_fingerprints(GtStrArray *md5_fingerprints,
                              GtStr  *fingerprints_filename,
                              unsigned long num_of_seqs)
{
  bool reading_succeeded = true;
  FILE *fingerprint_file = NULL;
  gt_assert(md5_fingerprints && fingerprints_filename);
  /* open file */
  if (gt_file_exists(gt_str_get(fingerprints_filename)))
    fingerprint_file = gt_fa_xfopen(gt_str_get(fingerprints_filename), "r");
  else
    reading_succeeded = false;
  /* reading file (each line contains a single MD5 sum) */
  if (reading_succeeded) {
    GtStr *line = gt_str_new();
    while (gt_str_read_next_line(line, fingerprint_file) != EOF) {
      gt_str_array_add(md5_fingerprints, line);
      gt_str_reset(line);
    }
    gt_str_delete(line);
    if (gt_str_array_size(md5_fingerprints) < num_of_seqs) {
      /* premature end of file (e.g., due to aborted construction) */
      reading_succeeded = false;
      gt_str_array_set_size(md5_fingerprints, 0);
    }
    else
      gt_assert(gt_str_array_size(md5_fingerprints) == num_of_seqs);
  }
  gt_fa_xfclose(fingerprint_file);
  return reading_succeeded;
}

static void add_fingerprints(GtStrArray *md5_fingerprints, GtBioseq *bs)
{
  unsigned long i;
  gt_assert(md5_fingerprints && bs);
  for (i = 0; i < gt_bioseq_number_of_sequences(bs); i++) {
    char *md5 = gt_md5_fingerprint(gt_bioseq_get_sequence(bs, i),
                                   gt_bioseq_get_sequence_length(bs, i));
    gt_str_array_add_cstr(md5_fingerprints, md5);
    gt_free(md5);
  }
}

static void strarray_dump_to_file(GtStrArray *sa, FILE *outfp)
{
  unsigned long i;
  gt_assert(sa && outfp);
  for (i = 0; i < gt_str_array_size(sa); i++) {
    gt_xfputs(gt_str_array_get(sa, i), outfp);
    gt_xfputc('\n', outfp);
  }
}

static void write_fingerprints(GtStrArray *md5_fingerprints,
                               GtStr *fingerprints_filename)
{
  FILE *fingerprints_file;
  gt_assert(md5_fingerprints && fingerprints_filename);
  fingerprints_file = gt_fa_xfopen(gt_str_get(fingerprints_filename), "w");
  strarray_dump_to_file(md5_fingerprints, fingerprints_file);
  gt_fa_xfclose(fingerprints_file);
}

static GtBioseqFingerprints* gt_bioseq_fingerprints_new(GtBioseq *bs)
{
  GtBioseqFingerprints *bsf;
  bool reading_succeeded = false;
  GtStr *fingerprints_filename;
  gt_assert(bs);
  bsf = gt_calloc(1, sizeof *bsf);
  bsf->md5_fingerprints = gt_str_array_new();
  fingerprints_filename = gt_str_clone(bs->sequence_file);
  gt_str_append_cstr(fingerprints_filename, GT_BIOSEQ_FINGERPRINTS);
  if (!bs->use_stdin && gt_file_exists(gt_str_get(fingerprints_filename)) &&
      !gt_file_is_newer(gt_str_get(bs->sequence_file),
                        gt_str_get(fingerprints_filename))) {
    /* only try to read the fingerprint file if the sequence file was not
       modified in the meantime */
    reading_succeeded = read_fingerprints(bsf->md5_fingerprints,
                                          fingerprints_filename,
                                          gt_bioseq_number_of_sequences(bs));
  }
  if (!reading_succeeded) {
    add_fingerprints(bsf->md5_fingerprints, bs);
    if (!bs->use_stdin)
      write_fingerprints(bsf->md5_fingerprints, fingerprints_filename);
  }
  gt_str_delete(fingerprints_filename);
  return bsf;
}

static void gt_bioseq_fingerprints_delete(GtBioseqFingerprints *bsf)
{
  if (!bsf) return;
  gt_str_array_delete(bsf->md5_fingerprints);
  gt_free(bsf);
}

static const char* gt_bioseq_fingerprints_get(GtBioseqFingerprints *bsf,
                                           unsigned long idx)
{
  gt_assert(bsf);
  return gt_str_array_get(bsf->md5_fingerprints, idx);
}

typedef struct {
  FILE *gt_bioseq_index,
       *gt_bioseq_raw;
  unsigned long offset;
  GtBioseq *bs;
} ConstructBioseqFilesInfo;

static int proc_description(const char *description, unsigned long length,
                            void *data, GT_UNUSED GtError *err)
{
  ConstructBioseqFilesInfo *info = (ConstructBioseqFilesInfo*) data;
  char *description_cstr;
  gt_error_check(err);
  if (info->bs->use_stdin) {
    description_cstr = gt_cstr_dup(description);
    gt_array_add(info->bs->descriptions, description_cstr);
  }
  else {
    if (length)
      gt_xfputs(description, info->gt_bioseq_index);
    gt_xfputc('\n', info->gt_bioseq_index);
  }
  return 0;
}

static int proc_sequence_part(const char *seqpart, unsigned long length,
                              void *data, GT_UNUSED GtError *err)
{
  ConstructBioseqFilesInfo *info = (ConstructBioseqFilesInfo*) data;
  gt_error_check(err);
  gt_assert(seqpart);
  if (info->bs->use_stdin) {
    info->bs->raw_sequence = gt_dynalloc(info->bs->raw_sequence,
                                         &info->bs->allocated,
                                         info->bs->raw_sequence_length +
                                         length);
    memcpy(info->bs->raw_sequence + info->bs->raw_sequence_length, seqpart,
           length);
    info->bs->raw_sequence_length += length;
  }
  else
    gt_xfputs(seqpart, info->gt_bioseq_raw);
  return 0;
}

static int proc_sequence_length(unsigned long sequence_length, void *data,
                                GT_UNUSED GtError *err)
{
  ConstructBioseqFilesInfo *info = (ConstructBioseqFilesInfo*) data;
  GtRange range;
  gt_error_check(err);
  if (info->bs->use_stdin) {
    range.start = info->offset;
    range.end = info->offset + sequence_length - 1;
    gt_array_add(info->bs->sequence_ranges, range);
  }
  else {
    fprintf(info->gt_bioseq_index, "%lu\n", info->offset);
    gt_assert(sequence_length);
    fprintf(info->gt_bioseq_index, "%lu\n", info->offset + sequence_length - 1);
  }
  info->offset += sequence_length;
  return 0;
}

/* this global variables are necessary for the signal handler below */
static ConstructBioseqFilesInfo gt_bioseq_files_info;
static const char *gt_bioseq_index_filename,
                  *gt_bioseq_raw_filename;

/* removes the incomplete bioseq files */
static void remove_bioseq_files(int sigraised)
{
  /* we don't care if fclose() succeeds, xunlink() will take care of it */
  (void) fclose(gt_bioseq_files_info.gt_bioseq_index);
  (void) fclose(gt_bioseq_files_info.gt_bioseq_raw);
  gt_xunlink(gt_bioseq_index_filename);
  gt_xunlink(gt_bioseq_raw_filename);
  (void) gt_xsignal(sigraised, SIG_DFL);
  gt_xraise(sigraised);
}

static int fill_bioseq(GtBioseq *bs, const char *index_filename,
                       const char *raw_filename, GtError *err)
{
  FILE *index_file;
  GtStr *index_line;
  unsigned long line_number = 1;
  char *description;
  GtRange range;
  int had_err = 0;

  gt_error_check(err);

  /* parse the index file and fill the sequence to index mapping */
  index_line = gt_str_new();
  index_file = gt_fa_xfopen(index_filename, "r");

  while (!had_err && gt_str_read_next_line(index_line, index_file) != EOF) {
    switch (line_number % 3) {
      case 1:
        /* process description */
        description = gt_cstr_dup(gt_str_get(index_line));
        gt_array_add(bs->descriptions, description);
        break;
      case 2:
        /* process sequence start */
        if (gt_parse_ulong(&range.start, gt_str_get(index_line))) {
          gt_error_set(err, "could not parse bioseq start in line %lu of file "
                         "\"%s\"", line_number, index_filename);
          had_err = -1;
        }
        break;
      case 0:
        /* process sequence end */
        if (gt_parse_ulong(&range.end, gt_str_get(index_line))) {
          gt_error_set(err,
                    "could not parse bioseq end in line %lu of file \"%s\"",
                    line_number, index_filename);
          had_err = -1;
        }
        else {
          gt_assert(range.start <= range.end); /* XXX */
          gt_array_add(bs->sequence_ranges, range);
        }
        break;
    }
    line_number++;
    gt_str_reset(index_line);
  }

  if (!had_err) {
    /* the number of descriptions equals the number of sequence ranges */
    gt_assert(gt_array_size(bs->descriptions) ==
           gt_array_size(bs->sequence_ranges));
    /* map the raw file */
    bs->raw_sequence = gt_fa_xmmap_read(raw_filename, &bs->raw_sequence_length);
  }

  gt_fa_xfclose(index_file);
  gt_str_delete(index_line);

  return had_err;
}

static int construct_bioseq_files(GtBioseq *bs, GtStr *gt_bioseq_index_file,
                                  GtStr *gt_bioseq_raw_file,
                                  GtFastaReaderType gt_fasta_reader_type,
                                  GtError *err)
{
  GtFastaReader *fasta_reader = NULL;
  GtStr *sequence_filename;
  int had_err;

  gt_error_check(err);

  /* open files & init */
  if (!bs->use_stdin) {
    gt_bioseq_files_info.gt_bioseq_index =
      gt_fa_xfopen((const char *) gt_str_get(gt_bioseq_index_file), "w");
    gt_bioseq_files_info.gt_bioseq_raw =
      gt_fa_xfopen(gt_str_get(gt_bioseq_raw_file), "w");
  }
  gt_bioseq_files_info.offset = 0;
  gt_bioseq_files_info.bs = bs;

  /* register the signal handler to remove incomplete files upon termination */
  if (!bs->use_stdin) {
    gt_bioseq_index_filename = gt_str_get(gt_bioseq_index_file);
    gt_bioseq_raw_filename = gt_str_get(gt_bioseq_raw_file);
    gt_sig_register_all(remove_bioseq_files);
  }

  /* read fasta file */
  sequence_filename = bs->use_stdin ? NULL : bs->sequence_file;
  switch (gt_fasta_reader_type) {
    case GT_FASTA_READER_REC:
      fasta_reader = gt_fasta_reader_rec_new(sequence_filename);
      break;
    case GT_FASTA_READER_FSM:
      fasta_reader = gt_fasta_reader_fsm_new(sequence_filename);
      break;
    case GT_FASTA_READER_SEQIT:
      fasta_reader = gt_fasta_reader_seqit_new(sequence_filename);
      break;
    default: gt_assert(0);
  }
  had_err = gt_fasta_reader_run(fasta_reader, proc_description,
                                proc_sequence_part, proc_sequence_length,
                                &gt_bioseq_files_info, err);
  gt_fasta_reader_delete(fasta_reader);

  /* unregister the signal handler */
  if (!bs->use_stdin)
    gt_sig_unregister_all();

  /* close files */
  if (!bs->use_stdin) {
    gt_fa_xfclose(gt_bioseq_files_info.gt_bioseq_index);
    gt_fa_xfclose(gt_bioseq_files_info.gt_bioseq_raw);
    if (had_err) {
      gt_xunlink(gt_bioseq_index_filename);
      gt_xunlink(gt_bioseq_raw_filename);
    }
  }

  return had_err;
}

static int gt_bioseq_fill(GtBioseq *bs, bool recreate,
                       GtFastaReaderType gt_fasta_reader_type, GtError *err)
{
  GtStr *gt_bioseq_index_file = NULL,
         *gt_bioseq_raw_file = NULL;
  int had_err = 0;

  gt_assert(!bs->raw_sequence);

  /* construct file names */
  if (!bs->use_stdin) {
    gt_bioseq_index_file = gt_str_clone(bs->sequence_file);
    gt_str_append_cstr(gt_bioseq_index_file, GT_BIOSEQ_INDEX);
    gt_bioseq_raw_file = gt_str_clone(bs->sequence_file);
    gt_str_append_cstr(gt_bioseq_raw_file, GT_BIOSEQ_RAW);
  }

  /* construct the bioseq files if necessary */
  if (recreate || bs->use_stdin ||
      !gt_file_exists(gt_str_get(gt_bioseq_index_file)) ||
      !gt_file_exists(gt_str_get(gt_bioseq_raw_file)) ||
      gt_file_is_newer(gt_str_get(bs->sequence_file),
                       gt_str_get(gt_bioseq_index_file)) ||
      gt_file_is_newer(gt_str_get(bs->sequence_file),
                       gt_str_get(gt_bioseq_raw_file))) {
    had_err = construct_bioseq_files(bs, gt_bioseq_index_file,
                                     gt_bioseq_raw_file, gt_fasta_reader_type,
                                     err);
  }

  if (!had_err && !bs->use_stdin) {
    /* fill the bioseq */
    had_err = fill_bioseq(bs, gt_str_get(gt_bioseq_index_file),
                          gt_str_get(gt_bioseq_raw_file), err);
  }

  /* free */
  gt_str_delete(gt_bioseq_index_file);
  gt_str_delete(gt_bioseq_raw_file);

  return had_err;
}

static GtBioseq* gt_bioseq_new_with_recreate_and_type(GtStr *sequence_file,
                                                 bool recreate,
                                                 GtFastaReaderType
                                                 gt_fasta_reader_type,
                                                 GtError *err)
{
  GtBioseq *bs;
  int had_err = 0;
  gt_error_check(err);
  bs = gt_calloc(1, sizeof (GtBioseq));
  if (!strcmp(gt_str_get(sequence_file), "-"))
    bs->use_stdin = true;
  if (!bs->use_stdin && !gt_file_exists(gt_str_get(sequence_file))) {
    gt_error_set(err, "sequence file \"%s\" does not exist or is not readable",
              gt_str_get(sequence_file));
    had_err = -1;
  }
  if (!had_err) {
    bs->sequence_file = gt_str_ref(sequence_file);
    bs->descriptions = gt_array_new(sizeof (char*));
    bs->sequence_ranges = gt_array_new(sizeof (GtRange));
    had_err = gt_bioseq_fill(bs, recreate, gt_fasta_reader_type, err);
  }
  if (had_err) {
    gt_bioseq_delete(bs);
    return NULL;
  }
  return bs;
}

GtBioseq* gt_bioseq_new(const char *sequence_file, GtError *err)
{
  GtBioseq *bs;
  GtStr *seqfile;
  gt_error_check(err);
  seqfile = gt_str_new_cstr(sequence_file);
  bs = gt_bioseq_new_with_recreate_and_type(seqfile, false, GT_FASTA_READER_REC,
                                            err);
  gt_str_delete(seqfile);
  return bs;
}

GtBioseq* gt_bioseq_new_recreate(const char *sequence_file, GtError *err)
{
  GtBioseq *bs;
  GtStr *seqfile;
  gt_error_check(err);
  seqfile = gt_str_new_cstr(sequence_file);
  bs = gt_bioseq_new_with_recreate_and_type(seqfile, true, GT_FASTA_READER_REC,
                                            err);
  gt_str_delete(seqfile);
  return bs;
}

GtBioseq* gt_bioseq_new_str(GtStr *sequence_file, GtError *err)
{
  return gt_bioseq_new_with_recreate_and_type(sequence_file, false,
                                           GT_FASTA_READER_REC, err);
}

GtBioseq* gt_bioseq_new_with_fasta_reader(const char *sequence_file,
                                           GtFastaReaderType fasta_reader,
                                           GtError *err)
{
  GtBioseq *bs;
  GtStr *seqfile;
  gt_error_check(err);
  seqfile = gt_str_new_cstr(sequence_file);
  bs = gt_bioseq_new_with_recreate_and_type(seqfile, true, fasta_reader, err);
  gt_str_delete(seqfile);
  return bs;
}

void gt_bioseq_delete(GtBioseq *bs)
{
  unsigned long i;
  if (!bs) return;
  gt_bioseq_fingerprints_delete(bs->fingerprints);
  gt_str_delete(bs->sequence_file);
  if (bs->seqs) {
    for (i = 0; i < gt_array_size(bs->descriptions); i++)
      gt_seq_delete(bs->seqs[i]);
    gt_free(bs->seqs);
  }
  for (i = 0; i < gt_array_size(bs->descriptions); i++)
    gt_free(*(char**) gt_array_get(bs->descriptions, i));
  gt_array_delete(bs->descriptions);
  gt_array_delete(bs->sequence_ranges);
  if (bs->use_stdin)
    gt_free(bs->raw_sequence);
  else
    gt_fa_xmunmap(bs->raw_sequence);
  gt_alphabet_delete(bs->alphabet);
  gt_free(bs);
}

static void determine_alphabet_if_necessary(GtBioseq *bs)
{
  gt_assert(bs);
  if (!bs->alphabet) {
    bs->alphabet = gt_alphabet_guess(gt_bioseq_get_raw_sequence(bs),
                                     gt_bioseq_get_raw_sequence_length(bs));
  }
}

GtAlphabet* gt_bioseq_get_alphabet(GtBioseq *bs)
{
  gt_assert(bs);
  determine_alphabet_if_necessary(bs);
  gt_assert(bs->alphabet);
  return bs->alphabet;
}

GtSeq* gt_bioseq_get_seq(GtBioseq *bs, unsigned long idx)
{
  gt_assert(bs);
  gt_assert(idx < gt_array_size(bs->descriptions));
  if (!bs->seqs)
    bs->seqs = gt_calloc(gt_array_size(bs->descriptions), sizeof (GtSeq*));
  determine_alphabet_if_necessary(bs);
  if (!bs->seqs[idx]) {
    bs->seqs[idx] = gt_seq_new(gt_bioseq_get_sequence(bs, idx),
                               gt_bioseq_get_sequence_length(bs, idx),
                               bs->alphabet);
    gt_seq_set_description(bs->seqs[idx], gt_bioseq_get_description(bs, idx));
  }
  return bs->seqs[idx];
}

const char* gt_bioseq_get_description(GtBioseq *bs, unsigned long idx)
{
  gt_assert(bs);
  return *(char**) gt_array_get(bs->descriptions, idx);
}

const char* gt_bioseq_get_sequence(GtBioseq *bs, unsigned long idx)
{
  GtRange sequence_range;
  gt_assert(bs);
  sequence_range = *(GtRange*) gt_array_get(bs->sequence_ranges, idx);
  return bs->raw_sequence + sequence_range.start;
}

const char* gt_bioseq_get_raw_sequence(GtBioseq *bs)
{
  gt_assert(bs);
  return bs->raw_sequence;
}

const char* gt_bioseq_get_md5_fingerprint(GtBioseq *bs, unsigned long idx)
{
  gt_assert(bs && idx < gt_bioseq_number_of_sequences(bs));
  if (!bs->fingerprints)
    bs->fingerprints = gt_bioseq_fingerprints_new(bs);
  gt_assert(gt_bioseq_fingerprints_get(bs->fingerprints, idx));
  return gt_bioseq_fingerprints_get(bs->fingerprints, idx);
}

unsigned long gt_bioseq_get_sequence_length(GtBioseq *bs, unsigned long idx)
{
  gt_assert(bs);
  return gt_range_length(gt_array_get(bs->sequence_ranges, idx));
}

unsigned long gt_bioseq_get_raw_sequence_length(GtBioseq *bs)
{
  gt_assert(bs);
  return bs->raw_sequence_length;
}

unsigned long gt_bioseq_number_of_sequences(GtBioseq *bs)
{
  gt_assert(bs);
  return gt_array_size(bs->descriptions);
}

void gt_bioseq_show_as_fasta(GtBioseq *bs, unsigned long width)
{
  unsigned long i;

  gt_assert(bs);

  for (i = 0; i < gt_bioseq_number_of_sequences(bs); i++) {
    gt_fasta_show_entry(gt_bioseq_get_description(bs, i),
                        gt_bioseq_get_sequence(bs, i),
                        gt_bioseq_get_sequence_length(bs, i), width);
  }
}

void gt_bioseq_show_sequence_as_fasta(GtBioseq *bs, unsigned long seqnum,
                                   unsigned long width)
{
  gt_assert(bs);
  gt_assert(seqnum < gt_bioseq_number_of_sequences(bs));

  gt_fasta_show_entry(gt_bioseq_get_description(bs, seqnum),
                      gt_bioseq_get_sequence(bs, seqnum),
                      gt_bioseq_get_sequence_length(bs, seqnum), width);

}

void gt_bioseq_show_gc_content(GtBioseq *bs)
{
  gt_assert(bs);
  determine_alphabet_if_necessary(bs);
  if (gt_alphabet_is_dna(bs->alphabet)) {
    printf("showing GC-content for sequence file \"%s\"\n",
           gt_str_get(bs->sequence_file));
    gt_gc_content_show(gt_bioseq_get_raw_sequence(bs),
                       gt_bioseq_get_raw_sequence_length(bs), bs->alphabet);
  }
}

void gt_bioseq_show_stat(GtBioseq *bs)
{
  unsigned long i, num_of_seqs;
  gt_assert(bs);
  num_of_seqs = gt_bioseq_number_of_sequences(bs);
  printf("showing statistics for sequence file \"%s\"\n",
         gt_str_get(bs->sequence_file));
  printf("number of sequences: %lu\n", num_of_seqs);
  printf("total length: %lu\n", gt_bioseq_get_raw_sequence_length(bs));
  for (i = 0; i < num_of_seqs; i++) {
    printf("sequence #%lu length: %lu\n", i+1,
           gt_bioseq_get_sequence_length(bs, i));
  }
}

void gt_bioseq_show_seqlengthdistri(GtBioseq *bs)
{
  GtDiscDistri *d;
  unsigned long i;
  gt_assert(bs);
  d = gt_disc_distri_new();
  for (i = 0; i < gt_bioseq_number_of_sequences(bs); i++)
    gt_disc_distri_add(d, gt_bioseq_get_sequence_length(bs, i));
  printf("sequence length distribution:\n");
  gt_disc_distri_show(d);
  gt_disc_distri_delete(d);
}
