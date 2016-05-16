/*
  Copyright (c) 2014 Joerg Winkler <j.winkler@posteo.de>
  Copyright (c) 2014 Center for Bioinformatics, University of Hamburg

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

#include <stdlib.h>
#include "core/array_api.h"
#include "core/assert_api.h"
#include "core/bittab_api.h"
#include "core/fasta_api.h"
#include "core/fastq.h"
#include "core/ma.h"
#include "core/mathsupport.h"
#include "core/seq_iterator_api.h"
#include "core/seq_iterator_fastq_api.h"
#include "core/str_array_api.h"
#include "core/types_api.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "tools/gt_fastq_sample.h"

typedef struct {
  GtUword length;
} GtFastqSampleArguments;

static void* gt_fastq_sample_arguments_new(void)
{
  GtFastqSampleArguments *arguments = gt_calloc((size_t) 1, sizeof *arguments);
  return arguments;
}

static void gt_fastq_sample_arguments_delete(void *tool_arguments)
{
  GtFastqSampleArguments *arguments = tool_arguments;
  if (arguments != NULL) {
    gt_free(arguments);
  }
}

static GtOptionParser* gt_fastq_sample_option_parser_new(void *tool_arguments)
{
  GtFastqSampleArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *optionlength;
  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("[option ...] -length <n> <fastq_file> "
                            "[<fastq_file>...]",
                            "Print samples by random choice from given FASTQ "\
                            "files using at least n sequence-chars.\n"\
                            "Output is fastq/fasta format depending on "\
                            "whether qualities are available.");

  /* -length */
  optionlength = gt_option_new_uword("length",
                                     "minimum number of chars to be chosen",
                                     &arguments->length,
                                     GT_UNDEF_UWORD);
  gt_option_is_mandatory(optionlength);
  gt_option_parser_add_option(op, optionlength);

  gt_option_parser_set_min_args(op, 1);

  return op;
}

static int gt_fastq_sample_arguments_check(GT_UNUSED int rest_argc,
                                           void *tool_arguments,
                                           GtError *err)
{
  GtFastqSampleArguments *arguments = tool_arguments;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(arguments);

  if (arguments->length < 1) {
    gt_error_set(err, "length must be a positive integer");
    had_err = -1;
  }
  return had_err;
}

static int gt_fastq_sample_print(GtStrArray *filenames,
                                 GtFastqSampleArguments *args, GtError *err)
{
  GtSeqIterator* fastq_iterator;
  GtUword sum, num_seq, len_count, seq_count, pos;
  GtBittab *vector;
  GtArray *lengths;
  int had_err = 0;
  GtUword len = 0;
  const GtUchar *qualities = NULL;
  const GtUchar *sequence = NULL;
  char *description = NULL;

  /* iterate through seq's: sum up lengths and count number of seq's */
  sum = 0;
  num_seq = 0;
  lengths = gt_array_new(sizeof(GtUword));
  fastq_iterator = gt_seq_iterator_fastq_new(filenames, err);
  if (gt_seq_iterator_has_qualities(fastq_iterator)) {
    gt_seq_iterator_set_quality_buffer(fastq_iterator, &qualities);
  }
  while ((had_err = gt_seq_iterator_next(fastq_iterator, &sequence, &len,
        &description, err)) == 1) {
    sum += len;
    num_seq += 1;
    gt_array_add(lengths, len);
  }
  gt_seq_iterator_delete(fastq_iterator);

  /* check input data for errors */
  if (had_err) return had_err;
  if (num_seq == 0) {
    gt_array_delete(lengths);
    gt_error_set(err, "file does not contain any sequence data");
    return -1;
  }
  if (sum < args->length) {
    gt_array_delete(lengths);
    gt_error_set(err, "requested length "GT_WU" exceeds length of sequences"
                 " ("GT_WU")", args->length, sum);
    return -1;
  }

  /* fill bit vector randomly until args->length is reached */
  vector = gt_bittab_new(num_seq);
  len_count = 0;
  seq_count = 0;
  pos = (num_seq != 1) ? gt_rand_max(num_seq-1): 0; /* random start position */
  while (len_count < args->length) {
    /* probability for setting a bit = args->length / sum    */
    if (gt_rand_max(sum-1) < args->length &&
        !gt_bittab_bit_is_set(vector, pos)) {
      gt_bittab_set_bit(vector, pos);
      len_count += *((GtUword *) gt_array_get(lengths, pos));
      seq_count += 1;
    }
    pos = (pos+1) % num_seq; /* iterate through bit vector cyclically */
  }
  gt_array_delete(lengths);

  /* print fastq entries according to bitvector */
  printf("total length "GT_WU" from "GT_WU" entries\n", len_count, seq_count);
  fastq_iterator = gt_seq_iterator_fastq_new(filenames, err);
  if (gt_seq_iterator_has_qualities(fastq_iterator)) {
    gt_seq_iterator_set_quality_buffer(fastq_iterator, &qualities);
  }
  pos = 0;
  while ((had_err = gt_seq_iterator_next(fastq_iterator, &sequence, &len,
        &description, err)) == 1) {
    gt_assert(pos < num_seq);
    if (gt_bittab_bit_is_set(vector, pos)) {
      if (gt_seq_iterator_has_qualities(fastq_iterator)) {
        gt_fastq_show_entry(description, (const char *) sequence,
                            (const char *) qualities, len, 0, false, NULL);
      } else {
        gt_fasta_show_entry(description, (const char *) sequence, len, 0, NULL);
      }
    }
    pos += 1;
  }
  gt_seq_iterator_delete(fastq_iterator);
  gt_bittab_delete(vector);
  return had_err;
}

static int gt_fastq_sample_runner(int argc, const char **argv, int parsed_args,
                              void *tool_arguments, GT_UNUSED GtError *err)
{
  GtFastqSampleArguments *arguments = tool_arguments;
  int had_err = 0;
  GtStrArray* filenames;
  int i;

  gt_error_check(err);
  gt_assert(arguments);

  filenames = gt_str_array_new();
  for (i = parsed_args; i < argc; i++) {
    gt_str_array_add_cstr(filenames, argv[i]);
  }
  gt_assert(gt_str_array_size(filenames) > 0);
  had_err = gt_fastq_sample_print(filenames, arguments, err);
  gt_str_array_delete(filenames);
  return had_err;
}

GtTool* gt_fastq_sample(void)
{
  return gt_tool_new(gt_fastq_sample_arguments_new,
                     gt_fastq_sample_arguments_delete,
                     gt_fastq_sample_option_parser_new,
                     gt_fastq_sample_arguments_check,
                     gt_fastq_sample_runner);
}
