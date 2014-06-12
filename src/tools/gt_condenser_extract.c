/*
  Copyright (c) 2013 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
  Copyright (c) 2013 Center for Bioinformatics, University of Hamburg

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

#include "core/logger.h"
#include "core/ma.h"
#include "core/str_api.h"
#include "core/types_api.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "core/xansi_api.h"
#include "extended/n_r_encseq.h"
#include "tools/gt_condenser_extract.h"

typedef struct {
  GtRange  range;
  GtStr   *original;
  bool     verbose;
} GtCondenserExtractArguments;

static void* gt_condenser_extract_arguments_new(void)
{
  GtCondenserExtractArguments *arguments = gt_calloc((size_t) 1,
                                                     sizeof *arguments);
  arguments->range.start =
    arguments->range.end = GT_UNDEF_UWORD;
  arguments->original = gt_str_new();
  return arguments;
}

static void gt_condenser_extract_arguments_delete(void *tool_arguments)
{
  GtCondenserExtractArguments *arguments = tool_arguments;
  if (arguments != NULL) {
    gt_str_delete(arguments->original);
    gt_free(arguments);
  }
}

static GtOptionParser*
gt_condenser_extract_option_parser_new(void *tool_arguments)
{
  GtCondenserExtractArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option;
  gt_assert(arguments);

  /* init */
  /*TODO soll das nur zu fasta sein oder kann man auch raw sequence
    extracten?*/
  op = gt_option_parser_new("[option ...] [archive]",
                            "Decompresses a condenser archive to fasta.");

  /* -original */
  option = gt_option_new_filename("original",
                                  "uncompressed encseq, needs to be present "
                                  "for development reasons.",
                                  arguments->original);
  gt_option_is_mandatory(option);
  gt_option_parser_add_option(op, option);

  /* -range */
  option = gt_option_new_range("range",
                               "Range of positions to extract"
                               ". If no "
                               "range is given, whole sequence "
                               "collection is extracted.",
                               &arguments->range, NULL);

  gt_option_parser_add_option(op, option);

  /* -verbose */
  option = gt_option_new_bool("verbose", "Print out verbose output to stderr.",
                              &arguments->verbose, false);
  gt_option_parser_add_option(op, option);

  return op;
}

static int gt_condenser_extract_arguments_check(int rest_argc,
                                           void *tool_arguments,
                                           GtError *err)
{
  GtCondenserExtractArguments *arguments = tool_arguments;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(arguments);

  if (rest_argc != 1) {
    had_err = -1;
    gt_error_set(err, "no file to extract from given use -help for usage");
  }
  else if (arguments->range.start != GT_UNDEF_UWORD) {
    if (arguments->range.start > arguments->range.end) {
      had_err = -1;
      gt_error_set(err, "give range of positions in the form 'a b' with a "
                   "<= b");
    }
  }
  return had_err;
}

static int gt_condenser_extract_runner(GT_UNUSED int argc,
                                  const char **argv,
                                  int parsed_args,
                                  void *tool_arguments,
                                  GtError *err)
{
  int had_err = 0;
  GtCondenserExtractArguments *arguments = tool_arguments;
  GtNREncseq *nre = NULL;
  GtEncseq *orig_encseq = NULL;
  GtEncseqLoader *esl;

  gt_error_check(err);
  gt_assert(arguments);

  /*load original encseq*/
  esl = gt_encseq_loader_new();
  orig_encseq = gt_encseq_loader_load(esl,
                                      gt_str_get(arguments->original),
                                      err);
  if (!orig_encseq) {
    had_err = -1;
  }
  gt_encseq_loader_delete(esl);
  if (!had_err) {
    nre = gt_n_r_encseq_new_from_file(argv[parsed_args], orig_encseq, err);
    if (nre == NULL) {
      had_err = -1;
    }
  }

  /*TODO get sequences by sequence ids: not yet implemented in n_r_encseq*/
  /*if (!had_err && arguments->range.start == GT_UNDEF_UWORD &&
       uedb != NULL) {
    GtUword idx,
            start = arguments->seqrange.start == GT_UNDEF_UWORD ?
                      0 :
                      arguments->seqrange.start,
            end = arguments->seqrange.end == GT_UNDEF_UWORD ?
                      uedb->nseq - 1 :
                      arguments->seqrange.end;
    for (idx = start; idx <= end && !had_err; idx++) {
      had_err = gt_unique_encseq_get_sequence_from_idx(idx, unique_encseq, uedb,
                                                       stdout, err);
    }
  }
  else if (!had_err) {*/
  if (!had_err) {
    GtNREncseqDecompressor *nred = gt_n_r_encseq_decompressor_new(nre);
    if (arguments->range.start == GT_UNDEF_ULONG &&
        arguments->range.end == GT_UNDEF_ULONG) {
      had_err = gt_n_r_encseq_decompressor_extract_origin_complete(stdout,
                                                                   nred,
                                                                   true,
                                                                   err);
    } else {
      had_err = gt_n_r_encseq_decompressor_extract_originrange(stdout,
                                                             nred,
                                                             &arguments->range,
                                                             false,
                                                             err);
    }
    gt_xfwrite_one("\n",stdout); /*TODO should better be in n_r_encseq.c?*/
    gt_n_r_encseq_decompressor_delete(nred);
  }
  gt_n_r_encseq_delete(nre);
  gt_encseq_delete(orig_encseq);
  return had_err;
}

GtTool* gt_condenser_extract(void)
{
  return gt_tool_new(gt_condenser_extract_arguments_new,
                     gt_condenser_extract_arguments_delete,
                     gt_condenser_extract_option_parser_new,
                     gt_condenser_extract_arguments_check,
                     gt_condenser_extract_runner);
}
