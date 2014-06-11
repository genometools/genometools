/*
  Copyright (c) 2013 Ole Eigenbrod <ole.eigenbrod@gmx.de>
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

#include "core/ma.h"
#include "core/unused_api.h"
#include "core/range_api.h"
#include "tools/gt_unique_encseq_extract.h"
#include "extended/unique_encseq.h"
#include "core/undef_api.h"
#include "core/logger.h"

typedef struct {
  GtRange range_option;
  GtRange seqrange_option;
  bool all_option,
       dbstats_coarse_option,
       dbstats_fine_option,
       dbstats_option,
       debug_logger_option;
} GtUniqueEncseqExtractArguments;

static void* gt_unique_encseq_extract_arguments_new(void)
{
  GtUniqueEncseqExtractArguments *arguments = gt_calloc((size_t) 1,
      sizeof *arguments);
  arguments->seqrange_option.start =
      arguments->seqrange_option.end = GT_UNDEF_ULONG;
  arguments->range_option.start =
      arguments->range_option.end = GT_UNDEF_ULONG;
  return arguments;
}

static void gt_unique_encseq_extract_arguments_delete(void *tool_arguments)
{
  GtUniqueEncseqExtractArguments *arguments = tool_arguments;
  if (arguments != NULL) {
    gt_free(arguments);
  }
}

static GtOptionParser* gt_unique_encseq_extract_option_parser_new(
    void *tool_arguments)
{
  GtUniqueEncseqExtractArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *all,
           *seqrange,
           *range,
           *debug,
           *dbstats,
           *dbstats_coarse,
           *dbstats_fine;
  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("[options] INPUTFILE",
                            "Decompresses a UniqueEncseq to FASTA format "
                            "or prints statistics about a given UniqueEncseq, "
                            "output is written to stdout.");

  /* -all */
  all = gt_option_new_bool("all",
      "decompression of the complete database",
       &arguments->all_option,
                              false);
  gt_option_parser_add_option(op, all);

  /* -seqrange */
  seqrange = gt_option_new_range("seqrange",
      "range of sequence indices to decompress",
      &arguments->seqrange_option,
      NULL);
  gt_option_parser_add_option(op, seqrange);

  /* -range */
  range = gt_option_new_range("range",
      "range of sequence positions to decompress",
      &arguments->range_option,
      NULL);
  gt_option_parser_add_option(op, range);

  /* -dbstats */
  dbstats = gt_option_new_bool("dbstats",
      "returns UniqueEncseq statistics, please specify further if "
      "fine (-dbstats_fine) or coarse (-dbstats_coarse) should be returned",
       &arguments->dbstats_option,
                              true);
  gt_option_parser_add_option(op, dbstats);

  /* -dbstats_coarse */
  dbstats_coarse = gt_option_new_bool("dbstats_coarse",
      "returns coarse stats of the UniqueEncseq (number of links, number of "
      "unique sequences, ...), option -dbstats is also required",
       &arguments->dbstats_coarse_option,
                              false);
  gt_option_parser_add_option(op, dbstats_coarse);

  /* -dbstats_fine */
  dbstats_fine = gt_option_new_bool("dbstats_fine",
      "returns fine stats of the UniqueEncseq (number of matches,mismatches,etc"
      "per linked sequence, to which sequence is it linked, ...), "
      "option -dbstats is also required",
       &arguments->dbstats_fine_option,
                              false);
  gt_option_parser_add_option(op, dbstats_fine);

  /* -debug */
  debug = gt_option_new_bool(
      "debug",
      "debug messages are printed",
      &arguments->debug_logger_option,
      false);
  gt_option_parser_add_option(op, debug);

 gt_option_is_mandatory_either_4(all, seqrange, range, dbstats);
 gt_option_imply_either_2(dbstats, dbstats_coarse, dbstats_fine);

  return op;
}

static int gt_unique_encseq_extract_arguments_check(int rest_argc,
                                                    void *tool_arguments,
                                                    GtError *err)
{
  GtUniqueEncseqExtractArguments GT_UNUSED *arguments = tool_arguments;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(arguments);
  if (!had_err && rest_argc != 1) {
    gt_error_set(err, "please specify exactly one input data structure!");
    had_err = -1;
  }

  return had_err;
}

static int gt_unique_encseq_extract_runner(GT_UNUSED int argc,
                                           const char **argv,
                                           int parsed_args,
                                           void *tool_arguments,
                                           GtError *err)
{
  GtUniqueEncseqExtractArguments *arguments = tool_arguments;
  GtUniqueEncseqDB *uedb;
  GtEncseqLoader *el;
  GtEncseq *unique_encseq = NULL;
  GtLogger *logger = gt_logger_new(
      arguments->debug_logger_option,
      GT_LOGGER_DEFLT_PREFIX, stdout);
  GtUword idx;
  int had_err = 0;

  gt_error_check(err);
  gt_assert(arguments);

  uedb = gt_unique_encseq_uedb_read(argv[parsed_args], err);
  if (uedb == NULL) {
    had_err = -1;
    gt_assert(gt_error_is_set(err));
  }
  if (!had_err) {
    el = gt_encseq_loader_new();
    unique_encseq = gt_encseq_loader_load(el, argv[parsed_args], err);
    if (unique_encseq == NULL)
      had_err = -1;
    gt_encseq_loader_delete(el);
  }
  if (!had_err)
    had_err = gt_unique_encseq_check_db(uedb, logger, err);

  if (!had_err && arguments->all_option) {
    for (idx = 0; idx <uedb->nseq && !had_err; idx++) {
      had_err = gt_unique_encseq_get_sequence_from_idx(idx, unique_encseq, uedb,
                                                       stdout, err);
    }
  }
  else if (!had_err && arguments->seqrange_option.start !=
      GT_UNDEF_ULONG && arguments->seqrange_option.end
      != GT_UNDEF_ULONG) {
    for (idx = arguments->seqrange_option.start;
        idx <= arguments->seqrange_option.end && !had_err;
        idx++) {
      had_err = gt_unique_encseq_get_sequence_from_idx(idx, unique_encseq, uedb,
          stdout, err);
    }
  }
  else if (!had_err && arguments->range_option.start !=
      GT_UNDEF_ULONG && arguments->range_option.end
      != GT_UNDEF_ULONG) {
    had_err = gt_unique_encseq_get_sequence_from_range(
        &arguments->range_option, unique_encseq, uedb,
        stdout, err);
  }
  else if (!had_err && arguments->dbstats_option &&
           arguments->dbstats_coarse_option) {
    gt_unique_encseq_database_stats_coarse(uedb, NULL);
  }
  else if (!had_err && arguments->dbstats_option &&
           arguments->dbstats_fine_option) {
    gt_unique_encseq_database_stats_fine(uedb, NULL);
  }

  gt_unique_encseq_delete_db(uedb);
  gt_encseq_delete(unique_encseq);
  gt_logger_delete(logger);
  return had_err;
}

GtTool* gt_unique_encseq_extract(void)
{
  return gt_tool_new(gt_unique_encseq_extract_arguments_new,
                     gt_unique_encseq_extract_arguments_delete,
                     gt_unique_encseq_extract_option_parser_new,
                     gt_unique_encseq_extract_arguments_check,
                     gt_unique_encseq_extract_runner);
}
