/*
  Copyright (c) 2011-2012 Joachim Bonnet <joachim.bonnet@studium.uni-hamburg.de>
  Copyright (c) 2012 Dirk Willrodt <willrodt@studium.uni-hamburg.de>
  Copyright (c) 2011 Center for Bioinformatics, University of Hamburg

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

#include <string.h>

#include "core/alphabet_api.h"
#include "core/basename_api.h"
#include "core/fa.h"
#include "core/fileutils_api.h"
#include "core/log_api.h"
#include "core/ma.h"
#include "core/mathsupport.h"
#include "core/range_api.h"
#include "core/showtime.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "core/xansi_api.h"
#include "extended/hcr.h"
#include "tools/gt_compreads_decompress.h"

typedef struct {
  bool descs,
       verbose;
  unsigned long bench;
  GtStr  *file,
         *smap,
         *alphabet,
         *name;
  GtRange rng;
} GtCsrHcrDecodeArguments;

static void* gt_compreads_decompress_arguments_new(void)
{
  GtCsrHcrDecodeArguments *arguments = gt_calloc((size_t) 1, sizeof *arguments);
  arguments->file = gt_str_new();
  arguments->smap = gt_str_new();
  arguments->name = gt_str_new();
  arguments->rng.start = GT_UNDEF_ULONG;
  arguments->rng.end = GT_UNDEF_ULONG;
  arguments->bench = 0;

  return arguments;
}

static void gt_compreads_decompress_arguments_delete(void *tool_arguments)
{
  GtCsrHcrDecodeArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_str_delete(arguments->file);
  gt_str_delete(arguments->name);
  gt_str_delete(arguments->smap);
  gt_free(arguments);
}

static GtOptionParser*
gt_compreads_decompress_option_parser_new(void *tool_arguments)
{
  GtCsrHcrDecodeArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option;
  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("[option ...] (-file file)",
                         "Decodes a file of compressed reads.");

  option = gt_option_new_bool("v", "be verbose",
                              &arguments->verbose, false);
  gt_option_parser_add_option(op, option);

  option = gt_option_new_bool("descs","enable description decoding",
                              &arguments->descs, false);
  gt_option_parser_add_option(op, option);

  option = gt_option_new_string("file", "specify base name of files containing"
                                " HCR.",
                                arguments->file, NULL);
  gt_option_is_mandatory(option);
  gt_option_parser_add_option(op, option);

  option = gt_option_new_string("name", "specify base name for decoded hcr"
                                " (suffix will be \".fastq\")",
                                arguments->name, NULL);
  gt_option_parser_add_option(op, option);

  option = gt_option_new_string("smap", "specify file containing alphabet"
                                "description (file must be an .al1 file)."
                                " If \"-smap\" is not set, dna alphabet is"
                                " used.",
                                arguments->smap, NULL);
  gt_option_parser_add_option(op, option);

  option = gt_option_new_range("range", "decode multiple consecutive reads."
                               " If range is not specified, the "
                               " entire file will be decoded.",
                               &arguments->rng, NULL);
  gt_option_parser_add_option(op, option);

  option = gt_option_new_ulong("benchmark", "decode given number random reads "
                               "and report the time to do this",
                               &arguments->bench, 0);
  gt_option_is_development_option(option);
  gt_option_parser_add_option(op, option);

  gt_option_parser_set_min_max_args(op, 0U, 0U);
  return op;
}

static int gt_compreads_decompress_benchmark(GtHcrDecoder *hcrd,
                                             unsigned long amount,
                                             GtTimer *timer,
                                             GtError *err) {
  char qual[BUFSIZ] = {0},
       seq[BUFSIZ] = {0};
  int had_err = 0;
  unsigned long rand,
                max_rand = gt_hcr_decoder_num_of_reads(hcrd) - 1,
                count;

  GtStr *timer_comment = gt_str_new_cstr("extracting ");
  GtStr *desc = gt_str_new();

  gt_str_append_ulong(timer_comment, amount);
  gt_str_append_cstr(timer_comment, " reads of ");
  gt_str_append_ulong(timer_comment, max_rand + 1);
  gt_str_append_cstr(timer_comment, "!");

  if (timer == NULL) {
    timer = gt_timer_new_with_progress_description("extract random reads");
    gt_timer_start(timer);
  }
  else {
    gt_timer_show_progress(timer, "extract random reads", stdout);
  }

  gt_log_log("%s",gt_str_get(timer_comment));
  for (count = 0; count < amount; count++) {
    if (!had_err) {
      rand = gt_rand_max(max_rand);
      gt_log_log("get read: %lu", rand);
      had_err = gt_hcr_decoder_decode(hcrd, rand, seq, qual, desc, err);
      gt_log_log("%s",gt_str_get(desc));
      gt_log_log("%s",seq);
      gt_log_log("%s",qual);
    }
  }
  gt_str_delete(timer_comment);
  gt_str_delete(desc);
  if (!gt_showtime_enabled())
    gt_timer_delete(timer);
  return had_err;
}

static int gt_compreads_decompress_runner(GT_UNUSED int argc,
                                    GT_UNUSED const char **argv,
                                    GT_UNUSED int parsed_args,
                                    void *tool_arguments, GtError *err)
{
  GtCsrHcrDecodeArguments *arguments = tool_arguments;
  int had_err = 0;
  GtAlphabet *alpha = NULL;
  GtHcrDecoder *hcrd = NULL;
  GtTimer *timer = NULL;
  unsigned long start,
                end;

  gt_error_check(err);
  gt_assert(arguments);

  if (gt_showtime_enabled()) {
    timer = gt_timer_new_with_progress_description("start");
    gt_timer_start(timer);
    gt_assert(timer);
  }

  if (gt_str_length(arguments->smap) > 0) {
    alpha = gt_alphabet_new_from_file_no_suffix(gt_str_get(arguments->smap),
                                                err);
    if (!alpha)
      had_err = -1;
  }
  else {
    alpha = gt_alphabet_new_dna();
    if (!alpha)
      had_err = -1;
  }

  if (!had_err) {
    if (timer != NULL)
      gt_timer_show_progress(timer, "decoding", stdout);

    if (gt_str_length(arguments->name) == 0) {
      char *basenameptr;
      basenameptr = gt_basename(gt_str_get(arguments->file));
      gt_str_set(arguments->name, basenameptr);
      gt_free(basenameptr);
    }
    hcrd = gt_hcr_decoder_new(gt_str_get(arguments->file), alpha,
                              arguments->descs, timer, err);
    if (hcrd == NULL)
      had_err = -1;
    else {
      if (arguments->bench != 0) {
        had_err = gt_compreads_decompress_benchmark(hcrd,
                                                    arguments->bench,
                                                    timer, err);
      }
      else {
        if (arguments->rng.start != GT_UNDEF_ULONG
            && arguments->rng.end != GT_UNDEF_ULONG) {
          if (arguments->rng.start >= gt_hcr_decoder_num_of_reads(hcrd)
                || arguments->rng.end >= gt_hcr_decoder_num_of_reads(hcrd)) {
            gt_error_set(err, "range %lu-%lu includes a read number exceeding "
                              "the total number of reads (%lu)",
                              arguments->rng.start,
                              arguments->rng.end,
                              gt_hcr_decoder_num_of_reads(hcrd));
            had_err = -1;
          }
          start = arguments->rng.start;
          end = arguments->rng.end;
        }
        else {
          start = 0;
          end = gt_hcr_decoder_num_of_reads(hcrd) - 1;
        }
        if (!had_err) {
          gt_log_log("filebasename: %s", gt_str_get(arguments->name));
          if (gt_hcr_decoder_decode_range(hcrd, gt_str_get(arguments->name),
                                          start, end, timer, err)
            != 0)
            had_err = -1;
        }
      }
    }
    gt_hcr_decoder_delete(hcrd);
  }
  gt_alphabet_delete(alpha);
  if (timer != NULL) {
    gt_timer_show_progress_final(timer, stdout);
    gt_timer_delete(timer);
  }
  if (had_err)
    gt_assert(gt_error_is_set(err));
  return had_err;
}

GtTool* gt_compreads_decompress(void)
{
  return gt_tool_new(gt_compreads_decompress_arguments_new,
                     gt_compreads_decompress_arguments_delete,
                     gt_compreads_decompress_option_parser_new,
                     NULL,
                     gt_compreads_decompress_runner);
}
