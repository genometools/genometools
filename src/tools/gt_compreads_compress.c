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
#include "core/ma.h"
#include "core/mathsupport.h"
#include "core/safearith.h"
#include "core/showtime.h"
#include "core/splitter_api.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "core/xansi_api.h"
#include "extended/hcr.h"
#include "extended/sampling.h"
#include "tools/gt_compreads_compress.h"

typedef struct {
  bool descs,
       pagewise,
       regular;
  GtStr  *smap,
         *method,
         *name;
  GtStrArray *files;
  unsigned long srate;
  GtQualRange qrng;
  GtRange arg_range;
} GtCsrHcrEncodeArguments;

static void* gt_compreads_compress_arguments_new(void)
{
  GtCsrHcrEncodeArguments *arguments = gt_calloc((size_t) 1, sizeof *arguments);
  arguments->smap = gt_str_new();
  arguments->name = gt_str_new();
  arguments->files = gt_str_array_new();
  arguments->qrng.start = GT_UNDEF_UINT;
  arguments->qrng.end = GT_UNDEF_UINT;
  arguments->arg_range.start = GT_UNDEF_ULONG;
  arguments->arg_range.end = GT_UNDEF_ULONG;
  arguments->method = gt_str_new();
  arguments->pagewise = false;
  arguments->regular = false;
  return arguments;
}

static void gt_compreads_compress_arguments_delete(void *tool_arguments)
{
  GtCsrHcrEncodeArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_str_delete(arguments->name);
  gt_str_array_delete(arguments->files);
  gt_str_delete(arguments->smap);
  gt_str_delete(arguments->method);
  gt_free(arguments);
}

static GtOptionParser*
gt_compreads_compress_option_parser_new(void *tool_arguments)
{
  GtCsrHcrEncodeArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option;
  static const char *methods[] = { "page", "regular", "none" };
  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("[option ...] (-files file [...])",
                         "Generates compact encoding for fastq data.");

  option = gt_option_new_bool("descs","encode descriptions",
                              &arguments->descs, false);
  gt_option_parser_add_option(op, option);

  option = gt_option_new_filename_array("files", "File(s) containing reads.",
                                        arguments->files);
  gt_option_parser_add_option(op, option);

  option = gt_option_new_string("name", "specify base name for HCR to be"
                                " generated. Only mandatory, if more than one"
                                " file was given.", arguments->name, NULL);
  gt_option_parser_add_option(op, option);

  option = gt_option_new_string("smap", "file containing alphabet description."
                                " If \"-smap\" is not set, dna alphabet is"
                                " used.",
                                arguments->smap, NULL);
  gt_option_parser_add_option(op, option);

  option = gt_option_new_range("qrange", "specify range of quality values."
                               " All values smaller or equal to the lower bound"
                               " will be converted to the lower bound. All"
                               " values equal or larger than the upper bound"
                               " will be converted to the upper bound.",
                               &arguments->arg_range, NULL);
  gt_option_parser_add_option(op, option);

  option = gt_option_new_ulong("srate", "sampling rate, set to sensible default"
                               " depending on sampling method",
                               &arguments->srate, GT_UNDEF_ULONG);
  gt_option_parser_add_option(op, option);

  option = gt_option_new_choice("stype", "type of sampling\n"
                                "one of regular - page - none",
                                arguments->method, "page",
                                methods);
  gt_option_parser_add_option(op, option);

  gt_option_parser_set_min_max_args(op, 0U, 0U);
  return op;
}

static int gt_compreads_compress_arguments_check(GT_UNUSED int rest_argc,
                                       void *tool_arguments,
                                       GtError *err)
{
  int had_err = 0;
  GtCsrHcrEncodeArguments *arguments = tool_arguments;
  GtSplitter *splitter = NULL;
  GtStr *buffer;
  gt_error_check(err);
  gt_assert(arguments);

  if (gt_str_length(arguments->name) == 0) {
    if (gt_str_array_size(arguments->files) > 1UL) {
      gt_error_set(err, "option \"-name\" needs to be specified if more"
                        " than one file is given");
      had_err = -1;
    }
    else {
      unsigned long i;
      char *basename;
      splitter = gt_splitter_new();
      basename = gt_basename(gt_str_array_get(arguments->files, 0));
      buffer = gt_str_new_cstr(basename);
      gt_splitter_split(splitter, gt_str_get(buffer), gt_str_length(buffer),
                        '.');
      for (i = 0; i < gt_splitter_size(splitter) - 1; i++) {
        gt_str_append_cstr(arguments->name,
                           gt_splitter_get_token(splitter, i));
        if (i < gt_splitter_size(splitter) - 2)
          gt_str_append_char(arguments->name, '.');
      }
      gt_free(basename);
      gt_splitter_delete(splitter);
      gt_str_delete(buffer);
    }
  }

  if (!had_err) {
    char *sampling_type = gt_str_get(arguments->method);
    static const char *methods[] = { "page", "regular", "none" };

    if (!strcmp(methods[0], sampling_type)) {
      arguments->pagewise = true;
      if (arguments->srate == GT_UNDEF_ULONG)
        arguments->srate = GT_SAMPLING_DEFAULT_PAGE_RATE;
      else if (arguments->srate == 0) {
        gt_error_set(err, "page sampling was chosen, but sampling rate was"
                          " set to %lu! this seems wrong.", arguments->srate);
        had_err = -1;
      }
    }
    else if (!strcmp(methods[1], sampling_type)) {
      arguments->regular = true;
      if (arguments->srate == GT_UNDEF_ULONG)
        arguments->srate = GT_SAMPLING_DEFAULT_REGULAR_RATE;
      else if (arguments->srate == 0) {
        gt_error_set(err, "regular sampling was chosen, but sampling rate "
                          " was set to %lu! this seems wrong.",
                     arguments->srate);
        had_err = -1;
      }
    }
    else if (!strcmp(methods[2], sampling_type)) {
      if (arguments->srate == GT_UNDEF_ULONG)
        arguments->srate = 0;
      else if (arguments->srate != 0) {
        gt_error_set(err, "no sampling was chosen, but sampling rate was"
                          " set to %lu! this seems wrong.", arguments->srate);
        had_err = -1;
      }
    }
    else {
      gt_error_set(err, "somethings wrong with the stype option");
      had_err = -1;
    }
  }

  if (!had_err) {
    if (arguments->arg_range.start != GT_UNDEF_ULONG) {
      if (arguments->arg_range.start <= (unsigned long) UINT_MAX) {
        gt_safe_assign(arguments->qrng.start, arguments->arg_range.start);
        if (arguments->arg_range.end <= (unsigned long) UINT_MAX)
          gt_safe_assign(arguments->qrng.end, arguments->arg_range.end);
        else
          had_err = -1;
      }
      else
        had_err = -1;
    }
    if (had_err)
      gt_error_set(err, "Range for qualities: value to large! larger than %u",
                   UINT_MAX);
  }
  return had_err;
}

static int gt_compreads_compress_runner(GT_UNUSED int argc,
                                    GT_UNUSED const char **argv,
                                    GT_UNUSED int parsed_args,
                                    void *tool_arguments, GtError *err)
{
  GtCsrHcrEncodeArguments *arguments = tool_arguments;
  int had_err = 0;
  GtAlphabet *alpha = NULL;
  GtHcrEncoder *hcre = NULL;
  GtTimer *timer = NULL;

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
      had_err = 1;
  }
  else {
    alpha = gt_alphabet_new_dna();
    if (!alpha)
      had_err = 1;
  }
  if (!had_err) {
    if (timer != NULL)
      gt_timer_show_progress(timer, "encoding", stdout);
    hcre = gt_hcr_encoder_new(arguments->files, alpha, arguments->descs,
                              arguments->qrng, timer, err);
    if (!hcre)
      had_err = 1;
    else {
      if (arguments->pagewise)
        gt_hcr_encoder_set_sampling_page(hcre);
      else if (arguments->regular)
        gt_hcr_encoder_set_sampling_regular(hcre);

      gt_hcr_encoder_set_sampling_rate(hcre, arguments->srate);

      if (gt_hcr_encoder_encode(hcre, gt_str_get(arguments->name),
                                timer, err) != 0)
        had_err = 1;
    }
    gt_hcr_encoder_delete(hcre);
  }
  gt_alphabet_delete(alpha);
  if (timer != NULL) {
    gt_timer_show_progress_final(timer, stdout);
    gt_timer_delete(timer);
  }
  return had_err;
}

GtTool* gt_compreads_compress(void)
{
  return gt_tool_new(gt_compreads_compress_arguments_new,
                  gt_compreads_compress_arguments_delete,
                  gt_compreads_compress_option_parser_new,
                  gt_compreads_compress_arguments_check,
                  gt_compreads_compress_runner);
}
