/*
  Copyright (c) 2012 Giorgio Gonnella <gonnella@zbh.uni-hamburg.de>
  Copyright (c) 2012 Center for Bioinformatics, University of Hamburg

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

#include "core/encseq_api.h"
#include "core/log_api.h"
#include "core/ma.h"
#include "core/unused_api.h"
#include "core/showtime.h"
#include "match/randomsamples.h"
#include "tools/gt_encseq2kmers.h"

typedef struct {
  unsigned int kmersize, samplingfactor;
  GtStr  *encseqinput;
  bool verbose;
} GtEncseq2kmersArguments;

static void* gt_encseq2kmers_arguments_new(void)
{
  GtEncseq2kmersArguments *arguments = gt_calloc(1, sizeof *arguments);
  arguments->encseqinput= gt_str_new();
  return arguments;
}

static void gt_encseq2kmers_arguments_delete(void *tool_arguments)
{
  GtEncseq2kmersArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_str_delete(arguments->encseqinput);
  gt_free(arguments);
}

static GtOptionParser* gt_encseq2kmers_option_parser_new(void *tool_arguments)
{
  GtEncseq2kmersArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option;
  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("-ii <indexname> -k <kmersize> [option ...]",
                         "Collect and process kmers of an encseq.");

  /* -k */
  option = gt_option_new_uint_min_max("k", "specify the kmer size",
                                  &arguments->kmersize, 31U, 1U, 32U);
  gt_option_parser_add_option(op, option);
  gt_option_is_mandatory(option);

  /* -sf */
  option = gt_option_new_uint_min("sf", "specify the sampling factor",
                                  &arguments->samplingfactor, 50U, 1U);
  gt_option_parser_add_option(op, option);
  gt_option_is_development_option(option);

  /* -ii */
  option = gt_option_new_string("ii", "specify the input sequence",
                                arguments->encseqinput, NULL);
  gt_option_parser_add_option(op, option);
  gt_option_is_mandatory(option);

  /* -v */
  option = gt_option_new_verbose(&arguments->verbose);
  gt_option_parser_add_option(op, option);

  gt_option_parser_set_max_args(op, 0);

  return op;
}

static int gt_encseq2kmers_arguments_check(GT_UNUSED int rest_argc,
                                       void *tool_arguments,
                                       GT_UNUSED GtError *err)
{
  GT_UNUSED GtEncseq2kmersArguments *arguments = tool_arguments;
  int had_err = 0;
  return had_err;
}

static int gt_encseq2kmers_runner(GT_UNUSED int argc,
    GT_UNUSED const char **argv, GT_UNUSED int parsed_args,
    void *tool_arguments, GtError *err)
{
  GtEncseq2kmersArguments *arguments = tool_arguments;
  GtEncseqLoader *el = NULL;
  GtEncseq *encseq = NULL;
  GtRandomSamples *codes;
  GtTimer *timer = NULL;
  bool haserr = false;

  gt_error_check(err);
  gt_assert(arguments);
  gt_log_log("indexname = %s", gt_str_get(arguments->encseqinput));
  gt_log_log("kmersize = %u", arguments->kmersize);
  gt_log_log("samplingfactor = %u", arguments->samplingfactor);
  el = gt_encseq_loader_new();
  gt_encseq_loader_drop_description_support(el);
  gt_encseq_loader_disable_autosupport(el);
  encseq = gt_encseq_loader_load(el, gt_str_get(arguments->encseqinput), err);
  if (encseq == NULL)
  {
    haserr = true;
  }
  if (!haserr)
  {
    if (gt_showtime_enabled())
    {
      timer = gt_timer_new_with_progress_description("for initialization");
      gt_timer_show_cpu_time_by_progress(timer);
      gt_timer_start(timer);
    }
    codes = gt_randomsamples_new(encseq, arguments->kmersize, timer);
    gt_randomsamples_sample(codes, arguments->samplingfactor);
    gt_randomsamples_delete(codes);
  }
  gt_encseq_delete(encseq);
  gt_encseq_loader_delete(el);
  if (timer != NULL)
  {
    gt_timer_show_progress_final(timer, stdout);
    gt_timer_delete(timer);
  }

  return haserr ? -1 : 0;
}

GtTool* gt_encseq2kmers(void)
{
  return gt_tool_new(gt_encseq2kmers_arguments_new,
                  gt_encseq2kmers_arguments_delete,
                  gt_encseq2kmers_option_parser_new,
                  gt_encseq2kmers_arguments_check,
                  gt_encseq2kmers_runner);
}
