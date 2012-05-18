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
#include "core/logger.h"
#include "core/ma.h"
#include "core/unused_api.h"
#include "core/showtime.h"
#include "core/thread.h"
#include "match/randomsamples.h"
#include "tools/gt_encseq2kmers.h"

typedef struct {
  unsigned int kmersize,
               samplingfactor,
               numofparts;
  unsigned long maximumspace;
  GtStr  *encseqinput,
         *memlimitarg;
  GtOption *refoptionmemlimit;
  bool verbose, quiet;
} GtEncseq2kmersArguments;

static void* gt_encseq2kmers_arguments_new(void)
{
  GtEncseq2kmersArguments *arguments = gt_calloc(1, sizeof *arguments);
  arguments->encseqinput= gt_str_new();
  arguments->numofparts = 0;
  arguments->memlimitarg = gt_str_new();
  arguments->maximumspace = 0UL; /* in bytes */
  return arguments;
}

static void gt_encseq2kmers_arguments_delete(void *tool_arguments)
{
  GtEncseq2kmersArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_str_delete(arguments->encseqinput);
  gt_option_delete(arguments->refoptionmemlimit);
  gt_str_delete(arguments->memlimitarg);
  gt_free(arguments);
}

static GtOptionParser* gt_encseq2kmers_option_parser_new(void *tool_arguments)
{
  GtEncseq2kmersArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option, *optionparts, *optionmemlimit, *q_option, *v_option;
  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("-ii <indexname> -k <kmersize> [option ...]",
                         "Collect and process kmers of an encseq.");

  /* -parts */
  optionparts = gt_option_new_uint("parts", "specify the number of parts",
                                  &arguments->numofparts, 0U);
  gt_option_parser_add_option(op, optionparts);

  /* -memlimit */
  optionmemlimit = gt_option_new_string("memlimit",
                       "specify maximal amount of memory to be used during "
                       "index construction (in bytes, the keywords 'MB' "
                       "and 'GB' are allowed)",
                       arguments->memlimitarg, NULL);
  gt_option_parser_add_option(op, optionmemlimit);
  gt_option_exclude(optionmemlimit, optionparts);
  arguments->refoptionmemlimit = gt_option_ref(optionmemlimit);

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
  v_option = gt_option_new_verbose(&arguments->verbose);
  gt_option_parser_add_option(op, v_option);

  /* -q */
  q_option = gt_option_new_bool("q",
      "suppress standard output messages",
      &arguments->quiet, false);
  gt_option_exclude(q_option, v_option);
  gt_option_parser_add_option(op, q_option);

  gt_option_parser_set_max_args(op, 0);

  return op;
}

static int gt_encseq2kmers_arguments_check(GT_UNUSED int rest_argc,
                                       void *tool_arguments,
                                       GT_UNUSED GtError *err)
{
  GT_UNUSED GtEncseq2kmersArguments *arguments = tool_arguments;
  bool haserr = false;
  if (gt_option_is_set(arguments->refoptionmemlimit))
  {
    if (gt_option_parse_spacespec(&arguments->maximumspace,
                                  "memlimit",
                                  arguments->memlimitarg,
                                  err) != 0)
    {
      haserr = true;
    }
    if (!haserr && !gt_ma_bookkeeping_enabled()) {
      gt_error_set(err, "option '-memlimit' requires "
                        "GT_MEM_BOOKKEEPING=on");
      haserr = true;
    }
  }
#ifdef GT_THREADS_ENABLED
  if (!haserr) {
    if (gt_jobs > 1 && gt_ma_bookkeeping_enabled()) {
      gt_error_set(err, "gt option '-j' and GT_MEM_BOOKKEEPING=on "
                        "are incompatible");
      haserr = true;
    }
  }
#endif
  return haserr ? -1 : 0;
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
  GtLogger *default_logger, *verbose_logger;
  bool haserr = false;

  gt_error_check(err);
  gt_assert(arguments);

  default_logger =
    gt_logger_new(!arguments->quiet, GT_LOGGER_DEFLT_PREFIX, stdout);
  verbose_logger =
    gt_logger_new(arguments->verbose, GT_LOGGER_DEFLT_PREFIX, stdout);

  gt_logger_log(default_logger, "gt encseq2kmers");
  gt_log_log("indexname = %s", gt_str_get(arguments->encseqinput));
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
    if ((codes = gt_randomsamples_new(encseq, timer, arguments->kmersize,
            default_logger, verbose_logger, err))
          == NULL)
    {
      haserr = true;
    }
    else
    {
      if (gt_randomsamples_run(codes, arguments->samplingfactor,
            arguments->numofparts, arguments->maximumspace) != 0)
      {
        haserr = true;
      }
      gt_randomsamples_delete(codes);
    }
  }
  gt_encseq_delete(encseq);
  gt_encseq_loader_delete(el);
  if (timer != NULL)
  {
    gt_timer_show_progress_final(timer, stdout);
    gt_timer_delete(timer);
  }

  gt_logger_delete(verbose_logger);
  gt_logger_delete(default_logger);
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
