/*
  Copyright (c) 2011 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2011 Giorgio Gonnella <gonnella@zbh.uni-hamburg.de>
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
#include "core/ma.h"
#include "core/unused_api.h"
#include "core/option_api.h"
#include "core/encseq.h"
#include "core/intbits.h"
#include "core/logger.h"
#include "core/minmax.h"
#ifdef GT_THREADS_ENABLED
#include "core/thread_api.h"
#endif
#include "match/rdj-spmfind.h"
#include "match/rdj-version.h"
#include "match/firstcodes.h"
#include "tools/gt_readjoiner_overlap.h"

typedef struct {
  bool singlestrand,
       elimtrans,
       verbose,
       quiet,
       showspm;
  unsigned int minmatchlength,
               numofparts,
               w_maxsize;
  unsigned long maximumspace,
                phase2extra;
  GtStr *encseqinput,
        *memlimitarg,
        *phase2extraarg;
  GtOption *refoptionmemlimit,
           *refoptionphase2extra;
  bool radixsmall;
  unsigned int radixparts;
  bool onlyallfirstcodes;
} GtReadjoinerOverlapArguments;

static void* gt_readjoiner_overlap_arguments_new(void)
{
  GtReadjoinerOverlapArguments *arguments = gt_calloc((size_t) 1,
      sizeof *arguments);
  arguments->encseqinput = gt_str_new();
  arguments->numofparts = 0;
  arguments->memlimitarg = gt_str_new();
  arguments->phase2extraarg = gt_str_new();
  arguments->phase2extra = 0UL; /* in bytes */
  arguments->maximumspace = 0UL; /* in bytes */
  return arguments;
}

static void gt_readjoiner_overlap_arguments_delete(void *tool_arguments)
{
  GtReadjoinerOverlapArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_str_delete(arguments->encseqinput);
  gt_option_delete(arguments->refoptionmemlimit);
  gt_option_delete(arguments->refoptionphase2extra);
  gt_str_delete(arguments->memlimitarg);
  gt_str_delete(arguments->phase2extraarg);
  gt_free(arguments);
}

static GtOptionParser* gt_readjoiner_overlap_option_parser_new(
    void *tool_arguments)
{
  GtReadjoinerOverlapArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option, *optionparts, *optionmemlimit, *q_option, *v_option;

  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("[option ...] [file]",
                            "Compute suffix prefix matches "
                            "from encoded sequence.");

  /* -readset */
  option = gt_option_new_string("readset", "specify the readset name",
                                arguments->encseqinput, NULL);
  gt_option_parser_add_option(op, option);
  gt_option_is_mandatory(option);

  /* -l */
  option = gt_option_new_uint_min("l", "specify the minimum SPM length",
                                  &arguments->minmatchlength, 0, 2U);
  gt_option_parser_add_option(op, option);
  gt_option_is_mandatory(option);

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

  /* -singlestrand */
  option = gt_option_new_bool("singlestrand", "do not use reverse complements "
      "of the reads", &arguments->singlestrand, false);
  gt_option_is_development_option(option);
  gt_option_parser_add_option(op, option);

  /* -elimtrans */
  option = gt_option_new_bool("elimtrans", "output only irreducible SPMs",
      &arguments->elimtrans, true);
  gt_option_is_development_option(option);
  gt_option_parser_add_option(op, option);

  /* -wmax */
  option = gt_option_new_uint("wmax", "specify the maximum width of w set;\n"
                              "use 0 to disable w set partitioning",
                              &arguments->w_maxsize, 32U);
  gt_option_is_development_option(option);
  gt_option_parser_add_option(op, option);

  /* -showspm */
  option = gt_option_new_bool("showspm", "show textual SPMs list on stdout",
      &arguments->showspm, false);
  gt_option_is_development_option(option);
  gt_option_parser_add_option(op, option);

  /* -phase2extra */
  option = gt_option_new_string("phase2extra",
                       "specify amount of additional space required for "
                       "the second phase of the computation involving the "
                       "processing of the intervals (in bytes, "
                       "the keywords 'MB' and 'GB' are allowed)",
                       arguments->phase2extraarg, NULL);
  gt_option_parser_add_option(op, option);
  arguments->refoptionphase2extra = gt_option_ref(option);
  gt_option_is_development_option(option);

  /* -v */
  v_option = gt_option_new_verbose(&arguments->verbose);
  gt_option_parser_add_option(op, v_option);

  /* -q */
  q_option = gt_option_new_bool("q", "suppress standard output messages",
                             &arguments->quiet, false);
  gt_option_parser_add_option(op, q_option);
  gt_option_exclude(q_option, v_option);

  /* -radixparts */
  option = gt_option_new_uint("radixparts", "specify the radixpart parameter",
      &arguments->radixparts, 1U);
  gt_option_is_development_option(option);
  gt_option_parser_add_option(op, option);

  /* -radixsmall */
  option = gt_option_new_bool("radixsmall", "specify the radixsmall parameter",
      &arguments->radixsmall, false);
  gt_option_is_development_option(option);
  gt_option_parser_add_option(op, option);

  /* -onlyallfirstcodes */
  option = gt_option_new_bool("onlyallfirstcodes", "only determines allcodes",
                              &arguments->onlyallfirstcodes, false);
  gt_option_is_development_option(option);
  gt_option_parser_add_option(op, option);

  gt_option_parser_set_version_func(op, gt_readjoiner_show_version);
  return op;
}

static int gt_readjoiner_overlap_arguments_check(int rest_argc,
                                         void *tool_arguments,
                                         GtError *err)
{
  GtReadjoinerOverlapArguments *arguments = tool_arguments;
  bool haserr = false;

  gt_error_check(err);
  if (rest_argc != 0)
  {
    gt_error_set(err,"unnecessary arguments");
    haserr = true;
  }
  if (!haserr && gt_option_is_set(arguments->refoptionmemlimit))
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
  if (!haserr && gt_option_is_set(arguments->refoptionphase2extra))
  {
    if (gt_option_parse_spacespec(&arguments->phase2extra,
          "phase2extra",
          arguments->phase2extraarg,
          err) != 0)
    {
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

static int gt_readjoiner_overlap_runner(GT_UNUSED int argc,
                                GT_UNUSED const char **argv,
                                GT_UNUSED int parsed_args,
                                void *tool_arguments,
                                GtError *err)
{
  GtReadjoinerOverlapArguments *arguments = tool_arguments;
  GtEncseqLoader *el = NULL;
  GtEncseq *encseq = NULL;
  GtLogger *default_logger, *verbose_logger;
  unsigned int kmersize;
  bool haserr = false;
  bool eqlen;
  unsigned long total_nof_irr_spm = 0, total_nof_trans_spm = 0;

  gt_error_check(err);
  gt_assert(arguments);

  default_logger = gt_logger_new(!arguments->quiet, GT_LOGGER_DEFLT_PREFIX,
      stdout);
  gt_logger_log(default_logger,
      "gt readjoiner overlap (version "GT_READJOINER_VERSION")");
  verbose_logger = gt_logger_new(arguments->verbose, GT_LOGGER_DEFLT_PREFIX,
      stdout);
  gt_logger_log(verbose_logger, "verbose output activated");
  kmersize = MIN((unsigned int) GT_UNITSIN2BITENC, arguments->minmatchlength);

  el = gt_encseq_loader_new();
  gt_encseq_loader_drop_description_support(el);
  gt_encseq_loader_disable_autosupport(el);
  encseq = gt_encseq_loader_load(el, gt_str_get(arguments->encseqinput),
                                 err);
  eqlen = gt_encseq_accesstype_get(encseq) == GT_ACCESS_TYPE_EQUALLENGTH;

  if (encseq == NULL)
    haserr = true;
  if (!haserr && !arguments->singlestrand)
  {
    if (gt_encseq_mirror(encseq, err) != 0)
      haserr = true;
  }
  if (!haserr)
  {
    unsigned int threadcount;
#ifdef GT_THREADS_ENABLED
    const unsigned int threads = gt_jobs;
#else
    const unsigned int threads = 1U;
#endif
    if (eqlen)
    {
      GtBUstate_spmeq **state_table
        = gt_malloc(sizeof (*state_table) * threads);

      for (threadcount = 0; threadcount < threads; threadcount++)
      {
        state_table[threadcount]
          = gt_spmfind_eqlen_state_new(encseq,
                (unsigned long)arguments->minmatchlength,
                (unsigned long)arguments->w_maxsize, arguments->elimtrans,
                arguments->showspm, gt_str_get(arguments->encseqinput),
                threadcount, default_logger, verbose_logger, err);
      }
      if (storefirstcodes_getencseqkmers_twobitencoding(encseq, kmersize,
            arguments->numofparts, arguments->maximumspace,
            arguments->minmatchlength, false, false,
            arguments->onlyallfirstcodes, 5U, arguments->phase2extra,
            arguments->radixsmall, arguments->radixparts,
            gt_spmfind_eqlen_process, gt_spmfind_eqlen_process_end,
            state_table, verbose_logger, err)
            != 0)
      {
        haserr = true;
      }
      for (threadcount = 0; threadcount < threads; threadcount++)
      {
        total_nof_irr_spm +=
          gt_spmfind_eqlen_nof_irr_spm(state_table[threadcount]);
        total_nof_trans_spm +=
          gt_spmfind_eqlen_nof_irr_spm(state_table[threadcount]);
        gt_spmfind_eqlen_state_delete(state_table[threadcount]);
      }
      gt_free(state_table);
    }
    else
    {
      GtBUstate_spmvar **state_table
        = gt_malloc(sizeof (*state_table) * threads);

      for (threadcount = 0; threadcount < threads; threadcount++)
      {
        state_table[threadcount]
           = gt_spmfind_varlen_state_new(encseq,
                  (unsigned long)arguments->minmatchlength,
                  (unsigned long)arguments->w_maxsize, arguments->elimtrans,
                  arguments->showspm, gt_str_get(arguments->encseqinput),
                  threadcount, default_logger, verbose_logger, err);
      }
      if (storefirstcodes_getencseqkmers_twobitencoding(encseq, kmersize,
            arguments->numofparts, arguments->maximumspace,
            arguments->minmatchlength, false, false,
            arguments->onlyallfirstcodes, 5U, arguments->phase2extra,
            arguments->radixsmall, arguments->radixparts,
            gt_spmfind_varlen_process, gt_spmfind_varlen_process_end,
            state_table, verbose_logger, err)
          != 0)
      {
        haserr = true;
      }
      for (threadcount = 0; threadcount < threads; threadcount++)
      {
        gt_spmfind_varlen_state_delete(state_table[threadcount]);
        total_nof_irr_spm +=
          gt_spmfind_eqlen_nof_irr_spm(state_table[threadcount]);
        total_nof_trans_spm +=
          gt_spmfind_eqlen_nof_irr_spm(state_table[threadcount]);
      }
      gt_free(state_table);
    }
  }
  if (!haserr)
  {
    gt_logger_log(default_logger, "number of %ssuffix-prefix matches = %lu",
        arguments->elimtrans ? "irreducible " : "", total_nof_irr_spm);
    if (arguments->elimtrans)
      gt_logger_log(default_logger, "number of transitive suffix-prefix "
          "matches = %lu", total_nof_trans_spm);
  }
  gt_logger_delete(default_logger);
  gt_logger_delete(verbose_logger);
  gt_encseq_delete(encseq);
  gt_encseq_loader_delete(el);
  return haserr ? -1 : 0;
}

GtTool* gt_readjoiner_overlap(void)
{
  return gt_tool_new(gt_readjoiner_overlap_arguments_new,
                     gt_readjoiner_overlap_arguments_delete,
                     gt_readjoiner_overlap_option_parser_new,
                     gt_readjoiner_overlap_arguments_check,
                     gt_readjoiner_overlap_runner);
}
