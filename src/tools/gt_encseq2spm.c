/*
  Copyright (c) 2011 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
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
#include "core/encseq_api.h"
#include "core/logger.h"
#include "core/intbits.h"
#include "core/minmax.h"
#include "core/showtime.h"
#ifdef GT_THREADS_ENABLED
#include "core/thread_api.h"
#endif
#include "tools/gt_encseq2spm.h"
#include "match/firstcodes.h"
#include "match/firstcodes-scan.h"
#include "match/esa-spmsk.h"

typedef struct
{
  bool checksuftab,
       singlestrand,
       verbose,
       outputspms,
       onlyaccum,
       onlyallfirstcodes,
       radixlarge,
       countspms;
  unsigned int minmatchlength,
               numofparts,
               radixparts,
               singlescan,
               addbscache_depth,
               forcek;
  unsigned long maximumspace,
                phase2extra;
  GtStr *encseqinput,
        *spmspec,
        *memlimitarg,
        *phase2extraarg;
  GtOption *refoptionmemlimit,
           *refoptionphase2extra;
} GtEncseq2spmArguments;

static void* gt_encseq2spm_arguments_new(void)
{
  GtEncseq2spmArguments *arguments = gt_calloc((size_t) 1, sizeof *arguments);
  arguments->outputspms = false;
  arguments->countspms = false;
  arguments->radixlarge = false;
  arguments->singlescan = 0;
  arguments->numofparts = 0;
  arguments->radixparts = 1U;
  arguments->encseqinput = gt_str_new();
  arguments->spmspec = gt_str_new();
  arguments->memlimitarg = gt_str_new();
  arguments->phase2extraarg = gt_str_new();
  arguments->phase2extra = 0UL; /* in bytes */
  arguments->maximumspace = 0UL; /* in bytes */
  return arguments;
}

static void gt_encseq2spm_arguments_delete(void *tool_arguments)
{
  GtEncseq2spmArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_str_delete(arguments->encseqinput);
  gt_str_delete(arguments->spmspec);
  gt_option_delete(arguments->refoptionmemlimit);
  gt_option_delete(arguments->refoptionphase2extra);
  gt_str_delete(arguments->memlimitarg);
  gt_str_delete(arguments->phase2extraarg);
  gt_free(arguments);
}

static GtOptionParser* gt_encseq2spm_option_parser_new(void *tool_arguments)
{
  GtEncseq2spmArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option, *optionparts, *optionmemlimit;

  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("[option ...] [file]",
                            "Compute suffix prefix matches "
                            "from encoded sequence.");

  /* -l */
  option = gt_option_new_uint_min("l", "specify the minimum length",
                                  &arguments->minmatchlength, 0, 1U);
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

  /* -checksuftab */
  option = gt_option_new_bool("checksuftab", "check the suffix table",
                             &arguments->checksuftab, false);
  gt_option_parser_add_option(op, option);
  gt_option_is_development_option(option);

  /* -singlestrand */
  option = gt_option_new_bool("singlestrand", "use only the forward strand "
                              "of the sequence",
                              &arguments->singlestrand, false);
  gt_option_parser_add_option(op, option);

  /* -spm */
  option = gt_option_new_string("spm", "specify output for spms",
                                arguments->spmspec, NULL);
  gt_option_parser_add_option(op, option);

  /* -ii */
  option = gt_option_new_string("ii", "specify the input sequence",
                                arguments->encseqinput, NULL);
  gt_option_parser_add_option(op, option);
  gt_option_is_mandatory(option);

  /* -onlyaccum */
  option = gt_option_new_bool("onlyaccum", "only accumulate codes",
                             &arguments->onlyaccum, false);
  gt_option_parser_add_option(op, option);
  gt_option_is_development_option(option);

  /* -onlyallfirstcodes */
  option = gt_option_new_bool("onlyallfirstcodes", "only determines allcodes",
                              &arguments->onlyallfirstcodes, false);
  gt_option_parser_add_option(op, option);
  gt_option_is_development_option(option);

  /* -addbscachedepth */
  option = gt_option_new_uint("addbscachedepth", "only determines allcodes",
                              &arguments->addbscache_depth, 5U);
  gt_option_parser_add_option(op, option);
  gt_option_is_development_option(option);

  /* -phase2extra */
  option = gt_option_new_string("phase2extra",
                       "specify  amount of additional space required for "
                       "the second phase of the computation involving the "
                       "processing of the intervals (in bytes, "
                       "the keywords 'MB' and 'GB' are allowed)",
                       arguments->phase2extraarg, NULL);
  gt_option_parser_add_option(op, option);
  arguments->refoptionphase2extra = gt_option_ref(option);
  gt_option_is_development_option(option);

  /* -radixlarge */
  option = gt_option_new_bool("radixlarge", "use large tables for radixsort",
                              &arguments->radixlarge, false);
  gt_option_parser_add_option(op, option);
  gt_option_is_development_option(option);

  /* -radixparts */
  option = gt_option_new_uint("radixparts", "specify the number of parts "
                              "for radixsort",
                              &arguments->radixparts, 1U);
  gt_option_parser_add_option(op, option);
  gt_option_is_development_option(option);

  /* -singlescan */
  option = gt_option_new_uint("singlescan", "run a single scan: 1=fast; "
                              "2=fast with check; 3=fast with output; "
                              "4=sfx-mapped4-version",
                              &arguments->singlescan, 0);
  gt_option_parser_add_option(op, option);
  gt_option_is_development_option(option);

  /* -forcek */
  option = gt_option_new_uint("forcek", "specify the value of k",
                              &arguments->forcek, 0);
  gt_option_is_development_option(option);
  gt_option_parser_add_option(op, option);

  option = gt_option_new_verbose(&arguments->verbose);
  gt_option_parser_add_option(op, option);

  return op;
}

static bool gt_encseq2spm_kmersize(GtEncseq2spmArguments *arguments,
    unsigned int *kmersize, GtError *err)
{
  bool haserr = false;
  gt_assert(kmersize != NULL);
  if (arguments->forcek > 0)
  {
    *kmersize = arguments->forcek;
    if (*kmersize > arguments->minmatchlength)
    {
      gt_error_set(err,"argument %u to option -forcek > l",
          *kmersize);
      haserr = true;
    }
    else if (*kmersize > (unsigned int)GT_UNITSIN2BITENC)
    {
      gt_error_set(err,
          "argument %u to option -forcek > %u (machine word size/2)",
          *kmersize, (unsigned int)GT_UNITSIN2BITENC);
      haserr = true;
    }
  }
  else
  {
    *kmersize = MIN((unsigned int) GT_UNITSIN2BITENC,
        arguments->minmatchlength);
  }
  gt_log_log("kmersize=%u", *kmersize);
  gt_assert(*kmersize > 0);
  return haserr;
}

static int gt_encseq2spm_arguments_check(int rest_argc,
                                         void *tool_arguments,
                                         GtError *err)
{
  GtEncseq2spmArguments *arguments = tool_arguments;
  bool haserr = false;

  gt_error_check(err);
  if (rest_argc != 0)
  {
    gt_error_set(err,"unnecessary arguments");
    haserr = true;
  }
  if (!haserr && gt_str_length(arguments->spmspec) > 0)
  {
    const char *spmspecstring = gt_str_get(arguments->spmspec);
    if (strcmp(spmspecstring,"show") == 0)
    {
      arguments->outputspms = true;
    } else
    {
      if (strcmp(spmspecstring,"count") == 0)
      {
        arguments->countspms = true;
      } else
      {
        gt_error_set(err,"illegal argument \"%s\" to option -spm",
                     spmspecstring);
        haserr = true;
      }
    }
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

static int gt_encseq2spm_runner(GT_UNUSED int argc,
                                GT_UNUSED const char **argv,
                                GT_UNUSED int parsed_args,
                                void *tool_arguments,
                                GtError *err)
{
  GtEncseq2spmArguments *arguments = tool_arguments;
  GtEncseqLoader *el = NULL;
  GtEncseq *encseq = NULL;
  bool haserr = false;

  gt_error_check(err);
  gt_assert(arguments);
  el = gt_encseq_loader_new();
  gt_encseq_loader_drop_description_support(el);
  gt_encseq_loader_disable_autosupport(el);
  encseq = gt_encseq_loader_load(el, gt_str_get(arguments->encseqinput),
                                 err);
  if (encseq == NULL)
  {
    haserr = true;
  }
  if (!haserr)
  {
    if (arguments->singlestrand)
    {
      gt_error_set(err,"option -singlestand is not implemented");
      haserr = true;
    } else
    {
      if (gt_encseq_mirror(encseq, err) != 0)
      {
        haserr = true;
      }
    }
  }

  if (!haserr && arguments->singlescan > 0)
  {
    GtTimer *timer = NULL;

    if (gt_showtime_enabled())
    {
      char *outmsg;

      switch (arguments->singlescan)
      {
        case 1:
          outmsg = "to run fast scanning";
          break;
        case 2:
          outmsg = "to run fast scanning with check";
          break;
        case 3:
          outmsg = "to run fast scanning with output";
          break;
        case 4:
          outmsg = "to run old scanning code";
          break;
        default:
          gt_error_set(err,"argument %u to option -singlescan not allowed",
                       arguments->singlescan);
          haserr = true;
      }
      if (!haserr)
      {
        timer = gt_timer_new_with_progress_description(outmsg);
        gt_timer_start(timer);
      }
    }
    if (!haserr)
    {
      unsigned int kmersize = 0;
      haserr = gt_encseq2spm_kmersize(arguments, &kmersize, err);
      if (!haserr)
      {
        if (arguments->singlescan == 4U)
        {
          gt_rungetencseqkmers(encseq,kmersize);
        } else
        {
          if (arguments->singlescan > 0)
          {
            gt_firstcode_runkmerscan(encseq,arguments->singlescan - 1,kmersize,
                arguments->minmatchlength);
          }
        }
      }
    }
    if (timer != NULL)
    {
      gt_timer_show_progress_final(timer, stdout);
      gt_timer_delete(timer);
    }
  }
  if (!haserr && arguments->singlescan == 0)
  {
    GtLogger *logger;
    const GtReadmode readmode = GT_READMODE_FORWARD;
    GtBUstate_spmsk **spmsk_states = NULL;
    unsigned int kmersize, threadcount;

#ifdef GT_THREADS_ENABLED
    const unsigned int threads = gt_jobs;
#else
    const unsigned int threads = 1U;
#endif

    if (arguments->countspms || arguments->outputspms)
    {
      spmsk_states = gt_malloc(sizeof (*spmsk_states) * threads);
      for (threadcount = 0; threadcount < threads; threadcount++)
      {
        spmsk_states[threadcount]
          = gt_spmsk_inl_new(encseq,
                             readmode,
                             (unsigned long) arguments->minmatchlength,
                             arguments->countspms,
                             arguments->outputspms,
                             gt_str_get(arguments->encseqinput));
      }
    }
    logger = gt_logger_new(arguments->verbose,GT_LOGGER_DEFLT_PREFIX, stdout);
    haserr = gt_encseq2spm_kmersize(arguments, &kmersize, err);
    if (!haserr)
    {
      if (storefirstcodes_getencseqkmers_twobitencoding(encseq,
                                                      kmersize,
                                                      arguments->numofparts,
                                                      arguments->maximumspace,
                                                      arguments->minmatchlength,
                                     /* use false */  arguments->checksuftab,
                                     /* use false */  arguments->onlyaccum,
                                     /* use false */  arguments->
                                                            onlyallfirstcodes,
                                     /* use 5U */     arguments->
                                                            addbscache_depth,
                                     /* specify the extra space needed for
                                        the function processing the interval */
                                                      arguments->phase2extra,
                                     /* use true */   arguments->radixlarge ?
                                                        false : true,
                                     /* use 2 without threads and
                                        use 1 with threads */
                                                      arguments->radixparts,
                                                      spmsk_states != NULL
                                                        ? gt_spmsk_inl_process
                                                        : NULL,
                                                      gt_spmsk_inl_process_end,
                                                      spmsk_states,
                                                      logger,
                                                      err) != 0)
      {
        haserr = true;
      }
    }
    if (spmsk_states != NULL)
    {
      unsigned long countmatches = 0;

      for (threadcount = 0; threadcount < threads; threadcount++)
      {
        countmatches += gt_spmsk_inl_delete(spmsk_states[threadcount]);
      }
      if (arguments->countspms)
      {
        printf("number of suffix-prefix matches=%lu\n",countmatches);
      }
      gt_free(spmsk_states);
    }
    gt_logger_delete(logger);
  }
  gt_encseq_delete(encseq);
  gt_encseq_loader_delete(el);
  return haserr ? -1 : 0;
}

GtTool* gt_encseq2spm(void)
{
  return gt_tool_new(gt_encseq2spm_arguments_new,
                     gt_encseq2spm_arguments_delete,
                     gt_encseq2spm_option_parser_new,
                     gt_encseq2spm_arguments_check,
                     gt_encseq2spm_runner);
}
