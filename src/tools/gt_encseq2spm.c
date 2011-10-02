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
#include "tools/gt_encseq2spm.h"
#include "match/firstcodes.h"
#include "match/esa-spmsk.h"

typedef struct
{
  bool checksuftab,
       mirrored,
       verbose,
       outputspms,
       onlyaccum,
       forceoverflow,
       countspms;
  unsigned int minmatchlength,
               numofparts;
  unsigned long maximumspace;
  GtStr *encseqinput,
        *spmspec,
        *memlimitarg;
  GtOption *refoptionmemlimit;
} GtEncseq2spmArguments;

static void* gt_encseq2spm_arguments_new(void)
{
  GtEncseq2spmArguments *arguments = gt_calloc((size_t) 1, sizeof *arguments);
  arguments->outputspms = false;
  arguments->countspms = false;
  arguments->numofparts = 0;
  arguments->encseqinput = gt_str_new();
  arguments->spmspec = gt_str_new();
  arguments->memlimitarg = gt_str_new();
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
  gt_str_delete(arguments->memlimitarg);
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

  /* -mirrored */
  option = gt_option_new_bool("mirrored", "use sequence with its mirror",
                             &arguments->mirrored, false);
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

  /* -foreceoverflow */
  option = gt_option_new_bool("forceoverflow", "force overflow of countocc",
                             &arguments->forceoverflow, false);
  gt_option_parser_add_option(op, option);
  gt_option_is_development_option(option);

  option = gt_option_new_verbose(&arguments->verbose);
  gt_option_parser_add_option(op, option);

  return op;
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
    if (gt_option_parse_memlimit(&arguments->maximumspace,
                                 arguments->memlimitarg,
                                 err) != 0)
    {
      haserr = true;
    }
  }
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
  if (!haserr && arguments->mirrored)
  {
    if (gt_encseq_mirror(encseq, err) != 0)
    {
      haserr = true;
    }
  }
  if (!haserr)
  {
    GtLogger *logger;
    const GtReadmode readmode = GT_READMODE_FORWARD;
    GtBUstate_spmsk *spmsk_state = NULL;

    if (arguments->countspms || arguments->outputspms)
    {
      spmsk_state = gt_spmsk_inl_new(encseq,
                                     readmode,
                                     (unsigned long) arguments->minmatchlength,
                                     arguments->countspms,
                                     arguments->outputspms,
                                     gt_str_get(arguments->encseqinput));
    }
    logger = gt_logger_new(arguments->verbose,GT_LOGGER_DEFLT_PREFIX, stdout);
    if (storefirstcodes_getencseqkmers_twobitencoding(encseq,
                                                      32U,
                                                      arguments->numofparts,
                                                      arguments->maximumspace,
                                                      arguments->minmatchlength,
                                     /* use false */  arguments->checksuftab,
                                     /* use false */  arguments->onlyaccum,
                                     /* use false */  arguments->forceoverflow,
                                     /* specify the extra space needed for
                                        the function processing the interval */
                                                      0,
                                                      spmsk_state != NULL
                                                        ? gt_spmsk_inl_process
                                                        : NULL,
                                                      spmsk_state,
                                                      logger,err) != 0)
    {
      haserr = true;
    }
    gt_spmsk_inl_delete(spmsk_state);
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
