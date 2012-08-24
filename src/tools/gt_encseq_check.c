/*
  Copyright (c) 2011 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
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

#include "core/encseq.h"
#include "core/ma.h"
#include "core/unused_api.h"
#include "core/logger.h"
#include "match/intcode-def.h"  /* XXX */
#include "match/test-mappedstr.pr"  /* XXX */
#include "tools/gt_encseq_check.h"

typedef struct {
  unsigned long scantrials,
                multicharcmptrials,
                prefixlength;
  bool verbose,
       nocheckunit;
} GtEncseqCheckArguments;

static void* gt_encseq_check_arguments_new(void)
{
  GtEncseqCheckArguments *arguments = gt_calloc(1, sizeof *arguments);
  return arguments;
}

static void gt_encseq_check_arguments_delete(void *tool_arguments)
{
  GtEncseqCheckArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_free(arguments);
}

static GtOptionParser* gt_encseq_check_option_parser_new(void *tool_arguments)
{
  GtEncseqCheckArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GT_UNUSED GtOption *option;
  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("[option ...] [file]",
                            "Check the consistency of an encoded "
                            "sequence file.");

  option = gt_option_new_ulong("scantrials", "specify number of scan trials",
                               &arguments->scantrials, 0);
  gt_option_parser_add_option(op, option);

  option = gt_option_new_ulong("multicharcmptrials",
                               "specify number of multicharacter trials",
                               &arguments->multicharcmptrials, 0);
  gt_option_parser_add_option(op, option);

  option = gt_option_new_ulong_min_max("prefixlength",
                                       "prefix length",
                                       &arguments->prefixlength, 0,
                                       0, MAXPREFIXLENGTH);
  gt_option_parser_add_option(op, option);

  option = gt_option_new_bool("nocheckunit","do not run checkextractunitatpos",
                              &arguments->nocheckunit, false);
  gt_option_parser_add_option(op, option);

  option = gt_option_new_verbose(&arguments->verbose);
  gt_option_parser_add_option(op, option);

  gt_option_parser_set_min_args(op, 1);

  return op;
}

static int gt_encseq_check_runner(GT_UNUSED int argc, const char **argv,
                                  int parsed_args, void *tool_arguments,
                                  GtError *err)
{
  GtEncseqCheckArguments *arguments = tool_arguments;
  int had_err = 0;
  GtEncseqLoader *encseq_loader;
  GtEncseq *encseq;
  GtLogger *logger = NULL;

  gt_error_check(err);
  gt_assert(arguments);

  logger = gt_logger_new(arguments->verbose, GT_LOGGER_DEFLT_PREFIX, stdout);
  encseq_loader = gt_encseq_loader_new();
  if (!(encseq = gt_encseq_loader_load(encseq_loader, argv[parsed_args], err)))
    had_err = -1;
  if (!had_err) {
    int readmode;

    gt_encseq_check_startpositions(encseq,logger);
    for (readmode = 0; readmode < 4; readmode++)
    {
      if (gt_alphabet_is_dna(gt_encseq_alphabet(encseq)) ||
           ((GtReadmode) readmode) == GT_READMODE_FORWARD ||
           ((GtReadmode) readmode) == GT_READMODE_REVERSE)
      {
        gt_logger_log(logger,"check consistency for readmode %s",
                      gt_readmode_show((GtReadmode) readmode));
        if (gt_encseq_check_consistency(encseq,
                           gt_encseq_filenames(encseq),
                           (GtReadmode) readmode,
                           arguments->scantrials,
                           arguments->multicharcmptrials,
                           gt_encseq_has_multiseq_support(encseq),
                           !arguments->nocheckunit,
                           logger,
                           err) != 0)
        {
          had_err = -1;
          break;
        }
      }
    }
    if (!had_err)
    {
      gt_encseq_check_specialranges(encseq);
    }
    if (!had_err)
    {
      gt_encseq_check_markpos(encseq);
    }
    if (!had_err)
    {
      had_err = gt_encseq_check_minmax(encseq, err);
    }
    if (!had_err &&
           arguments->prefixlength > 0)
    {
      if (gt_verifymappedstr(encseq,
                             arguments->prefixlength,
                             err) != 0)
      {
        had_err = -1;
      }
    }
  }
  gt_encseq_delete(encseq);
  gt_encseq_loader_delete(encseq_loader);
  gt_logger_delete(logger);
  return had_err;
}

GtTool* gt_encseq_check(void)
{
  return gt_tool_new(gt_encseq_check_arguments_new,
                     gt_encseq_check_arguments_delete,
                     gt_encseq_check_option_parser_new,
                     NULL,
                     gt_encseq_check_runner);
}
