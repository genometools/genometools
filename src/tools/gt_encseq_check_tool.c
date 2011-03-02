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
#include "match/test-mappedstr.pr"  /* XXX */
#include "tools/gt_encseq_check_tool.h"

typedef struct {
  bool bool_option_encseq_check_tool;
  GtStr  *str_option_encseq_check_tool;
  unsigned long scantrials,
                multicharcmptrials;
} GtEncseqCheckToolArguments;

static void* gt_encseq_check_tool_arguments_new(void)
{
  GtEncseqCheckToolArguments *arguments = gt_calloc(1, sizeof *arguments);
  arguments->str_option_encseq_check_tool = gt_str_new();
  return arguments;
}

static void gt_encseq_check_tool_arguments_delete(void *tool_arguments)
{
  GtEncseqCheckToolArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_str_delete(arguments->str_option_encseq_check_tool);
  gt_free(arguments);
}

static GtOptionParser* gt_encseq_check_tool_option_parser_new(void
                                                                *tool_arguments)
{
  GtEncseqCheckToolArguments *arguments = tool_arguments;
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

  gt_option_parser_set_min_args(op, 1);

  return op;
}

static int gt_encseq_check_tool_runner(GT_UNUSED int argc, const char **argv,
                                       int parsed_args, void *tool_arguments,
                                       GtError *err)
{
  GtEncseqCheckToolArguments *arguments = tool_arguments;
  int had_err = 0;
  GtEncseqLoader *encseq_loader;
  GtEncseq *encseq;
  gt_error_check(err);
  gt_assert(arguments);

  encseq_loader = gt_encseq_loader_new();
  if (!(encseq = gt_encseq_loader_load(encseq_loader, argv[parsed_args], err)))
    had_err = -1;
  if (!had_err) {
    int readmode;
    for (readmode = 0; readmode < 4; readmode++)
    {
      if (gt_alphabet_is_dna(gt_encseq_alphabet(encseq)) ||
           ((GtReadmode) readmode) == GT_READMODE_FORWARD ||
           ((GtReadmode) readmode) == GT_READMODE_REVERSE)
      {
        if (gt_encseq_check_consistency(encseq,
                           gt_encseq_filenames(encseq),
                           (GtReadmode) readmode,
                           arguments->scantrials,
                           arguments->multicharcmptrials,
                           gt_encseq_has_multiseq_support(encseq),
                           err) != 0)
        {
          had_err = true;
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
    /* if (!had_err &&
        readmode == GT_READMODE_FORWARD &&
        suffixarray.prefixlength > 0)
    {
      if (gt_verifymappedstr(encseq,
                             prefixlength,
                             err) != 0)
      {
        had_err = true;
      }
    } */
  }
  gt_encseq_delete(encseq);
  gt_encseq_loader_delete(encseq_loader);
  return had_err;
}

GtTool* gt_encseq_check_tool(void)
{
  return gt_tool_new(gt_encseq_check_tool_arguments_new,
                  gt_encseq_check_tool_arguments_delete,
                  gt_encseq_check_tool_option_parser_new,
                  NULL,
                  gt_encseq_check_tool_runner);
}
