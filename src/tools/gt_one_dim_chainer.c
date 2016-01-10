/*
  Copyright (c) 2015 Fabian Sobanski <0sobansk@informatik.uni-hamburg.de>
  Copyright (c) 2015 Christopher Keil <christopher.keil@studium.uni-hamburg.de>
  Copyright (c) 2015 Center for Bioinformatics, University of Hamburg

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

#include <math.h>

#include "core/ma.h"
#include "core/unused_api.h"
#include "extended/one_dim_chainer_match.h"
#include "tools/gt_one_dim_chainer.h"

typedef struct {
  GtUword overlap;
} GtOneDimChainerArguments;

/* Allocates memory and returns a pointer for <arguments> needed for GtOptions.
*/
static void* gt_one_dim_chainer_arguments_new(void)
{
  GtOneDimChainerArguments *arguments =
    gt_calloc((size_t) 1, sizeof *arguments);

  return arguments;
}

/* Frees memory which has been allocated for <arguments> needed for GtOptions.
*/
static void gt_one_dim_chainer_arguments_delete(void *tool_arguments)
{
  GtOneDimChainerArguments *arguments = tool_arguments;
  if (arguments != NULL) {
    gt_free(arguments);
  }
}

/* Uses the <tool_arguments> to create a GtOptionParser which is used to process
passed options for running the Gt_one_dim_chainer tool.  */
static GtOptionParser* gt_one_dim_chainer_option_parser_new(
    void *tool_arguments)
{
  GtOneDimChainerArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option;
  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("[option ...] [file]", /* XXX */
                            "DESCRIBE YOUR TOOL IN ONE LINE HERE."); /* XXX */

  /* -overlap */
  option = gt_option_new_uword("overlap", "number of bases that are allowed"
                               "to overlap in a chain",
                               &arguments->overlap, 0);
  gt_option_parser_add_option(op, option);

  return op;
}

/* Checks the <tool_arguments> for errors and handles errors by passing a
message to the GtError object <err>. Returns -1 if errors have been detected. */
static int gt_one_dim_chainer_arguments_check(GT_UNUSED int rest_argc,
                                       void *tool_arguments,
                                       GT_UNUSED GtError *err)
{
  GtOneDimChainerArguments *arguments = tool_arguments;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(arguments);

  return had_err;
}

/* Executes the chaining algorithm and prints the results.
The command line arguments <argv> should contain a match file name from which
matches are passed to a GtSeedextendMatchIterator.
Options are passed via <tool_arguments>. If errors occur, the GtError
object <err> is updated with a message and the function returns -1. */
static int gt_one_dim_chainer_runner(int argc, const char **argv,
                              GT_UNUSED int parsed_args, void *tool_arguments,
                              GtError *err)
{
  GtOneDimChainerArguments *arguments = tool_arguments;
  int had_err = 0;

  gt_error_check(err);
  gt_assert(arguments);

  /* Read in match file */
  GtStr *matchfilename = gt_str_new_cstr(argv[argc-1]);
  GtOneDimChainerMatch *match,
                       *chainstart = gt_malloc(sizeof *chainstart);

  /* Set up match iterator */
  GtSeedextendMatchIterator *semi
       = gt_seedextend_match_iterator_new(matchfilename, err);
  gt_str_delete(matchfilename);
  if (semi == NULL)
  {
    return -1;
  }
  gt_seedextend_match_iterator_querymatchoutoptions_set(
      semi, true, false, 0, false, false);
  (void) gt_seedextend_match_iterator_all_sorted(semi, false);

  had_err = gt_1d_chainer_calc_chain(arguments->overlap,
                                     chainstart, semi, err);

  /* Print the chain of matches */
  for (match = chainstart; match != NULL; match = match->suc)
  {
    GtQuerymatch *querymatch = match->querymatch;
    printf("%" PRIu64 "\t" GT_WU "\t" GT_WU "\n",
           gt_querymatch_queryseqnum(querymatch),
           gt_querymatch_querystart(querymatch),
           gt_querymatch_querystart(querymatch) +
           gt_querymatch_querylen(querymatch));
  }

  /* Delete the chain of matches and the iterator from memory */
  gt_1d_chainer_match_delete(chainstart);
  gt_seedextend_match_iterator_delete(semi);
  return had_err;
}

/* A wrapper which creates and returns  new GtTool from passed arguments,
options and the defined runner method. */
GtTool* gt_one_dim_chainer(void)
{
  return gt_tool_new(gt_one_dim_chainer_arguments_new,
                     gt_one_dim_chainer_arguments_delete,
                     gt_one_dim_chainer_option_parser_new,
                     gt_one_dim_chainer_arguments_check,
                     gt_one_dim_chainer_runner);
}
