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

#include "core/ma.h"
#include "core/unused_api.h"
#include "extended/priority_queue.h"
#include "match/seed-extend-iter.h"
#include "match/querymatch.h"
#include "tools/gt_one_dim_chainer.h"

struct Match;   /* forward declaration */

typedef struct {
  bool bool_option_one_dim_chainer;
  GtStr  *str_option_one_dim_chainer;
} GtOneDimChainerArguments;

/* A struct that defines a match object with respective state vars */
typedef struct Match {

  /* predecessor index */
  struct Match *prec;

  uint64_t queryseqnum;

  GtUword querystart;
  GtUword queryend; 

  GtUword length;

} Match;

static void* gt_one_dim_chainer_arguments_new(void)
{
  GtOneDimChainerArguments *arguments = gt_calloc((size_t) 1, sizeof *arguments);
  arguments->str_option_one_dim_chainer = gt_str_new();
  return arguments;
}

static void gt_one_dim_chainer_arguments_delete(void *tool_arguments)
{
  GtOneDimChainerArguments *arguments = tool_arguments;
  if (arguments != NULL) {
    gt_str_delete(arguments->str_option_one_dim_chainer);
    gt_free(arguments);
  }
}

static GtOptionParser* gt_one_dim_chainer_option_parser_new(void *tool_arguments)
{
  GtOneDimChainerArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option;
  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("[option ...] [file]", /* XXX */
                            "DESCRIBE YOUR TOOL IN ONE LINE HERE."); /* XXX */

  /* -bool */
  option = gt_option_new_bool("bool", "bool option one_dim_chainer",
                              &arguments->bool_option_one_dim_chainer, false);
  gt_option_parser_add_option(op, option);

  /* -str */
  option = gt_option_new_string("str", "str option one_dim_chainer",
                                arguments->str_option_one_dim_chainer, NULL);
  gt_option_parser_add_option(op, option);

  return op;
}

static int gt_one_dim_chainer_arguments_check(GT_UNUSED int rest_argc,
                                       void *tool_arguments,
                                       GT_UNUSED GtError *err)
{
  GtOneDimChainerArguments *arguments = tool_arguments;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(arguments);

  /* XXX: do some checking after the option have been parsed (usally this is not
     necessary and this function can be removed completely). */
  if (gt_str_length(arguments->str_option_one_dim_chainer))
    printf("%s\n", gt_str_get(arguments->str_option_one_dim_chainer));

  return had_err;
}

static int gt_one_dim_chainer_runner(int argc, const char **argv, int parsed_args,
                              void *tool_arguments, GT_UNUSED GtError *err)
{
  GtOneDimChainerArguments *arguments = tool_arguments;
  int had_err = 0;

  gt_error_check(err);
  gt_assert(arguments);

  /* XXX */
  if (arguments->bool_option_one_dim_chainer)
    printf("argc=%d, parsed_args=%d\n", argc, parsed_args);
  printf("argv[0]=%s\n", argv[0]);

  GtSeedextendMatchIterator *semi
       = gt_seedextend_match_iterator_new(matchfile,err);

  if (semi == NULL)
  {
    had_err = -1;
  } else
  {
    while (true)
    {
      uint64_t queryseqnum;
      GtUword querystart, queryend, maxchainlen = 0;
      GtQuerymatch *querymatchptr = gt_seedextend_match_iterator_next(semi);
      if (querymatchptr == NULL)
      {
        break;
      }
      queryseqnum = gt_querymatch_queryseqnum(querymatchptr);
      querystart = gt_querymatch_querystart(querymatchptr);
      queryend = querystart + gt_querymatch_querylen(querymatchptr) - 1;
      
    }
  }

  return had_err;
}

GtTool* gt_one_dim_chainer(void)
{
  return gt_tool_new(gt_one_dim_chainer_arguments_new,
                     gt_one_dim_chainer_arguments_delete,
                     gt_one_dim_chainer_option_parser_new,
                     gt_one_dim_chainer_arguments_check,
                     gt_one_dim_chainer_runner);
}
