/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

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

#include "core/error.h"
#include "core/option_api.h"
#include "core/versionfunc.h"
#include "match/test-mtrieins.pr"
#include "tools/gt_trieins.h"

static GtOPrval parse_options(bool *onlyins,int *parsed_args,
                              int argc, const char **argv, GtError *err)
{
  GtOptionParser *op;
  GtOption *option;
  GtOPrval oprval;

  gt_error_check(err);
  op = gt_option_parser_new("[options] indexname",
                         "Perform trie insertions and check consistency.");
  gt_option_parser_set_mail_address(op, "<kurtz@zbh.uni-hamburg.de>");
  option= gt_option_new_bool("ins","perform only insertions",onlyins,false);
  gt_option_parser_add_option(op, option);
  gt_option_parser_set_min_max_args(op, 1U, 1U);
  oprval = gt_option_parser_parse(op, parsed_args, argc, argv, gt_versionfunc,
                                  err);
  gt_option_parser_delete(op);
  return oprval;
}

int gt_trieins(int argc, const char **argv, GtError *err)
{
  bool haserr = false;
  int parsed_args;
  bool onlyins = false;

  gt_error_check(err);

  switch (parse_options(&onlyins,&parsed_args, argc, argv, err)) {
    case GT_OPTION_PARSER_OK: break;
    case GT_OPTION_PARSER_ERROR: return -1;
    case GT_OPTION_PARSER_REQUESTS_EXIT: return 0;
  }
  gt_assert(parsed_args == 1);

  if (gt_test_trieins(onlyins,argv[parsed_args],err) != 0)
  {
    haserr = true;
  }
  return haserr ? -1 : 0;
}
