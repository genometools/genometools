/*
  Copyright (c) 2006-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2008 Center for Bioinformatics, University of Hamburg

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
#include "libgtcore/error.h"
#include "libgtcore/option.h"
#include "libgtcore/unused.h"
#include "libgtcore/versionfunc.h"
#include "libgtext/multiset_matching.h"
#include "tools/gt_multiset_matching.h"

static OPrval parse_options(int *parsed_args, int argc, const char **argv,
                            Error *err)
{
  OptionParser *op;
  OPrval oprval;
  error_check(err);
  op = option_parser_new("[option ...] multiset_string text",
                         "Match multiset defined by multiset_string against "
                         "text.");
  option_parser_set_min_max_args(op, 2, 2);
  oprval = option_parser_parse(op, parsed_args, argc, argv, versionfunc, err);
  option_parser_delete(op);
  return oprval;
}

static void show_match(unsigned long pos, UNUSED void *data)
{
  printf("%lu\n", pos + 1);
}

int gt_multiset_matching(int argc, const char **argv, Error *err)
{
  int parsed_args;
  error_check(err);

  /* option parsing */
  switch (parse_options(&parsed_args, argc, argv, err)) {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR: return -1;
    case OPTIONPARSER_REQUESTS_EXIT: return 0;
  }

  /* matching */
  assert(parsed_args + 1 < argc);
  multiset_matching((unsigned char*) argv[parsed_args],
                    strlen(argv[parsed_args]),
                    (unsigned char*) argv[parsed_args+1],
                    strlen(argv[parsed_args+1]), NULL, show_match);

  return 0;
}
