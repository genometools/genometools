/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg

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
#include "libgtcore/array2dim.h"
#include "libgtcore/option.h"
#include "libgtcore/versionfunc.h"
#include "libgtext/multilcp.h"

static OPrval parse_options(int *parsed_args, int argc, const char **argv,
                            Env *env)
{
  OptionParser *op;
  OPrval oprval;
  env_error_check(env);
  op = option_parser_new("[option ...] seq1 seq2",
                         "Compute lcp lengths of seq1 and seq2 in "
                         "O(|seq1|*|seq2|) time and show them.",env);
  oprval = option_parser_parse_min_max_args(op, parsed_args, argc, argv,
                                            versionfunc, 2, 2, env);
  option_parser_delete(op, env);
  return oprval;
}

int gt_multilcp(int argc, const char **argv, Env *env)
{
  const char *seq1, *seq2;
  int parsed_args, len1, len2, **multilcptab;
  env_error_check(env);

  /* option parsing */
  switch (parse_options(&parsed_args, argc, argv, env)) {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR: return -1;
    case OPTIONPARSER_REQUESTS_EXIT: return 0;
  }

  assert(parsed_args + 1 < argc);
  seq1 = argv[parsed_args];
  seq2 = argv[parsed_args + 1];
  len1 = strlen(seq1);
  len2 = strlen(seq2);
  if (len1 == 0 || len2 == 0) {
    env_error_set(env, "sequence of length 0 not allowed");
    return -1;
  }
  multilcptab = multilcp_compute(seq1, len1, seq2, len2, env);
  multilcp_show(multilcptab, len1, len2);
  array2dim_delete(multilcptab);

  return 0;
}
