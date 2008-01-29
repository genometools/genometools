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
#include "libgtcore/versionfunc.h"
#include "libgtext/matchcount.h"
#include "tools/gt_matchcount.h"

static OPrval parse_options(int *parsed_args, int argc, const char **argv,
                            Error *err)
{
  OptionParser *op;
  OPrval oprval;
  error_check(err);
  op = option_parser_new("[option ...] k seq1 seq2",
                         "Compute the match-count for each substring pair of "
                         "length k from seq1 and seq2.");
  option_parser_set_min_max_args(op, 3, 3);
  oprval = option_parser_parse(op, parsed_args, argc, argv, versionfunc, err);
  option_parser_delete(op);
  return oprval;
}

static void proc_match_count(int u_pos, int v_pos, int match_count)
{
  printf("mc(%2d, %2d)=%2d\n", u_pos, v_pos, match_count);
}

int gt_matchcount(int argc, const char **argv, Error *err)
{
  const char *seq1, *seq2;
  int k, len1, len2, parsed_args, had_err = 0;
  error_check(err);

  /* option parsing */
  switch (parse_options(&parsed_args, argc, argv, err)) {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR: return -1;
    case OPTIONPARSER_REQUESTS_EXIT: return 0;
  }
  assert(parsed_args + 2 < argc);

  if (sscanf(argv[parsed_args], "%d", &k) != 1 || k <= 0) {
    error_set(err, "first argument <k> must be positive integer");
    had_err = -1;
  }

  if (!had_err) {
    /* store pointer to sequences */
    seq1 = argv[parsed_args+1];
    seq2 = argv[parsed_args+2];

    /* determine sequence lengths */
    len1 = (int) strlen(seq1);
    len2 = (int) strlen(seq2);

    /* compute match count */
    printf("args=%d %s %s\n", k, seq1, seq2);
    matchcount(seq1, len1, seq2, len2, k, proc_match_count);
  }

  return had_err;
}
