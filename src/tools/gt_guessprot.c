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

#include <stdlib.h>
#include "core/fileutils_api.h"
#include "core/option_api.h"
#include "core/versionfunc.h"
#include "tools/gt_guessprot.h"

static GtOPrval parse_options(int *parsed_args, int argc, const char **argv,
                              GtError *err)
{
  GtOptionParser *op;
  GtOPrval oprval;
  gt_error_check(err);
  op = gt_option_parser_new("filenames",
                         "Guess if sequence in filenames is protein or DNA.");
  gt_option_parser_set_min_args(op, 1U);
  oprval = gt_option_parser_parse(op, parsed_args, argc, argv, gt_versionfunc,
                                  err);
  gt_option_parser_delete(op);
  return oprval;
}

int gt_guessprot(int argc, const char **argv, GtError *err)
{
  int i, parsed_args, retval;
  GtStrArray *filenametab;

  gt_error_check(err);

  switch (parse_options(&parsed_args, argc, argv, err)) {
    case GT_OPTION_PARSER_OK: break;
    case GT_OPTION_PARSER_ERROR: return -1;
    case GT_OPTION_PARSER_REQUESTS_EXIT: return 0;
  }

  filenametab = gt_str_array_new();
  for (i=parsed_args; i < argc; i++)
  {
    gt_str_array_add_cstr(filenametab,argv[i]);
  }
  retval = gt_files_guess_if_protein_sequences(filenametab,err);
  gt_str_array_delete(filenametab);
  if (retval < 0)
  {
    return -1;
  }
  if (retval == 1)
  {
    /*@ignore@*/
    exit(EXIT_FAILURE); /* guess it is a protein sequence */
    /*@end@*/

  } else
  {
    return 0;
  }
}
