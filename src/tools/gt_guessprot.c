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

#include "libgtcore/versionfunc.h"
#include "libgtcore/option.h"
#include "libgtmatch/guessprot.pr"
#include "tools/gt_guessprot.h"

static OPrval parse_options(int *parsed_args, int argc, const char **argv,
                            Error *err)
{
  OptionParser *op;
  OPrval oprval;
  error_check(err);
  op = option_parser_new("filenames",
                         "guess if sequence in filenames is protein or DNA.");
  option_parser_set_min_args(op, 1U);
  oprval = option_parser_parse(op, parsed_args, argc, argv, versionfunc, err);
  option_parser_delete(op);
  return oprval;
}

int gt_guessprot(int argc, const char **argv, Error *err)
{
  int i, parsed_args, retval;
  StrArray *filenametab;

  error_check(err);

  switch (parse_options(&parsed_args, argc, argv, err)) {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR: return -1;
    case OPTIONPARSER_REQUESTS_EXIT: return 0;
  }

  filenametab = strarray_new();
  for (i=parsed_args; i < argc; i++)
  {
    strarray_add_cstr(filenametab,argv[i]);
  }
  retval = guessifproteinsequencestream(filenametab,err);
  strarray_delete(filenametab);
  if (retval < 0)
  {
    return -1;
  }
  if (retval == 1)
  {
    /*@ignore@*/
    exit(1); /* XXX */
    /*@end@*/

  } else
  {
    return 0;
  }
}
