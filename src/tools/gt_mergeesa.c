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

#include "libgtcore/error.h"
#include "libgtcore/option.h"
#include "libgtcore/versionfunc.h"
#include "libgtmatch/verbose-def.h"
#include "libgtmatch/test-mergeesa.pr"
#include "tools/gt_mergeesa.h"

static OPrval parse_options(Str *indexname,StrArray *indexnametab,
                            int *parsed_args, int argc,
                            const char **argv, Error *err)
{
  OptionParser *op;
  OPrval oprval;
  Option *option;

  error_check(err);
  op = option_parser_new("storeindex <mkvindex1> <mkvindex2> ...",
                         "Merge indexes into one index.");
  option_parser_set_mailaddress(op,"<kurtz@zbh.uni-hamburg.de>");
  option = option_new_filenamearray("ii",
                                    "specify input index files (mandatory)",
                                    indexnametab);
  option_is_mandatory(option);
  option_parser_add_option(op, option);

  option = option_new_string("indexname",
                             "specify index to be created",
                             indexname, NULL);

  option_is_mandatory(option);
  option_parser_add_option(op, option);

  oprval = option_parser_parse(op, parsed_args, argc, argv, versionfunc, err);
  option_parser_delete(op);
  return oprval;
}

int gt_mergeesa(int argc, const char **argv, Error *err)
{
  Str *storeindex;
  StrArray *indexnametab;
  bool haserr = false;
  int parsed_args;

  error_check(err);

  storeindex = str_new();
  indexnametab = strarray_new();
  switch (parse_options(storeindex, indexnametab, &parsed_args, argc, argv,
                        err)) {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR:
         haserr = true; break;
    case OPTIONPARSER_REQUESTS_EXIT: return 0;
  }
  if (!haserr)
  {
    unsigned long i;
    Verboseinfo *verboseinfo;

    printf("# storeindex=%s\n",str_get(storeindex));
    for (i=0; i<strarray_size(indexnametab); i++)
    {
      printf("# input=%s\n",strarray_get(indexnametab,i));
    }
    verboseinfo = newverboseinfo(false);
    if (performtheindexmerging(storeindex,
                              indexnametab,
                              verboseinfo,
                              err) != 0)
    {
      haserr = true;
    }
    freeverboseinfo(&verboseinfo);
  }
  str_delete(storeindex);
  strarray_delete(indexnametab);
  return haserr ? -1 : 0;
}
