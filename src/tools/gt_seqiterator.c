/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#include <sys/types.h>
#include <sys/stat.h>
#include "libgtcore/fileutils.h"
#include "libgtcore/option.h"
#include "libgtcore/seqiterator.h"
#include "libgtcore/versionfunc.h"
#include "libgtcore/xposix.h"
#include "libgtcore/progressbar.h"

static OPrval parse_options(bool *verbose,int *parsed_args,
                            int argc, const char **argv,Env *env)
{
  OptionParser *op;
  Option *option;
  OPrval oprval;

  env_error_check(env);
  op = option_parser_new("[options] file [...]",
                         "Parse the supplied Fasta files.",
                         env);
  option_parser_set_mailaddress(op,"<kurtz@zbh.uni-hamburg.de>");
  option= option_new_bool("v","be verbose",verbose,false,env);
  option_parser_add_option(op, option, env);
  oprval = option_parser_parse_min_args(op, parsed_args, argc, argv,
                                        versionfunc, (unsigned int) 1,
                                        env);
  option_parser_delete(op, env);
  return oprval;
}

int gt_seqiterator(int argc, const char **argv, Env *env)
{
  StrArray *files;
  SeqIterator *seqit;
  const Uchar *sequence;
  char *desc;
  unsigned long len;
  int i, parsed_args, had_err;
  off_t totalsize;
  bool verbose = false;

  env_error_check(env);

  /* option parsing */
  switch (parse_options(&verbose,&parsed_args, argc, argv, env)) {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR: return -1;
    case OPTIONPARSER_REQUESTS_EXIT: return 0;
  }

  files = strarray_new(env);
  for (i = parsed_args; i < argc; i++)
  {
    strarray_add_cstr(files, argv[i], env);
  }

  totalsize = files_estimate_total_size(files);
  printf("# estimated total size is %llu\n",(unsigned long long) totalsize);
  seqit = seqiterator_new(files, NULL, true, env);
  if (verbose)
  {
    progressbar_start(seqiterator_getcurrentcounter(seqit, (unsigned long long)
                                                           totalsize),
                                                           (unsigned long long)
                                                           totalsize);
  }
  while (true)
  {
    had_err = seqiterator_next(seqit, &sequence, &len, &desc, env);
    if (had_err != 1)
    {
      break;
    }
    env_ma_free(desc, env);
  }
  seqiterator_delete(seqit, env);
  if (verbose)
  {
    progressbar_stop();
  }
  strarray_delete(files, env);
  return had_err;
}
