/*
  Copyright (c) 2007 Thomas Jahns <Thomas.Jahns@gmx.net>

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

#include "gt_packedindex_chk_integrity.h"
#include "libgtcore/ensure.h"
#include "libgtcore/env.h"
#include "libgtcore/option.h"
#include "libgtcore/versionfunc.h"

#include "libgtmatch/eis-encidxseq.h"

#define DEFAULT_PROGRESS_INTERVAL  100000UL

/***************************************************************************
 * routines to validate basic index integrity
 ***************************************************************************/

struct chkIndexOptions
{
  unsigned long skipCount;
  unsigned long progressInterval;
};

static OPrval
parseChkIndexOptions(int *parsed_args, int argc, const char *argv[],
                     struct chkIndexOptions *optOut, Env *env);

extern int
gt_packedindex_chk_integrity(int argc, const char *argv[], Env *env)
{
  struct encIdxSeq *seq;
  struct chkIndexOptions options;
  Str *inputProject;
  int parsedArgs;
  int had_err = 0;
  env_error_check(env);

  switch (parseChkIndexOptions(&parsedArgs, argc, argv, &options, env))
  {
    case OPTIONPARSER_OK:
      break;
    case OPTIONPARSER_ERROR:
      return -1;
    case OPTIONPARSER_REQUESTS_EXIT:
      return 0;
  }

  inputProject = str_new_cstr(argv[parsedArgs], env);
  env_error_check(env);
  seq = loadBlockEncIdxSeq(inputProject, EIS_FEATURE_REGION_SUMS, env);
  ensure(had_err, seq);
  if (had_err)
  {
    env_error_set(env, "Failed to load index: %s", str_get(inputProject));
  }
  else
  {
    fprintf(stderr, "# Using index over sequence "FormatSeqpos
            " symbols long.\n", EISLength(seq));
    {
      int corrupt;
      ensure(had_err,
             !(corrupt = verifyIntegrity(seq, inputProject, options.skipCount,
                                         options.progressInterval, stderr,
                                         env)));
      if (corrupt)
      {
        fputs("Integrity check failed for index.\n", stderr);
        fputs(EISintegrityCheckResultStrings[corrupt], stderr);
        fputs("\n", stderr);
      }
    }
  }
  if (seq)
    deleteEncIdxSeq(seq, env);
  if (inputProject)
    str_delete(inputProject, env);
  return had_err?-1:0;
}

static OPrval
parseChkIndexOptions(int *parsed_args, int argc, const char *argv[],
                     struct chkIndexOptions *optOut,
                     Env *env)
{
  OptionParser *op;
  Option *option;
  OPrval oprval;

  env_error_check(env);
  op = option_parser_new("indexname",
                         "Map <indexname> block composition index"
                         "and bwt and check index integrity.",
                         env);
  option = option_new_ulong("skip", "number of symbols to skip",
                            &optOut->skipCount, 0,
                            env);
  option_parser_add_option(op, option, env);

  option = option_new_ulong("ticks", "print dot after this many symbols"
                            " tested okay", &optOut->progressInterval,
                            DEFAULT_PROGRESS_INTERVAL, env);
  option_parser_add_option(op, option, env);

  oprval = option_parser_parse_min_max_args(op, parsed_args, argc,
                                            (const char **)argv,
                                            versionfunc, 1, 1, env);
  option_parser_delete(op, env);
  return oprval;
}
