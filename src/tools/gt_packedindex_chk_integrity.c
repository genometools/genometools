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

/**
 * \file gt_packedindex_chk_integrity
 * calls routines to validate basic index integrity
 */

#include "gt_packedindex_chk_integrity.h"
#include "libgtcore/ensure.h"
#include "libgtcore/error.h"
#include "libgtcore/option.h"
#include "libgtcore/versionfunc.h"
#include "libgtmatch/eis-encidxseq.h"
#include "tools/gt_packedindex_chk_integrity.h"

#define DEFAULT_PROGRESS_INTERVAL  100000UL

struct chkIndexOptions
{
  unsigned long skipCount;
  unsigned long progressInterval;
  int checkFlags;
};

static OPrval
parseChkIndexOptions(int *parsed_args, int argc, const char *argv[],
                     struct chkIndexOptions *param, Error *err);

extern int
gt_packedindex_chk_integrity(int argc, const char *argv[], Error *err)
{
  struct encIdxSeq *seq;
  struct chkIndexOptions options;
  Str *inputProject;
  int parsedArgs;
  int had_err = 0;
  error_check(err);

  switch (parseChkIndexOptions(&parsedArgs, argc, argv, &options, err))
  {
    case OPTIONPARSER_OK:
      break;
    case OPTIONPARSER_ERROR:
      return -1;
    case OPTIONPARSER_REQUESTS_EXIT:
      return 0;
  }

  inputProject = str_new_cstr(argv[parsedArgs]);
  seq = loadBlockEncIdxSeq(inputProject, EIS_FEATURE_REGION_SUMS, err);
  ensure(had_err, seq);
  if (had_err)
  {
    error_set(err, "Failed to load index: %s", str_get(inputProject));
  }
  else
  {
    fprintf(stderr, "# Using index over sequence "FormatSeqpos
            " symbols long.\n", EISLength(seq));
    {
      int corrupt;
      ensure(
        had_err,
        !(corrupt = EISVerifyIntegrity(seq, inputProject, options.skipCount,
                                       options.progressInterval, stderr,
                                       options.checkFlags, err)));
      if (corrupt)
      {
        fputs("Integrity check failed for index.\n", stderr);
        fputs(EISIntegrityCheckResultStrings[corrupt], stderr);
        fputs("\n", stderr);
      }
    }
  }
  if (seq)
    deleteEncIdxSeq(seq);
  if (inputProject)
    str_delete(inputProject);
  return had_err?-1:0;
}

static OPrval
parseChkIndexOptions(int *parsed_args, int argc, const char *argv[],
                     struct chkIndexOptions *param, Error *err)
{
  OptionParser *op;
  Option *option;
  OPrval oprval;
  bool extRankCheck;

  error_check(err);
  op = option_parser_new("indexname",
                         "Map <indexname> block composition index"
                         "and bwt and check index integrity.");
  option = option_new_ulong("skip", "number of symbols to skip",
                            &param->skipCount, 0);
  option_parser_add_option(op, option);

  option = option_new_ulong("ticks", "print dot after this many symbols"
                            " tested okay", &param->progressInterval,
                            DEFAULT_PROGRESS_INTERVAL);
  option_parser_add_option(op, option);

  option = option_new_bool("ext-rank-check",
                           "do additional checks of rank query results",
                           &extRankCheck, false);
  option_parser_add_option(op, option);

  option_parser_set_min_max_args(op, 1, 1);
  oprval = option_parser_parse(op, parsed_args, argc, (const char**) argv,
                               versionfunc, err);
  option_parser_delete(op);
  param->checkFlags = EIS_VERIFY_BASIC
    | (extRankCheck?EIS_VERIFY_EXT_RANK:0);
  return oprval;
}
