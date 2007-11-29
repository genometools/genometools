/*
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

#include "libgtcore/bioseq.h"
#include "libgtcore/fasta.h"
#include "libgtcore/grep.h"
#include "libgtcore/option.h"
#include "libgtcore/outputfile.h"
#include "libgtcore/versionfunc.h"

typedef struct {
  Str *pattern;
  unsigned long width;
  GenFile *outfp;
} ExtractSeqArguments;

static OPrval parse_options(int *parsed_args, ExtractSeqArguments *arguments,
                            int argc, const char **argv, Error *err)
{
  OptionParser *op;
  OutputFileInfo *ofi;
  Option *option;
  OPrval oprval;
  error_check(err);

  /* init */
  op = option_parser_new("[option ...] [sequence_file ...]",
                         "Extract sequences from given sequence file(s).");
  ofi = outputfileinfo_new();

  /* -match */
  option = option_new_string("match", "extract all sequences whose description "
                             "matches the given pattern.\nThe given pattern "
                             "must be a valid extended regular expression.",
                             arguments->pattern, NULL);
  option_is_mandatory(option);
  option_parser_add_option(op, option);

  /* -width */
  option = option_new_ulong("width", "set output width for showing of "
                            "sequences (0 disables formatting)",
                            &arguments->width, 0);
  option_parser_add_option(op, option);

  /* output file options */
  outputfile_register_options(op, &arguments->outfp, ofi);

  /* parse */
  oprval = option_parser_parse(op, parsed_args, argc, argv, versionfunc, err);

  /* free */
  outputfileinfo_delete(ofi);
  option_parser_delete(op);

  return oprval;
}

static int extractseq(GenFile *outfp, Bioseq *bs, const char *pattern,
                      unsigned long width, Error *err)
{
  const char *desc;
  unsigned long i;
  bool match;
  int had_err = 0;

  error_check(err);
  assert(bs && pattern);

  for (i = 0; !had_err && i < bioseq_number_of_sequences(bs); i++) {
    desc = bioseq_get_description(bs, i);
    assert(desc);
    had_err = grep(&match, pattern, desc, err);
    if (!had_err && match) {
      fasta_show_entry_generic(desc, bioseq_get_sequence(bs, i),
                               bioseq_get_sequence_length(bs, i), width, outfp);
    }
  }

  return had_err;
}

int gt_extractseq(int argc, const char **argv, Env *env)
{
  ExtractSeqArguments arguments;
  Bioseq *bs;
  int parsed_args, had_err = 0;
  env_error_check(env);

  /* option parsing */
  arguments.pattern = str_new();
  switch (parse_options(&parsed_args, &arguments, argc, argv, env_error(env))) {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR:
      str_delete(arguments.pattern);
      return -1;
    case OPTIONPARSER_REQUESTS_EXIT:
      str_delete(arguments.pattern);
      return 0;
  }

  if (parsed_args == argc) { /* no file given, use stdin */
    if (!(bs = bioseq_new("-", env_error(env))))
      had_err = -1;
    if (!had_err) {
      had_err = extractseq(arguments.outfp, bs, str_get(arguments.pattern),
                           arguments.width, env_error(env));
    }
    bioseq_delete(bs);
  }

  /* process all files */
  while (!had_err && parsed_args < argc) {
    if (!(bs = bioseq_new(argv[parsed_args], env_error(env))))
      had_err = -1;
    if (!had_err) {
      had_err = extractseq(arguments.outfp, bs, str_get(arguments.pattern),
                           arguments.width, env_error(env));
    }
    bioseq_delete(bs);
    parsed_args++;
  }

  /* free */
  str_delete(arguments.pattern);
  genfile_close(arguments.outfp);

  return had_err;
}
