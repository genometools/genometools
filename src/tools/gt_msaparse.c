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

#include "libgtcore/fileutils.h"
#include "libgtcore/option.h"
#include "libgtcore/versionfunc.h"
#include "libgtext/msa.h"

typedef struct {
  bool show,
       consensus,
       sumofpairs;
} MSAparse_arguments;

static OPrval parse_options(int *parsed_args, MSAparse_arguments *arguments,
                            int argc, const char **argv, Error *err)
{
  OptionParser *op;
  Option *o;
  OPrval oprval;
  error_check(err);
  op = option_parser_new("[option ...] MSA_file",
                         "Parse multiple sequence alignment (MSA) file and "
                         "optionally show score(s).");
  /* -show */
  o = option_new_bool("show", "show the parsed MSA on stdout", &arguments->show,
                      false);
  option_parser_add_option(op, o);
  /* -consensus */
  o = option_new_bool("consensus", "show consensus distance",
                      &arguments->consensus, false);
  option_parser_add_option(op, o);
  /* -sumofpairs */
  o = option_new_bool("sumofpairs", "show optimal sum of pairwise scores",
                      &arguments->sumofpairs, false);
  option_parser_add_option(op, o);
  /* parse */
  oprval = option_parser_parse_min_max_args(op, parsed_args, argc, argv,
                                            versionfunc, 1, 1, err);
  option_parser_delete(op);
  return oprval;
}

int gt_msaparse(int argc, const char **argv, Env *env)
{
  MSAparse_arguments arguments;
  int parsed_args, had_err = 0;
  MSA *msa = NULL;
  env_error_check(env);

  /* option parsing */
  switch (parse_options(&parsed_args, &arguments, argc, argv, env_error(env))) {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR: return -1;
    case OPTIONPARSER_REQUESTS_EXIT: return 0;
  }

  /* make sure sequence_file exists */
  assert(parsed_args < argc);
  if (!file_exists(argv[parsed_args])) {
    env_error_set(env, "MSA_file '%s' does not exist", argv[parsed_args]);
    had_err = -1;
  }

  if (!had_err) {
    /* multiple sequence alignment construction */
    msa = msa_new(argv[parsed_args], env_error(env));
    if (!msa)
      had_err = -1;
  }

  if (!had_err) {
    /* output */
    if (arguments.show)
      msa_show(msa);
    if (arguments.consensus)
      printf("consensus distance: %lu\n", msa_consensus_distance(msa, env));
    if (arguments.sumofpairs) {
      printf("sum of pairwise scores: %lu\n",
             msa_sum_of_pairwise_scores(msa, env));
    }
  }

  /* free */
  msa_delete(msa);

  return had_err;
}
