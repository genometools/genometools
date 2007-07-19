/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "libgtcore/bioseq.h"
#include "libgtcore/cstr.h"
#include "libgtcore/fileutils.h"
#include "libgtcore/option.h"
#include "libgtcore/versionfunc.h"
#include "libgtcore/xansi.h"
#include "libgtext/qgramdist.h"

static OPrval parse_options(int *parsed_args, unsigned int *q, int argc,
                            const char **argv, Env *env)
{
  OptionParser *op;
  Option *o;
  OPrval oprval;
  env_error_check(env);
  op = option_parser_new("[option ...] seq_file_1 seq_file_2",
                         "Compute q-gram distance for each sequence "
                         "combination.", env);
  o = option_new_uint_min("q", "set q", q, 3, 1, env);
  option_parser_add_option(op, o, env);
  oprval = option_parser_parse_min_max_args(op, parsed_args, argc, argv,
                                            versionfunc, 2, 2, env);
  option_parser_delete(op, env);
  return oprval;
}

int gt_qgramdist(int argc, const char **argv, Env *env)
{
  Bioseq *bioseq_1 = NULL, *bioseq_2 = NULL;
  unsigned long i, j, dist;
  Seq *seq_1, *seq_2;
  int parsed_args, had_err = 0;
  unsigned int q;
  env_error_check(env);

  /* option parsing */
  switch (parse_options(&parsed_args, &q, argc, argv, env)) {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR: return -1;
    case OPTIONPARSER_REQUESTS_EXIT: return 0;
  }
  assert(parsed_args+1 < argc);

  /* make sure seq_file_1 exists */
  if (!file_exists(argv[parsed_args])) {
    env_error_set(env, "seq_file_1 \"%s\" does not exist", argv[parsed_args]);
    had_err = -1;
  }

  /* make sure seq_file_2 exists */
  if (!had_err && !file_exists(argv[parsed_args+1])) {
    env_error_set(env, "seq_file_2 \"%s\" does not exist", argv[parsed_args+1]);
    had_err = -1;
  }

  /* init */
  if (!had_err) {
    bioseq_1 = bioseq_new(argv[parsed_args], env);
    if (!bioseq_1)
      had_err = -1;
    if (!had_err) {
      bioseq_2 = bioseq_new(argv[parsed_args+1], env);
      if (!bioseq_2)
        had_err = -1;
    }

    /* compute q-gram distance for all sequence combinations */
    for (i = 0; i < bioseq_number_of_sequences(bioseq_1); i++) {
      for (j = 0; j < bioseq_number_of_sequences(bioseq_2); j++) {
        seq_1 = bioseq_get_seq(bioseq_1, i, env);
        seq_2 = bioseq_get_seq(bioseq_2, j, env);
        dist = qgramdist(seq_1, seq_2, q, env);
        printf("qgramdist_%u_(", q);
        cstr_show(seq_get_orig(seq_1), seq_length(seq_1), stdout);
        xputchar(',');
        cstr_show(seq_get_orig(seq_2), seq_length(seq_2), stdout);
        printf(")=%lu\n", dist);
      }
    }
  }

  /* free */
  bioseq_delete(bioseq_2, env);
  bioseq_delete(bioseq_1, env);

  return had_err;
}
