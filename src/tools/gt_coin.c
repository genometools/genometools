/*
  Copyright (c) 2005-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2005-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <ctype.h>
#include "gt.h"

static OPrval parse_options(int *parsed_args, int argc, const char **argv,
                            Env *env)
{
  OptionParser *op;
  OPrval oprval;
  env_error_check(env);
  op = option_parser_new("sequence_of_coin_tosses", "Decode "
                         "'sequence_of_coin_tosses' and show the result on "
                         "stdout.", env);
  oprval = option_parser_parse_min_max_args(op, parsed_args, argc, argv,
                                            versionfunc, 1, 1, env);
  option_parser_delete(op, env);
  return oprval;
}

int gt_coin(int argc, const char **argv, Env *env)
{
  unsigned int i, *emissions, *state_sequence = NULL, num_of_emissions;
  int parsed_args, has_err = 0;
  HMM *hmm = NULL;
  env_error_check(env);

  /* option parsing */
  switch (parse_options(&parsed_args, argc, argv, env)) {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR: return -1;
    case OPTIONPARSER_REQUESTS_EXIT: return 0;
  }
  assert(parsed_args == 1);

  /* save sequence */
  num_of_emissions = strlen(argv[1]);
  emissions = env_ma_malloc(env, sizeof (unsigned int) * num_of_emissions);
  for (i = 0; i < num_of_emissions; i++) {
    emissions[i] = toupper(argv[1][i]);
    switch (emissions[i]) {
      case 'H':
        emissions[i] = HEAD;
        break;
      case 'T':
        emissions[i] = TAIL;
        break;
      default:
        env_error_set(env, "emissions[%u]=%c is not a valid character (only "
                      "`H' and `T' allowed)", i, (char) emissions[i]);
        has_err = -1;
    }
  }

  if (!has_err) {
    /* create the HMM */
    hmm = coin_hmm_loaded(env);

    /* decoding */
    state_sequence = env_ma_malloc(env,
                                   sizeof (unsigned int) * num_of_emissions);
    hmm_decode(hmm, state_sequence, emissions, num_of_emissions, env);

    /* print most probable state sequence state sequence */
    for (i = 0 ; i < num_of_emissions; i++) {
      switch (state_sequence[i]) {
        case COIN_FAIR:
          xputchar('F');
          break;
        case COIN_LOADED:
          xputchar('L');
          break;
        default: assert(0);
      }
    }
    xputchar('\n');
  }

  /* free */
  hmm_delete(hmm, env);
  env_ma_free(emissions, env);
  env_ma_free(state_sequence, env);

  return has_err;
}
