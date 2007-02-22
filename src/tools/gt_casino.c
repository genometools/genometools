/*
  Copyright (c) 2005-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2005-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "gt.h"

static OPrval parse_options(int *parsed_args, int argc, char **argv, Env *env)
{
  OptionParser *op;
  OPrval oprval;
  env_error_check(env);
  op = option_parser_new("sequence_of_die_rolls", "Decode "
                         "'sequence_of_die_rolls' and show the result on "
                         "stdout.", env);
  oprval = option_parser_parse_min_max_args(op, parsed_args, argc, argv,
                                            versionfunc, 1, 1, env);
  option_parser_delete(op, env);
  return oprval;
}

int gt_casino(int argc, char *argv[], Env *env)
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
    emissions[i] = argv[1][i];
    switch (emissions[i]) {
      case '1':
        emissions[i] = ONE;
        break;
      case '2':
        emissions[i] = TWO;
        break;
      case '3':
        emissions[i] = THREE;
        break;
      case '4':
        emissions[i] = FOUR;
        break;
      case '5':
        emissions[i] = FIVE;
        break;
      case '6':
        emissions[i] = SIX;
        break;
      default:
        env_error_set(env, "emissions[%u]=%c is not a valid character (only "
                      "`1' to `6' allowed)", i, emissions[i]);
        has_err = -1;
    }
  }

  if (!has_err) {
    /* create the HMM */
    hmm = dice_hmm_loaded(env);

    /* decoding */
    state_sequence = env_ma_malloc(env,
                                   sizeof (unsigned int) * num_of_emissions);
    hmm_decode(hmm, state_sequence, emissions, num_of_emissions, env);

    /* print most probable state sequence state sequence */
    for (i = 0 ; i < num_of_emissions; i++) {
      switch (state_sequence[i]) {
        case DICE_FAIR:
          putchar('F');
          break;
        case DICE_LOADED:
          putchar('L');
          break;
        default: assert(0);
      }
    }
    putchar('\n');
  }

  /* free */
  hmm_delete(hmm, env);
  env_ma_free(emissions, env);
  env_ma_free(state_sequence, env);

  return has_err;
}
