/*
  Copyright (c) 2005-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2005-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <gt.h>

int gt_casino(int argc, char *argv[])
{
  unsigned int i, *emissions, *state_sequence, num_of_emissions;
  HMM *hmm;

  /* argument checking */
  if (argc != 2) {
    fprintf(stderr, "Usage: %s sequence_of_die_rolls\n", argv[0]);
    return EXIT_FAILURE;
  }

  /* save sequence */
  num_of_emissions = strlen(argv[1]);
  emissions = xmalloc(sizeof(unsigned int) * num_of_emissions);
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
        fprintf(stderr , "emissions[%u]=%c is not a valid character (only `1' "
                         "to `6' allowed)", i, emissions[i]);
        return EXIT_FAILURE;
    }
  }

  /* create the HMM */
  hmm = dice_hmm_loaded();

  /* decoding */
  state_sequence = xmalloc(sizeof(unsigned int) * num_of_emissions);
  hmm_decode(hmm, state_sequence, emissions, num_of_emissions);

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

  /* free */
  hmm_free(hmm);
  free(emissions);
  free(state_sequence);

  return (EXIT_SUCCESS);
}
