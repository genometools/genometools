/*
  Copyright (c) 2005-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2005-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <ctype.h>
#include <gt.h>

int gt_coin(int argc, char *argv[])
{
  unsigned int i, *emissions, *state_sequence, num_of_emissions;
  HMM *hmm;

  /* argument checking */
  if (argc != 2) {
    fprintf(stderr, "Usage: %s sequence_of_coin_tosses\n", argv[0]);
    return EXIT_FAILURE;
  }

  /* save sequence */
  num_of_emissions = strlen(argv[1]);
  emissions = xmalloc(sizeof(unsigned int) * num_of_emissions);
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
        fprintf(stderr , "emissions[%u]=%c is not a valid character (only `H' "
                         "and `T' allowed)", i, emissions[i]);
        return EXIT_FAILURE;
    }
  }

  /* create the HMM */
  hmm = coin_hmm_loaded();

  /* decoding */
  state_sequence = xmalloc(sizeof(unsigned int) * num_of_emissions);
  hmm_decode(hmm, state_sequence, emissions, num_of_emissions);

  /* print most probable state sequence state sequence */
  for (i = 0 ; i < num_of_emissions; i++) {
    switch (state_sequence[i]) {
      case COIN_FAIR:
        putchar('F');
        break;
      case COIN_LOADED:
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

  return EXIT_SUCCESS;
}
