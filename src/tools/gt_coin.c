/*
  Copyright (c) 2005-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2005-2008 Center for Bioinformatics, University of Hamburg

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

#include <ctype.h>
#include "libgtcore/ma.h"
#include "libgtcore/option.h"
#include "libgtcore/versionfunc.h"
#include "libgtcore/xansi.h"
#include "libgtext/coin_hmm.h"
#include "tools/gt_coin.h"

static OPrval parse_options(int *parsed_args, int argc, const char **argv,
                            Error *err)
{
  OptionParser *op;
  OPrval oprval;
  error_check(err);
  op = option_parser_new("sequence_of_coin_tosses", "Decode "
                         "'sequence_of_coin_tosses' and show the result on "
                         "stdout.");
  option_parser_set_min_max_args(op, 1, 1);
  oprval = option_parser_parse(op, parsed_args, argc, argv, versionfunc, err);
  option_parser_delete(op);
  return oprval;
}

int gt_coin(int argc, const char **argv, Error *err)
{
  unsigned int i, *emissions, *state_sequence = NULL, num_of_emissions;
  int parsed_args, had_err = 0;
  HMM *hmm = NULL;
  error_check(err);

  /* option parsing */
  switch (parse_options(&parsed_args, argc, argv, err)) {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR: return -1;
    case OPTIONPARSER_REQUESTS_EXIT: return 0;
  }
  assert(parsed_args == 1);

  /* save sequence */
  num_of_emissions = strlen(argv[1]);
  emissions = ma_malloc(sizeof (unsigned int) * num_of_emissions);
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
        error_set(err, "emissions[%u]=%c is not a valid character (only "
                       "`H' and `T' allowed)", i, (char) emissions[i]);
        had_err = -1;
    }
  }

  if (!had_err) {
    /* create the HMM */
    hmm = coin_hmm_loaded();

    /* decoding */
    state_sequence = ma_malloc(sizeof (unsigned int) * num_of_emissions);
    hmm_decode(hmm, state_sequence, emissions, num_of_emissions);

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
  hmm_delete(hmm);
  ma_free(emissions);
  ma_free(state_sequence);

  return had_err;
}
