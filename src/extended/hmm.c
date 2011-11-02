/*
  Copyright (c) 2005-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2005-2007 Center for Bioinformatics, University of Hamburg

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

/*
  The algorithms in this file are described in chapter 3 of the book

  [DEKM98] R. Durbin, S.R. Eddy, A. Krogh, and G. Mitchison.
  Biological sequence analysis: probabilistic models of proteins and nucleic
  acids. Cambridge University Press, 1998.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "core/array2dim_api.h"
#include "core/assert_api.h"
#include "core/chardef.h"
#include "core/ensure.h"
#include "core/mathsupport.h"
#include "core/xansi_api.h"
#include "extended/coin_hmm.h"
#include "extended/dice_hmm.h"
#include "extended/hmm.h"

#define MINUSINFINITY   -99999.0
#define PSEUDOCOUNT     1

struct GtHMM {
  unsigned int num_of_states,
               num_of_symbols;
  double *initial_state_prob, /* log values */
         **transition_prob,   /* log values */
         **emission_prob;     /* log values */
};

GtHMM* gt_hmm_new(unsigned int num_of_states, unsigned int num_of_symbols)
{
  GtHMM *hmm;
  unsigned int i, j;

  gt_assert(num_of_states && num_of_symbols);

  /* alloc */
  hmm = gt_malloc(sizeof *hmm);
  hmm->initial_state_prob = gt_malloc(sizeof (double) * num_of_states);
  gt_array2dim_malloc(hmm->transition_prob, num_of_states, num_of_states);
  gt_array2dim_malloc(hmm->emission_prob, num_of_states, num_of_symbols);

  /* init */
  hmm->num_of_states = num_of_states;
  hmm->num_of_symbols = num_of_symbols;

  for (i = 0; i < num_of_states; i++)
    gt_hmm_set_initial_state_probability(hmm, i, 1.0 / num_of_states);

  for (i = 0; i < num_of_states; i++) {
    for (j = 0; j < num_of_states; j++)
      gt_hmm_set_transition_probability(hmm, i, j, 0.0);
  }

  for (i = 0; i < num_of_states; i++) {
    for (j = 0; j < num_of_symbols; j++)
      gt_hmm_set_emission_probability(hmm, i, j, 0.0);
  }

  return hmm;
}

void gt_hmm_set_initial_state_probability(GtHMM *hmm, unsigned int state_num,
                                       double probability)
{
  gt_assert(hmm);
  gt_assert(probability >= 0.0 && probability <= 1.0);
  gt_assert(state_num < hmm->num_of_states);
  if (probability == 0.0)
    hmm->initial_state_prob[state_num] = MINUSINFINITY;
  else
    hmm->initial_state_prob[state_num] = log(probability);
}

double gt_hmm_get_initial_state_probability(const GtHMM *hmm,
                                            unsigned int state_num)
{
  gt_assert(hmm);
  gt_assert(state_num < hmm->num_of_states);
  if (hmm->initial_state_prob[state_num] == MINUSINFINITY)
    return 0.0;
  return exp(hmm->initial_state_prob[state_num]);
}

void gt_hmm_set_transition_probability(GtHMM *hmm,
                                    unsigned int from_state_num,
                                    unsigned int to_state_num,
                                    double probability)
{
  gt_assert(hmm);
  gt_assert(probability >= 0.0 && probability <= 1.0);
  gt_assert(from_state_num < hmm->num_of_states);
  gt_assert(to_state_num < hmm->num_of_states);
  if (probability == 0.0)
    hmm->transition_prob[from_state_num][to_state_num] = MINUSINFINITY;
  else
    hmm->transition_prob[from_state_num][to_state_num] = log(probability);
}

double gt_hmm_get_transition_probability(const GtHMM *hmm,
                                      unsigned int from_state_num,
                                      unsigned int to_state_num)
{
  gt_assert(hmm);
  gt_assert(from_state_num < hmm->num_of_states);
  gt_assert(to_state_num < hmm->num_of_states);
  if (hmm->transition_prob[from_state_num][to_state_num] == MINUSINFINITY)
    return 0.0;
  return exp(hmm->transition_prob[from_state_num][to_state_num]);
}

void gt_hmm_set_missing_transition_probabilities(GtHMM *hmm)
{
  unsigned int row, column, num_of_missing_entries;
  double prob, sum_of_probabilities;

  gt_assert(hmm);

  for (row = 0; row < hmm->num_of_states; row++) {
    /* determine sum of probabilities */
    sum_of_probabilities = 0.0;
    num_of_missing_entries = 0;
    for (column = 0; column < hmm->num_of_states; column++) {
      prob = gt_hmm_get_transition_probability(hmm, row, column);
      if (prob == 0.0)
        num_of_missing_entries++;
      else
        sum_of_probabilities += prob;
    }

    /* set missing probabilites (equal distribution) */
    if (num_of_missing_entries) {
      for (column = 0; column < hmm->num_of_states; column++) {
        prob = gt_hmm_get_transition_probability(hmm, row, column);
        if (prob == 0.0) {
          gt_hmm_set_transition_probability(hmm, row, column,
                                         (1.0 - sum_of_probabilities) /
                                         num_of_missing_entries);
        }
      }
    }
  }
}

void gt_hmm_set_emission_probability(GtHMM *hmm,
                                  unsigned int state_num,
                                  unsigned int symbol_num,
                                  double probability)
{
  gt_assert(hmm);
  gt_assert(probability >= 0.0 && probability <= 1.0);
  gt_assert(state_num < hmm->num_of_states);
  symbol_num = (symbol_num == WILDCARD) ? hmm->num_of_symbols - 1 : symbol_num;
  gt_assert(symbol_num < hmm->num_of_symbols);
  if (probability == 0.0)
    hmm->emission_prob[state_num][symbol_num] = MINUSINFINITY;
  else
    hmm->emission_prob[state_num][symbol_num] = log(probability);
}

double gt_hmm_get_emission_probability(const GtHMM *hmm,
                                    unsigned int state_num,
                                    unsigned int symbol_num)
{
  gt_assert(hmm);
  gt_assert(state_num < hmm->num_of_states);
  symbol_num = (symbol_num == WILDCARD) ? hmm->num_of_symbols - 1 : symbol_num;
  gt_assert(symbol_num < hmm->num_of_symbols);
  if (hmm->emission_prob[state_num][symbol_num] == MINUSINFINITY)
    return 0.0;
  return exp(hmm->emission_prob[state_num][symbol_num]);
}

static bool gt_hmm_has_valid_initial_state_probs(const GtHMM *hmm)
{
  unsigned int state;
  double sum_of_probabilities = 0.0;

  gt_assert(hmm);

  for (state = 0; state < hmm->num_of_states; state++)
    sum_of_probabilities += gt_hmm_get_initial_state_probability(hmm, state);

  return gt_double_equals_one(sum_of_probabilities);
}

static bool gt_hmm_has_valid_emissions(const GtHMM *hmm)
{
  unsigned int state, symbol;
  double sum_of_probabilities;

  gt_assert(hmm);

  for (state = 0; state < hmm->num_of_states; state++) {
    sum_of_probabilities = 0.0;
    for (symbol = 0; symbol < hmm->num_of_symbols; symbol++)
      sum_of_probabilities += gt_hmm_get_emission_probability(hmm,
                                                              state, symbol);
    if (!gt_double_equals_one(sum_of_probabilities))
      return false;
  }
  return true;
}

static bool gt_hmm_has_valid_states(const GtHMM *hmm)
{
  unsigned int state_1, state_2;
  double sum_of_probabilities;

  gt_assert(hmm);

  for (state_1 = 0; state_1 < hmm->num_of_states; state_1++) {
    sum_of_probabilities = 0.0;
    for (state_2 = 0; state_2 < hmm->num_of_states; state_2++) {
      sum_of_probabilities += gt_hmm_get_transition_probability(hmm,state_1,
                                                             state_2);

    }
    if (!gt_double_equals_one(sum_of_probabilities))
      return false;
  }

  return true;
}

bool gt_hmm_is_valid(const GtHMM *hmm)
{
  gt_assert(hmm);
  return (gt_hmm_has_valid_initial_state_probs(hmm) &&
          gt_hmm_has_valid_emissions(hmm) &&
          gt_hmm_has_valid_states(hmm));

}

void gt_hmm_init_random(GtHMM *hmm)
{
  double random_value, cumulative_prob;
  unsigned int i, j;
  gt_assert(hmm);

  /* initialize initial state probabilities in random fashion */
  cumulative_prob = 0.0;
  for (i = 0; i < hmm->num_of_states - 1; i++) {
    random_value = gt_rand_max_double(1.0 - cumulative_prob);
    gt_hmm_set_initial_state_probability(hmm, i, random_value);
    cumulative_prob += random_value;
  }
  gt_assert(cumulative_prob <= 1.0);
  gt_hmm_set_initial_state_probability(hmm, i,  1.0 - cumulative_prob);

  /* initialize transition probabilities in random fashion */
  for (i = 0; i < hmm->num_of_states; i++) {
    cumulative_prob = 0.0;
    for (j = 0; j < hmm->num_of_states - 1; j++) {
      random_value = gt_rand_max_double(1.0 - cumulative_prob);
      gt_hmm_set_transition_probability(hmm, i, j, random_value);
      cumulative_prob += random_value;
    }
    gt_assert(cumulative_prob <= 1.0);
    gt_hmm_set_transition_probability(hmm, i, j, 1.0 - cumulative_prob);
  }

  /* initialize emission probabilities in random fashion */
  for (i = 0; i < hmm->num_of_states; i++) {
    cumulative_prob = 0.0;
    for (j = 0; j < hmm->num_of_symbols - 1; j++) {
      random_value = gt_rand_max_double(1.0 - cumulative_prob);
      gt_hmm_set_emission_probability(hmm, i, j, random_value);
      cumulative_prob += random_value;
    }
    gt_assert(cumulative_prob <= 1.0);
    gt_hmm_set_emission_probability(hmm, i, j, 1.0 - cumulative_prob);
  }

  gt_assert(gt_hmm_is_valid(hmm));
}

/* [DEKM98, p. 56] */
void gt_hmm_decode(const GtHMM *hmm,
                unsigned int *state_sequence,
                const unsigned int *emissions,
                unsigned int num_of_emissions)
{
  double **max_probabilities, tmp_prob;
  unsigned int **backtrace, colidx, precolidx, emission;
  int row, column, num_of_rows, num_of_columns, previous_row;

  gt_assert(hmm);
  gt_assert(gt_hmm_is_valid(hmm));
  gt_assert(num_of_emissions);

  /* alloc tables */
  num_of_rows = hmm->num_of_states;
  num_of_columns = num_of_emissions;
  gt_array2dim_malloc(max_probabilities, num_of_rows, 2);
  gt_array2dim_malloc(backtrace, num_of_rows, num_of_columns);

  /* fill DP table */
  for (row = 0; row < num_of_rows; row++) { /* first column */
    if (emissions[0] == WILDCARD)
      emission = hmm->num_of_symbols - 1;
    else
      emission = emissions[0];
    gt_assert(emission < hmm->num_of_symbols);
    max_probabilities[row][0] = hmm->initial_state_prob[row] +
                                hmm->emission_prob[row][emission];
    backtrace[row][0] = row;
  }

  for (column = 1; column < num_of_columns; column++) { /* other columns */
    if (emissions[column] == WILDCARD)
      emission = hmm->num_of_symbols - 1;
    else
      emission = emissions[column];
    gt_assert(emission < hmm->num_of_symbols);
    colidx = column & 1;
    precolidx = (column - 1) & 1;
    for (row = 0; row < num_of_rows; row++) {
      max_probabilities[row][colidx] = max_probabilities[0][precolidx] +
                                     hmm->transition_prob[0][row] +
                                     hmm->emission_prob[row][emission];
      backtrace[row][column] = 0;
      for (previous_row = 1; previous_row < num_of_rows; previous_row++) {
        tmp_prob = max_probabilities[previous_row][precolidx] +
                   hmm->transition_prob[previous_row][row] +
                   hmm->emission_prob[row][emission];
        if (tmp_prob > max_probabilities[row][colidx]) {
          max_probabilities[row][colidx] = tmp_prob;
          backtrace[row][column] = previous_row;
        }
      }
    }
  }

  /* backtracing, determine end state with maximum probability */
  colidx = (num_of_columns - 1) & 1;
  tmp_prob = max_probabilities[0][colidx];
  state_sequence[num_of_columns - 1] = 0;
  for (row = 1; row < num_of_rows; row++) {
    if (max_probabilities[row][colidx] > tmp_prob)
      state_sequence[num_of_columns - 1] = row;
  }

  /* backtracing, follow the links */
  for (column = num_of_columns - 2; column >= 0; column--)
    state_sequence[column] = backtrace[state_sequence[column + 1]][column + 1];

  /* free tables */
  gt_array2dim_delete(backtrace);
  gt_array2dim_delete(max_probabilities);
}

/* [DEKM98, p. 58] */
static void compute_forward_table(double **f, const GtHMM *hmm,
                                  const unsigned int *emissions,
                                  unsigned long num_of_emissions)
{
  int row, column, previous_row;
  double tmp_prob;

  gt_assert(f && hmm && emissions && num_of_emissions);

  /* fill DP table 'f' (hmm->num_of_states * num_of_emissions) */
  for (row = 0; row < hmm->num_of_states; row++) { /* first column */
    gt_assert(emissions[0] < hmm->num_of_symbols);
    f[row][0] = hmm->initial_state_prob[row] +
                hmm->emission_prob[row][emissions[0]];
  }

  for (column = 1; column < num_of_emissions; column++) { /* other columns */
    gt_assert(emissions[column] < hmm->num_of_symbols);
    for (row = 0; row < hmm->num_of_states; row++) {
      f[row][column] = hmm->emission_prob[row][emissions[column]];
      tmp_prob = f[0][column-1] + hmm->transition_prob[0][row];
      for (previous_row = 1; previous_row < hmm->num_of_states;
           previous_row++) {
        /* XXX: replace the logsum() call with a tabulated version */
        tmp_prob = gt_logsum(tmp_prob, f[previous_row][column-1] +
                                       hmm->transition_prob[previous_row][row]);
      }
      f[row][column] += tmp_prob;
    }
  }
}

/* [DEKM98, p. 58] */
double gt_hmm_forward(const GtHMM* hmm, const unsigned int *emissions,
                   unsigned int num_of_emissions)
{
  unsigned int i;
  double **f, P;

  gt_assert(hmm && emissions && num_of_emissions);
  gt_array2dim_malloc(f, hmm->num_of_states, num_of_emissions);

  /* XXX: we do not need the full table here, the last column would suffice */
  compute_forward_table(f, hmm, emissions, num_of_emissions);

  /* compute P(x) */
  P = f[0][num_of_emissions-1];
  for (i = 1; i < hmm->num_of_states; i++) {
    /* XXX: replace the logsum() call with a tabulated version */
    P = gt_logsum(P, f[i][num_of_emissions-1]);
  }

  gt_array2dim_delete(f);
  return P;
}

/* [DEKM98, p. 59] */
static void compute_backward_table(double **b, const GtHMM *hmm,
                                   const unsigned int *emissions,
                                   unsigned long num_of_emissions)
{
  int row, column, next_row;
  double tmp_prob;

  gt_assert(b && hmm && emissions && num_of_emissions);

  /* fill DP table 'b' (hmm->num_of_states * num_of_emissions) */
  for (row = 0; row < hmm->num_of_states; row++) { /* last column */
    b[row][num_of_emissions-1] = 0.0; /* probability = 1.0 */
  }

  for (column = num_of_emissions - 2; column >= 0; column--) { /* other cols */
    gt_assert(emissions[column] < hmm->num_of_symbols);
    for (row = 0; row < hmm->num_of_states; row++) {
      tmp_prob = hmm->transition_prob[row][0] +
                 hmm->emission_prob[0][emissions[column+1]] +
                 b[0][column+1];
      for (next_row = 1; next_row < hmm->num_of_states; next_row++) {
        /* XXX: replace the logsum() call with a tabulated version */
        tmp_prob = gt_logsum(tmp_prob, hmm->transition_prob[row][next_row] +
                             hmm->emission_prob[next_row][emissions[column+1]] +
                             b[next_row][column+1]);
      }
      b[row][column] = tmp_prob;
    }
  }
}

/* [DEKM98, p. 59] */
double gt_hmm_backward(const GtHMM* hmm, const unsigned int *emissions,
                    unsigned int num_of_emissions)
{
  unsigned int i;
  double **b, P;

  gt_assert(hmm && emissions && num_of_emissions);
  gt_array2dim_malloc(b, hmm->num_of_states, num_of_emissions);

  /* XXX: we do not need the full table here, the last column would suffice */
  compute_backward_table(b, hmm, emissions, num_of_emissions);

  /* compute P(x) */
  P = hmm->initial_state_prob[0] + hmm->emission_prob[0][emissions[0]] +
      b[0][0];
  for (i = 1; i < hmm->num_of_states; i++) {
    /* XXX: replace the logsum() call with a tabulated version */
    P = gt_logsum(P, hmm->initial_state_prob[i] +
                     hmm->emission_prob[i][emissions[0]] + b[i][0]);
  }

  gt_array2dim_delete(b);
  return P;
}

void gt_hmm_emit(GtHMM *hmm, unsigned long num_of_emissions,
              void (*proc_emission)(unsigned int symbol, void *data),
              void *data)
{
  unsigned long i;
  unsigned int state, symbol, next_state;
  double random_value, cumulative_prob;
  gt_assert(hmm);

  /* determine initial state */
  random_value = gt_rand_0_to_1();
  cumulative_prob = 0.0;
  for (state = 0; state < hmm->num_of_states - 1; state++) {
    cumulative_prob += gt_hmm_get_initial_state_probability(hmm, state);
    if (cumulative_prob > random_value)
      break;
  }

  /* emit loop */
  for (i = 0; i < num_of_emissions; i++) {
    /* emit character */
    random_value = gt_rand_0_to_1();
    cumulative_prob = 0.0;
    for (symbol = 0; symbol < hmm->num_of_symbols - 1; symbol++) {
      cumulative_prob += gt_hmm_get_emission_probability(hmm, state, symbol);
      if (cumulative_prob > random_value)
        break;
    }
    if (proc_emission)
      proc_emission(symbol, data);
    /* go to next state */
    random_value = gt_rand_0_to_1();
    cumulative_prob = 0.0;
    for (next_state = 0; next_state < hmm->num_of_states - 1; next_state++) {
      cumulative_prob += gt_hmm_get_transition_probability(hmm,
                                                           state, next_state);
      if (cumulative_prob > random_value)
        break;
    }
    state = next_state;
  }
}

double gt_hmm_rmsd(const GtHMM *hmm_a, const GtHMM* hmm_b)
{
  unsigned int i, j;
  double difference = 0.0, rmsd = 0.0;
  gt_assert(hmm_a && hmm_b);
  gt_assert(hmm_a->num_of_states == hmm_b->num_of_states);
  gt_assert(hmm_a->num_of_symbols == hmm_b->num_of_symbols);
  gt_assert(gt_hmm_is_valid(hmm_a));
  gt_assert(gt_hmm_is_valid(hmm_b));
  /* add transitions probabilities */
  for (i = 0; i < hmm_a->num_of_states; i++) {
    for (j = 0; j < hmm_a->num_of_states; j++) {
      difference = gt_hmm_get_transition_probability(hmm_a, i, j) -
                   gt_hmm_get_transition_probability(hmm_b, i, j);
      rmsd += difference * difference;
    }
  }
  /* add emission probabilities */
  for (i = 0; i < hmm_a->num_of_states; i++) {
    for (j = 0; j < hmm_a->num_of_symbols; j++) {
      difference = gt_hmm_get_emission_probability(hmm_a, i, j) -
                   gt_hmm_get_emission_probability(hmm_b, i, j);
      rmsd += difference * difference;
    }
  }
  return sqrt(rmsd);
}

void gt_hmm_show(const GtHMM *hmm, FILE *fp)
{
  unsigned int i, j;
  gt_assert(hmm && fp);
  fprintf(fp, "# of states: %u\n", hmm->num_of_states);
  fprintf(fp, "# of symbols: %u\n", hmm->num_of_symbols);
  fprintf(fp, "initial state probabilities:\n");
  for (i = 0; i < hmm->num_of_states; i++) {
    fprintf(fp, "%2u: %f", i, gt_hmm_get_initial_state_probability(hmm, i));
  }
  gt_xfputc('\n', fp);
  fprintf(fp, "transition probabilities:\n");
  for (i = 0; i < hmm->num_of_states; i++) {
    fprintf(fp, "%2u:", i);
    for (j = 0; j < hmm->num_of_states; j++) {
      fprintf(fp, " %.2f", gt_hmm_get_transition_probability(hmm, i, j));
    }
    gt_xfputc('\n', fp);
  }
  fprintf(fp, "emission probabilities:\n");
  for (i = 0; i < hmm->num_of_states; i++) {
    fprintf(fp, "%2u:", i);
    for (j = 0; j < hmm->num_of_symbols; j++) {
      fprintf(fp, " %.2f", gt_hmm_get_emission_probability(hmm, i, j));
    }
    gt_xfputc('\n', fp);
  }
}

int gt_hmm_unit_test(GtError *err)
{
  /* the last coin string must be the longest */
  static char *coin_tosses[] = { "H", "T", "HH", "HT", "TH", "TT", "HTHT",
                                 "HHHHHTTTTT", "HTTHTHTHHTHTHHHTHTHTHTHTHHHTH",
                                 "HHHHHHHHHHHHHHHHHHHTTTTTTTTTTTTTTT", "HTTHTH",
                                 "HHTHHTTHTTHTHTHTHTHTTTTTTTHHHTTHHHHHHTHT",
                                 "HHTTTHHTTTHHTHTHTTTHTHHHTHTHHHTHTHHHTHHHTH"
                                 "HHTHTHHHHTTTHTTHHHTTTHTTTHHTHTHTHHTHHTHTHH" };
  /* the last dice string must be the longest */
  static char *dice_rolls[] = { "1", "2", "3", "4", "5", "6", "11", "12", "56",
                                "156246", "165565254154", "66614566161",
                                "12345654321551515144561456131641135452134",
                                "66666666666666666666666666666666666666666666",
                                "3151162464466442453113216311641521336251445436"
                                "3165662656666665116645313265124563666463163666"
                                "3162326455236266666625151631222555441666566563"
                                "5643243641315134651463534111264146262533563661"
                                "6366646623252441366166116325256246225526525226"
                                "6435353336233121625364414432335163243633665562"
                                "466662632666612355245242" };
  unsigned int *encoded_seq;
  GtAlphabet *alpha;
  size_t i, j, len, size;
  GtHMM *fair_hmm, *loaded_hmm;
  int had_err = 0;
  gt_error_check(err);

  /* test the GtHMM class with the coin GtHMMs */
  fair_hmm = gt_coin_hmm_fair();
  loaded_hmm = gt_coin_hmm_loaded();
  alpha = gt_coin_hmm_alphabet();
  size = sizeof (coin_tosses) / sizeof (coin_tosses[0]);
  encoded_seq = gt_malloc(sizeof (int) * strlen(coin_tosses[size-1]));

  for (i = 0; i < size && !had_err; i++) {
    len = strlen(coin_tosses[i]);
    for (j = 0; j < len; j++)
      encoded_seq[j] = gt_alphabet_encode(alpha, coin_tosses[i][j]);
    /* XXX: remove exp() calls */
    gt_ensure(had_err,
           gt_double_equals_double(exp(gt_hmm_forward(fair_hmm, encoded_seq,
                                                      len)),
                                   exp(gt_hmm_backward(fair_hmm, encoded_seq,
                                                    len))));
    gt_ensure(had_err,
           gt_double_equals_double(exp(gt_hmm_forward(loaded_hmm, encoded_seq,
                                                   len)),
                                   exp(gt_hmm_backward(loaded_hmm, encoded_seq,
                                                    len))));
  }

  gt_free(encoded_seq);
  gt_alphabet_delete(alpha);
  gt_ensure(had_err,
            gt_double_equals_double(gt_hmm_rmsd(fair_hmm, fair_hmm), 0.0));
  gt_ensure(had_err,
            gt_double_equals_double(gt_hmm_rmsd(loaded_hmm, loaded_hmm), 0.0));
  gt_hmm_delete(loaded_hmm);
  gt_hmm_delete(fair_hmm);

  /* test the GtHMM class with the dice GtHMMs */
  fair_hmm = gt_dice_hmm_fair();
  loaded_hmm = gt_dice_hmm_loaded();
  alpha = gt_dice_hmm_alphabet();
  size = sizeof (dice_rolls) / sizeof (dice_rolls[0]);
  encoded_seq = gt_malloc(sizeof (int) * strlen(dice_rolls[size-1]));

  for (i = 0; i < size && !had_err; i++) {
    len = strlen(dice_rolls[i]);
    for (j = 0; j < len; j++) {
      encoded_seq[j] = gt_alphabet_encode(alpha, dice_rolls[i][j]);
    }
    /* XXX: remove exp() calls */
    gt_ensure(had_err,
           gt_double_equals_double(exp(gt_hmm_forward(fair_hmm, encoded_seq,
                                                      len)),
                                   exp(gt_hmm_backward(fair_hmm, encoded_seq,
                                                       len))));
    gt_ensure(had_err,
           gt_double_equals_double(exp(gt_hmm_forward(loaded_hmm, encoded_seq,
                                                      len)),
                                   exp(gt_hmm_backward(loaded_hmm, encoded_seq,
                                                       len))));
  }

  gt_free(encoded_seq);
  gt_alphabet_delete(alpha);
  gt_ensure(had_err,
            gt_double_equals_double(gt_hmm_rmsd(fair_hmm, fair_hmm), 0.0));
  gt_ensure(had_err,
            gt_double_equals_double(gt_hmm_rmsd(loaded_hmm, loaded_hmm), 0.0));
  gt_hmm_delete(loaded_hmm);
  gt_hmm_delete(fair_hmm);

  return had_err;
}

void gt_hmm_delete(GtHMM *hmm)
{
  if (!hmm) return;
  gt_free(hmm->initial_state_prob);
  gt_array2dim_delete(hmm->transition_prob);
  gt_array2dim_delete(hmm->emission_prob);
  gt_free(hmm);
}
