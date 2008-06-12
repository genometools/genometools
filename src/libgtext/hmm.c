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
#include <assert.h>
#include <math.h>
#include "libgtcore/array2dim.h"
#include "libgtcore/ensure.h"
#include "libgtcore/mathsupport.h"
#include "libgtcore/xansi.h"
#include "libgtext/coin_hmm.h"
#include "libgtext/dice_hmm.h"
#include "libgtext/hmm.h"

#define MINUSINFINITY   -99999.0
#define PSEUDOCOUNT     1

struct HMM {
  unsigned int num_of_states,
               num_of_symbols;
  double *initial_state_prob, /* log values */
         **transition_prob,   /* log values */
         **emission_prob;     /* log values */
};

HMM* hmm_new(unsigned int num_of_states, unsigned int num_of_symbols)
{
  HMM *hmm;
  unsigned int i, j;

  assert(num_of_states && num_of_symbols);

  /* alloc */
  hmm = ma_malloc(sizeof *hmm);
  hmm->initial_state_prob = ma_malloc(sizeof (double) * num_of_states);
  array2dim_malloc(hmm->transition_prob, num_of_states, num_of_states);
  array2dim_malloc(hmm->emission_prob, num_of_states, num_of_symbols);

  /* init */
  hmm->num_of_states = num_of_states;
  hmm->num_of_symbols = num_of_symbols;

  for (i = 0; i < num_of_states; i++)
    hmm_set_initial_state_probability(hmm, i, 1.0 / num_of_states);

  for (i = 0; i < num_of_states; i++) {
    for (j = 0; j < num_of_states; j++)
      hmm_set_transition_probability(hmm, i, j, 0.0);
  }

  for (i = 0; i < num_of_states; i++) {
    for (j = 0; j < num_of_symbols; j++)
      hmm_set_emission_probability(hmm, i, j, 0.0);
  }

  return hmm;
}

void hmm_set_initial_state_probability(HMM *hmm, unsigned int state_num,
                                       double probability)
{
  assert(hmm);
  assert(probability >= 0.0 && probability <= 1.0);
  assert(state_num < hmm->num_of_states);
  if (probability == 0.0)
    hmm->initial_state_prob[state_num] = MINUSINFINITY;
  else
    hmm->initial_state_prob[state_num] = log(probability);
}

double hmm_get_initial_state_probability(const HMM *hmm, unsigned int state_num)
{
  assert(hmm);
  assert(state_num < hmm->num_of_states);
  if (hmm->initial_state_prob[state_num] == MINUSINFINITY)
    return 0.0;
  return exp(hmm->initial_state_prob[state_num]);
}

void hmm_set_transition_probability(HMM *hmm,
                                    unsigned int from_state_num,
                                    unsigned int to_state_num,
                                    double probability)
{
  assert(hmm);
  assert(probability >= 0.0 && probability <= 1.0);
  assert(from_state_num < hmm->num_of_states);
  assert(to_state_num < hmm->num_of_states);
  if (probability == 0.0)
    hmm->transition_prob[from_state_num][to_state_num] = MINUSINFINITY;
  else
    hmm->transition_prob[from_state_num][to_state_num] = log(probability);
}

double hmm_get_transition_probability(const HMM *hmm,
                                      unsigned int from_state_num,
                                      unsigned int to_state_num)
{
  assert(hmm);
  assert(from_state_num < hmm->num_of_states);
  assert(to_state_num < hmm->num_of_states);
  if (hmm->transition_prob[from_state_num][to_state_num] == MINUSINFINITY)
    return 0.0;
  return exp(hmm->transition_prob[from_state_num][to_state_num]);
}

void hmm_set_missing_transition_probabilities(HMM *hmm)
{
  unsigned int row, column, num_of_missing_entries;
  double prob, sum_of_probabilities;

  assert(hmm);

  for (row = 0; row < hmm->num_of_states; row++) {
    /* determine sum of probabilities */
    sum_of_probabilities = 0.0;
    num_of_missing_entries = 0;
    for (column = 0; column < hmm->num_of_states; column++) {
      prob = hmm_get_transition_probability(hmm, row, column);
      if (prob == 0.0)
        num_of_missing_entries++;
      else
        sum_of_probabilities += prob;
    }

    /* set missing probabilites (equal distribution) */
    if (num_of_missing_entries) {
      for (column = 0; column < hmm->num_of_states; column++) {
        prob = hmm_get_transition_probability(hmm, row, column);
        if (prob == 0.0) {
          hmm_set_transition_probability(hmm, row, column,
                                         (1.0 - sum_of_probabilities) /
                                         num_of_missing_entries);
        }
      }
    }
  }
}

void hmm_set_emission_probability(HMM *hmm,
                                  unsigned int state_num,
                                  unsigned int symbol_num,
                                  double probability)
{
  assert(hmm);
  assert(probability >= 0.0 && probability <= 1.0);
  assert(state_num < hmm->num_of_states);
  assert(symbol_num < hmm->num_of_symbols);
  if (probability == 0.0)
    hmm->emission_prob[state_num][symbol_num] = MINUSINFINITY;
  else
    hmm->emission_prob[state_num][symbol_num] = log(probability);
}

double hmm_get_emission_probability(const HMM *hmm,
                                    unsigned int state_num,
                                    unsigned int symbol_num)
{
  assert(hmm);
  assert(state_num < hmm->num_of_states);
  assert(symbol_num < hmm->num_of_symbols);
  if (hmm->emission_prob[state_num][symbol_num] == MINUSINFINITY)
    return 0.0;
  return exp(hmm->emission_prob[state_num][symbol_num]);
}

static bool hmm_has_valid_initial_state_probs(const HMM *hmm)
{
  unsigned int state;
  double sum_of_probabilities = 0.0;

  assert(hmm);

  for (state = 0; state < hmm->num_of_states; state++)
    sum_of_probabilities += hmm_get_initial_state_probability(hmm, state);

  return double_equals_one(sum_of_probabilities);
}

static bool hmm_has_valid_emissions(const HMM *hmm)
{
  unsigned int state, symbol;
  double sum_of_probabilities;

  assert(hmm);

  for (state = 0; state < hmm->num_of_states; state++) {
    sum_of_probabilities = 0.0;
    for (symbol = 0; symbol < hmm->num_of_symbols; symbol++)
      sum_of_probabilities += hmm_get_emission_probability(hmm, state, symbol);
    if (!double_equals_one(sum_of_probabilities))
      return false;
  }
  return true;
}

static bool hmm_has_valid_states(const HMM *hmm)
{
  unsigned int state_1, state_2;
  double sum_of_probabilities;

  assert(hmm);

  for (state_1 = 0; state_1 < hmm->num_of_states; state_1++) {
    sum_of_probabilities = 0.0;
    for (state_2 = 0; state_2 < hmm->num_of_states; state_2++) {
      sum_of_probabilities += hmm_get_transition_probability(hmm,state_1,
                                                             state_2);

    }
    if (!double_equals_one(sum_of_probabilities))
      return false;
  }

  return true;
}

bool hmm_is_valid(const HMM *hmm)
{
  assert(hmm);
  return (hmm_has_valid_initial_state_probs(hmm) &&
          hmm_has_valid_emissions(hmm) &&
          hmm_has_valid_states(hmm));

}

void hmm_init_random(HMM *hmm)
{
  double random_value, cumulative_prob;
  unsigned int i, j;
  assert(hmm);

  /* initialize initial state probabilities in random fashion */
  cumulative_prob = 0.0;
  for (i = 0; i < hmm->num_of_states - 1; i++) {
    random_value = rand_max_double(1.0 - cumulative_prob);
    hmm_set_initial_state_probability(hmm, i, random_value);
    cumulative_prob += random_value;
  }
  assert(cumulative_prob <= 1.0);
  hmm_set_initial_state_probability(hmm, i,  1.0 - cumulative_prob);

  /* initialize transition probabilities in random fashion */
  for (i = 0; i < hmm->num_of_states; i++) {
    cumulative_prob = 0.0;
    for (j = 0; j < hmm->num_of_states - 1; j++) {
      random_value = rand_max_double(1.0 - cumulative_prob);
      hmm_set_transition_probability(hmm, i, j, random_value);
      cumulative_prob += random_value;
    }
    assert(cumulative_prob <= 1.0);
    hmm_set_transition_probability(hmm, i, j, 1.0 - cumulative_prob);
  }

  /* initialize emission probabilities in random fashion */
  for (i = 0; i < hmm->num_of_states; i++) {
    cumulative_prob = 0.0;
    for (j = 0; j < hmm->num_of_symbols - 1; j++) {
      random_value = rand_max_double(1.0 - cumulative_prob);
      hmm_set_emission_probability(hmm, i, j, random_value);
      cumulative_prob += random_value;
    }
    assert(cumulative_prob <= 1.0);
    hmm_set_emission_probability(hmm, i, j, 1.0 - cumulative_prob);
  }

  assert(hmm_is_valid(hmm));
}

/* [DEKM98, p. 56] */
void hmm_decode(const HMM *hmm,
                unsigned int *state_sequence,
                const unsigned int *emissions,
                unsigned int num_of_emissions)
{
  double **max_probabilities, tmp_prob;
  unsigned int **backtrace, colidx, precolidx;
  int row, column, num_of_rows, num_of_columns, previous_row;

  assert(hmm);
  assert(hmm_is_valid(hmm));
  assert(num_of_emissions);

  /* alloc tables */
  num_of_rows = hmm->num_of_states;
  num_of_columns = num_of_emissions;
  array2dim_malloc(max_probabilities, num_of_rows, 2);
  array2dim_malloc(backtrace, num_of_rows, num_of_columns);

  /* fill DP table */
  for (row = 0; row < num_of_rows; row++) { /* first column */
    assert(emissions[0] < hmm->num_of_symbols);
    max_probabilities[row][0] = hmm->initial_state_prob[row] +
                                hmm->emission_prob[row][emissions[0]];
    backtrace[row][0] = row;
  }

  for (column = 1; column < num_of_columns; column++) { /* other columns */
    assert(emissions[column] < hmm->num_of_symbols);
    colidx = column & 1;
    precolidx = (column - 1) & 1;
    for (row = 0; row < num_of_rows; row++) {
      max_probabilities[row][colidx] = max_probabilities[0][precolidx] +
                                     hmm->transition_prob[0][row] +
                                     hmm->emission_prob[row][emissions[column]];
      backtrace[row][column] = 0;
      for (previous_row = 1; previous_row < num_of_rows; previous_row++) {
        tmp_prob = max_probabilities[previous_row][precolidx] +
                   hmm->transition_prob[previous_row][row] +
                   hmm->emission_prob[row][emissions[column]];
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
  array2dim_delete(backtrace);
  array2dim_delete(max_probabilities);
}

/* [DEKM98, p. 58] */
static void compute_forward_table(double **f, const HMM *hmm,
                                  const unsigned int *emissions,
                                  unsigned long num_of_emissions)
{
  int row, column, previous_row;
  double tmp_prob;

  assert(f && hmm && emissions && num_of_emissions);

  /* fill DP table 'f' (hmm->num_of_states * num_of_emissions) */
  for (row = 0; row < hmm->num_of_states; row++) { /* first column */
    assert(emissions[0] < hmm->num_of_symbols);
    f[row][0] = hmm->initial_state_prob[row] +
                hmm->emission_prob[row][emissions[0]];
  }

  for (column = 1; column < num_of_emissions; column++) { /* other columns */
    assert(emissions[column] < hmm->num_of_symbols);
    for (row = 0; row < hmm->num_of_states; row++) {
      f[row][column] = hmm->emission_prob[row][emissions[column]];
      tmp_prob = f[0][column-1] + hmm->transition_prob[0][row];
      for (previous_row = 1; previous_row < hmm->num_of_states;
           previous_row++) {
        /* XXX: replace the logsum() call with a tabulated version */
        tmp_prob = logsum(tmp_prob, f[previous_row][column-1] +
                                    hmm->transition_prob[previous_row][row]);
      }
      f[row][column] += tmp_prob;
    }
  }
}

/* [DEKM98, p. 58] */
double hmm_forward(const HMM* hmm, const unsigned int *emissions,
                   unsigned int num_of_emissions)
{
  unsigned int i;
  double **f, P;

  assert(hmm && emissions && num_of_emissions);
  array2dim_malloc(f, hmm->num_of_states, num_of_emissions);

  /* XXX: we do not need the full table here, the last column would suffice */
  compute_forward_table(f, hmm, emissions, num_of_emissions);

  /* compute P(x) */
  P = f[0][num_of_emissions-1];
  for (i = 1; i < hmm->num_of_states; i++) {
    /* XXX: replace the logsum() call with a tabulated version */
    P = logsum(P, f[i][num_of_emissions-1]);
  }

  array2dim_delete(f);
  return P;
}

/* [DEKM98, p. 59] */
static void compute_backward_table(double **b, const HMM *hmm,
                                   const unsigned int *emissions,
                                   unsigned long num_of_emissions)
{
  int row, column, next_row;
  double tmp_prob;

  assert(b && hmm && emissions && num_of_emissions);

  /* fill DP table 'b' (hmm->num_of_states * num_of_emissions) */
  for (row = 0; row < hmm->num_of_states; row++) { /* last column */
    b[row][num_of_emissions-1] = 0.0; /* probability = 1.0 */
  }

  for (column = num_of_emissions - 2; column >= 0; column--) { /* other cols */
    assert(emissions[column] < hmm->num_of_symbols);
    for (row = 0; row < hmm->num_of_states; row++) {
      tmp_prob = hmm->transition_prob[row][0] +
                 hmm->emission_prob[0][emissions[column+1]] +
                 b[0][column+1];
      for (next_row = 1; next_row < hmm->num_of_states; next_row++) {
        /* XXX: replace the logsum() call with a tabulated version */
        tmp_prob = logsum(tmp_prob, hmm->transition_prob[row][next_row] +
                          hmm->emission_prob[next_row][emissions[column+1]] +
                          b[next_row][column+1]);
      }
      b[row][column] = tmp_prob;
    }
  }
}

/* [DEKM98, p. 59] */
double hmm_backward(const HMM* hmm, const unsigned int *emissions,
                    unsigned int num_of_emissions)
{
  unsigned int i;
  double **b, P;

  assert(hmm && emissions && num_of_emissions);
  array2dim_malloc(b, hmm->num_of_states, num_of_emissions);

  /* XXX: we do not need the full table here, the last column would suffice */
  compute_backward_table(b, hmm, emissions, num_of_emissions);

  /* compute P(x) */
  P = hmm->initial_state_prob[0] + hmm->emission_prob[0][emissions[0]] +
      b[0][0];
  for (i = 1; i < hmm->num_of_states; i++) {
    /* XXX: replace the logsum() call with a tabulated version */
    P = logsum(P, hmm->initial_state_prob[i] +
                  hmm->emission_prob[i][emissions[0]] + b[i][0]);
  }

  array2dim_delete(b);
  return P;
}

void hmm_emit(HMM *hmm, unsigned long num_of_emissions,
              void (*proc_emission)(unsigned int symbol, void *data),
              void *data)
{
  unsigned long i;
  unsigned int state, symbol, next_state;
  double random_value, cumulative_prob;
  assert(hmm);

  /* determine initial state */
  random_value = rand_0_to_1();
  cumulative_prob = 0.0;
  for (state = 0; state < hmm->num_of_states - 1; state++) {
    cumulative_prob += hmm_get_initial_state_probability(hmm, state);
    if (cumulative_prob > random_value)
      break;
  }

  /* emit loop */
  for (i = 0; i < num_of_emissions; i++) {
    /* emit character */
    random_value = rand_0_to_1();
    cumulative_prob = 0.0;
    for (symbol = 0; symbol < hmm->num_of_symbols - 1; symbol++) {
      cumulative_prob += hmm_get_emission_probability(hmm, state, symbol);
      if (cumulative_prob > random_value)
        break;
    }
    if (proc_emission)
      proc_emission(symbol, data);
    /* go to next state */
    random_value = rand_0_to_1();
    cumulative_prob = 0.0;
    for (next_state = 0; next_state < hmm->num_of_states - 1; next_state++) {
      cumulative_prob += hmm_get_transition_probability(hmm, state, next_state);
      if (cumulative_prob > random_value)
        break;
    }
    state = next_state;
  }
}

double hmm_rmsd(const HMM *hmm_a, const HMM* hmm_b)
{
  unsigned int i, j;
  double difference = 0.0, rmsd = 0.0;
  assert(hmm_a && hmm_b);
  assert(hmm_a->num_of_states == hmm_b->num_of_states);
  assert(hmm_a->num_of_symbols == hmm_b->num_of_symbols);
  assert(hmm_is_valid(hmm_a));
  assert(hmm_is_valid(hmm_b));
  /* add transitions probabilities */
  for (i = 0; i < hmm_a->num_of_states; i++) {
    for (j = 0; j < hmm_a->num_of_states; j++) {
      difference = hmm_get_transition_probability(hmm_a, i, j) -
                   hmm_get_transition_probability(hmm_b, i, j);
      rmsd += difference * difference;
    }
  }
  /* add emission probabilities */
  for (i = 0; i < hmm_a->num_of_states; i++) {
    for (j = 0; j < hmm_a->num_of_symbols; j++) {
      difference = hmm_get_emission_probability(hmm_a, i, j) -
                   hmm_get_emission_probability(hmm_b, i, j);
      rmsd += difference * difference;
    }
  }
  return sqrt(rmsd);
}

void hmm_show(const HMM *hmm, FILE *fp)
{
  unsigned int i, j;
  assert(hmm && fp);
  fprintf(fp, "# of states: %u\n", hmm->num_of_states);
  fprintf(fp, "# of symbols: %u\n", hmm->num_of_symbols);
  fprintf(fp, "initial state probabilities:\n");
  for (i = 0; i < hmm->num_of_states; i++) {
    fprintf(fp, "%2u: %f", i, hmm_get_initial_state_probability(hmm, i));
  }
  xfputc('\n', fp);
  fprintf(fp, "transition probabilities:\n");
  for (i = 0; i < hmm->num_of_states; i++) {
    fprintf(fp, "%2u:", i);
    for (j = 0; j < hmm->num_of_states; j++) {
      fprintf(fp, " %.2f", hmm_get_transition_probability(hmm, i, j));
    }
    xfputc('\n', fp);
  }
  fprintf(fp, "emission probabilities:\n");
  for (i = 0; i < hmm->num_of_states; i++) {
    fprintf(fp, "%2u:", i);
    for (j = 0; j < hmm->num_of_symbols; j++) {
      fprintf(fp, " %.2f", hmm_get_emission_probability(hmm, i, j));
    }
    xfputc('\n', fp);
  }
}

int hmm_unit_test(Error *err)
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
  Alpha *alpha;
  size_t i, j, len, size;
  HMM *fair_hmm, *loaded_hmm;
  int had_err = 0;
  error_check(err);

  /* test the HMM class with the coin HMMs */
  fair_hmm = coin_hmm_fair();
  loaded_hmm = coin_hmm_loaded();
  alpha = coin_hmm_alpha();
  size = sizeof (coin_tosses) / sizeof (coin_tosses[0]);
  encoded_seq = ma_malloc(sizeof (int) * strlen(coin_tosses[size-1]));

  for (i = 0; i < size && !had_err; i++) {
    len = strlen(coin_tosses[i]);
    for (j = 0; j < len; j++)
      encoded_seq[j] = alpha_encode(alpha, coin_tosses[i][j]);
    /* XXX: remove exp() calls */
    ensure(had_err,
           double_equals_double(exp(hmm_forward(fair_hmm, encoded_seq, len)),
                                exp(hmm_backward(fair_hmm, encoded_seq, len))));
    ensure(had_err,
           double_equals_double(exp(hmm_forward(loaded_hmm, encoded_seq, len)),
                                exp(hmm_backward(loaded_hmm, encoded_seq, len)))
                               );
  }

  ma_free(encoded_seq);
  alpha_delete(alpha);
  ensure(had_err, double_equals_double(hmm_rmsd(fair_hmm, fair_hmm), 0.0));
  ensure(had_err, double_equals_double(hmm_rmsd(loaded_hmm, loaded_hmm), 0.0));
  hmm_delete(loaded_hmm);
  hmm_delete(fair_hmm);

  /* test the HMM class with the dice HMMs */
  fair_hmm = dice_hmm_fair();
  loaded_hmm = dice_hmm_loaded();
  alpha = dice_hmm_alpha();
  size = sizeof (dice_rolls) / sizeof (dice_rolls[0]);
  encoded_seq = ma_malloc(sizeof (int) * strlen(dice_rolls[size-1]));

  for (i = 0; i < size && !had_err; i++) {
    len = strlen(dice_rolls[i]);
    for (j = 0; j < len; j++) {
      encoded_seq[j] = alpha_encode(alpha, dice_rolls[i][j]);
    }
    /* XXX: remove exp() calls */
    ensure(had_err,
           double_equals_double(exp(hmm_forward(fair_hmm, encoded_seq, len)),
                                exp(hmm_backward(fair_hmm, encoded_seq, len))));
    ensure(had_err,
           double_equals_double(exp(hmm_forward(loaded_hmm, encoded_seq, len)),
                                exp(hmm_backward(loaded_hmm, encoded_seq, len)))
                               );
  }

  ma_free(encoded_seq);
  alpha_delete(alpha);
  ensure(had_err, double_equals_double(hmm_rmsd(fair_hmm, fair_hmm), 0.0));
  ensure(had_err, double_equals_double(hmm_rmsd(loaded_hmm, loaded_hmm), 0.0));
  hmm_delete(loaded_hmm);
  hmm_delete(fair_hmm);

  return had_err;
}

void hmm_delete(HMM *hmm)
{
  if (!hmm) return;
  ma_free(hmm->initial_state_prob);
  array2dim_delete(hmm->transition_prob);
  array2dim_delete(hmm->emission_prob);
  ma_free(hmm);
}
