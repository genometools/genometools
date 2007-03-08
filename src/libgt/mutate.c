/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include <string.h>
#include <libgt/mathsupport.h>
#include <libgt/mutate.h>

#define MUTATED_DESC_PRIMER " [mutated with rate "

static char* mutate_description(const char *description, unsigned int rate,
                                Env *env)
{
  unsigned long mutated_description_len;
  char *mutated_description;
  int rval;
  env_error_check(env);
  assert(description);
  assert(rate >= 0 && rate <= 100);
  mutated_description_len = strlen(description) + strlen(MUTATED_DESC_PRIMER)
                            + 3  /* for the rate */
                            + 1  /* terminal ']' */
                            + 1; /* terminal '\n' */
  mutated_description = env_ma_malloc(env,
                                      sizeof (char) * mutated_description_len);
  rval = snprintf(mutated_description, mutated_description_len, "%s%s%u]",
                  description, MUTATED_DESC_PRIMER, rate);
  assert(rval < mutated_description_len);
  return mutated_description;
}

static char random_character(Alpha *alpha)
{
  /* we do not want to get wildcard characters */
  return alpha_decode(alpha, rand_max(alpha_size(alpha) - 1 - 1));
}

static char* mutate_seq(const char *seq, unsigned long len, Alpha *alpha,
                        unsigned int rate, Env *env)
{
  unsigned long i, j, allocated, substitution_events = 0, insertion_events = 0,
                deletion_events = 0, total_events = 0;
  unsigned int cc;
  double rand_prob, mutate_prob;
  char *mutated_seq;
  env_error_check(env);
  assert(seq && alpha);
  assert(rate >= 0 && rate <= 100);
  mutate_prob = (double) rate / 100.0;
  allocated = len * 2; /* XXX: possibly reduce this memory consumption */
  mutated_seq = env_ma_malloc(env, sizeof (char) * allocated);
  for (i = 0, j = 0; i < len; i++) {
    cc = alpha_encode(alpha, seq[i]);
    if (rand_0_to_1() <= mutate_prob) {
      /* mutate */
      rand_prob = rand_0_to_1();
      if (rand_prob <= 0.8) {
        /* substitution (80% probability) */
        mutated_seq[j++] = random_character(alpha); /* add random character */
        substitution_events++;
      }
      else if (rand_prob <= 0.9) {
        /* insertion (10% probability) */
        mutated_seq[j++] = alpha_decode(alpha, cc); /* keep orig. character */
        mutated_seq[j++] = random_character(alpha); /* add random character */
        insertion_events++;
      }
      else {
      /* deletion (10% probability) */
        deletion_events++;
      }
      total_events++;
    }
    else
      mutated_seq[j++] = alpha_decode(alpha, cc); /* keep original character */
  }
  mutated_seq[j] = '\0'; /* terminate */
  env_log_log(env, "total number of mutation events: %lu", total_events);
  env_log_log(env, "number of substitution events: %lu", substitution_events);
  env_log_log(env, "number of insertion events: %lu", insertion_events);
  env_log_log(env, "number of deletion events: %lu", deletion_events);
  return mutated_seq;
}

Seq* mutate(const char *description, const char *orig_seq, unsigned long len,
            Alpha *alpha, unsigned int rate, Env *env)
{
  char *mutated_description, *mutated_seq;
  Seq *seq;
  env_error_check(env);
  assert(description && orig_seq && alpha);
  assert(rate >= 0 && rate <= 100);
  mutated_description = mutate_description(description, rate, env);
  mutated_seq = mutate_seq(orig_seq, len, alpha, rate, env);
  seq = seq_new_own(mutated_seq, strlen(mutated_seq), alpha, env);
  seq_set_description_own(seq, mutated_description, env);
  return seq;
}
